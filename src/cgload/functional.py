"""Attaching functional annotation to loaded features (milestone 7).

The whole difficulty is the join. eggNOG and InterProScan are given a protein
FASTA, and whatever was in those headers is what comes back — which may be the
transcript ID, the protein ID, the locus tag, or any of those with a suffix the
user's pipeline added. Nothing in the functional file says which.

So the match is a lookup across several candidate keys, and the **match rate is
reported and enforced** rather than assumed. A file where 40% of queries matched
nothing has told you something important, and the default is to stop.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy import Connection, Engine, insert, select
from sqlalchemy.exc import IntegrityError

from cgload.db.init import assert_vocabulary_current
from cgload.db.schema import (
    annotation_run,
    assembly,
    feature,
    feature_attribute,
    protein_annotation,
)
from cgload.parsers.functional import FunctionalFile, parse_functional
from cgload.storage import sha256_of

#: Attribute keys that may carry the identifier a protein FASTA was headed with,
#: in the order they are tried. `source_id` is tried first and is not in this
#: list because it is a column rather than an attribute.
#:
#: Order matters only for reporting which key matched; a query matching two
#: different features through two different keys is an error, not a preference.
LOOKUP_ATTRIBUTES: tuple[str, ...] = (
    "protein_id",
    "transcript_id",
    "orig_protein_id",
    "orig_transcript_id",
    "locus_tag",
    "Name",
    "gene",
)

#: Feature types a functional annotation may attach to. A protein annotation
#: belongs on the coding sequence or its transcript, never on the gene: a gene
#: with several isoforms has several proteins, and attaching to the gene would
#: make them indistinguishable.
ANNOTATABLE_TYPES = frozenset({"CDS", "mRNA"})

#: Rows per INSERT. Bounded so that a genome-scale functional file does not build
#: its entire row set in memory inside the load transaction (D-044).
INSERT_BATCH_SIZE = 5_000


class FunctionalLoadError(RuntimeError):
    """A functional file that cannot be attached without guessing."""


@dataclass
class FunctionalReport:
    path: Path
    tool: str
    tool_version: str | None
    query_count: int
    matched_count: int
    unmatched_count: int
    hit_count: int
    matched_by_key: dict[str, int] = field(default_factory=dict)
    unmatched_examples: list[str] = field(default_factory=list)
    analyses: dict[str, int] = field(default_factory=dict)

    @property
    def match_rate(self) -> float:
        return self.matched_count / self.query_count if self.query_count else 0.0


def load_functional(
    engine: Engine,
    assembly_id: int,
    path: Path,
    *,
    tool: str | None = None,
    minimum_match_rate: float = 1.0,
) -> FunctionalReport:
    """Attach a functional file to an assembly already in the database.

    Used to add annotation after the fact. During a fresh ``load`` the work
    happens in ``attach_functional`` inside the load's own transaction, so a
    failure here cannot leave a genome loaded with its annotation half-attached.
    """
    assert_vocabulary_current(engine)
    with engine.begin() as connection:
        if connection.execute(
            select(assembly.c.assembly_id).where(assembly.c.assembly_id == assembly_id)
        ).scalar_one_or_none() is None:
            raise FunctionalLoadError(
                f"no assembly with id {assembly_id}. Load the genome first, then "
                f"attach its functional annotation."
            )
        return _attach(
            connection,
            assembly_id,
            Path(path),
            tool=tool,
            minimum_match_rate=minimum_match_rate,
        )


def attach_functional(connection: Connection, assembly_id: int, entry) -> FunctionalReport:
    """Attach one file inside an open transaction. Called by the loader."""
    return _attach(
        connection,
        assembly_id,
        entry.path,
        tool=entry.tool,
        minimum_match_rate=entry.minimum_match_rate,
    )


def _attach(
    connection: Connection,
    assembly_id: int,
    path: Path,
    *,
    tool: str | None,
    minimum_match_rate: float,
) -> FunctionalReport:
    """``minimum_match_rate`` defaults to 1.0 at every caller: every query must
    match a feature. Deliberately strict, because a functional file that
    half-matches is the classic silent partial success -- the load reports
    success, the database looks annotated, and the missing half is unrecoverable
    (D-040).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"functional annotation not found: {path}")

    parsed = parse_functional(path, tool=tool)
    index = _build_index(connection, assembly_id)
    resolved, report = _resolve(parsed, index, path)

    if report.match_rate < minimum_match_rate:
        examples = ", ".join(report.unmatched_examples[:5])
        raise FunctionalLoadError(
            f"{path.name}: {report.unmatched_count} of {report.query_count} query "
            f"proteins ({100 * (1 - report.match_rate):.1f}%) match no feature in "
            f"assembly {assembly_id}. First unmatched: {examples}. Either the file "
            f"was produced from a different annotation, or the FASTA headers were "
            f"renamed; pass a lower minimum match rate to load the matching subset "
            f"deliberately."
        )

    try:
        run_id = connection.execute(
            insert(annotation_run).values(
                assembly_id=assembly_id,
                tool=parsed.tool,
                tool_version=parsed.tool_version,
                path=str(path),
                sha256=sha256_of(path),
                query_count=report.query_count,
                matched_count=report.matched_count,
                unmatched_count=report.unmatched_count,
                loaded_at=datetime.now(UTC),
            )
        ).inserted_primary_key[0]
    except IntegrityError as exc:
        raise FunctionalLoadError(
            f"{path.name}: this file is already attached to assembly "
            f"{assembly_id}. Underlying error: {exc.orig}"
        ) from exc

    # Deduplicated and inserted in batches rather than materialising every row
    # at once. Measured on real emapper-2.1.4 output: 57 stored findings per
    # protein, 90 of them GO terms. A 13,000-protein fungal genome therefore
    # yields around three-quarters of a million rows from one file, and building
    # that as a single list plus a single dictionary inside an already-open
    # transaction holding the whole genome was the load's memory ceiling (D-044).
    seen: set[tuple[int, str, str, int | None]] = set()
    batch: list[dict[str, object]] = []
    stored = 0

    def flush() -> None:
        nonlocal stored, batch
        if batch:
            connection.execute(insert(protein_annotation), batch)
            stored += len(batch)
            batch = []

    for hit, feature_id in resolved:
        # The same accession may legitimately appear twice for one protein --
        # InterProScan reports repeated domains at different positions, and both
        # are real. Exact duplicates are collapsed because eggNOG lists a term
        # once per source column, so one fact arrives from two of them.
        key = (feature_id, hit.analysis, hit.accession, hit.protein_start)
        if key in seen:
            continue
        seen.add(key)
        batch.append(
            {
                "annotation_run_id": run_id,
                "feature_id": feature_id,
                "query_id": hit.query_id,
                "analysis": hit.analysis,
                "accession": hit.accession,
                "description": hit.description,
                "protein_start": hit.protein_start,
                "protein_end": hit.protein_end,
                "score": hit.score,
                "evalue": hit.evalue,
            }
        )
        if len(batch) >= INSERT_BATCH_SIZE:
            flush()
    flush()
    report.hit_count = stored
    return report


def _build_index(connection: Connection, assembly_id: int) -> dict[str, set[int]]:
    """Every identifier that could name a protein, mapped to feature rows.

    Built once per load rather than queried per hit: a eukaryotic genome has tens
    of thousands of proteins and an InterProScan run has hundreds of thousands of
    hits, so a per-hit query would dominate the load.
    """
    index: dict[str, set[int]] = defaultdict(set)

    rows = connection.execute(
        select(feature.c.feature_id, feature.c.source_id, feature.c.feature_type).where(
            feature.c.assembly_id == assembly_id,
            feature.c.segment_index == 0,
            feature.c.feature_type.in_(ANNOTATABLE_TYPES),
        )
    ).all()
    identifiers = {row.feature_id: row.source_id for row in rows}
    for row in rows:
        if row.source_id:
            index[row.source_id].add(row.feature_id)

    if not identifiers:
        return index

    attributes = connection.execute(
        select(
            feature_attribute.c.feature_id,
            feature_attribute.c.key,
            feature_attribute.c.value,
        ).where(
            feature_attribute.c.feature_id.in_(identifiers),
            feature_attribute.c.key.in_(LOOKUP_ATTRIBUTES),
        )
    ).all()
    for row in attributes:
        value = row.value.strip()
        if not value:
            continue
        index[value].add(row.feature_id)
        # RefSeq writes `gnl|WGS:AACM|FGSG_11579T0` in orig_protein_id, and a
        # pipeline that fed the bare id to eggNOG would otherwise never match.
        if "|" in value:
            index[value.rsplit("|", 1)[-1]].add(row.feature_id)

    return index


def _resolve(
    parsed: FunctionalFile, index: dict[str, set[int]], path: Path
) -> tuple[list[tuple[object, int]], FunctionalReport]:
    """Match each query to exactly one feature, or record it as unmatched."""
    report = FunctionalReport(
        path=path,
        tool=parsed.tool,
        tool_version=parsed.tool_version,
        query_count=len(parsed.queries),
        matched_count=0,
        unmatched_count=0,
        hit_count=0,
    )

    matches: dict[str, int] = {}
    for query in sorted(parsed.queries):
        candidates = _candidates(query, index)
        if not candidates:
            report.unmatched_count += 1
            if len(report.unmatched_examples) < 20:
                report.unmatched_examples.append(query)
            continue
        if len(candidates) > 1:
            raise FunctionalLoadError(
                f"{path.name}: query {query!r} matches {len(candidates)} features. "
                f"cgload will not choose one, because attaching a protein's "
                f"annotation to the wrong isoform is unrecoverable and invisible. "
                f"Disambiguate the FASTA headers the tool was given."
            )
        matches[query] = next(iter(candidates))
        report.matched_count += 1

    resolved: list[tuple[object, int]] = []
    analyses: dict[str, int] = defaultdict(int)
    for hit in parsed.hits:
        feature_id = matches.get(hit.query_id)
        if feature_id is None:
            continue
        resolved.append((hit, feature_id))
        analyses[hit.analysis] += 1
    report.analyses = dict(analyses)

    return resolved, report


def _candidates(query: str, index: dict[str, set[int]]) -> set[int]:
    """Try the query as given, then with common pipeline decorations removed.

    Only exact keys are used. Nothing here does prefix matching or similarity.
    Pairing proteins to transcripts by how much their product names resemble each
    other silently mispairs some and drops the rest, and a load that drops
    annotation while reporting success is the failure this project exists to
    prevent.
    """
    if query in index:
        return index[query]

    # A protein FASTA header is often `ID description`, and some tools keep the
    # whole line. Splitting on whitespace is not a similarity match: it takes the
    # first field, which is the identifier by FASTA convention.
    first_field = query.split()[0] if query.split() else query
    if first_field != query and first_field in index:
        return index[first_field]

    if "|" in first_field:
        tail = first_field.rsplit("|", 1)[-1]
        if tail in index:
            return index[tail]

    return set()
