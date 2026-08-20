"""The loader: assembly plus annotation, in one transaction (milestone 4).

Everything the parsers produced becomes rows here, and nothing else in the
project writes to the ``feature`` table. The transaction spans the database and
the data directory (D-019), so a failure at any point — a bad file, a
constraint violation, a full disk — leaves both exactly as they were.

Read this alongside ``cgload.parsers.normalise``: that module decides what is
true about a file, this one decides what is stored. Keeping the two apart is
what lets the verification of D-008 be written against the file rather than
against the loader's own beliefs.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy import Connection, Engine, bindparam, delete, insert, select
from sqlalchemy.exc import IntegrityError

from cgload import __version__
from cgload.assembly import (
    ASSEMBLY_FILENAME,
    AssemblyExistsError,
    _upsert_organism,
    _upsert_strain,
)
from cgload.config import AnnotationFile, OrganismConfig
from cgload.db import vocabulary as vocab
from cgload.db.init import assert_vocabulary_current
from cgload.db.schema import (
    assembly,
    feature,
    feature_attribute,
    sequence_region,
    source_file,
)
from cgload.db.verification import compare_with_database, count_database_for_file
from cgload.envelope import apply_policy, format_report
from cgload.fasta import index_and_read
from cgload.functional import FunctionalReport, attach_functional
from cgload.parsers import profiles as dialects
from cgload.parsers.genbank import GenBankError, tokenize_genbank
from cgload.parsers.normalise import (
    NormalisedFile,
    NormaliseError,
    check_declared_proteins,
    counts_by_source_type,
    normalise,
)
from cgload.parsers.tokenizer import tokenize
from cgload.storage import sha256_of, staged_load
from cgload.verify import Outcome, TierResult, count_source_naively

ANNOTATION_DIRECTORY = "annotation"
PROTEIN_FILENAME = "proteins.faa"

#: Attribute keys not worth a row. Bulky, wholly redundant with a stored column,
#: or pure prediction diagnostics (D-025). Everything else is kept, because
#: `explain` must be able to show what the input said.
SKIPPED_ATTRIBUTES = frozenset({"ID", "Parent"})


class LoadError(RuntimeError):
    """A load that cannot proceed without guessing. Always names the file."""


@dataclass
class FileReport:
    path: Path
    dialect: str
    confidence: str
    evidence: str
    feature_count: int
    counts_by_type: dict[str, int]
    unexpected_types: set[str] = field(default_factory=set)
    missing_types: set[str] = field(default_factory=set)
    protein_count: int = 0


@dataclass
class LoadReport:
    assembly_id: int
    data_directory: Path
    region_count: int
    total_length: int
    files: list[FileReport] = field(default_factory=list)
    functional: list[FunctionalReport] = field(default_factory=list)
    verification: list[TierResult] = field(default_factory=list)
    envelope_repairs: list = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def feature_count(self) -> int:
        return sum(report.feature_count for report in self.files)


def load_organism(
    engine: Engine,
    config: OrganismConfig,
    data_root: Path,
    *,
    copy: bool = True,
    force: bool = False,
    repair_envelopes: bool = False,
) -> LoadReport:
    """Load one organism: assembly, then every annotation file.

    Ordering matters and is deliberate. Every file is tokenised, detected and
    normalised *before* anything is inserted, so a dialect that cannot be
    identified or a hierarchy that cannot be resolved stops the load before the
    database is touched at all. Parent resolution across files happens once,
    against the union of features (D-014).
    """
    assert_vocabulary_current(engine)

    if not config.fasta.exists():
        raise FileNotFoundError(f"assembly not found: {config.fasta}")
    for annotation in config.annotation:
        if not annotation.path.exists():
            raise FileNotFoundError(f"annotation not found: {annotation.path}")
    for entry in config.functional:
        if not entry.path.exists():
            raise FileNotFoundError(f"functional annotation not found: {entry.path}")

    parsed = [_parse(annotation, repair_envelopes=repair_envelopes)
              for annotation in config.annotation]
    _reject_cross_file_collisions(parsed)

    with staged_load(
        engine,
        data_root,
        organism=f"{config.genus}_{config.species}",
        strain=config.strain,
        assembly=config.assembly_name,
        force=force,
    ) as load:
        connection = load.connection

        staged_fasta, fasta_hash = load.place(config.fasta, ASSEMBLY_FILENAME, copy=copy)
        index_path, regions = index_and_read(staged_fasta)

        organism_id = _upsert_organism(connection, config)
        strain_id = _upsert_strain(connection, organism_id, config.strain)
        _reject_duplicate_assembly(connection, strain_id, config, force=force)

        assembly_id = connection.execute(
            insert(assembly).values(
                strain_id=strain_id,
                name=config.assembly_name,
                version=config.assembly_version,
                accession=config.accession,
                fasta_path=load.relative(staged_fasta),
                fasta_sha256=fasta_hash,
                vocabulary_version=vocab.VOCABULARY_VERSION,
                vocabulary_hash=vocab.content_hash(),
                cgload_version=__version__,
                loaded_at=datetime.now(UTC),
            )
        ).inserted_primary_key[0]

        connection.execute(
            insert(source_file).values(
                assembly_id=assembly_id,
                path=str(config.fasta),
                sha256=fasta_hash,
                file_format="fasta",
            )
        )

        connection.execute(
            insert(sequence_region),
            [
                {
                    "assembly_id": assembly_id,
                    "source_id": region.source_id,
                    "length": region.length,
                }
                for region in regions
            ],
        )
        region_ids = dict(
            connection.execute(
                select(sequence_region.c.source_id, sequence_region.c.sequence_region_id).where(
                    sequence_region.c.assembly_id == assembly_id
                )
            ).all()
        )

        file_ids: dict[Path, int] = {}
        report = LoadReport(
            assembly_id=assembly_id,
            data_directory=load.target,
            region_count=len(regions),
            total_length=sum(region.length for region in regions),
        )

        for parse in parsed:
            report.warnings.extend(
                _check_declared_regions(parse, regions, config)
            )
            try:
                file_id = connection.execute(
                    insert(source_file).values(
                        assembly_id=assembly_id,
                        path=str(parse.path),
                        sha256=sha256_of(parse.path),
                        file_format="genbank" if parse.is_genbank else "gff3",
                        detected_dialect=parse.detection.profile.name,
                        detection_confidence=parse.detection.confidence.value,
                        id_prefix=parse.id_prefix,
                    )
                ).inserted_primary_key[0]
            except IntegrityError as exc:
                # Every constraint violation becomes a domain error naming the
                # file. A SQLAlchemy traceback is technically loud but tells the
                # user nothing they can act on.
                raise LoadError(
                    f"{parse.path.name}: the database rejected this file's provenance "
                    f"row, which usually means the same path is being loaded twice into "
                    f"one assembly. Underlying error: {exc.orig}"
                ) from exc

            file_ids[parse.path] = file_id
            _insert_features(connection, assembly_id, file_id, region_ids, parse)
            report.files.append(parse.as_report())
            report.envelope_repairs.extend(parse.envelope_repairs)
            report.warnings.extend(parse.warnings)

        proteins = {
            gene: sequence
            for parse in parsed
            for gene, sequence in parse.normalised.proteins.items()
        }
        if proteins:
            _write_proteins(load, proteins)

        load.manifest.update(
            {
                "cgload_version": __version__,
                "vocabulary_version": vocab.VOCABULARY_VERSION,
                "vocabulary_hash": vocab.content_hash(),
                "organism": f"{config.genus} {config.species}",
                "strain": config.strain,
                "assembly": config.assembly_name,
                "accession": config.accession,
                "ncbi_taxon_id": config.ncbi_taxon_id,
                "loaded_at": datetime.now(UTC).isoformat(),
                "regions": len(regions),
                "total_length": report.total_length,
                "files": [
                    {
                        "role": "assembly",
                        "source": str(config.fasta),
                        "stored": load.relative(staged_fasta),
                        "sha256": fasta_hash,
                        "mode": "copy" if copy else "link",
                    },
                    {"role": "index", "stored": load.relative(index_path)},
                ]
                + [
                    {
                        "role": "annotation",
                        "source": str(parse.path),
                        "sha256": sha256_of(parse.path),
                        "dialect": parse.detection.profile.name,
                        "confidence": parse.detection.confidence.value,
                        "evidence": parse.detection.evidence,
                        "id_prefix": parse.id_prefix,
                        "features": parse.normalised.feature_count,
                        "counts_by_source_type": counts_by_source_type(parse.normalised),
                    }
                    for parse in parsed
                ],
                "proteins": len(proteins),
                "warnings": report.warnings,
            }
        )

        # Functional annotation attaches to features, so it must come after
        # them -- and inside the same transaction, or a failure here would leave
        # a genome loaded and its annotation half-attached.
        for entry in config.functional:
            report.functional.append(
                attach_functional(connection, assembly_id, entry)
            )

        # The primary completeness check, run automatically rather than left for
        # the user to remember (D-046). A load that silently lost features is the
        # failure this project exists to prevent, so the check that catches it is
        # not optional and not a separate step.
        for parse in parsed:
            source_file_id = file_ids[parse.path]
            naive = count_source_naively(parse.path)
            stored = count_database_for_file(connection, assembly_id, source_file_id)
            result = compare_with_database(
                naive, stored, label=f"source count [{parse.path.name}]"
            )
            report.verification.append(result)
            if result.outcome is Outcome.FAIL:
                detail = "; ".join(
                    f"{name}: file has {expected:,}, database has {actual:,}"
                    for name, (expected, actual) in result.mismatches.items()
                )
                raise LoadError(
                    f"{parse.path.name}: the load did not store what the file "
                    f"contains, so it is being rolled back. {detail}. This is a bug "
                    f"in cgload, not in your file; please report it with the file "
                    f"that produced it."
                )

        # A repaired load is declared, not merely successful: the receipt records
        # every coordinate cgload changed and what it changed it from.
        load.manifest["envelope_repairs"] = [
            {
                "feature": item.feature_id,
                "type": item.feature_type,
                "seqid": item.seqid,
                "from": [item.old_start, item.old_end],
                "to": [item.new_start, item.new_end],
            }
            for item in report.envelope_repairs
        ]

        load.manifest["verification"] = [
            {"check": item.name, "outcome": item.outcome.value, "detail": item.detail}
            for item in report.verification
        ]

        load.manifest["functional"] = [
            {
                "path": str(item.path),
                "tool": item.tool,
                "tool_version": item.tool_version,
                "queries": item.query_count,
                "matched": item.matched_count,
                "unmatched": item.unmatched_count,
                "hits": item.hit_count,
                "analyses": item.analyses,
            }
            for item in report.functional
        ]

        return report


# --------------------------------------------------------------------------
# Parsing, before any row is written
# --------------------------------------------------------------------------


@dataclass
class ParsedFile:
    path: Path
    id_prefix: str | None
    detection: dialects.Detection
    normalised: NormalisedFile
    is_genbank: bool = False
    #: `##sequence-region` declarations, cross-checked against the FASTA.
    declared_regions: dict[str, tuple[int, int]] = field(default_factory=dict)
    #: Parent envelopes recomputed from their contents (D-055). Empty unless
    #: --repair-envelopes was passed, and always reported and manifested.
    envelope_repairs: list = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def as_report(self) -> FileReport:
        return FileReport(
            path=self.path,
            dialect=self.detection.profile.name,
            confidence=self.detection.confidence.value,
            evidence=self.detection.evidence,
            feature_count=self.normalised.feature_count,
            counts_by_type=counts_by_source_type(self.normalised),
            unexpected_types=self.normalised.unexpected_types,
            missing_types=self.normalised.missing_types,
            protein_count=len(self.normalised.proteins),
        )


def _is_genbank(path: Path) -> bool:
    """Settle the format by content.

    A GenBank flat file opens with a LOCUS line; a GFF3 opens with `##gff-version`
    or a tab-delimited row. Reading the first line is cheaper than trusting the
    name, and the name lies: the *M. aeruginosa* record downloaded as `.gbff`
    while actually being gzip (D-035).
    """
    with path.open("rb") as handle:
        head = handle.read(2)
        if head == b"\x1f\x8b":  # gzip
            import gzip

            with gzip.open(path, "rt", errors="replace") as unzipped:
                first = unzipped.readline()
        else:
            handle.seek(0)
            first = handle.readline().decode("utf-8", errors="replace")
    return first.startswith("LOCUS")


def _parse(annotation: AnnotationFile, *, repair_envelopes: bool = False) -> ParsedFile:
    """Tokenise, detect and normalise one file. No database involvement."""
    is_genbank = _is_genbank(annotation.path)
    if is_genbank:
        try:
            tokenized = tokenize_genbank(annotation.path)
        except GenBankError as exc:
            # Already names the file and line; re-raised so `load` reports it as
            # an input error rather than a traceback.
            raise LoadError(str(exc)) from exc
    else:
        if annotation.looks_like_genbank:
            raise LoadError(
                f"{annotation.path.name} is named like a GenBank file but does not "
                f"begin with a LOCUS line. Check the download is complete."
            )
        # --repair-envelopes has to reach the tokenizer, not only the envelope
        # stage: an inverted row is refused at read time, long before any
        # envelope logic runs, so without this the flag cannot help the one
        # maize row that needs it (D-056).
        tokenized = tokenize(annotation.path, tolerate_inverted=repair_envelopes)

    try:
        detection = dialects.detect(tokenized, override=annotation.dialect)
    except dialects.DetectionError as exc:
        raise LoadError(f"{annotation.path.name}: {exc}") from exc

    try:
        normalised = normalise(
            tokenized.features,
            detection.profile,
            id_prefix=annotation.id_prefix,
            comments=tokenized.comments,
            unexpected=dialects.unexpected_types(detection.profile, tokenized.features),
            missing=dialects.missing_types(detection.profile, tokenized.features),
        )
    except NormaliseError as exc:
        raise LoadError(f"{annotation.path.name}: {exc}") from exc

    # A parent must be at least as wide as what it contains (D-055). Every other
    # rule cgload has looks at one row in isolation and cannot see this.
    kept, repairs, remaining = apply_policy(normalised.features, repair=repair_envelopes)
    if remaining:
        backwards = sum(1 for item in remaining if item.is_self_contradiction)
        contained = len(remaining) - backwards
        parts = []
        if contained:
            parts.append(
                f"{contained} feature(s) lie outside the parent that claims them"
            )
        if backwards:
            parts.append(f"{backwards} row(s) are written back to front")
        raise LoadError(
            f"{annotation.path.name}: {' and '.join(parts)}, so at least one "
            f"coordinate in each case is wrong.\n{format_report([], remaining)}\n"
            f"If these are the known NCBI trans-splicing envelopes, "
            f"--repair-envelopes recomputes each parent from what it contains."
        )
    normalised.features = kept

    parsed = ParsedFile(
        path=annotation.path,
        id_prefix=annotation.id_prefix,
        detection=detection,
        normalised=normalised,
        is_genbank=is_genbank,
        declared_regions=tokenized.declared_regions,
        envelope_repairs=repairs,
    )

    problems = check_declared_proteins(normalised)
    if problems:
        listed = "\n  ".join(problems[:5])
        raise LoadError(
            f"{annotation.path.name}: the declared protein lengths do not follow from "
            f"the CDS spans ({len(problems)} gene(s)). This means the coordinates were "
            f"read wrongly or a segment was dropped, so the load is stopping:\n  {listed}"
        )

    if normalised.unexpected_types:
        parsed.warnings.append(
            f"{annotation.path.name}: loaded {len(normalised.unexpected_types)} feature "
            f"type(s) not declared by the {detection.profile.name} profile: "
            f"{', '.join(sorted(normalised.unexpected_types))}. They are stored verbatim "
            f"and linked only where the file gave an explicit parent."
        )
    if normalised.missing_types:
        parsed.warnings.append(
            f"{annotation.path.name}: the {detection.profile.name} profile expects "
            f"{', '.join(sorted(normalised.missing_types))}, which this file does not "
            f"contain. Check the file is complete."
        )
    if detection.confidence is dialects.Confidence.LOW:
        parsed.warnings.append(
            f"{annotation.path.name}: dialect {detection.profile.name} matched on "
            f"{detection.evidence}. Pass dialect = \"{detection.profile.name}\" in the "
            f"config to make the choice explicit."
        )
    return parsed


def _reject_cross_file_collisions(parsed: list[ParsedFile]) -> None:
    """The same ID in two files is an error naming both (D-014).

    Not a merge, not a silent namespace. Two files describing the same assembly
    may legitimately both be loaded — RefSeq annotation plus a funannotate file
    of novel genes — but if they reuse an ID, one gene would swallow the other.
    """
    seen_paths: dict[Path, int] = {}
    for index, parse in enumerate(parsed):
        if parse.path in seen_paths:
            raise LoadError(
                f"{parse.path.name} is listed twice in the config (entries "
                f"{seen_paths[parse.path] + 1} and {index + 1}). Loading it twice would "
                f"double every count; remove the duplicate entry."
            )
        seen_paths[parse.path] = index

    owners: dict[str, Path] = {}
    for parse in parsed:
        for item in parse.normalised.features:
            if item.source_id is None:
                continue
            previous = owners.get(item.source_id)
            if previous is not None and previous != parse.path:
                raise LoadError(
                    f"ID {item.source_id!r} is defined in both {previous.name} and "
                    f"{parse.path.name}. cgload will not merge or overwrite them; give "
                    f"one file an id_prefix in the config to keep them distinct."
                )
            owners[item.source_id] = parse.path


def _reject_duplicate_assembly(
    connection: Connection, strain_id: int, config: OrganismConfig, *, force: bool
) -> None:
    version = assembly.c.version
    condition = (
        version.is_(None) if config.assembly_version is None else version == config.assembly_version
    )
    existing = connection.execute(
        select(assembly.c.assembly_id).where(
            assembly.c.strain_id == strain_id,
            assembly.c.name == config.assembly_name,
            condition,
        )
    ).scalar_one_or_none()
    if existing is None:
        return
    if not force:
        raise AssemblyExistsError(
            f"assembly {config.assembly_name!r} for strain {config.strain!r} is already "
            f"loaded (assembly_id={existing}). Load it under a different version, or "
            f"pass --force to replace it."
        )

    # Under --force the superseded rows are deleted inside the same transaction
    # as the replacement.
    #
    # Returning silently here left both assemblies in the database sharing one
    # data directory, since --force also replaces the directory. One of the two
    # rows then described coordinates against a FASTA that no longer existed,
    # and nothing re-checks `fasta_sha256`, so the inconsistency was
    # undetectable (D-073).
    #
    # The delete is transactional: if the load fails, the rollback restores the
    # old rows, and storage.staged_load restores the old directory.
    connection.execute(delete(assembly).where(assembly.c.assembly_id == existing))


def _check_declared_regions(parse: ParsedFile, regions, config: OrganismConfig) -> list[str]:
    """Cross-check, never overwrite (finding 11).

    RefSeq's `##sequence-region` and its `region` row carry a length and a taxon
    id. The config is authoritative — the file is the thing being validated — so
    a disagreement is a warning naming both values.
    """
    warnings: list[str] = []
    lengths = {region.source_id: region.length for region in regions}

    for item in parse.normalised.features:
        if item.source_type != "region":
            continue
        for value in item.attributes.get("Dbxref", []):
            if value.startswith("taxon:") and config.ncbi_taxon_id is not None:
                declared = int(value.removeprefix("taxon:"))
                if declared != config.ncbi_taxon_id:
                    warnings.append(
                        f"{parse.path.name}: the file declares taxon {declared} but the "
                        f"config says {config.ncbi_taxon_id}. The config is being used."
                    )

    for source_id, (start, end) in parse.declared_regions.items():
        actual = lengths.get(source_id)
        if actual is not None and end - start + 1 != actual:
            warnings.append(
                f"{parse.path.name}: ##sequence-region says {source_id} is "
                f"{end - start + 1} bp but the assembly has {actual} bp."
            )
    return warnings


# --------------------------------------------------------------------------
# Insertion
# --------------------------------------------------------------------------


def _insert_features(
    connection: Connection,
    assembly_id: int,
    file_id: int,
    region_ids: dict[str, int],
    parse: ParsedFile,
) -> None:
    """Write one row per segment, then resolve parents, then attributes.

    Parents are set in a second pass because a child may appear before its
    parent in the file, and because ``feature.parent_id`` points at a row that
    may not exist yet on the first pass.
    """
    unknown = {
        item.seqid for item in parse.normalised.features if item.seqid not in region_ids
    }
    if unknown:
        raise LoadError(
            f"{parse.path.name}: {len(unknown)} sequence name(s) in the annotation are "
            f"absent from the assembly, first is {sorted(unknown)[0]!r}. The annotation "
            f"and the FASTA describe different assemblies."
        )

    rows = []
    for item in parse.normalised.features:
        for segment in item.segments:
            rows.append(
                {
                    "assembly_id": assembly_id,
                    "sequence_region_id": region_ids[item.seqid],
                    "source_file_id": file_id,
                    "source_id": item.source_id,
                    "source_type": item.source_type,
                    "feature_type": item.feature_type,
                    "start": segment.start,
                    "end": segment.end,
                    "strand": item.strand,
                    "phase": segment.phase,
                    "score": segment.score,
                    "segment_index": segment.segment_index,
                    "parent_id": None,
                    "linkage_method": item.linkage_method,
                    "source_program": item.source_program,
                }
            )

    try:
        connection.execute(insert(feature), rows)
    except IntegrityError as exc:
        raise LoadError(
            f"{parse.path.name}: the database rejected a feature row. This usually means "
            f"two rows share an ID and a span, which would be a duplicate rather than a "
            f"discontinuous feature. Underlying error: {exc.orig}"
        ) from exc

    _resolve_parent_ids(connection, assembly_id, file_id, parse)
    _insert_attributes(connection, assembly_id, file_id, parse)


def _resolve_parent_ids(
    connection: Connection, assembly_id: int, file_id: int, parse: ParsedFile
) -> None:
    """Point each segment at its parent's first segment.

    A child attaches to the parent's lowest-coordinate segment, so
    ``parent_id`` is single-valued and a walk up the graph is unambiguous. The
    full parent span is reachable through ``source_id``, which is why the
    choice of segment costs nothing.
    """
    stored = connection.execute(
        select(
            feature.c.feature_id,
            feature.c.source_id,
            feature.c.segment_index,
        ).where(feature.c.assembly_id == assembly_id, feature.c.source_file_id == file_id)
    ).all()

    first_segment: dict[str, int] = {}
    for row in stored:
        if row.source_id is None:
            continue
        if row.segment_index == 0:
            first_segment[row.source_id] = row.feature_id

    by_source_id: dict[str, list[int]] = defaultdict(list)
    for row in stored:
        if row.source_id is not None:
            by_source_id[row.source_id].append(row.feature_id)

    updates = []
    for item in parse.normalised.features:
        if item.parent_source_id is None or item.source_id is None:
            continue
        parent_id = first_segment.get(item.parent_source_id)
        if parent_id is None:
            raise LoadError(
                f"{parse.path.name}: parent {item.parent_source_id!r} of "
                f"{item.source_id!r} was not stored. This is an internal inconsistency, "
                f"not a problem with the file."
            )
        for feature_id in by_source_id[item.source_id]:
            updates.append({"row_id": feature_id, "parent": parent_id})

    if updates:
        connection.execute(
            feature.update()
            .where(feature.c.feature_id == bindparam("row_id"))
            .values(parent_id=bindparam("parent")),
            updates,
        )


def _insert_attributes(
    connection: Connection, assembly_id: int, file_id: int, parse: ParsedFile
) -> None:
    """Store column-9 attributes against the feature's first segment.

    Attributes describe the feature, not the segment, and RefSeq repeats them on
    every row — so storing them per segment would multiply the table that
    ``explain`` reads by the segment count for no gain.
    """
    first_segment = dict(
        connection.execute(
            select(feature.c.source_id, feature.c.feature_id).where(
                feature.c.assembly_id == assembly_id,
                feature.c.source_file_id == file_id,
                feature.c.segment_index == 0,
            )
        ).all()
    )

    rows = []
    for item in parse.normalised.features:
        if item.source_id is None:
            continue
        feature_id = first_segment.get(item.source_id)
        if feature_id is None:
            continue
        for key, values in item.attributes.items():
            if key in SKIPPED_ATTRIBUTES:
                continue
            for value in values:
                rows.append({"feature_id": feature_id, "key": key, "value": value})

    if rows:
        connection.execute(insert(feature_attribute), rows)


def _write_proteins(load, proteins: dict[str, str]) -> None:
    """Declared translations to disk, not to the database (D-025).

    A protein is sequence, and D-004 puts sequence on disk. Some are 2 kB, and
    ``feature_attribute`` is read on every ``explain`` call.
    """
    directory = load.staging / ANNOTATION_DIRECTORY
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / PROTEIN_FILENAME
    with path.open("w") as handle:
        for gene, sequence in sorted(proteins.items()):
            handle.write(f">{gene}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")
