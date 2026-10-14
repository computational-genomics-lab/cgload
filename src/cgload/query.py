"""``cgload query``: read-only questions about one loaded assembly (D-085).

Three kinds of question, each answered with whole features:

* **by identifier**: the feature whose ID is given, and with ``descendants``
  everything under it (a gene with its transcripts, exons and CDS);
* **by region**: every feature with at least one row overlapping a region, in
  samtools notation (``chr``, ``chr:start-end``, ``chr:start``);
* **by functional term**: every feature carrying a GO, KEGG, Pfam, EC or other
  accession that eggNOG-mapper or InterProScan assigned to it.

Every statement is built with SQLAlchemy Core, so every value the user typed is
a bound parameter. Sets of rows are passed between steps as subqueries, never
as IN lists built in Python, because those fail past the driver's parameter
limit on a large genome (AUDIT-8). The parent chain is walked with recursive
CTEs, which SQLite, PostgreSQL and MySQL 8 all support. A recursive step uses
UNION rather than UNION ALL, so a cycle in the parent links, which the loader
should never store, cannot make a query run forever.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from sqlalchemy import Connection, and_, func, or_, select

from cgload.db.export_source import StoredFeature, features_where
from cgload.db.schema import annotation_run, assembly, feature, protein_annotation, sequence_region
from cgload.parsers.functional import canonical_kegg_pathway


class QueryError(ValueError):
    """The question cannot be answered as asked: an unknown assembly, sequence
    or identifier, or a malformed region."""


@dataclass(frozen=True)
class Region:
    """1-based, inclusive, as in GFF3 and samtools."""

    seqid: str
    start: int
    end: int
    region_id: int
    length: int

    def __str__(self) -> str:
        return f"{self.seqid}:{self.start:,}-{self.end:,}"


@dataclass(frozen=True)
class TermHit:
    """One functional assignment that matched, on the row it is attached to."""

    row_id: int
    analysis: str
    accession: str
    description: str | None
    protein_start: int | None
    protein_end: int | None
    evalue: float | None


@dataclass
class QueryResult:
    features: list[StoredFeature]
    #: Row id of each feature's first segment -> identifier of the top of its
    #: parent chain (the feature itself when it has no parent).
    top_level: dict[int, str | None]
    hits: list[TermHit]
    description: str


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------


def require_assembly(connection: Connection, assembly_id: int) -> None:
    found = connection.execute(
        select(assembly.c.assembly_id).where(assembly.c.assembly_id == assembly_id)
    ).first()
    if found is None:
        raise QueryError(f"no assembly with id {assembly_id}")


_RANGE = re.compile(r"^(?P<start>[\d,]+)(?:-(?P<end>[\d,]+))?$")


def resolve_region(connection: Connection, assembly_id: int, text: str) -> Region:
    """Read ``chr``, ``chr:start-end`` or ``chr:start`` (to the sequence's end).

    A sequence name may itself contain a colon (human alt contigs do), so the
    whole text is tried as a name first, as samtools does. Commas in numbers
    are ignored. An end past the sequence is cut to its length; a start past
    it is an error.
    """
    text = text.strip()
    if not text:
        raise QueryError("the region is empty")

    def lookup(name: str):
        return connection.execute(
            select(sequence_region.c.sequence_region_id, sequence_region.c.length).where(
                sequence_region.c.assembly_id == assembly_id,
                sequence_region.c.source_id == name,
            )
        ).first()

    whole = lookup(text)
    if whole is not None:
        return Region(text, 1, whole.length, whole.sequence_region_id, whole.length)

    name, _, span = text.rpartition(":")
    found = lookup(name) if name else None
    if found is None:
        examples = connection.execute(
            select(sequence_region.c.source_id)
            .where(sequence_region.c.assembly_id == assembly_id)
            .order_by(sequence_region.c.sequence_region_id)
            .limit(3)
        ).scalars().all()
        shown = ", ".join(examples) if examples else "none are stored"
        raise QueryError(
            f"no sequence named {(name or text)!r} in assembly {assembly_id} "
            f"(for example: {shown})"
        )
    match = _RANGE.match(span.strip())
    if match is None:
        raise QueryError(f"cannot read {span!r} as start-end; write it like {name}:1000-2000")
    start = int(match["start"].replace(",", ""))
    end = int(match["end"].replace(",", "")) if match["end"] else found.length
    if start < 1 or end < start:
        raise QueryError(f"{text!r}: the start must be at least 1 and not after the end")
    if start > found.length:
        raise QueryError(f"{text!r}: {name} is only {found.length:,} bp long")
    return Region(name, start, min(end, found.length), found.sequence_region_id, found.length)


def term_variants(term: str) -> list[str]:
    """The stored spellings a typed accession may have (D-085).

    eggNOG writes KEGG orthologs as ``ko:K00799`` and people type ``K00799``;
    KEGG pathways are stored in their ``ko`` form and people paste the ``map``
    form (D-032); GO and EC identifiers arrive with or without their prefix in
    either case. Nothing else is rewritten, and nothing is matched by prefix or
    similarity: a Pfam name from eggNOG (``Pkinase``) and a Pfam accession from
    InterProScan (``PF00069``) are different strings, and cgload has no table
    to map one to the other.
    """
    stripped = term.strip()
    if not stripped:
        raise QueryError("the term is empty")
    variants = {stripped}
    lower = stripped.lower()
    if lower.startswith("go:"):
        variants.add("GO:" + stripped[3:])
    elif lower.startswith("ec:"):
        variants.add(stripped[3:])
    elif lower.startswith("ko:") and stripped[3:4].upper() == "K":
        variants.update({stripped[3:], "ko:" + stripped[3:]})
    elif re.fullmatch(r"[Kk]\d{5}", stripped):
        variants.update({stripped.upper(), "ko:" + stripped.upper()})
    variants.add(canonical_kegg_pathway(stripped))
    return sorted(variants)


# --------------------------------------------------------------------------
# Row selections. Each returns a SELECT of feature.feature_id, used as a
# subquery by the next step.
# --------------------------------------------------------------------------


def _rows_of_id(assembly_id: int, source_id: str):
    return select(feature.c.feature_id).where(
        feature.c.assembly_id == assembly_id, feature.c.source_id == source_id
    )


def _with_descendants(assembly_id: int, rows):
    """``rows``, every other row of their features, and everything below them.

    The walk starts from every row of each matched feature, not only the rows
    that matched. A child's ``parent_id`` names its parent's first row, and a
    region can overlap only a later row of a parent written on several rows,
    as NCBI writes a trans-spliced mRNA. Starting from the matched rows alone
    lost the children that hang from the first row.
    """
    start = select(feature.c.feature_id).where(_whole_features(assembly_id, rows))
    tree = start.cte("descendants", recursive=True)
    walked = tree.alias()
    child = feature.alias("child")
    tree = tree.union(select(child.c.feature_id).where(child.c.parent_id == walked.c.feature_id))
    return select(tree.c.feature_id)


def _with_ancestors(rows):
    start = select(feature.c.feature_id, feature.c.parent_id).where(
        feature.c.feature_id.in_(rows)
    )
    chain = start.cte("ancestors", recursive=True)
    walked = chain.alias()
    above = feature.alias("above")
    chain = chain.union(
        select(above.c.feature_id, above.c.parent_id).where(
            above.c.feature_id == walked.c.parent_id
        )
    )
    return select(chain.c.feature_id)


def _rows_in_region(region: Region):
    """Rows overlapping the region, including rows that wrap round a circle.

    A row that crosses the origin of a circular sequence is stored as GFF3
    writes it, ending past the sequence's length, and NCBI writes some rows
    wholly past it (D-081). Such a row covers position p at p + length, so the
    region is tried a second time, one length further on. On a linear sequence
    the second test can never hold: the loader refuses any row that ends past
    the end of a linear sequence.
    """
    shifted_start = region.start + region.length
    shifted_end = region.end + region.length
    return select(feature.c.feature_id).where(
        feature.c.sequence_region_id == region.region_id,
        or_(
            and_(feature.c.start <= region.end, feature.c.end >= region.start),
            and_(feature.c.start <= shifted_end, feature.c.end >= shifted_start),
        ),
    )


def _of_types(rows, types: tuple[str, ...]):
    if not types:
        return rows
    return select(feature.c.feature_id).where(
        feature.c.feature_id.in_(rows), feature.c.feature_type.in_(sorted(set(types)))
    )


def _whole_features(assembly_id: int, rows):
    """Every row of every feature that has a row in ``rows``."""
    named = select(feature.c.source_id).where(
        feature.c.feature_id.in_(rows), feature.c.source_id.is_not(None)
    )
    return and_(
        feature.c.assembly_id == assembly_id,
        or_(
            feature.c.source_id.in_(named),
            and_(feature.c.source_id.is_(None), feature.c.feature_id.in_(rows)),
        ),
    )


def _top_level(connection: Connection, rows) -> dict[int, str | None]:
    """First-segment row id -> identifier at the top of its parent chain."""
    start = select(
        feature.c.feature_id.label("origin"),
        feature.c.feature_id.label("node"),
        feature.c.parent_id.label("up"),
    ).where(feature.c.feature_id.in_(rows), feature.c.segment_index == 0)
    chain = start.cte("lineage", recursive=True)
    walked = chain.alias()
    above = feature.alias("above")
    chain = chain.union(
        select(walked.c.origin, above.c.feature_id, above.c.parent_id).where(
            above.c.feature_id == walked.c.up
        )
    )
    tops = connection.execute(
        select(chain.c.origin, feature.c.source_id)
        .join(feature, feature.c.feature_id == chain.c.node)
        .where(chain.c.up.is_(None))
    ).all()
    return {row.origin: row.source_id for row in tops}


# --------------------------------------------------------------------------
# The three questions
# --------------------------------------------------------------------------


def _finish(
    connection: Connection,
    assembly_id: int,
    rows,
    *,
    ancestors: bool,
    hits: list[TermHit],
    description: str,
) -> QueryResult:
    selected = _with_ancestors(rows) if ancestors else rows
    features = features_where(connection, _whole_features(assembly_id, selected))
    top = _top_level(connection, select(feature.c.feature_id).where(
        _whole_features(assembly_id, selected)
    ))
    return QueryResult(features, top, hits, description)


def by_id(
    connection: Connection,
    assembly_id: int,
    source_id: str,
    *,
    descendants: bool = False,
    types: tuple[str, ...] = (),
    ancestors: bool = False,
) -> QueryResult:
    require_assembly(connection, assembly_id)
    rows = _rows_of_id(assembly_id, source_id)
    if connection.execute(rows.limit(1)).first() is None:
        raise LookupError(f"no feature with ID {source_id!r} in assembly {assembly_id}")
    if descendants:
        rows = _with_descendants(assembly_id, rows)
    rows = _of_types(rows, types)
    label = f"ID {source_id}" + (" and its descendants" if descendants else "")
    return _finish(connection, assembly_id, rows, ancestors=ancestors, hits=[],
                   description=label)


def in_region(
    connection: Connection,
    assembly_id: int,
    region_text: str,
    *,
    descendants: bool = False,
    types: tuple[str, ...] = (),
    ancestors: bool = False,
) -> QueryResult:
    require_assembly(connection, assembly_id)
    region = resolve_region(connection, assembly_id, region_text)
    rows = _rows_in_region(region)
    if descendants:
        rows = _with_descendants(assembly_id, rows)
    rows = _of_types(rows, types)
    return _finish(connection, assembly_id, rows, ancestors=ancestors, hits=[],
                   description=f"features overlapping {region}")


def with_term(
    connection: Connection,
    assembly_id: int,
    term: str,
    *,
    analysis: str | None = None,
    descendants: bool = False,
    ancestors: bool = False,
) -> QueryResult:
    require_assembly(connection, assembly_id)
    variants = term_variants(term)
    conditions = [
        annotation_run.c.assembly_id == assembly_id,
        protein_annotation.c.accession.in_(variants),
    ]
    if analysis:
        conditions.append(func.lower(protein_annotation.c.analysis) == analysis.lower())
    matched = (
        select(
            protein_annotation.c.feature_id,
            protein_annotation.c.analysis,
            protein_annotation.c.accession,
            protein_annotation.c.description,
            protein_annotation.c.protein_start,
            protein_annotation.c.protein_end,
            protein_annotation.c.evalue,
        )
        .join(annotation_run,
              annotation_run.c.annotation_run_id == protein_annotation.c.annotation_run_id)
        .where(*conditions)
        .order_by(protein_annotation.c.feature_id, protein_annotation.c.analysis,
                  protein_annotation.c.accession, protein_annotation.c.protein_start)
    )
    hits = [TermHit(*row) for row in connection.execute(matched).all()]
    rows = select(protein_annotation.c.feature_id).join(
        annotation_run,
        annotation_run.c.annotation_run_id == protein_annotation.c.annotation_run_id,
    ).where(*conditions)
    if descendants:
        rows = _with_descendants(assembly_id, rows)
    shown = term if not analysis else f"{analysis} {term}"
    return _finish(connection, assembly_id, rows, ancestors=ancestors, hits=hits,
                   description=f"features annotated with {shown}")


def stored_analyses(connection: Connection, assembly_id: int) -> list[str]:
    """Analysis names present for an assembly, for an empty term query's hint."""
    return list(
        connection.execute(
            select(protein_annotation.c.analysis)
            .join(annotation_run,
                  annotation_run.c.annotation_run_id == protein_annotation.c.annotation_run_id)
            .where(annotation_run.c.assembly_id == assembly_id)
            .distinct()
            .order_by(protein_annotation.c.analysis)
        ).scalars()
    )
