"""Reading stored features back out, for export.

``cgload.export`` serialises feature objects and knows nothing about the
database. That is the right split — it makes the serialiser testable against
parser output directly — but it left `export` unable to do the one thing the
command exists for: write out what is *in the database*, after the input files
are gone.

This module is that half. It reconstructs the same shape ``normalise`` produces,
so the exporter cannot tell whether it was handed parser output or database rows,
and the round-trip test therefore covers the path a user actually takes.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from sqlalchemy import Connection, and_, select

from cgload.db.schema import (
    assembly,
    feature,
    feature_attribute,
    organism,
    sequence_region,
    source_file,
    strain,
)


@dataclass
class StoredSegment:
    """One row. Mirrors ``normalise.Segment``."""

    line_number: int
    start: int
    end: int
    phase: int | None
    score: float | None
    segment_index: int = 0
    #: The row's own strand (D-096).
    strand: str | None = None
    #: The row's own attributes, when it was stored with attributes that differ
    #: from the first row's (D-096); otherwise None.
    attributes: dict[str, list[str]] | None = None


@dataclass
class StoredFeature:
    """One feature, shaped like ``normalise.NormalisedFeature``.

    Deliberately duck-typed rather than reusing the parser's class: importing it
    here would give the exporter a path back into parser code, and the round-trip
    check is only meaningful while the two halves stay separable.
    """

    source_id: str | None
    source_type: str
    feature_type: str
    seqid: str
    strand: str
    source_program: str
    segments: list[StoredSegment]
    parent_source_id: str | None = None
    linkage_method: str | None = None
    attributes: dict[str, list[str]] = field(default_factory=dict)
    #: The database row of the first segment: the row that attributes, parent
    #: links and functional annotation point at.
    row_id: int | None = None
    #: The parent's first-segment row, as stored in ``parent_id``.
    parent_row_id: int | None = None

    @property
    def start(self) -> int:
        return min(segment.start for segment in self.segments)

    @property
    def end(self) -> int:
        return max(segment.end for segment in self.segments)

    @property
    def is_anonymous(self) -> bool:
        return self.source_id is None


@dataclass(frozen=True)
class AssemblyContext:
    """What the export header needs besides the features themselves."""

    label: str
    accession: str | None
    vocabulary_hash: str
    cgload_version: str
    sequence_regions: dict[str, int]
    dialects: tuple[str, ...]


def load_context(connection: Connection, assembly_id: int) -> AssemblyContext:
    row = connection.execute(
        select(
            organism.c.genus,
            organism.c.species,
            strain.c.name.label("strain"),
            assembly.c.name.label("assembly"),
            assembly.c.accession,
            assembly.c.vocabulary_hash,
            assembly.c.cgload_version,
        )
        .select_from(assembly.join(strain).join(organism))
        .where(assembly.c.assembly_id == assembly_id)
    ).one_or_none()
    if row is None:
        raise LookupError(f"no assembly with id {assembly_id}")

    regions = dict(
        connection.execute(
            select(sequence_region.c.source_id, sequence_region.c.length).where(
                sequence_region.c.assembly_id == assembly_id
            )
        ).all()
    )
    dialects = tuple(
        sorted(
            value
            for value in connection.execute(
                select(source_file.c.detected_dialect)
                .where(
                    source_file.c.assembly_id == assembly_id,
                    source_file.c.detected_dialect.is_not(None),
                )
                .distinct()
            ).scalars()
        )
    )
    return AssemblyContext(
        label=f"{row.genus} {row.species} / {row.strain} / {row.assembly}",
        accession=row.accession,
        vocabulary_hash=row.vocabulary_hash,
        cgload_version=row.cgload_version,
        sequence_regions=regions,
        dialects=dialects,
    )


def load_features(
    connection: Connection, assembly_id: int, *, source_file_id: int | None = None
) -> list[StoredFeature]:
    """Every feature of an assembly, or of one of its source files."""
    conditions = [feature.c.assembly_id == assembly_id]
    if source_file_id is not None:
        conditions.append(feature.c.source_file_id == source_file_id)
    return features_where(connection, and_(*conditions))


def features_where(connection: Connection, condition) -> list[StoredFeature]:
    """Rebuild features from the rows matching ``condition``, segments regrouped.

    ``condition`` is a SQL expression over the ``feature`` table. ``query``
    passes one that selects whole features; a condition that selects only some
    segments of a feature yields those segments only.

    Anonymous rows (D-024) cannot be grouped by identifier and are each their own
    feature, which is exactly how they were stored — so the count that comes out
    equals the count that went in without any special case.

    The parent's identifier comes from a join on ``parent_id``, so it is right
    even when the parent itself is not among the rows selected.
    """
    parent = feature.alias("parent")
    rows = connection.execute(
        select(
            feature,
            sequence_region.c.source_id.label("seqid"),
            source_file.c.detected_dialect.label("dialect"),
            parent.c.source_id.label("parent_source_id"),
        )
        .join(sequence_region, sequence_region.c.sequence_region_id == feature.c.sequence_region_id)
        .join(source_file, source_file.c.source_file_id == feature.c.source_file_id)
        .outerjoin(parent, parent.c.feature_id == feature.c.parent_id)
        .where(condition)
        .order_by(feature.c.feature_id, feature.c.segment_index)
    ).all()
    if not rows:
        return []

    attributes: dict[int, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    # A join, not an IN list of every row id: one bound parameter per row broke
    # export past the driver's limit -- SQLite's (999 / 32,766 / 250,000 by
    # build) and PostgreSQL's 65,535 (AUDIT-8). Ordered, because multi-valued
    # attributes such as start_range=.,3250 are order-sensitive and a server
    # database returns unordered rows in any order it likes.
    for row in connection.execute(
        select(
            feature_attribute.c.feature_id,
            feature_attribute.c.key,
            feature_attribute.c.value,
        )
        .join(feature, feature.c.feature_id == feature_attribute.c.feature_id)
        .where(condition)
        .order_by(feature_attribute.c.feature_attribute_id)
    ).all():
        attributes[row.feature_id][row.key].append(row.value)

    named: dict[str, StoredFeature] = {}
    ordered: list[StoredFeature] = []

    for row in rows:
        own = {key: list(values) for key, values in attributes.get(row.feature_id, {}).items()}
        segment = StoredSegment(
            line_number=row.feature_id,
            start=row.start,
            end=row.end,
            phase=row.phase,
            score=row.score,
            segment_index=row.segment_index,
            strand=row.strand,
        )
        if row.source_id is not None and row.source_id in named:
            item = named[row.source_id]
            item.segments.append(segment)
            if row.segment_index == 0:
                item.row_id = row.feature_id
                item.parent_row_id = row.parent_id
                if own:
                    item.attributes = own
            elif own:
                # A later row stored with attributes of its own (D-096).
                segment.attributes = own
            continue

        item = StoredFeature(
            source_id=row.source_id,
            source_type=row.source_type,
            feature_type=row.feature_type,
            seqid=row.seqid,
            strand=row.strand,
            source_program=row.source_program or row.dialect or "cgload",
            segments=[segment],
            parent_source_id=row.parent_source_id,
            linkage_method=row.linkage_method,
            attributes={
                key: list(values) for key, values in attributes.get(row.feature_id, {}).items()
            },
            row_id=row.feature_id,
            parent_row_id=row.parent_id,
        )
        if row.source_id is not None:
            named[row.source_id] = item
        ordered.append(item)

    for item in ordered:
        item.segments.sort(key=lambda segment: (segment.start, segment.end))
        # The feature's strand: its rows' strand, or `?` when they lie on both
        # strands, as a trans-spliced feature's do (D-068, D-096).
        strands = {segment.strand for segment in item.segments}
        item.strand = strands.pop() if len(strands) == 1 else "?"
    ordered.sort(key=lambda item: (item.seqid, item.start, item.end))
    return ordered
