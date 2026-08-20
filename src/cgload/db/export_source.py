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

from sqlalchemy import Connection, select

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
    """Rebuild features from rows, grouping segments back together.

    Anonymous rows (D-024) cannot be grouped by identifier and are each their own
    feature, which is exactly how they were stored — so the count that comes out
    equals the count that went in without any special case.
    """
    conditions = [feature.c.assembly_id == assembly_id]
    if source_file_id is not None:
        conditions.append(feature.c.source_file_id == source_file_id)

    regions = dict(
        connection.execute(
            select(sequence_region.c.sequence_region_id, sequence_region.c.source_id)
        ).all()
    )
    programs_by_file = dict(
        connection.execute(
            select(source_file.c.source_file_id, source_file.c.detected_dialect)
        ).all()
    )

    rows = connection.execute(
        select(feature)
        .where(*conditions)
        .order_by(feature.c.feature_id, feature.c.segment_index)
    ).all()
    if not rows:
        return []

    parent_of = {row.feature_id: row.parent_id for row in rows}
    identifier_of = {row.feature_id: row.source_id for row in rows}

    attributes: dict[int, dict[str, list[str]]] = defaultdict(lambda: defaultdict(list))
    for row in connection.execute(
        select(
            feature_attribute.c.feature_id,
            feature_attribute.c.key,
            feature_attribute.c.value,
        ).where(
            feature_attribute.c.feature_id.in_([row.feature_id for row in rows])
        )
    ).all():
        attributes[row.feature_id][row.key].append(row.value)

    named: dict[str, StoredFeature] = {}
    ordered: list[StoredFeature] = []

    for row in rows:
        segment = StoredSegment(
            line_number=row.feature_id,
            start=row.start,
            end=row.end,
            phase=row.phase,
            score=row.score,
            segment_index=row.segment_index,
        )
        if row.source_id is not None and row.source_id in named:
            named[row.source_id].segments.append(segment)
            continue

        parent_id = parent_of.get(row.feature_id)
        item = StoredFeature(
            source_id=row.source_id,
            source_type=row.source_type,
            feature_type=row.feature_type,
            seqid=regions[row.sequence_region_id],
            strand=row.strand,
            source_program=row.source_program
            or programs_by_file.get(row.source_file_id)
            or "cgload",
            segments=[segment],
            parent_source_id=identifier_of.get(parent_id) if parent_id else None,
            linkage_method=row.linkage_method,
            attributes={key: list(values) for key, values in attributes[row.feature_id].items()},
        )
        if row.source_id is not None:
            named[row.source_id] = item
        ordered.append(item)

    for item in ordered:
        item.segments.sort(key=lambda segment: (segment.start, segment.end))
    ordered.sort(key=lambda item: (item.seqid, item.start, item.end))
    return ordered
