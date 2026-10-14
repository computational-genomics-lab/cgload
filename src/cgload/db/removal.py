"""Deleting every row one assembly owns.

Two callers need this: ``load --force``, which deletes the old assembly inside
the transaction that loads its replacement, and ``remove``. Both used to be
missing. ``load --force`` deleted only the ``assembly`` row, and since foreign
keys are enforced (on SQLite too, by ``build_engine``), every replacement of an
assembly that had any features stopped with "FOREIGN KEY constraint failed".
It failed safely, because the transaction rolled back, but it never worked.

Rows go in dependency order, children before the rows they reference
(D-083):

1. functional annotation, then its runs;
2. feature attributes;
3. the parent links between the assembly's own features, then the features;
4. sequence regions and source files;
5. the assembly row;
6. with ``prune``, the strain and organism if nothing else refers to them.

Step 3 clears ``parent_id`` before deleting. SQLite and PostgreSQL check a
foreign key at the end of the statement, so one ``DELETE`` of parents and
children together would pass on both, but MySQL's InnoDB checks row by row and
refuses to delete a parent before its children.

Only rows this assembly owns are touched. Nothing here deletes a row of another
assembly that happens to point at this one; by construction none should exist,
and if one did, the foreign key error is the right outcome rather than a quiet
cascade into someone else's data.

Every statement selects its rows with a subquery, never an IN list built in
Python, for the reason given in AUDIT-8.
"""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import Connection, delete, func, select, update

from cgload.db.counts import feature_count_expression
from cgload.db.schema import (
    annotation_run,
    assembly,
    feature,
    feature_attribute,
    organism,
    protein_annotation,
    sequence_region,
    source_file,
    strain,
)


@dataclass(frozen=True)
class Footprint:
    """What an assembly holds, counted before it is removed."""

    features: int
    rows: int
    attributes: int
    functional_assignments: int
    annotation_runs: int
    sequence_regions: int
    source_files: int


def _scalar(connection: Connection, statement) -> int:
    return int(connection.execute(statement).scalar_one() or 0)


def footprint(connection: Connection, assembly_id: int) -> Footprint:
    owned = feature.c.assembly_id == assembly_id
    runs = select(annotation_run.c.annotation_run_id).where(
        annotation_run.c.assembly_id == assembly_id
    )
    return Footprint(
        features=_scalar(connection, select(feature_count_expression()).where(owned)),
        rows=_scalar(connection, select(func.count()).select_from(feature).where(owned)),
        attributes=_scalar(
            connection,
            select(func.count())
            .select_from(feature_attribute.join(feature))
            .where(owned),
        ),
        functional_assignments=_scalar(
            connection,
            select(func.count())
            .select_from(protein_annotation)
            .where(protein_annotation.c.annotation_run_id.in_(runs)),
        ),
        annotation_runs=_scalar(
            connection,
            select(func.count())
            .select_from(annotation_run)
            .where(annotation_run.c.assembly_id == assembly_id),
        ),
        sequence_regions=_scalar(
            connection,
            select(func.count())
            .select_from(sequence_region)
            .where(sequence_region.c.assembly_id == assembly_id),
        ),
        source_files=_scalar(
            connection,
            select(func.count())
            .select_from(source_file)
            .where(source_file.c.assembly_id == assembly_id),
        ),
    )


def delete_assembly_rows(
    connection: Connection, assembly_id: int, *, prune: bool = False
) -> None:
    """Delete everything ``assembly_id`` owns, inside the caller's transaction.

    ``prune`` also deletes the strain, and then the organism, when no other
    assembly or strain refers to them. ``remove`` prunes, so that removing the
    only assembly of an organism leaves the database as it was before the load.
    ``load --force`` does not, because the replacement reuses the strain.
    """
    strain_id = connection.execute(
        select(assembly.c.strain_id).where(assembly.c.assembly_id == assembly_id)
    ).scalar_one_or_none()
    if strain_id is None:
        raise LookupError(f"no assembly with id {assembly_id}")

    runs = select(annotation_run.c.annotation_run_id).where(
        annotation_run.c.assembly_id == assembly_id
    )
    owned_rows = select(feature.c.feature_id).where(feature.c.assembly_id == assembly_id)

    connection.execute(
        delete(protein_annotation).where(protein_annotation.c.annotation_run_id.in_(runs))
    )
    connection.execute(delete(annotation_run).where(annotation_run.c.assembly_id == assembly_id))
    connection.execute(
        delete(feature_attribute).where(feature_attribute.c.feature_id.in_(owned_rows))
    )
    connection.execute(
        update(feature)
        .where(feature.c.assembly_id == assembly_id, feature.c.parent_id.is_not(None))
        .values(parent_id=None)
    )
    connection.execute(delete(feature).where(feature.c.assembly_id == assembly_id))
    connection.execute(
        delete(sequence_region).where(sequence_region.c.assembly_id == assembly_id)
    )
    connection.execute(delete(source_file).where(source_file.c.assembly_id == assembly_id))
    connection.execute(delete(assembly).where(assembly.c.assembly_id == assembly_id))

    if not prune:
        return
    still_used = select(assembly.c.assembly_id).where(assembly.c.strain_id == strain_id)
    if connection.execute(still_used.limit(1)).first() is not None:
        return
    organism_id = connection.execute(
        select(strain.c.organism_id).where(strain.c.strain_id == strain_id)
    ).scalar_one()
    connection.execute(delete(strain).where(strain.c.strain_id == strain_id))
    other_strains = select(strain.c.strain_id).where(strain.c.organism_id == organism_id)
    if connection.execute(other_strains.limit(1)).first() is None:
        connection.execute(delete(organism).where(organism.c.organism_id == organism_id))
