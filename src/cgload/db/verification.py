"""The database half of the completeness check.

Kept out of ``cgload.verify`` on purpose. That module must not import parser
code, and a test enforces it by scanning the source text for the strings
``tokenize`` and ``normalise`` — so anything needing the schema lives here
instead, and the naive counting pass stays provably free of cgload's own
reading logic.

The division is not cosmetic. ``verify.py`` answers *what does the file say*.
This module answers *what did we store*. Neither can see the other's method,
which is the property D-008 is built on.
"""

from __future__ import annotations

from collections import Counter

from sqlalchemy import Connection, func, select

from cgload.db.counts import feature_count_expression
from cgload.db.schema import feature, sequence_region, source_file
from cgload.verify import Outcome, TierResult

#: Types the naive pass counts but the loader deliberately does not store as
#: features. Empty today and kept as the place to record any such divergence,
#: because an undocumented one would look identical to a lost feature.
NOT_STORED: frozenset[str] = frozenset()


def count_database(connection: Connection, assembly_id: int) -> Counter:
    """Features per verbatim source type, for one assembly.

    Groups by ``source_type`` rather than ``feature_type`` (D-017): the naive
    pass reads column 3 as written, so comparing against the canonical column
    would make every synonym look like a discrepancy. Counts features rather
    than rows, so a six-segment CDS is one (D-008, D-024).
    """
    rows = connection.execute(
        select(feature.c.source_type, feature_count_expression())
        .where(feature.c.assembly_id == assembly_id)
        .group_by(feature.c.source_type)
    ).all()
    return Counter({row[0]: row[1] for row in rows})


def count_database_for_file(
    connection: Connection, assembly_id: int, source_file_id: int
) -> Counter:
    """As above, restricted to one input file.

    Needed because an assembly may be built from several annotation files
    (D-014), and comparing one file's naive count against the whole assembly's
    stored count would report every other file's features as surplus.
    """
    rows = connection.execute(
        select(feature.c.source_type, feature_count_expression())
        .where(
            feature.c.assembly_id == assembly_id,
            feature.c.source_file_id == source_file_id,
        )
        .group_by(feature.c.source_type)
    ).all()
    return Counter({row[0]: row[1] for row in rows})


def compare_with_database(
    source_counts: Counter, database_counts: Counter, *, label: str = "source count"
) -> TierResult:
    """The primary check: what the file says against what was stored.

    Disagreement is always a failure. Unlike the reference comparison there is
    no innocent explanation: both numbers describe the same file, so a
    difference means something was lost or invented between reading and saving.
    """
    expected = {
        name: count for name, count in source_counts.items() if name not in NOT_STORED
    }
    actual = dict(database_counts)
    differing = {
        name
        for name in set(expected) | set(actual)
        if expected.get(name, 0) != actual.get(name, 0)
    }

    total_source = sum(
        count for name, count in source_counts.items() if name not in NOT_STORED
    )
    total_database = sum(database_counts.values())

    if differing:
        return TierResult(
            name=label,
            outcome=Outcome.FAIL,
            detail=(
                f"{len(differing)} feature type(s) differ between the file and the "
                f"database; features were lost or invented during loading"
            ),
            expected=expected,
            actual=actual,
        )

    return TierResult(
        name=label,
        outcome=Outcome.PASS,
        detail=(
            f"{total_source:,} features in the file, {total_database:,} in the "
            f"database, across {len(source_counts)} type(s)"
        ),
    )


def annotation_files(connection: Connection, assembly_id: int) -> list[tuple[int, str]]:
    """Every annotation input for an assembly, as (source_file_id, path)."""
    rows = connection.execute(
        select(source_file.c.source_file_id, source_file.c.path).where(
            source_file.c.assembly_id == assembly_id,
            source_file.c.file_format != "fasta",
        )
    ).all()
    return [(row[0], row[1]) for row in rows]


def sequences_loaded(connection: Connection, assembly_id: int) -> int:
    return connection.execute(
        select(func.count()).select_from(sequence_region).where(
            sequence_region.c.assembly_id == assembly_id
        )
    ).scalar_one()
