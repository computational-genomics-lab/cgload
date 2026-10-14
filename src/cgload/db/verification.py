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
from pathlib import Path

from sqlalchemy import Connection, func, select

from cgload.db.counts import feature_count_expression
from cgload.db.schema import assembly, feature, sequence_region, source_file
from cgload.storage import sha256_of
from cgload.verify import LINKED_TYPES, Outcome, TierResult, check_parent_edges

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


def annotation_files_with_format(
    connection: Connection, assembly_id: int
) -> list[tuple[int, str, str, str | None]]:
    """As ``annotation_files``, plus each file's format and detected dialect."""
    rows = connection.execute(
        select(
            source_file.c.source_file_id,
            source_file.c.path,
            source_file.c.file_format,
            source_file.c.detected_dialect,
        ).where(
            source_file.c.assembly_id == assembly_id,
            source_file.c.file_format != "fasta",
        )
    ).all()
    return [(row[0], row[1], row[2], row[3]) for row in rows]


def count_parent_edges(
    connection: Connection,
    assembly_id: int,
    source_file_id: int,
    *,
    computed_edges: bool,
) -> tuple[dict[str | None, int], int]:
    """Edges per linkage method for one file, and the features stating none.

    Counts features, not rows: each feature has exactly one row with
    ``segment_index`` 0.

    Where the file states its own parents (GFF3), every feature the file gave a
    parent must be linked, whatever its type; the loader records that statement
    as ``linkage_method``. That includes the anonymous intron and codon rows
    whose parents were silently dropped before AUDIT-2. A feature the file left
    parentless is top-level by the file's own statement (D-069); those of a type
    that usually takes a parent are returned separately, for the report.

    Where every edge is computed (GenBank), every feature of a type in
    ``LINKED_TYPES`` is expected to have one (D-074).

    A feature that should be linked and is not is counted under ``None``; one
    that is linked although the file named no parent, under ``"unrecorded"``.
    """
    conditions = [
        feature.c.assembly_id == assembly_id,
        feature.c.source_file_id == source_file_id,
        feature.c.segment_index == 0,
    ]
    if computed_edges:
        conditions.append(feature.c.feature_type.in_(sorted(LINKED_TYPES)))
    rows = connection.execute(
        select(feature.c.parent_id, feature.c.linkage_method, feature.c.feature_type)
        .where(*conditions)
    ).all()
    per_method: Counter = Counter()
    top_level = 0
    for parent_id, method, feature_type in rows:
        if not computed_edges and method is None:
            if parent_id is not None:
                per_method["unrecorded"] += 1
            elif feature_type in LINKED_TYPES:
                top_level += 1
            continue
        per_method[method if parent_id is not None else None] += 1
    return dict(per_method), top_level


def parent_edge_tier(
    connection: Connection,
    assembly_id: int,
    source_file_id: int,
    *,
    file_format: str,
    dialect: str | None,
    label: str = "parent edges",
) -> TierResult:
    """The parent-edge tier for one stored file (AUDIT-1)."""
    computed = file_format == "genbank"
    per_method, top_level = count_parent_edges(
        connection, assembly_id, source_file_id, computed_edges=computed
    )
    result = check_parent_edges(
        per_method, dialect="genbank" if computed else (dialect or file_format)
    )
    if per_method.get("unrecorded"):
        result = TierResult(
            label,
            Outcome.FAIL,
            f"{per_method['unrecorded']:,} feature(s) are linked although the file "
            "named no parent for them; an edge was invented during loading",
            {"unrecorded": 0},
            {"unrecorded": per_method["unrecorded"]},
        )
    detail = result.detail
    if top_level:
        detail += f"; {top_level:,} more named no parent and are top-level (D-069)"
    return TierResult(label, result.outcome, detail, result.expected, result.actual)


def stored_assembly_tier(
    connection: Connection, assembly_id: int, *, data_dir: Path | None = None
) -> TierResult:
    """Is the stored assembly still there, and still the file that was loaded?

    Every coordinate in the database refers to it. Its SHA-256 was recorded at
    load time and never checked again, so a deleted or replaced FASTA left a
    database that verified PASS (AUDIT-4, D-079). This tier needs the data
    directory: the one given with ``--data-dir``, or the one recorded at load.
    """
    name = "stored assembly"
    row = connection.execute(
        select(assembly.c.fasta_path, assembly.c.fasta_sha256, assembly.c.data_root).where(
            assembly.c.assembly_id == assembly_id
        )
    ).one_or_none()
    if row is None or row.fasta_path is None:
        return TierResult(name, Outcome.SKIPPED, "no stored assembly is recorded")
    root = data_dir if data_dir is not None else (Path(row.data_root) if row.data_root else None)
    if root is None:
        return TierResult(
            name,
            Outcome.SKIPPED,
            "the database does not record where the data directory is; pass --data-dir",
        )
    stored = Path(root) / row.fasta_path
    if not stored.is_file():
        return TierResult(
            name,
            Outcome.FAIL,
            f"the stored assembly is missing: {stored}. Every coordinate in the "
            f"database refers to it. If the data directory moved, pass --data-dir.",
        )
    index = stored.with_name(stored.name + ".fai")
    if not index.is_file():
        return TierResult(name, Outcome.FAIL, f"the assembly index is missing: {index}")
    digest = sha256_of(stored)
    if row.fasta_sha256 and digest != row.fasta_sha256:
        return TierResult(
            name,
            Outcome.FAIL,
            f"{stored} has changed since it was loaded: SHA-256 {digest[:12]}..., "
            f"recorded {row.fasta_sha256[:12]}...",
        )
    return TierResult(
        name, Outcome.PASS, f"{stored} and its index are present; SHA-256 matches the load"
    )
