"""Feature counts, defined once.

The feature-level count is not ``COUNT(*)`` and not plain
``COUNT(DISTINCT source_id)``. Real input contains rows with no ID at all --
AUGUSTUS emits intron, start_codon and stop_codon rows carrying only a Parent --
and each of those is its own feature, because with no ID there is nothing to
group segments by.

So::

    features = COUNT(DISTINCT source_id)      # named, segments collapsed
             + COUNT(*) WHERE source_id IS NULL   # anonymous, one each

Defined here rather than inline so that ``stats``, ``verify`` and ``query``
cannot drift apart on the one number the project's headline claim rests on.
"""

from __future__ import annotations

from sqlalchemy import Integer, case, func, select

from cgload.db.schema import feature


def feature_count_expression():
    """SQL expression for the number of features in the selected rows."""
    return func.count(func.distinct(feature.c.source_id)) + func.coalesce(
        func.sum(case((feature.c.source_id.is_(None), 1), else_=0)), 0
    ).cast(Integer)


def counts_by(column, assembly_id: int):
    """Feature counts grouped by ``source_type`` or ``feature_type``."""
    return (
        select(column, feature_count_expression().label("features"))
        .where(feature.c.assembly_id == assembly_id)
        .group_by(column)
        .order_by(feature_count_expression().desc())
    )
