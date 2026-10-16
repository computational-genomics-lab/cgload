"""Genomic bins, so that a region query reads only the rows near the region (D-094).

The hierarchical binning scheme of the UCSC browser and of the SAM/BAM index,
in the generalised form of the CSI index: each level divides a sequence into
bins eight times smaller than the level above, the smallest 16 kb. A row is
stored with the smallest bin that holds it whole, so the rows that can overlap
a region are exactly those whose bin is one of the few dozen bins that overlap
the region at any level. With ``min_shift`` 14 and depth 8 the scheme covers
coordinates up to 2^38, beyond any sequence and beyond a circular sequence's
rows that run past its origin.

Without bins the index on (sequence, start, end) bounded a region query only
from above: every row starting before the region's end was read, half a
chromosome for a region in its middle, and PostgreSQL, unable to estimate how
few of them overlapped, planned the query as scans of the whole table.

The two functions are those of the CSI specification, written for 1-based
inclusive coordinates.
"""

from __future__ import annotations

MIN_SHIFT = 14
DEPTH = 8
#: Beyond this many bins a region covers so much of its sequence that reading
#: the sequence's rows is cheaper than matching a long list of bins.
MAX_QUERY_BINS = 500


def bin_of(start: int, end: int, *, min_shift: int = MIN_SHIFT, depth: int = DEPTH) -> int:
    """The bin of a row spanning ``start``..``end`` (1-based, inclusive)."""
    begin, last = start - 1, end - 1
    shift = min_shift
    offset = ((1 << depth * 3) - 1) // 7
    for level in range(depth, 0, -1):
        if begin >> shift == last >> shift:
            return offset + (begin >> shift)
        shift += 3
        offset -= 1 << ((level - 1) * 3)
    return 0


def bins_overlapping(
    start: int, end: int, *, min_shift: int = MIN_SHIFT, depth: int = DEPTH
) -> list[int]:
    """Every bin that can hold a row overlapping ``start``..``end``."""
    begin, last = start - 1, end - 1
    found: list[int] = []
    shift = min_shift + depth * 3
    offset = 0
    for level in range(depth + 1):
        found.extend(range(offset + (begin >> shift), offset + (last >> shift) + 1))
        shift -= 3
        offset += 1 << (level * 3)
    return found
