"""D-094: the bins that let a region query read only the rows near the region."""

from __future__ import annotations

import random

from cgload.db.bins import DEPTH, MIN_SHIFT, bin_of, bins_overlapping

BAI = {"min_shift": 14, "depth": 5}


def test_the_bam_index_parameters_give_the_bins_of_the_sam_specification():
    # SAM v1 section 5.3: reg2bin over zero-based [beg, end) with 16 kb leaves.
    assert bin_of(1, 1, **BAI) == 4681
    assert bin_of(16_385, 32_768, **BAI) == 4682
    assert bin_of(1, 16_385, **BAI) == 585  # crosses a 16 kb boundary: 128 kb level
    assert bin_of(1, 1 << 29, **BAI) == 0
    assert bin_of(1, 1 << 26, **BAI) == 1


def test_a_row_lies_in_one_of_the_bins_of_every_region_it_overlaps():
    rnd = random.Random(94)
    limit = 1 << (MIN_SHIFT + 3 * DEPTH)
    for _ in range(20_000):
        # Rows of any size, regions up to a few megabases, near each other.
        size = rnd.choice([10, 1_000, 100_000, 10_000_000, 1 << 34])
        start = rnd.randint(1, limit - size)
        end = start + rnd.randint(0, size)
        reach = rnd.choice([10, 1_000, 100_000, 3_000_000])
        q_start = rnd.randint(max(1, start - reach), end + reach)
        q_end = q_start + rnd.randint(0, reach)
        overlaps = start <= q_end and end >= q_start
        if overlaps:
            found = bins_overlapping(q_start, q_end)
            assert bin_of(start, end) in found, (start, end, q_start, q_end)


def test_a_megabase_region_needs_a_few_dozen_bins():
    assert len(bins_overlapping(100_000_001, 101_000_000)) < 100
    assert bins_overlapping(1, 1)[0] == 0
