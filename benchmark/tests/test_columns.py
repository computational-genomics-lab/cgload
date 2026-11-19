"""columns.py: the row-by-row comparison of GFF3 columns 2, 6 and 8."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import columns  # noqa: E402

ROWS = [
    "chr1\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=g1",
    "chr1\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=m1;Parent=g1",
    "chr1\tRefSeq\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=m1",
    "chr1\tGnomon\tCDS\t500\t900\t.\t+\t2\tID=c1;Parent=m1",
    "chr1\tAUGUSTUS\tintron\t301\t499\t0.75\t+\t.\tParent=m1",
    "chr1\tBestRefSeq%2CGnomon\tgene\t1000\t1900\t12.5\t-\t.\tID=g2",
]


def _write(path: Path, rows: list[str]) -> Path:
    path.write_text("##gff-version 3\n" + "\n".join(rows) + "\n")
    return path


def test_identical_files_have_no_difference(tmp_path):
    counts = columns.compare(_write(tmp_path / "a.gff3", ROWS),
                             _write(tmp_path / "b.gff3", ROWS))["counts"]
    assert counts["rows_compared"] == len(ROWS)
    assert counts["source_differs"] == counts["score_differs"] == counts["phase_differs"] == 0


def test_each_column_is_compared_row_by_row(tmp_path):
    changed = list(ROWS)
    changed[3] = changed[3].replace("Gnomon", "RefSeq")      # a later row's column 2
    changed[4] = changed[4].replace("0.75", "0.8")          # an anonymous row's score
    changed[2] = changed[2].replace("\t0\t", "\t1\t")        # a phase
    counts = columns.compare(_write(tmp_path / "a.gff3", ROWS),
                             _write(tmp_path / "b.gff3", changed))["counts"]
    assert (counts["source_differs"], counts["score_differs"], counts["phase_differs"]) == (1, 1, 1)


def test_encoding_and_number_format_alone_are_kept_apart(tmp_path):
    same = list(ROWS)
    same[5] = same[5].replace("BestRefSeq%2CGnomon", "BestRefSeq,Gnomon").replace("12.5", "12.50")
    counts = columns.compare(_write(tmp_path / "a.gff3", ROWS),
                             _write(tmp_path / "b.gff3", same))["counts"]
    assert counts["source_differs"] == counts["score_differs"] == 0
    assert counts["source_encoding_only"] == counts["score_text_only"] == 1


def test_a_feature_with_rows_missing_is_counted_not_compared(tmp_path):
    counts = columns.compare(_write(tmp_path / "a.gff3", ROWS),
                             _write(tmp_path / "b.gff3", ROWS[:3] + ROWS[4:]))["counts"]
    assert counts["features_row_count_differs"] == 1


def test_the_values_that_cgload_0_2_1_wrote_differently_are_counted_from_the_input(tmp_path):
    rows = list(ROWS)
    rows[4] = rows[4].replace("0.75", "1234567.5")         # %g keeps six digits
    counts = columns.compare(_write(tmp_path / "a.gff3", rows),
                             _write(tmp_path / "b.gff3", rows))["counts"]
    assert counts["input_score_beyond_six_digits"] == 1
    assert counts["input_source_differs_from_first_row"] == 1   # c1's Gnomon row
