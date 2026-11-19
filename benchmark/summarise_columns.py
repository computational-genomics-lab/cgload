#!/usr/bin/env python3
"""Summarise run_columns.py's table as the paragraph the paper reports.

    python3 benchmark/summarise_columns.py RESULTS/columns.tsv

Checks that every genome of genomes.tsv is in the table and loaded, adds the
counts up, and prints the sentences for the paper's Results with the numbers
filled in. It exits 1, and says why, when a genome is missing or a column of
any row came back different.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

SUMMED = ("rows_input", "rows_compared", "source_differs", "source_encoding_only",
          "score_differs", "score_text_only", "phase_differs", "features_row_count_differs",
          "features_missing", "features_extra", "rows_span_differs",
          "input_score_beyond_six_digits", "input_source_differs_from_first_row")


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1:
        raise SystemExit(__doc__)
    with open(argv[0], newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    wanted = {genome["accession"] for genome in bc.genomes()}
    present = {row["accession"] for row in rows}
    problems = []
    if wanted - present:
        problems.append(f"missing genomes: {sorted(wanted - present)}")
    loaded = [row for row in rows if row["outcome"] in ("loaded", "loaded_with_repair")]
    if len(loaded) != len(rows):
        problems.append("not loaded: " + ", ".join(
            f"{row['accession']} ({row['outcome']})" for row in rows if row not in loaded))
    versions = {row["cgload_version"] for row in rows}
    total = {name: sum(int(row[name]) for row in loaded if row.get(name) not in (None, "", "NA"))
             for name in SUMMED}
    differs = total["source_differs"] + total["score_differs"] + total["phase_differs"]
    if differs:
        problems.append(f"{differs:,} values came back different; see columns/*.json")
    if total["features_row_count_differs"] or total["features_missing"]:
        problems.append("some features did not come back with all their rows")

    def n(name: str) -> str:
        return f"{total[name]:,}"

    print(f"cgload version(s): {', '.join(sorted(versions))}; genomes: {len(loaded)}")
    print()
    print(
        f"cgload {', '.join(sorted(versions))} wrote back columns 2, 6 and 8 of every row "
        f"in all {len(loaded)} genomes as the files had them: {n('rows_compared')} rows, "
        f"with {n('source_differs')} differences in column 2, {n('score_differs')} in "
        f"column 6 and {n('phase_differs')} in column 8. The input files hold "
        f"{n('input_score_beyond_six_digits')} scores with more than six significant "
        f"digits and {n('input_source_differs_from_first_row')} rows whose column 2 "
        f"differs from the first row of their feature, the only values that version 0.2.1 "
        f"would have written differently.")
    print()
    print({name: total[name] for name in SUMMED})
    for problem in problems:
        print("PROBLEM:", problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
