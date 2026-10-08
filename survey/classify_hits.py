#!/usr/bin/env python3
"""Assign every flagged row from a cgsurvey run to one class and count them.

    python3 classify_hits.py results/hits.tsv results/survey.tsv > class_counts.csv

Classes, as defined in the manuscript:
  A  convention: strand '?', or a feature crossing the origin of a circular sequence
  B  documented exception: a child outside its parent, where the child or the
     parent carries an exception=, part= or partial= attribute
  C  coordinate rule broken: start greater than end, or a row past the sequence
     end that does not fit the origin rule
  D  unexplained: a child outside its parent with none of the class B attributes

A row with several flags is counted once, in its most serious class (D > C > B > A).
Origin-crossing rows are counted in survey.tsv, not listed in hits.tsv, so they are
added to class A from there. Uses only the Python standard library.
"""
import collections
import csv
import sys

RANK = {"A": 1, "B": 2, "C": 3, "D": 4}


def row_class(hit: dict) -> str:
    kind = hit["kind"]
    if kind in ("inverted", "out-of-range"):
        return "C"
    if kind.startswith("escapes-parent"):
        return "B" if hit["excused_by"] else "D"
    if kind == "strand-?":
        return "A"
    raise SystemExit(f"unexpected flag kind: {kind}")


def main(hits_path: str, survey_path: str) -> None:
    survey = {r["accession"]: r for r in csv.DictReader(open(survey_path), delimiter="\t")}
    worst: dict[tuple, str] = {}
    for hit in csv.DictReader(open(hits_path), delimiter="\t"):
        key = (hit["accession"], hit["line"])
        c = row_class(hit)
        if key not in worst or RANK[c] > RANK[worst[key]]:
            worst[key] = c
    counts = collections.defaultdict(collections.Counter)
    for (acc, _line), c in worst.items():
        counts[acc][c] += 1
    for acc, r in survey.items():
        counts[acc]["A"] += int(r["origin_spanning_convention"])
    out = csv.writer(sys.stdout, lineterminator="\n")
    out.writerow(["species", "accession", "A_convention", "B_exception_attribute",
                  "C_breaks_coordinate_rule", "D_no_explanation", "total"])
    for acc, r in survey.items():
        c = counts[acc]
        out.writerow([r["species"], acc, c["A"], c["B"], c["C"], c["D"], sum(c.values())])
    totals = {k: sum(counts[a][k] for a in survey) for k in "ABCD"}
    flagged = sum(1 for a in survey if sum(counts[a].values()))
    print(f"# {len(survey)} genomes, {flagged} with flagged rows; class totals "
          f"A={totals['A']} B={totals['B']} C={totals['C']} D={totals['D']}, "
          f"all={sum(totals.values())}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
