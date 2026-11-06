#!/usr/bin/env python3
"""Compare a loader's result again, without running the loader again.

    python3 benchmark/recompare.py ACCESSION --tool chado --results RESULTS --work WORK

For a row of compare_loaders.tsv whose comparison did not finish: the second
run's human Chado load ran into its 24-hour limit, left nothing, and the
comparison of 3.2 million lost features with its 15 MB log then took longer
than the load (compare.py, ``_Names``). The tool's own measurements in the row
stay as they are; only the comparison's columns are filled in, its hard-case
rows added and its JSON kept, as run_compare.py does. The table's old line is
kept in RESULTS/superseded/.

Only Chado, and only where its read-back found nothing: then the store it
read is empty by definition, and an empty one stands in for the one removed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402
import gffmodel  # noqa: E402
import run_compare  # noqa: E402


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--tool", required=True, choices=["chado"])
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args(argv)
    args.truth = args.results / "truth" / f"{args.accession}.sqlite"
    table = args.results / "compare_loaders.tsv"
    lines = table.read_text().splitlines(keepends=True)
    header = lines[0].rstrip("\n").split("\t")
    at = [i for i, line in enumerate(lines[1:], 1)
          if line.rstrip("\n").split("\t")[:2] == [args.accession, args.tool]]
    if len(at) != 1:
        print(f"{args.accession} {args.tool}: {len(at)} rows in {table}, need one", file=sys.stderr)
        return 1
    row = dict(zip(header, lines[at[0]].rstrip("\n").split("\t"), strict=True))
    if row["truth"] != "NA":
        print(f"{args.accession} {args.tool}: already compared", file=sys.stderr)
        return 1
    steps = args.work / "compare" / args.tool / args.accession / "steps"
    stored = steps.parent / "stored.sqlite"
    read_back = json.loads((steps / "extract.out").read_text())
    if read_back.get("features") != 0:
        print(f"{args.accession} {args.tool}: the read-back found features; run the step "
              "again instead", file=sys.stderr)
        return 1
    gffmodel.create(stored).close()
    run_compare.compare_features(args, row, steps, stored, args.tool,
                                 "seqid,type,strand,segs,parents",
                                 [steps / "load.out", steps / "load.err"])
    stored.unlink(missing_ok=True)
    if row["truth"] == "NA":
        print(f"{args.accession} {args.tool}: the comparison failed again; see {steps}",
              file=sys.stderr)
        return 1
    aside = args.results / "superseded" / datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    aside.mkdir(parents=True)
    shutil.copy2(table, aside / table.name)
    lines[at[0]] = "\t".join(bc._cell(row[name]) for name in header) + "\n"
    table.write_text("".join(lines))
    print(json.dumps(row))
    return 0


if __name__ == "__main__":
    sys.exit(main())
