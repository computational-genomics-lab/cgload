#!/usr/bin/env python3
"""Set aside a finished step's results so that run.sh measures it again.

    python3 benchmark/forget.py RESULTS --tool chado --tool txdbmaker --timing

run.sh skips whatever its tables already hold, so a step whose results are void
is repeated only once they are gone. Nothing is deleted: every table this
changes is first copied whole to RESULTS/superseded/STAMP/, and each tool's
comparison files are moved there, so what was replaced, and why, stays with
the results (PROTOCOL.md, section 12).

* ``--tool NAME``: that loader's rows in compare_loaders.tsv and
  compare_cases.tsv, and its files in compare/.
* ``--timing``: the whole of query_timing.tsv.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from datetime import UTC, datetime
from pathlib import Path

TOOL_TABLES = ("compare_loaders.tsv", "compare_cases.tsv")


def drop_rows(path: Path, tools: set[str]) -> int:
    """Remove the tools' rows, leaving every other line exactly as it was."""
    lines = path.read_text().splitlines(keepends=True)
    column = lines[0].rstrip("\n").split("\t").index("tool")
    kept = [lines[0]] + [line for line in lines[1:]
                         if line.rstrip("\n").split("\t")[column] not in tools]
    path.write_text("".join(kept))
    return len(lines) - len(kept)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", type=Path)
    parser.add_argument("--tool", action="append", default=[])
    parser.add_argument("--timing", action="store_true")
    args = parser.parse_args(argv)
    results: Path = args.results
    if not (results / "cgload_scale.tsv").exists():
        print(f"{results} does not look like a benchmark results folder", file=sys.stderr)
        return 1
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")
    aside = results / "superseded" / stamp
    aside.mkdir(parents=True)
    tools = set(args.tool)
    if tools:
        for name in TOOL_TABLES:
            table = results / name
            if table.exists():
                shutil.copy2(table, aside / name)
                dropped = drop_rows(table, tools)
                print(f"{name}: {dropped} rows of {', '.join(sorted(tools))} set aside")
        moved = 0
        for tool in sorted(tools):
            for path in sorted((results / "compare").glob(f"*.{tool}.json")):
                (aside / "compare").mkdir(exist_ok=True)
                shutil.move(str(path), aside / "compare" / path.name)
                moved += 1
        print(f"compare/: {moved} files set aside")
    if args.timing and (results / "query_timing.tsv").exists():
        shutil.move(str(results / "query_timing.tsv"), aside / "query_timing.tsv")
        print("query_timing.tsv: set aside")
    print(f"kept in {aside}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
