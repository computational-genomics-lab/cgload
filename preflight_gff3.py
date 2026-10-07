#!/usr/bin/env python3
"""Report every structural violation in a GFF3 file, in one pass.

Loading a 1.2-million-feature genome to discover one bad row, fixing it, and
loading again is a poor way to find out what a file contains. This walks the
file once, checks four things, and prints everything wrong.

It shares no code with cgload's parser, deliberately (D-008): a pre-flight
check that used the parser could only find what the parser already refuses.

    python preflight_gff3.py genomic.gff

Checks, in order of how quietly they fail:

1. **A child must lie inside its parent's envelope.** This is the one that
   matters. A row whose start exceeds its end fails loudly and you fix it; a
   row whose end is merely *wrong* loads cleanly and corrupts the database.
   Only the parent-child relationship catches the second kind.
2. Start must not exceed end.
3. Coordinates must lie within the declared ##sequence-region.
4. Strand must be one of + - . ?
"""

from __future__ import annotations

import gzip
import re
import sys
from collections import defaultdict


def opener(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    return gzip.open(path, "rt", errors="replace") if magic == b"\x1f\x8b" else open(
        path, errors="replace"
    )


ID = re.compile(r"(?:^|;)ID=([^;]*)")
PARENT = re.compile(r"(?:^|;)Parent=([^;]*)")


def main(path: str) -> int:
    regions: dict[str, int] = {}
    # id -> (line, seqid, type, start, end, strand)
    features: dict[str, tuple] = {}
    # (parent_id) -> list of (line, type, start, end)
    children: dict[str, list] = defaultdict(list)
    problems: list[tuple[int, str, str]] = []

    with opener(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    regions[parts[1]] = int(parts[3])
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue
            seqid, _source, ftype, raw_start, raw_end, _score, strand = columns[:7]
            try:
                start, end = int(raw_start), int(raw_end)
            except ValueError:
                problems.append((number, "coordinate", f"{raw_start}..{raw_end} is not numeric"))
                continue

            if start > end:
                problems.append(
                    (number, "inverted", f"{ftype} start {start:,} exceeds end {end:,}")
                )
            if strand not in ("+", "-", ".", "?"):
                problems.append((number, "strand", f"{strand!r} is not + - . or ?"))
            length = regions.get(seqid)
            if length is not None and (start > length or end > length):
                problems.append(
                    (
                        number,
                        "out-of-range",
                        f"{ftype} {start:,}..{end:,} exceeds {seqid} length {length:,}",
                    )
                )

            found = ID.search(columns[8])
            if found:
                fid = found.group(1)
                prev = features.get(fid)
                if prev is None:
                    features[fid] = (number, seqid, ftype, start, end, strand)
                else:
                    # GFF3 lets one feature span several rows that share an ID
                    # (a trans-spliced gene, a CDS in pieces). Its span is the
                    # union of those rows. Keeping only the last row made the
                    # children of the earlier rows look like they escaped.
                    pline, pseq, ptype, pstart, pend, pstrand = prev
                    features[fid] = (pline, pseq, ptype, min(pstart, start),
                                     max(pend, end), pstrand)
            parent = PARENT.search(columns[8])
            if parent:
                for one in parent.group(1).split(","):
                    children[one.strip()].append((number, ftype, start, end))

    # The check that catches silent corruption.
    escapes = 0
    for parent_id, kids in children.items():
        if parent_id not in features:
            continue
        pline, _seqid, ptype, pstart, pend, pstrand = features[parent_id]
        for kline, ktype, kstart, kend in kids:
            if kstart < pstart or kend > pend:
                escapes += 1
                problems.append(
                    (
                        kline,
                        "escapes-parent",
                        f"{ktype} {kstart:,}..{kend:,} lies outside {ptype} "
                        f"{parent_id} {pstart:,}..{pend:,} (line {pline}, strand {pstrand!r})",
                    )
                )

    if not problems:
        print(f"{path}: no structural problems found")
        return 0

    problems.sort()
    by_kind: dict[str, int] = defaultdict(int)
    for _line, kind, _detail in problems:
        by_kind[kind] += 1

    print(f"{path}: {len(problems):,} structural problems\n")
    for kind, count in sorted(by_kind.items(), key=lambda x: -x[1]):
        print(f"  {kind:16s} {count:,}")
    print()
    for line, kind, detail in problems[:40]:
        print(f"  line {line:>9,}  {kind:16s} {detail}")
    if len(problems) > 40:
        print(f"  ... and {len(problems) - 40:,} more")

    if escapes:
        print(
            "\nNOTE: 'escapes-parent' is the dangerous one. Those rows pass every\n"
            "per-row check and load cleanly, with a parent envelope that is wrong.\n"
            "Nothing downstream can detect it afterwards."
        )
    return 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
