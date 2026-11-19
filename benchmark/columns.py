#!/usr/bin/env python3
"""Compare GFF3 columns 2, 6 and 8 (source, score, phase) row by row.

    python3 benchmark/columns.py INPUT.gff3[.gz] OUTPUT.gff3 [--json OUT.json]
        [--ignore-attribute KEY]

The round-trip comparison (truth.py, compare.py) reduces each feature to its
sequence, type, strand, row spans, parents and attributes; it never reads
columns 2, 6 and 8. This program closes that gap. It reads the input and a
loader's output the same way and compares, for every row, the program that
made it (column 2), its score (column 6) and its phase (column 8).

Rows are matched by feature and by position within the feature: the rows of
one feature, sorted by span, are compared in order. A feature is its ID, or,
for a row without one, the key ``gffmodel.anon_key`` gives it, as in truth.py.
A feature whose number of rows differs between the two files is counted
apart, and its rows are not compared.

What counts as a difference:

* column 2: different text after percent-decoding (``BestRefSeq%2CGnomon``
  and ``BestRefSeq,Gnomon`` are the same program; a difference in encoding
  alone is counted separately);
* column 6: a different number (``.`` is no score); ``1.0`` against ``1`` is
  counted separately, as text only;
* column 8: different text.

It also counts, from the input alone, the values that cgload 0.2.1 wrote
differently (decisions D-109 and D-110 in docs/decisions.md): scores with more
than six significant digits, which its export rounded, and rows whose column 2
differs from the first row of their feature, which its export replaced by the
first row's.

Standard library and gffmodel only: it shares no code with any loader.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gffmodel  # noqa: E402

#: How many examples of each kind of difference the JSON keeps.
EXAMPLES = 10


def read_rows(path: Path, ignore: frozenset[str] = frozenset()) -> dict[str, list[tuple]]:
    """Feature key -> its rows as (start, end, column 2, column 6, column 8), in
    the order of the file."""
    rows: dict[str, list[tuple]] = defaultdict(list)
    anonymous_seen: dict[str, int] = {}
    for _number, columns in gffmodel.iter_rows(path):
        seqid, source, ftype, start, end, score, strand, phase, column9 = columns
        try:
            start_n, end_n = int(start), int(end)
        except ValueError:
            continue
        attributes = gffmodel.parse_attributes(column9)
        for name in ignore:
            attributes.pop(name, None)
        identifier = attributes.get("ID", [None])[0]
        if identifier is None:
            parents = ",".join(sorted(set(attributes.get("Parent", []))))
            key = gffmodel.unique_anon(
                gffmodel.anon_key(seqid, ftype, start_n, end_n, strand, parents),
                anonymous_seen)
        else:
            key = identifier
        rows[key].append((start_n, end_n, source, score, phase))
    return rows


def _score(text: str) -> float | None:
    if text == ".":
        return None
    try:
        return float(text)
    except ValueError:
        return float("nan")      # never equal to anything, so it is reported


def rounded_by_g(text: str) -> bool:
    """Would ``f"{score:g}"`` (cgload 0.2.1's export) change this score's value?"""
    value = _score(text)
    return value is not None and value == value and float(f"{value:g}") != value


def compare(input_path: Path, output_path: Path,
            ignore: frozenset[str] = frozenset()) -> dict:
    before = read_rows(input_path)
    after = read_rows(output_path, ignore)
    counts = defaultdict(int)
    examples: dict[str, list] = defaultdict(list)

    def note(kind: str, key: str, row_in: tuple, row_out: tuple | None = None) -> None:
        counts[kind] += 1
        if len(examples[kind]) < EXAMPLES:
            examples[kind].append({"feature": key, "input": row_in, "output": row_out})

    counts["rows_input"] = sum(len(r) for r in before.values())
    for key, rows_in in before.items():
        first_source = unquote(rows_in[0][2])
        for row in rows_in:
            if rounded_by_g(row[3]):
                note("input_score_beyond_six_digits", key, row)
            if unquote(row[2]) != first_source:
                note("input_source_differs_from_first_row", key, row)
    counts["rows_output"] = sum(len(r) for r in after.values())
    for key, rows_in in before.items():
        rows_out = after.get(key)
        if rows_out is None:
            counts["features_missing"] += 1
            counts["rows_missing"] += len(rows_in)
            continue
        if len(rows_out) != len(rows_in):
            note("features_row_count_differs", key, (len(rows_in),), (len(rows_out),))
            continue
        for row_in, row_out in zip(sorted(rows_in), sorted(rows_out), strict=True):
            counts["rows_compared"] += 1
            if row_in[:2] != row_out[:2]:
                counts["rows_span_differs"] += 1
            source_in, source_out = row_in[2], row_out[2]
            if unquote(source_in) != unquote(source_out):
                note("source_differs", key, row_in, row_out)
            elif source_in != source_out:
                note("source_encoding_only", key, row_in, row_out)
            score_in, score_out = _score(row_in[3]), _score(row_out[3])
            if score_in != score_out:
                note("score_differs", key, row_in, row_out)
            elif row_in[3] != row_out[3]:
                note("score_text_only", key, row_in, row_out)
            if row_in[4] != row_out[4]:
                note("phase_differs", key, row_in, row_out)
    counts["features_extra"] = sum(1 for key in after if key not in before)
    for kind in ("input_score_beyond_six_digits", "input_source_differs_from_first_row",
                 "source_differs", "source_encoding_only", "score_differs",
                 "score_text_only", "phase_differs", "features_row_count_differs",
                 "features_missing", "rows_missing", "rows_span_differs"):
        counts.setdefault(kind, 0)
    return {"input": str(input_path), "output": str(output_path),
            "counts": dict(sorted(counts.items())), "examples": dict(examples)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--json", type=Path, help="also write the result here")
    parser.add_argument("--ignore-attribute", action="append", default=[], metavar="KEY",
                        help="an attribute the loader adds on export (repeatable)")
    args = parser.parse_args(argv)
    result = compare(args.input, args.output, frozenset(args.ignore_attribute))
    text = json.dumps(result, indent=1)
    if args.json:
        args.json.write_text(text + "\n")
    print(json.dumps(result["counts"], indent=1))
    differs = sum(result["counts"][kind] for kind in
                  ("source_differs", "score_differs", "phase_differs"))
    return 1 if differs else 0


if __name__ == "__main__":
    sys.exit(main())
