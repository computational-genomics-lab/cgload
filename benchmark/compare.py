#!/usr/bin/env python3
"""Compare what a loader stored with the reference count of the same file.

    python3 benchmark/compare.py TRUTH.sqlite STORED.sqlite --tool NAME [options]

Both databases hold the table described in gffmodel.py: TRUTH.sqlite from
truth.py run on the input, STORED.sqlite from the loader's own store (or, for
cgload, from truth.py run on what cgload exports).

Every feature of the file ends in exactly one of three states:

* kept: the loader holds it once, and every compared field equals the file's;
* changed: the loader holds it, but a compared field differs, or it holds it as
  several separate records;
* lost: the loader does not hold it.

A change the loader itself reported, as cgload reports each envelope it
repairs, is counted as declared (``--declared``). A loss or change is counted as
announced when the loader's own output (``--log``) names the feature, or names
its type on a line that says the type was skipped, dropped, ignored, rejected
or unknown; otherwise it is silent. Features the loader holds that the file
does not are counted as extra, by type.

Standard library only; no loader code.
"""

from __future__ import annotations

import argparse
import json
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gffmodel  # noqa: E402

FIELDS = ("seqid", "type", "strand", "segs", "parents", "attrs")
ANNOUNCING = re.compile(r"skip|drop|ignor|reject|unknown|discard|remov|not load|could not",
                        re.IGNORECASE)


def _log_text(paths: list[Path]) -> tuple[str, list[str]]:
    text = []
    for path in paths:
        if path.exists():
            text.append(path.read_text(errors="replace"))
    joined = "\n".join(text)
    announcing = [line for line in joined.splitlines() if ANNOUNCING.search(line)]
    return joined, announcing


def compare(truth: Path, stored: Path, *, fields=FIELDS, types: set[str] | None = None,
            declared: set[str] = frozenset(), logs: list[Path] = (), examples: int = 5) -> dict:
    connection = sqlite3.connect(truth)
    connection.execute("ATTACH DATABASE ? AS s", (str(stored),))
    log_text, announcing_lines = _log_text(list(logs))
    announcing_text = "\n".join(announcing_lines)

    def announced(key: str, ftype: str) -> bool:
        if key and key in log_text:
            return True
        return bool(ftype) and re.search(rf"\b{re.escape(ftype)}\b", announcing_text) is not None

    totals = Counter()
    by_type: dict[str, Counter] = defaultdict(Counter)
    by_case: dict[str, Counter] = defaultdict(Counter)
    field_diff = Counter()
    case_field_diff: dict[str, Counter] = defaultdict(Counter)
    shown: dict[str, list[str]] = defaultdict(list)

    columns = ", ".join(f"t.{name}" if name != "end" else 't."end"' for name in (
        "key", "seqid", "type", "strand", "segs", "parents", "attrs", "flags"))
    stored_columns = ", ".join(f"x.{name}" for name in (
        "key", "seqid", "type", "strand", "segs", "parents", "attrs", "records"))
    query = (f"SELECT {columns}, {stored_columns} FROM main.feature t "  # noqa: S608 - fixed names
             "LEFT JOIN s.feature x ON x.key = t.key")
    for row in connection.execute(query):
        (key, seqid, ftype, strand, segs, parents, attrs, flags,
         s_key, s_seqid, s_type, s_strand, s_segs, s_parents, s_attrs, s_records) = row
        if types is not None and ftype not in types:
            totals["not_modelled"] += 1
            by_type[ftype]["not_modelled"] += 1
            continue
        cases = [name for name in flags.split(",") if name] if flags else []
        totals["truth"] += 1
        by_type[ftype]["truth"] += 1
        truth_parents = set(parents.split(",")) if parents else set()
        totals["truth_links"] += len(truth_parents)
        for name in cases:
            by_case[name]["truth"] += 1
        if s_key is None:
            state = "lost"
            loud = announced(key, ftype)
            diffs = []
        else:
            stored_parents = set(s_parents.split(",")) if s_parents else set()
            totals["kept_links"] += len(truth_parents & stored_parents)
            values = {"seqid": (seqid, s_seqid), "type": (ftype, s_type),
                      "strand": (strand, s_strand), "segs": (segs, s_segs),
                      "parents": (parents or "", s_parents or ""), "attrs": (attrs, s_attrs)}
            diffs = [name for name in fields if values[name][0] != values[name][1]]
            if (s_records or 1) > 1:
                diffs.append("split")
            if not diffs:
                state = "kept"
                loud = False
            elif key in declared:
                state = "changed_declared"
                loud = True
            else:
                state = "changed"
                loud = announced(key, ftype)
        totals[state] += 1
        by_type[ftype][state] += 1
        if state in ("lost", "changed"):
            totals[f"{state}_{'announced' if loud else 'silent'}"] += 1
        for name in diffs:
            field_diff[name] += 1
        for name in cases:
            by_case[name][state] += 1
            if state in ("lost", "changed"):
                by_case[name][f"{state}_{'announced' if loud else 'silent'}"] += 1
            for diff in diffs:
                case_field_diff[name][diff] += 1
        label = state if not diffs else f"{state}:{'+'.join(diffs)}"
        if state != "kept" and len(shown[label]) < examples:
            shown[label].append(key)

    stored_types = Counter(dict(connection.execute(
        "SELECT type, count(*) FROM s.feature GROUP BY type").fetchall()))
    extra_types = Counter(dict(connection.execute(
        "SELECT x.type, count(*) FROM s.feature x LEFT JOIN main.feature t ON t.key = x.key "
        "WHERE t.key IS NULL GROUP BY x.type").fetchall()))
    stored_links = connection.execute(
        "SELECT coalesce(sum(length(parents) - length(replace(parents, ',', '')) + 1), 0) "
        "FROM s.feature WHERE parents != ''").fetchone()[0]
    connection.close()
    return {
        "fields_compared": list(fields),
        "totals": dict(totals),
        "stored_features": sum(stored_types.values()),
        "stored_links": stored_links,
        "extra_by_type": dict(extra_types),
        "field_differences": dict(field_diff),
        "by_type": {name: dict(counts) for name, counts in sorted(by_type.items())},
        "by_case": {name: dict(by_case[name]) for name in gffmodel.HARD_CASES if name in by_case},
        "case_field_differences": {name: dict(c) for name, c in case_field_diff.items()},
        "examples": dict(shown),
    }


CASE_COLUMNS = ("truth", "kept", "changed", "changed_declared", "lost",
                "changed_silent", "lost_silent")


def case_rows(result: dict, genome: str, tool: str) -> list[list]:
    rows = []
    for name, counts in result["by_case"].items():
        rows.append([genome, tool, name] + [counts.get(c, 0) for c in CASE_COLUMNS]
                    + [",".join(f"{k}:{v}" for k, v in sorted(
                        result["case_field_differences"].get(name, {}).items()))])
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("truth", type=Path)
    parser.add_argument("stored", type=Path)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--genome", default="")
    parser.add_argument("--fields", default=",".join(FIELDS),
                        help="fields to compare (default: all of %(default)s)")
    parser.add_argument("--types", type=Path,
                        help="file with one feature type per line: only these are compared, "
                             "for a loader that models only some types")
    parser.add_argument("--declared", type=Path,
                        help="file with one feature key per line that the loader reported changing")
    parser.add_argument("--log", type=Path, action="append", default=[],
                        help="the loader's output, to tell announced from silent (repeatable)")
    parser.add_argument("--json", type=Path, help="write the full result here")
    parser.add_argument("--cases-tsv", type=Path, help="append one row per hard case here")
    args = parser.parse_args(argv)

    fields = tuple(name for name in args.fields.split(",") if name)
    unknown = set(fields) - set(FIELDS)
    if unknown:
        parser.error(f"unknown field(s): {', '.join(sorted(unknown))}")
    types = None
    if args.types:
        types = {line.strip() for line in args.types.read_text().splitlines() if line.strip()}
    declared = set()
    if args.declared and args.declared.exists():
        declared = {line.strip() for line in args.declared.read_text().splitlines() if line.strip()}
    result = compare(args.truth, args.stored, fields=fields, types=types, declared=declared,
                     logs=args.log)
    result.update({"genome": args.genome, "tool": args.tool})
    if args.json:
        args.json.write_text(json.dumps(result, indent=1) + "\n")
    if args.cases_tsv:
        gffmodel.locked_append(
            args.cases_tsv,
            "\t".join(("genome", "tool", "case", *CASE_COLUMNS, "fields_changed")),
            ("\t".join(map(str, row)) for row in case_rows(result, args.genome, args.tool)))
    totals = result["totals"]
    print(f"{args.tool} {args.genome}: {totals.get('truth', 0):,} features in the file; "
          f"kept {totals.get('kept', 0):,}, changed {totals.get('changed', 0):,} "
          f"(+{totals.get('changed_declared', 0):,} declared), lost {totals.get('lost', 0):,}; "
          f"extra {sum(result['extra_by_type'].values()):,}; "
          f"links {totals.get('kept_links', 0):,} of {totals.get('truth_links', 0):,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
