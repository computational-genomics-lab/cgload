#!/usr/bin/env python3
"""Turn the result tables into the tables the paper reports.

    python3 benchmark/summarise.py RESULTS_DIR

Reads the tables the runners appended to and writes RESULTS_DIR/tables/:

  scale.tsv        per genome: size, and for SQLite and PostgreSQL the load's
                   outcome, time, peak memory, database size, verify, round trip
  tools.tsv        per loader: genomes loaded or failed, features kept, changed
                   and lost over all genomes, parent links kept
  cases.tsv        per loader and hard case: kept, changed, lost, summed
  faults.tsv       per fault: how many genomes it was detected in, against
                   what was expected
  kills.tsv        per loader and genome: what interrupted loads left behind
  timing.tsv       per genome, tool, task and view: median first run and
                   median repeated run, over targets
  agreement.tsv    per genome and in total: cgload against the survey
  headline.json    the numbers the abstract and results quote

Every number here is computed from a file in RESULTS_DIR, so each one in the
paper can be traced to a row. Standard library only.
"""

from __future__ import annotations

import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402


def read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value, kind=float):
    try:
        return kind(value)
    except (TypeError, ValueError):
        return None


def write(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    columns = list(rows[0])
    for row in rows[1:]:
        columns.extend(key for key in row if key not in columns)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def scale(results: Path) -> list[dict]:
    genomes = {g["accession"]: g for g in bc.genomes()}
    runs = defaultdict(dict)
    for row in read(results / "cgload_scale.tsv"):
        runs[row["accession"]][row["backend"]] = row
    out = []
    for accession, genome in genomes.items():
        summary_path = results / "truth" / f"{accession}.json"
        truth = json.loads(summary_path.read_text()) if summary_path.exists() else {}
        row = {"accession": accession, "organism": genome["organism"],
               "clade": genome["clade"], "rows": truth.get("rows"),
               "features": truth.get("features"), "parent_links": truth.get("parent_links")}
        for backend in ("sqlite", "postgres"):
            run = runs[accession].get(backend, {})
            mb = number(run.get("db_bytes"))
            row.update({
                f"{backend}_outcome": run.get("outcome", "not run"),
                f"{backend}_load_s": run.get("load_wall_s"),
                f"{backend}_peak_mb": run.get("load_peak_mb"),
                f"{backend}_db_mb": round(mb / 1e6, 1) if mb else None,
                f"{backend}_verify": "pass" if run.get("verify_exit") == "0" else
                (run.get("verify_exit") and "fail") or None,
                f"{backend}_rt_kept": run.get("rt_kept"),
                f"{backend}_rt_changed_declared": run.get("rt_changed_declared"),
                f"{backend}_rt_changed": run.get("rt_changed"),
                f"{backend}_rt_lost": run.get("rt_lost"),
            })
        data_mb = number(runs[accession].get("sqlite", {}).get("data_bytes"))
        row["data_dir_mb"] = round(data_mb / 1e6, 1) if data_mb else None
        row["refusal"] = runs[accession].get("sqlite", {}).get("refusal", "")
        out.append(row)
    return out


def tools(results: Path) -> list[dict]:
    per_tool = defaultdict(Counter)
    units = {}
    for row in read(results / "cgload_scale.tsv"):
        tool = f"cgload-{row['backend']}"
        units[tool] = "features"
        per_tool[tool]["genomes"] += 1
        per_tool[tool][row.get("outcome") or "failed"] += 1
        for key in ("truth", "kept", "changed", "changed_declared", "lost", "lost_silent",
                    "changed_silent", "truth_links", "kept_links"):
            per_tool[tool][key] += number(row.get(f"rt_{key}"), int) or 0
    for row in read(results / "compare_loaders.tsv"):
        tool = row["tool"]
        if row.get("unit") not in (None, "", "NA"):
            units[tool] = row["unit"]
        per_tool[tool]["genomes"] += 1
        per_tool[tool][row.get("outcome") or "failed"] += 1
        for key in ("truth", "kept", "changed", "changed_declared", "lost", "extra",
                    "not_modelled", "parts_unattached", "truth_links", "kept_links",
                    "lost_silent", "changed_silent", "genes_truth", "genes_held"):
            per_tool[tool][key] += number(row.get(key), int) or 0
    out = []
    for tool, counts in sorted(per_tool.items()):
        out.append({"tool": tool, "unit": units.get(tool, ""), **{key: counts[key] for key in (
            "genomes", "loaded", "loaded_with_repair", "refused", "failed", "timed_out",
            "truth", "kept", "changed", "changed_declared", "lost", "lost_silent",
            "changed_silent", "extra", "not_modelled", "parts_unattached", "truth_links",
            "kept_links", "genes_truth", "genes_held")}})
    return out


def cases(results: Path) -> list[dict]:
    totals = defaultdict(Counter)
    fields = defaultdict(Counter)
    for name in ("cgload_fidelity_cases.tsv", "compare_cases.tsv"):
        for row in read(results / name):
            key = (row["tool"], row["case"])
            for column in ("truth", "kept", "changed", "changed_declared", "lost",
                           "changed_silent", "lost_silent"):
                totals[key][column] += number(row.get(column), int) or 0
            for item in (row.get("fields_changed") or "").split(","):
                field, _, count = item.partition(":")
                if field and field != "NA" and count:
                    fields[key][field] += int(count)
    return [{"tool": tool, "case": case, **dict(counts),
             "fields_changed": ",".join(f"{field}:{n}" for field, n in
                                        sorted(fields[(tool, case)].items()))}
            for (tool, case), counts in sorted(totals.items())]


def faults(results: Path) -> list[dict]:
    grouped = defaultdict(list)
    for row in read(results / "faults.tsv"):
        grouped[(row["fault_id"], row["fault"], row["layer"], row["expected"])].append(row)
    out = []
    for (fault_id, fault, layer, expected), rows in sorted(grouped.items()):
        out.append({
            "fault_id": fault_id, "fault": fault, "layer": layer, "expected": expected,
            "genomes": len(rows),
            "detected": sum(r["observed"] == "detected" for r in rows),
            "as_expected": sum(r["as_expected"] == "True" for r in rows),
            "how": "; ".join(sorted({r["message"][:120] for r in rows if r["message"] not in
                                     ("", "NA")}))[:400],
        })
    return out


KILL_STATES = ("empty", "partial", "complete", "intact", "damaged", "no_file", "unreadable",
               "setup_failed")


def kills(results: Path) -> list[dict]:
    grouped = defaultdict(list)
    for row in read(results / "kills.tsv"):
        grouped[(row["tool"], row["accession"])].append(row)
    out = []
    for (tool, accession), rows in sorted(grouped.items()):
        states = Counter(r["state"] for r in rows)
        partial = [r for r in rows if r["state"] == "partial"]
        fractions = [number(r["records_left"]) / number(r["records_complete"])
                     for r in partial if number(r["records_left"]) is not None
                     and number(r["records_complete"])]
        out.append({"tool": tool, "accession": accession, "trials": len(rows),
                    **{f"state_{name}": states.get(name, 0) for name in KILL_STATES},
                    "partial_opens": sum(r["opens"] == "True" for r in partial),
                    "file_left_does_not_open": sum(r["opens"] == "False" and r["state"] not in
                                                   ("no_file", "empty") for r in rows),
                    "partial_records_min_pct": round(100 * min(fractions), 1) if fractions
                    else None,
                    "partial_records_max_pct": round(100 * max(fractions), 1) if fractions
                    else None,
                    "finished_before_kill": sum(r["finished_first"] == "True" for r in rows),
                    "orphan_directory": sum(r.get("orphan_directory") == "True" for r in rows),
                    "reload_failed": sum(r["reload_ok"] == "False" for r in rows),
                    "verify_failed_after": sum(r["verify_after"] == "fail" for r in rows),
                    "max_leftover_mb": round(max((number(r["leftover_bytes"]) or 0)
                                                 for r in rows) / 1e6, 1)})
    return out


def timing(results: Path) -> list[dict]:
    grouped = defaultdict(lambda: defaultdict(list))
    answers = defaultdict(set)
    dropped = defaultdict(set)
    for row in read(results / "query_timing.tsv"):
        wall = number(row["wall_s"])
        if wall is None:
            continue
        key = (row["accession"], row["view"], row["tool"], row["task"])
        grouped[key][row["target"]].append((int(row["run"]), wall))
        answers[key].add(row["answer"])
        dropped[key].add(row["caches_dropped"])
    out = []
    for key, targets in sorted(grouped.items()):
        firsts, repeats = [], []
        for runs in targets.values():
            runs.sort()
            firsts.append(runs[0][1])
            later = [wall for run, wall in runs if run > 1]
            if later:
                repeats.append(statistics.median(later))
        accession, view, tool, task = key
        out.append({"accession": accession, "view": view, "tool": tool, "task": task,
                    "targets": len(targets),
                    "first_run_median_s": round(statistics.median(firsts), 4),
                    "repeat_median_s": round(statistics.median(repeats), 4) if repeats else None,
                    "caches_dropped": ",".join(sorted(dropped[key]))})
    return out


def agreement(results: Path) -> tuple[list[dict], dict]:
    rows = read(results / "agreement_summary.tsv")
    totals = Counter()
    for row in rows:
        for key in ("survey_rows", "cgload_features", "both_accept", "both_object",
                    "survey_documented_cgload_objects", "survey_error_cgload_accepts",
                    "cgload_only"):
            totals[key] += number(row.get(key), int) or 0
    totals["genomes_compared"] = sum(row.get("same_file") == "True" for row in rows)
    return rows, dict(totals)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    results = Path(argv[0])
    tables = results / "tables"
    scale_rows = scale(results)
    write(tables / "scale.tsv", scale_rows)
    tool_rows = tools(results)
    write(tables / "tools.tsv", tool_rows)
    write(tables / "cases.tsv", cases(results))
    fault_rows = faults(results)
    write(tables / "faults.tsv", fault_rows)
    write(tables / "kills.tsv", kills(results))
    write(tables / "timing.tsv", timing(results))
    agreement_rows, agreement_totals = agreement(results)
    write(tables / "agreement.tsv", agreement_rows)
    loaded = [r for r in scale_rows if str(r.get("sqlite_outcome", "")).startswith("loaded")]
    headline = {
        "genomes": len(scale_rows),
        "loaded_sqlite": sum(1 for r in scale_rows if r["sqlite_outcome"] == "loaded"),
        "loaded_sqlite_with_repair": sum(1 for r in scale_rows
                                         if r["sqlite_outcome"] == "loaded_with_repair"),
        "loaded_postgres": sum(1 for r in scale_rows if r["postgres_outcome"] == "loaded"),
        "loaded_postgres_with_repair": sum(1 for r in scale_rows
                                           if r["postgres_outcome"] == "loaded_with_repair"),
        "rows_loaded": sum(number(r["rows"], int) or 0 for r in loaded),
        "features_loaded": sum(number(r["features"], int) or 0 for r in loaded),
        "slowest_load_s": max((number(r["sqlite_load_s"]) or 0 for r in scale_rows), default=None),
        "largest_peak_mb": max((number(r["sqlite_peak_mb"]) or 0 for r in scale_rows),
                               default=None),
        "faults_as_expected": f"{sum(r['as_expected'] for r in fault_rows)} of "
                              f"{sum(r['genomes'] for r in fault_rows)}",
        "tools": {r["tool"]: {k: r[k] for k in ("genomes", "loaded", "loaded_with_repair",
                                                 "refused", "failed", "timed_out", "kept",
                                                 "changed", "changed_declared", "lost",
                                                 "lost_silent", "changed_silent")}
                  for r in tool_rows},
        "agreement": agreement_totals,
    }
    (tables / "headline.json").write_text(json.dumps(headline, indent=1) + "\n")
    print(json.dumps(headline, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
