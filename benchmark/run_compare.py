#!/usr/bin/env python3
"""Load one genome with one of the other loaders, and compare what it kept.

    python3 benchmark/run_compare.py ACCESSION --tool gffutils|gffutils-default|txdbmaker|chado \\
        --inputs DIR --work DIR --results DIR

The loaders and how each is run (PROTOCOL.md, section 6):

* ``gffutils``: ``create_db`` with ``merge_strategy="create_unique"``, the
  documented way to load a file whose rows repeat an ID without editing it.
* ``gffutils-default``: ``create_db`` with its default, ``"error"``, to record
  what a user who changes nothing gets. Only the outcome is recorded.
* ``txdbmaker``: ``makeTxDbFromGFF``, through run_txdb.R, with the sequences the
  file marks circular.
* ``chado``: python-chado's ``load_gff`` into a new Chado 1.31 database,
  through run_chado.py, in its own Python environment (BENCH_CHADO_PYTHON).

The load is the measured step, with a time limit (--timeout, default 12 h). A
finished load is read back and compared with the reference count
(RESULTS/truth/ACCESSION.sqlite). One row goes to RESULTS/compare_loaders.tsv
and the hard-case rows to RESULTS/compare_cases.tsv. A genome and tool already
in the table is skipped.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

COLUMNS = (
    "accession", "tool", "tool_version", "outcome", "exit_code", "wall_s", "cpu_s", "peak_mb",
    "error", "unit", "truth", "kept", "changed", "changed_declared", "lost", "extra",
    "not_modelled", "parts_unattached", "truth_links", "kept_links", "fields_changed",
    "lost_silent", "changed_silent", "genes_truth", "genes_held",
)


def outcome_of(measurement: dict) -> str:
    if measurement.get("timed_out"):
        return "timed_out"
    return "loaded" if measurement["exit_code"] == 0 else "failed"


def error_of(steps: Path, name: str) -> str:
    text = bc.text_of(steps, name)
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    for line in reversed(lines):
        if "Error" in line or "error" in line or "Exception" in line:
            return line[:400]
    return bc.first_lines("\n".join(lines[-3:]), 3)


def circular_seqids(truth: Path) -> list[str]:
    with sqlite3.connect(truth) as connection:
        return [row[0] for row in connection.execute(
            "SELECT seqid FROM sequence WHERE circular = 1 ORDER BY seqid")]


def fill_from_compare(row: dict, result: dict) -> None:
    totals = result["totals"]
    row.update({
        "unit": "features", "truth": totals.get("truth", 0), "kept": totals.get("kept", 0),
        "changed": totals.get("changed", 0),
        "changed_declared": totals.get("changed_declared", 0), "lost": totals.get("lost", 0),
        "extra": sum(result["extra_by_type"].values()),
        "not_modelled": totals.get("not_modelled", 0),
        "truth_links": totals.get("truth_links", 0), "kept_links": totals.get("kept_links", 0),
        "fields_changed": ",".join(f"{k}:{v}" for k, v in
                                   sorted(result["field_differences"].items())),
        "lost_silent": totals.get("lost_silent", 0),
        "changed_silent": totals.get("changed_silent", 0),
    })


def compare_features(args, row, steps, stored, tool, fields, logs) -> None:
    result_json = steps.parent / "compare.json"
    command = [sys.executable, str(bc.HERE / "compare.py"), str(args.truth), str(stored),
               "--tool", tool, "--genome", args.accession, "--fields", fields,
               "--json", str(result_json),
               "--cases-tsv", str(args.results / "compare_cases.tsv")]
    for log in logs:
        command += ["--log", str(log)]
    bc.measured("compare", command, steps)
    if result_json.exists():
        result = json.loads(result_json.read_text())
        fill_from_compare(row, result)
        keep = args.results / "compare" / f"{args.accession}.{tool}.json"
        keep.parent.mkdir(parents=True, exist_ok=True)
        keep.write_text(result_json.read_text())


def run_gffutils(args, row, work, steps, strategy: str) -> None:
    gff = args.inputs / args.accession / "genomic.gff"
    db = work / "gffutils.db"
    load = bc.measured("load", [sys.executable, str(bc.HERE / "run_gffutils.py"), "load",
                                str(gff), str(db), "--strategy", strategy],
                       steps, timeout=args.timeout)
    row.update(outcome=outcome_of(load), exit_code=load["exit_code"], wall_s=load["wall_s"],
               cpu_s=load.get("cpu_s"), peak_mb=load["peak_rss_mb"])
    try:
        import gffutils

        row["tool_version"] = f"gffutils {gffutils.__version__} merge_strategy={strategy}"
    except ImportError:
        row["tool_version"] = f"gffutils (not importable) merge_strategy={strategy}"
    if load["exit_code"] != 0:
        row["error"] = error_of(steps, "load")
        return
    if strategy == "error":
        return
    stored = work / "stored.sqlite"
    bc.measured("extract", [sys.executable, str(bc.HERE / "run_gffutils.py"), "extract",
                            str(db), str(stored)], steps)
    compare_features(args, row, steps, stored, "gffutils", "seqid,type,strand,segs,parents,attrs",
                     [steps / "load.out", steps / "load.err"])
    if args.accession not in args.keep_db:
        db.unlink(missing_ok=True)
    stored.unlink(missing_ok=True)


def run_txdbmaker(args, row, work, steps) -> None:
    gff = args.inputs / args.accession / "genomic.gff"
    out = work / "txdb"
    rscript = os.environ.get("BENCH_RSCRIPT", "Rscript")
    load = bc.measured("load", [rscript, str(bc.HERE / "run_txdb.R"), str(gff), str(out),
                                ",".join(circular_seqids(args.truth))],
                       steps, timeout=args.timeout)
    row.update(outcome=outcome_of(load), exit_code=load["exit_code"], wall_s=load["wall_s"],
               cpu_s=load.get("cpu_s"), peak_mb=load["peak_rss_mb"])
    versions = out / "versions.txt"
    if versions.exists():
        lines = versions.read_text().splitlines()
        row["tool_version"] = "; ".join(line for line in lines[1:] if not line.endswith(
            "not installed") and not line.startswith("rtracklayer"))
    if load["exit_code"] != 0:
        row["error"] = error_of(steps, "load")
        return
    result_json = work / "compare.json"
    bc.measured("compare", [sys.executable, str(bc.HERE / "compare_txdb.py"), str(args.truth),
                            str(out), "--gff", str(gff), "--log", str(steps / "load.out"),
                            "--log", str(steps / "load.err"), "--genome", args.accession,
                            "--json", str(result_json),
                            "--cases-tsv", str(args.results / "compare_cases.tsv")], steps)
    if result_json.exists():
        result = json.loads(result_json.read_text())
        totals = result["totals"]
        row.update({
            "unit": "transcripts", "truth": totals.get("truth", 0),
            "kept": totals.get("kept", 0), "changed": totals.get("changed", 0),
            "changed_declared": 0, "lost": totals.get("lost", 0),
            "extra": result.get("extra_transcripts", 0),
            "not_modelled": sum(result.get("not_modelled_by_type", {}).values()),
            "parts_unattached": sum(result.get("parts_under_no_transcript", {}).values()),
            "fields_changed": ",".join(f"{k}:{v}" for k, v in
                                       sorted(result["field_differences"].items())),
            "genes_truth": result["genes"].get("truth", 0),
            "genes_held": result["genes"].get("held", 0),
            "lost_silent": totals.get("lost_silent", 0),
            "changed_silent": totals.get("changed_silent", 0),
        })
        keep = args.results / "compare" / f"{args.accession}.txdbmaker.json"
        keep.parent.mkdir(parents=True, exist_ok=True)
        keep.write_text(result_json.read_text())
    if args.accession not in args.keep_db:
        shutil.rmtree(out, ignore_errors=True)


def run_chado(args, row, work, steps) -> None:
    import psycopg

    python = os.environ.get("BENCH_CHADO_PYTHON")
    conninfo = os.environ.get("BENCH_PG_CONNINFO")
    schema = os.environ.get("BENCH_CHADO_SCHEMA")
    ontology = os.environ.get("BENCH_SO_OBO")
    psql = os.environ.get("BENCH_PSQL", "psql")
    missing = [name for name, value in (("BENCH_CHADO_PYTHON", python),
                                        ("BENCH_PG_CONNINFO", conninfo),
                                        ("BENCH_CHADO_SCHEMA", schema),
                                        ("BENCH_SO_OBO", ontology)) if not value]
    if missing:
        row.update(outcome="not_run", error="not set: " + ", ".join(missing))
        return
    database = "chado_" + bc.safe_name(args.accession)
    with psycopg.connect(conninfo + " dbname=postgres", autocommit=True) as connection:
        connection.execute(f'DROP DATABASE IF EXISTS "{database}"')
        connection.execute(f'CREATE DATABASE "{database}"')
    dsn = f"{conninfo} dbname={database}"
    bc.measured("setup", [python, str(bc.HERE / "run_chado.py"), "setup", "--dsn", dsn,
                          "--schema", schema, "--so", ontology, "--psql", psql], steps)
    genome = bc.genome(args.accession)
    load = bc.measured("load", [python, "-W", "ignore", str(bc.HERE / "run_chado.py"), "load",
                                "--dsn", dsn, str(args.inputs / args.accession / "genomic.gff"),
                                "--genus", genome["genus"], "--species", genome["species"]],
                       steps, timeout=args.timeout)
    row.update(outcome=outcome_of(load), exit_code=load["exit_code"], wall_s=load["wall_s"],
               cpu_s=load.get("cpu_s"), peak_mb=load["peak_rss_mb"],
               tool_version="python-chado 2.3.9, Chado 1.31")
    if load["exit_code"] != 0:
        row["error"] = error_of(steps, "load")
    stored = work / "stored.sqlite"
    # Read back even after a failure: whatever a failed load left behind is a result.
    bc.measured("extract", [python, str(bc.HERE / "run_chado.py"), "extract", "--dsn", dsn,
                            str(stored)], steps)
    if stored.exists():
        compare_features(args, row, steps, stored, "chado", "seqid,type,strand,segs,parents",
                         [steps / "load.out", steps / "load.err"])
        stored.unlink(missing_ok=True)
    if args.accession not in args.keep_db:
        with psycopg.connect(conninfo + " dbname=postgres", autocommit=True) as connection:
            connection.execute(f'DROP DATABASE IF EXISTS "{database}"')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--tool", required=True,
                        choices=["gffutils", "gffutils-default", "txdbmaker", "chado"])
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=12 * 3600)
    parser.add_argument("--keep-db", default="",
                        help="comma-separated accessions whose databases are kept for timing")
    args = parser.parse_args(argv)
    args.keep_db = {name for name in args.keep_db.split(",") if name}
    args.truth = args.results / "truth" / f"{args.accession}.sqlite"
    table = args.results / "compare_loaders.tsv"
    if bc.done(table, accession=args.accession, tool=args.tool):
        print(f"{args.accession} {args.tool}: already in {table}, skipped")
        return 0
    work = args.work / "compare" / args.tool / args.accession
    shutil.rmtree(work, ignore_errors=True)
    steps = work / "steps"
    steps.mkdir(parents=True)
    row = dict.fromkeys(COLUMNS)
    row.update(accession=args.accession, tool=args.tool)
    if args.tool == "gffutils":
        run_gffutils(args, row, work, steps, "create_unique")
    elif args.tool == "gffutils-default":
        run_gffutils(args, row, work, steps, "error")
    elif args.tool == "txdbmaker":
        run_txdbmaker(args, row, work, steps)
    else:
        run_chado(args, row, work, steps)
    bc.append_row(table, row)
    print(json.dumps(row))
    return 0


if __name__ == "__main__":
    sys.exit(main())
