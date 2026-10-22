#!/usr/bin/env python3
"""Load one genome into cgload on one backend, and record what it cost and what it kept.

    python3 benchmark/run_cgload.py ACCESSION --backend sqlite|postgres --inputs DIR \\
        --work DIR --results DIR [--pg-conninfo "host=... port=... user=..."]

INPUTS/ACCESSION/ holds genomic.gff, assembly.fa and, for the fungi,
genomic.gbff; RESULTS holds the reference count (truth/ACCESSION.sqlite, from
truth.py). The steps, each measured with measure.py:

1. ``cgload init`` on a new database (an SQLite file, or a new PostgreSQL
   database named after the genome).
2. ``cgload load``. If cgload refuses the file, the reasons are kept, the
   database is checked to be still empty (a refused load must change nothing),
   and the load is repeated with ``--repair-envelopes``.
3. ``cgload verify``, with ``--against`` the GenBank file where there is one.
4. The size of the database and of the data directory.
5. ``cgload export``, then truth.py on what it wrote and compare.py against the
   reference count: the round trip. Envelopes cgload reports repairing (in
   manifest.json) are counted as declared changes.

Appends one row to RESULTS/cgload_scale.tsv and the round trip's hard-case rows
to RESULTS/cgload_fidelity_cases.tsv. Skips a genome and backend already in
the table, so a stopped run can be restarted.
"""

from __future__ import annotations

import argparse
import json
import os
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

TIMEOUT = 6 * 3600

#: The columns of cgload_scale.tsv, in order; a step that did not run leaves NA.
COLUMNS = (
    "accession", "backend", "cgload_commit", "init_exit", "plain_exit", "plain_wall_s",
    "refusal", "empty_after_refusal", "outcome", "refusal_after_repair", "load_wall_s",
    "load_cpu_s", "load_peak_mb", "warnings", "repairs", "envelope_exceptions", "verify_exit",
    "verify_wall_s", "verify_peak_mb", "tier_source_count", "tier_parent_edges",
    "tier_stored_assembly", "tier_cross_parser", "db_bytes", "data_bytes", "export_exit",
    "export_wall_s", "rt_truth", "rt_kept", "rt_changed", "rt_changed_declared", "rt_lost",
    "rt_extra", "rt_truth_links", "rt_kept_links", "rt_fields_changed",
)


def write_config(path: Path, genome: dict, inputs: Path) -> None:
    def quoted(text: str) -> str:
        return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'

    path.write_text(
        "[organism]\n"
        f"genus = {quoted(genome['genus'])}\n"
        f"species = {quoted(genome['species'])}\n"
        f"strain = {quoted(genome['strain'])}\n\n"
        "[assembly]\n"
        f"name = {quoted(genome['assembly_name'])}\n"
        f"accession = {quoted(genome['accession'])}\n"
        f"fasta = {quoted(str(inputs / 'assembly.fa'))}\n\n"
        "[[files.annotation]]\n"
        f"path = {quoted(str(inputs / 'genomic.gff'))}\n"
    )


def postgres(conninfo: str, database: str) -> str:
    """Make a new, empty database; return its SQLAlchemy URL."""
    import psycopg

    with psycopg.connect(conninfo + " dbname=postgres", autocommit=True) as connection:
        connection.execute(f'DROP DATABASE IF EXISTS "{database}"')
        connection.execute(f'CREATE DATABASE "{database}"')
    parts = dict(item.split("=", 1) for item in conninfo.split())
    user = parts.pop("user", "")
    query = "&".join(f"{key}={value}" for key, value in parts.items())
    return f"postgresql+psycopg://{user}@/{database}?{query}"


def database_bytes(backend: str, sqlite_path: Path, conninfo: str, database: str) -> int | None:
    if backend == "sqlite":
        return sqlite_path.stat().st_size if sqlite_path.exists() else None
    import psycopg

    with psycopg.connect(conninfo + f" dbname={database}") as connection:
        return connection.execute("SELECT pg_database_size(current_database())").fetchone()[0]


def assembly_ids(backend: str, sqlite_path: Path, conninfo: str, database: str) -> list[int]:
    query = "SELECT assembly_id FROM assembly ORDER BY assembly_id"
    if backend == "sqlite":
        with sqlite3.connect(sqlite_path) as connection:
            return [row[0] for row in connection.execute(query)]
    import psycopg

    with psycopg.connect(conninfo + f" dbname={database}") as connection:
        return [row[0] for row in connection.execute(query)]


def assemblies_in(backend: str, sqlite_path: Path, conninfo: str, database: str) -> int:
    return len(assembly_ids(backend, sqlite_path, conninfo, database))


def run(args) -> int:
    genome = bc.genome(args.accession)
    table = args.results / "cgload_scale.tsv"
    if bc.done(table, accession=args.accession, backend=args.backend):
        print(f"{args.accession} {args.backend}: already in {table}, skipped")
        return 0
    inputs = args.inputs / args.accession
    work = args.work / "cgload" / args.backend / args.accession
    if work.exists():
        import shutil

        shutil.rmtree(work)
    steps = work / "steps"
    steps.mkdir(parents=True)
    config = work / "organism.toml"
    write_config(config, genome, inputs)
    data = work / "data"
    sqlite_path = work / "cgload.db"
    database = "cgb_" + bc.safe_name(args.accession)
    if args.backend == "sqlite":
        url = f"sqlite:///{sqlite_path}"
    else:
        url = postgres(args.pg_conninfo, database)
    env = {"CGLOAD_DB_URL": url}
    row = dict.fromkeys(COLUMNS)
    row.update({"accession": args.accession, "backend": args.backend,
                "cgload_commit": bc.git_commit(bc.HERE.parent)})

    init = bc.measured("init", ["cgload", "init"], steps, env=env)
    row["init_exit"] = init["exit_code"]
    load_command = ["cgload", "load", "--config", str(config), "--data-dir", str(data)]
    plain = bc.measured("load", load_command, steps, timeout=TIMEOUT, env=env)
    row["plain_exit"] = plain["exit_code"]
    row["plain_wall_s"] = plain["wall_s"]
    row["refusal"] = ""
    final, final_name = plain, "load"
    if plain["exit_code"] != 0:
        row["refusal"] = bc.first_lines(bc.text_of(steps, "load"), 4)
        try:
            row["empty_after_refusal"] = assemblies_in(args.backend, sqlite_path,
                                                       args.pg_conninfo, database) == 0
        except Exception as error:  # noqa: BLE001 - recorded, not hidden
            row["empty_after_refusal"] = f"could not check: {error}"
        final = bc.measured("load_repair", [*load_command, "--repair-envelopes"], steps,
                            timeout=TIMEOUT, env=env)
        final_name = "load_repair"
    if final["exit_code"] == 0:
        row["outcome"] = "loaded" if final_name == "load" else "loaded_with_repair"
    elif final.get("timed_out"):
        row["outcome"] = "timed_out"
    else:
        row["outcome"] = "refused"
        row["refusal_after_repair"] = bc.first_lines(bc.text_of(steps, final_name), 4)
    row["load_wall_s"] = final["wall_s"]
    row["load_cpu_s"] = final.get("cpu_s")
    row["load_peak_mb"] = final["peak_rss_mb"]
    load_text = bc.text_of(steps, final_name)
    row["warnings"] = sum(1 for line in load_text.splitlines() if "warning:" in line)

    if final["exit_code"] == 0:
        manifest = next(data.glob("**/manifest.json"), None)
        repairs = []
        if manifest:
            content = json.loads(manifest.read_text())
            repairs = [item["feature"] for item in content.get("envelope_repairs", [])]
            row["envelope_exceptions"] = len(content.get("envelope_exceptions", []))
        row["repairs"] = len(repairs)
        (work / "declared.txt").write_text("\n".join(repairs) + "\n")

        # Not always 1: on PostgreSQL a refused load can use up an id, since a
        # sequence does not roll back.
        ids = assembly_ids(args.backend, sqlite_path, args.pg_conninfo, database)
        assembly = str(ids[0]) if len(ids) == 1 else "1"
        verify_command = ["cgload", "verify", "--assembly-id", assembly]
        gbff = inputs / "genomic.gbff"
        if gbff.exists():
            verify_command += ["--against", str(gbff)]
        verify = bc.measured("verify", verify_command, steps, timeout=TIMEOUT, env=env)
        row["verify_exit"] = verify["exit_code"]
        row["verify_wall_s"] = verify["wall_s"]
        row["verify_peak_mb"] = verify["peak_rss_mb"]
        report = bc.text_of(steps, "verify")
        # The whole report, not only each tier's verdict: a failing tier names
        # the types and counts that disagree.
        (args.results / "verify").mkdir(parents=True, exist_ok=True)
        (args.results / "verify" / f"{args.accession}.{args.backend}.txt").write_text(report)
        found = bc.tiers(report)
        for name in ("source count", "parent edges", "stored assembly", "cross-parser"):
            row["tier_" + bc.safe_name(name)] = found.get(name, "absent")

        row["db_bytes"] = database_bytes(args.backend, sqlite_path, args.pg_conninfo, database)
        row["data_bytes"] = bc.tree_bytes(data)

        exported = work / "export.gff3"
        export = bc.measured("export", ["cgload", "export", "--assembly-id", assembly,
                                        "--format", "gff3", "--output", str(exported)],
                             steps, timeout=TIMEOUT, env=env)
        row["export_exit"] = export["exit_code"]
        row["export_wall_s"] = export["wall_s"]
        if export["exit_code"] == 0:
            stored = work / "stored.sqlite"
            bc.measured("truth_export", [sys.executable, str(bc.HERE / "truth.py"), str(exported),
                                         str(stored), "--ignore-attribute",
                                         "cgload_linkage_method"], steps)
            result_json = work / "fidelity.json"
            bc.measured("compare", [
                sys.executable, str(bc.HERE / "compare.py"),
                str(args.results / "truth" / f"{args.accession}.sqlite"), str(stored),
                "--tool", f"cgload-{args.backend}", "--genome", args.accession,
                "--declared", str(work / "declared.txt"),
                "--log", str(steps / f"{final_name}.out"),
                "--log", str(steps / f"{final_name}.err"),
                "--json", str(result_json),
                "--cases-tsv", str(args.results / "cgload_fidelity_cases.tsv")], steps)
            if result_json.exists():
                result = json.loads(result_json.read_text())
                totals = result["totals"]
                for name in ("truth", "kept", "changed", "changed_declared", "lost",
                             "truth_links", "kept_links"):
                    row["rt_" + name] = totals.get(name, 0)
                row["rt_extra"] = sum(result["extra_by_type"].values())
                row["rt_fields_changed"] = ",".join(
                    f"{k}:{v}" for k, v in sorted(result["field_differences"].items()))
                (args.results / "fidelity").mkdir(parents=True, exist_ok=True)
                target = args.results / "fidelity" / f"{args.accession}.{args.backend}.json"
                target.write_text(result_json.read_text())
            if not args.keep_export:
                exported.unlink(missing_ok=True)
                stored.unlink(missing_ok=True)

    unexpected = set(row) - set(COLUMNS)
    if unexpected:
        raise SystemExit(f"columns missing from COLUMNS: {sorted(unexpected)}")
    bc.append_row(table, row)
    print(json.dumps(row))
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--backend", choices=["sqlite", "postgres"], required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--pg-conninfo", default=os.environ.get("BENCH_PG_CONNINFO", ""))
    parser.add_argument("--keep-export", action="store_true")
    args = parser.parse_args(argv)
    if args.backend == "postgres" and not args.pg_conninfo:
        parser.error("--pg-conninfo (or BENCH_PG_CONNINFO) is needed for postgres")
    return run(args)


if __name__ == "__main__":
    sys.exit(main())
