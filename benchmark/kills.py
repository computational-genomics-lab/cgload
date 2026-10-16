#!/usr/bin/env python3
"""Interrupted loads: what does a loader leave behind when it is killed?

    python3 benchmark/kills.py ACCESSION --tool cgload|cgload-force|gffutils|chado \\
        --inputs DIR --work DIR --results DIR [--trials 10]

Each trial starts a load and kills it, with SIGKILL to its whole process group,
at a time drawn uniformly between 5% and 95% of how long a full load of the
same genome took, alone: from the scale step for cgload, from the control load
(below) for gffutils and Chado. A seeded random generator picks the moments.

* ``cgload``: a load into a new database. Afterwards the database must hold no
  assembly and no feature, and the same load, run again, must succeed.
* ``cgload-force``: the genome is loaded once, then ``load --force`` replaces it
  and is killed. Afterwards the first load must still be there and still pass
  ``verify``.
* ``gffutils``: ``create_db`` with create_unique, killed. Is a database file left,
  does it open, and how many records and parent links does it hold?
* ``chado``: python-chado's ``load_gff`` into a new Chado database, killed. How
  many features and feature relationships does the database hold?

What a killed load left is counted in the loader's own units (rows of its
feature table, rows of its link table) and set against the same counts taken
the same way from a load that was not killed: the earlier load for
``cgload-force``, the repeated load for ``cgload``, and for gffutils and Chado
one uninterrupted control load per genome (RESULTS/kills_complete.tsv). The
state is then

* ``empty``: nothing of the genome is stored;
* ``complete``: everything a full load stores is there (the kill came after
  the data was written);
* ``partial``: some of it is there;
* ``intact`` / ``damaged`` (``cgload-force``): the earlier load still passes
  ``verify``, or not;
* ``no_file`` (gffutils): no database file was left.

One row per trial goes to RESULTS/kills.tsv.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import random
import shutil
import signal
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

COLUMNS = ("accession", "tool", "trial", "full_load_s", "kill_at_s", "finished_first",
           "state", "opens", "assemblies", "records_left", "records_complete", "links_left",
           "links_complete", "verify_after", "reload_ok", "leftover_bytes",
           "orphan_directory", "leftover_after_reload", "note")
CONTROL_COLUMNS = ("accession", "tool", "records", "links", "wall_s", "note")
#: What each loader's counts are counted in.
CGLOAD_COUNTS = ("SELECT count(*) FROM feature",
                 "SELECT count(*) FROM feature WHERE parent_id IS NOT NULL")
GFFUTILS_COUNTS = ("SELECT count(*) FROM features",
                   "SELECT count(*) FROM relations WHERE level = 1")
CHADO_COUNTS = ("SELECT count(*) FROM feature", "SELECT count(*) FROM feature_relationship")


def full_duration(args) -> float | None:
    """How long a full cgload load of this genome took, from the scale step,
    which runs each load alone."""
    table = args.results / "cgload_scale.tsv"
    if not table.exists():
        return None
    with table.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["accession"] == args.accession and row.get("backend") == "sqlite":
                try:
                    return float(row["load_wall_s"])
                except (TypeError, ValueError):
                    return None
    return None


def start_and_kill(command: list[str], delay: float, env: dict, log: Path) -> bool:
    """Start, wait ``delay`` seconds, kill the process group. True if it ended first."""
    with log.open("w") as handle:
        process = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT,
                                   start_new_session=True, env={**os.environ, **env})
        deadline = time.monotonic() + delay
        while time.monotonic() < deadline:
            if process.poll() is not None:
                return True
            time.sleep(0.05)
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    return False


def sqlite_count(db: Path, query: str):
    """Opened read-write, as the next user of the file would open it: SQLite
    rolls back a transaction that a killed process left half-written (its hot
    journal) only for a connection that may write."""
    try:
        with sqlite3.connect(db) as connection:
            return connection.execute(query).fetchone()[0]
    except sqlite3.Error as error:
        return f"unreadable: {error}"


def state_of(records, links, complete: tuple | None) -> str:
    if not isinstance(records, int) or not isinstance(links, int):
        return "unreadable"
    if records == 0:
        return "empty"
    if complete is not None and (records, links) == tuple(complete):
        return "complete"
    return "partial"


def leftover(data: Path) -> int:
    """Bytes in the data directory that belong to no stored assembly: staging
    directories and directories set aside by --force."""
    if not data.exists():
        return 0
    loose = [p for p in data.rglob("*") if p.is_file()
             and (".staging" in p.parts or any("superseded" in part for part in p.parts))]
    return sum(p.stat().st_size for p in loose)


def config_for(args) -> tuple[Path, bool]:
    config = args.work / "cgload" / "sqlite" / args.accession / "organism.toml"
    repair = bc.done(args.results / "cgload_scale.tsv", accession=args.accession,
                     backend="sqlite", outcome="loaded_with_repair")
    return config, repair


def trial_cgload(args, number: int, delay: float, force: bool) -> dict:
    work = args.work / "kills" / args.tool / args.accession / str(number)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    db, data = work / "cgload.db", work / "data"
    env = {"CGLOAD_DB_URL": f"sqlite:///{db}"}
    config, repair = config_for(args)
    load = ["cgload", "load", "--config", str(config), "--data-dir", str(data)]
    if repair:
        load.append("--repair-envelopes")

    def run(command: list[str]) -> subprocess.CompletedProcess:
        return subprocess.run(command, env={**os.environ, **env}, capture_output=True,
                              text=True, check=False)

    def counts() -> tuple:
        return tuple(sqlite_count(db, query) for query in CGLOAD_COUNTS)

    def verified() -> bool:
        return run(["cgload", "verify", "--assembly-id", "1"]).returncode == 0

    run(["cgload", "init"])
    row = {}
    if force:
        first = run(load)
        if first.returncode != 0:
            return {"state": "setup_failed", "note": bc.first_lines(first.stderr, 2)}
        row["records_complete"], row["links_complete"] = counts()
        load = [*load, "--force"]
    row["finished_first"] = start_and_kill(load, delay, env, work / "load.log")
    row["assemblies"] = sqlite_count(db, "SELECT count(*) FROM assembly")
    row["records_left"], row["links_left"] = counts()
    row["leftover_bytes"] = leftover(data)
    if force:
        ok = verified()
        row["verify_after"] = "pass" if ok else "fail"
        same = (row["records_left"], row["links_left"]) == (row["records_complete"],
                                                           row["links_complete"])
        row["state"] = "intact" if ok and row["assemblies"] == 1 and same else "damaged"
    elif row["assemblies"] == 1:
        # Killed after the load committed: the genome is stored; is it whole?
        ok = verified()
        row["verify_after"] = "pass" if ok else "fail"
        row["records_complete"], row["links_complete"] = counts()
        row["state"] = "complete" if ok else "damaged"
    else:
        row["state"] = "empty" if row["assemblies"] == 0 and row["records_left"] == 0 \
            else "partial"
        # A directory in place with no rows: the kill came between the rename
        # and the commit (D-095). The repeated load must replace it.
        row["orphan_directory"] = any(
            path.is_dir() for path in data.glob("*/*/*") if ".staging" not in path.parts)
        again = run(load)
        row["reload_ok"] = again.returncode == 0
        row["leftover_after_reload"] = leftover(data)
        if again.returncode == 0:
            row["records_complete"], row["links_complete"] = counts()
        else:
            row["note"] = bc.first_lines(again.stderr, 2)
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return row


def gffutils_load(args, db: Path) -> list[str]:
    return [sys.executable, str(bc.HERE / "run_gffutils.py"), "load",
            str(args.inputs / args.accession / "genomic.gff"), str(db),
            "--strategy", "create_unique"]


def gffutils_opens(db: Path) -> tuple[bool, str | None]:
    """Whether gffutils itself can open the file and list what it holds."""
    try:
        import gffutils

        database = gffutils.FeatureDB(str(db))
        sum(1 for _ in database.featuretypes())
    except Exception as error:  # noqa: BLE001 - what failed is the result
        return False, f"{type(error).__name__}: {error}"[:300]
    return True, None


def trial_gffutils(args, number: int, delay: float, complete: tuple) -> dict:
    work = args.work / "kills" / args.tool / args.accession / str(number)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    db = work / "gffutils.db"
    row = {"finished_first": start_and_kill(gffutils_load(args, db), delay, {},
                                            work / "load.log"),
           "records_complete": complete[0], "links_complete": complete[1]}
    if not db.exists():
        row.update(state="no_file", opens=False)
    else:
        row["records_left"], row["links_left"] = (sqlite_count(db, q) for q in GFFUTILS_COUNTS)
        row["opens"], row["note"] = gffutils_opens(db)
        row["state"] = state_of(row["records_left"], row["links_left"], complete)
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return row


def chado_database(args, name: str):
    """A new Chado database with the schema and the Sequence Ontology; yields its DSN
    and drops it afterwards."""
    import psycopg

    conninfo = os.environ["BENCH_PG_CONNINFO"]
    database = f"chado_kill_{bc.safe_name(args.accession)}_{name}"
    with psycopg.connect(conninfo + " dbname=postgres", autocommit=True) as connection:
        connection.execute(f'DROP DATABASE IF EXISTS "{database}"')
        connection.execute(f'CREATE DATABASE "{database}"')
    dsn = f"{conninfo} dbname={database}"
    subprocess.run([os.environ["BENCH_CHADO_PYTHON"], str(bc.HERE / "run_chado.py"), "setup",
                    "--dsn", dsn, "--schema", os.environ["BENCH_CHADO_SCHEMA"],
                    "--so", os.environ["BENCH_SO_OBO"],
                    "--psql", os.environ.get("BENCH_PSQL", "psql")],
                   capture_output=True, check=False)
    return database, dsn


def chado_drop(database: str) -> None:
    import psycopg

    with psycopg.connect(os.environ["BENCH_PG_CONNINFO"] + " dbname=postgres",
                         autocommit=True) as connection:
        connection.execute(f'DROP DATABASE IF EXISTS "{database}"')


def chado_load(args, dsn: str) -> list[str]:
    genome = bc.genome(args.accession)
    return [os.environ["BENCH_CHADO_PYTHON"], "-W", "ignore", str(bc.HERE / "run_chado.py"),
            "load", "--dsn", dsn, str(args.inputs / args.accession / "genomic.gff"),
            "--genus", genome["genus"], "--species", genome["species"]]


def chado_counts(dsn: str) -> tuple:
    import psycopg

    with psycopg.connect(dsn) as connection:
        return tuple(connection.execute(query).fetchone()[0] for query in CHADO_COUNTS)


def trial_chado(args, number: int, delay: float, complete: tuple) -> dict:
    database, dsn = chado_database(args, str(number))
    work = args.work / "kills" / args.tool / args.accession / str(number)
    work.mkdir(parents=True, exist_ok=True)
    row = {"finished_first": start_and_kill(chado_load(args, dsn), delay, {},
                                            work / "load.log"),
           "records_complete": complete[0], "links_complete": complete[1]}
    row["records_left"], row["links_left"] = chado_counts(dsn)
    row["state"] = state_of(row["records_left"], row["links_left"], complete)
    chado_drop(database)
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return row


def control(args) -> tuple[tuple[int, int], float]:
    """What an uninterrupted load of this genome stores, counted as the trials
    count it (gffutils and Chado), and how long the load took. Run once per
    genome, alone, and kept in RESULTS/kills_complete.tsv."""
    table = args.results / "kills_complete.tsv"
    if table.exists():
        with table.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["accession"] == args.accession and row["tool"] == args.tool:
                    return (int(row["records"]), int(row["links"])), float(row["wall_s"])
    work = args.work / "kills" / args.tool / args.accession / "complete"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    if args.tool == "gffutils":
        db = work / "gffutils.db"
        started = time.monotonic()
        done_ = subprocess.run(gffutils_load(args, db), capture_output=True, text=True,
                               check=False)
        wall = time.monotonic() - started
        counts = tuple(sqlite_count(db, query) for query in GFFUTILS_COUNTS)
    else:
        database, dsn = chado_database(args, "complete")
        started = time.monotonic()
        done_ = subprocess.run(chado_load(args, dsn), capture_output=True, text=True,
                               check=False)
        wall = time.monotonic() - started
        counts = chado_counts(dsn)
        chado_drop(database)
    if done_.returncode != 0 or not all(isinstance(n, int) for n in counts):
        raise SystemExit(f"the control load of {args.accession} with {args.tool} failed: "
                         f"{bc.first_lines(done_.stderr, 3)}")
    bc.append_row(table, dict(zip(CONTROL_COLUMNS, (args.accession, args.tool, *counts,
                                                    round(wall, 1), None), strict=True)))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    return counts, wall


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--tool", required=True,
                        choices=["cgload", "cgload-force", "gffutils", "chado"])
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args(argv)
    os.environ.pop("CGLOAD_DB_URL", None)
    complete = None
    if args.tool in ("gffutils", "chado"):
        complete, full = control(args)
    else:
        full = full_duration(args)
    if full is None:
        raise SystemExit(f"no full-load time for {args.accession} {args.tool}; run the scale "
                         "step first")
    seed = int(hashlib.sha256(f"{args.accession}:{args.tool}".encode()).hexdigest()[:12], 16)
    rng = random.Random(seed)
    delays = [rng.uniform(0.05, 0.95) * full for _ in range(args.trials)]
    table = args.results / "kills.tsv"
    for number, delay in enumerate(delays, start=1):
        if bc.done(table, accession=args.accession, tool=args.tool, trial=number):
            continue
        if args.tool == "gffutils":
            found = trial_gffutils(args, number, delay, complete)
        elif args.tool == "chado":
            found = trial_chado(args, number, delay, complete)
        else:
            found = trial_cgload(args, number, delay, force=args.tool == "cgload-force")
        row = dict.fromkeys(COLUMNS)
        row.update(accession=args.accession, tool=args.tool, trial=number, full_load_s=full,
                   kill_at_s=round(delay, 2))
        row.update({key: value for key, value in found.items() if key in COLUMNS})
        bc.append_row(table, row)
        print(f"{args.accession} {args.tool} trial {number} at {delay:.1f}s: {row['state']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
