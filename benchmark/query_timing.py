#!/usr/bin/env python3
"""Everyday questions, timed: how long a user waits for an answer.

    python3 benchmark/query_timing.py ACCESSION --work DIR --results DIR \\
        [--targets 20] [--repeats 5]

The tasks (PROTOCOL.md, section 6, step 8):

  T1 gene model   one gene and every feature below it, by the gene's ID
  T2 region       every feature overlapping a 1 Mb window
  T3 counts       how many features of each type
  T4 sequence     every feature on one sequence, written out as GFF3

The tools: cgload on SQLite and on PostgreSQL (the databases the scale step
built), gffutils (the database the comparison step built), and tabix on the
bgzipped annotation cgload writes beside its database. tabix has no index of
IDs, so it does T2 and T4 only; T3 for tabix is a full scan of the file, the
only way it can answer.

Two views of the same tasks:

* command: each answer is one new process, as a user typing a command gets it,
  interpreter start-up included. For gffutils the command is gffutils_query.py,
  a few lines of its documented API that print what was found as GFF3 lines:
  gffutils' own ``gffutils-cli`` cannot do these tasks (in 0.14 its ``region``
  raises NotImplementedError and its ``children`` lists grandchildren twice).
* library: one process opens the database once and asks T1 to T3 in turn, as a
  script would; only the call is timed. (T4 writes a file, a command's job.)

Targets are drawn with a seeded generator from the reference count: genes that
have children, 1 Mb windows on sequences at least that long, and the three
sequences for T4 (the one with most features and two more with at least 1,000).
Each target is asked --repeats times in a row. The first time is reported
separately; when the operating system's page cache can be dropped (root, or
sudo without a password), it is dropped before each first time, and the table
says whether it was. One row per run goes to RESULTS/query_timing.tsv.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import random
import shutil
import sqlite3
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

COLUMNS = ("accession", "view", "tool", "task", "target", "run", "first", "caches_dropped",
           "wall_s", "exit_code", "answer")
WINDOW = 1_000_000


def drop_caches() -> bool:
    commands = [["sh", "-c", "sync; echo 3 > /proc/sys/vm/drop_caches"],
                ["sudo", "-n", "sh", "-c", "sync; echo 3 > /proc/sys/vm/drop_caches"]]
    for command in commands:
        try:
            if subprocess.run(command, capture_output=True, check=False).returncode == 0:
                return True
        except OSError:
            continue
    return False


def pick_targets(truth: Path, accession: str, count: int) -> dict:
    rng = random.Random(int(hashlib.sha256(accession.encode()).hexdigest()[:12], 16))
    with sqlite3.connect(truth) as connection:
        parents = {row[0] for row in connection.execute(
            "SELECT DISTINCT parents FROM feature WHERE parents != '' AND parents NOT LIKE '%,%'")}
        genes = sorted(key for (key,) in connection.execute(
            "SELECT key FROM feature WHERE type = 'gene'") if key in parents)
        lengths = dict(connection.execute(
            "SELECT seqid, length FROM sequence WHERE length IS NOT NULL").fetchall())
        per_sequence = dict(connection.execute(
            "SELECT seqid, count(*) FROM feature GROUP BY seqid").fetchall())
    long_enough = sorted(s for s, n in lengths.items() if n >= WINDOW and per_sequence.get(s))
    windows = []
    if long_enough:
        weights = [lengths[s] for s in long_enough]
        for _ in range(count):
            seqid = rng.choices(long_enough, weights=weights)[0]
            start = rng.randint(1, lengths[seqid] - WINDOW + 1)
            windows.append((seqid, start, start + WINDOW - 1))
    busiest = max(per_sequence, key=per_sequence.get)
    others = sorted(s for s, n in per_sequence.items() if n >= 1000 and s != busiest)
    sequences = [busiest, *rng.sample(others, min(2, len(others)))]
    return {"genes": rng.sample(genes, min(count, len(genes))), "windows": windows,
            "sequences": sequences}


class Tools:
    def __init__(self, args):
        self.args = args
        base = args.work / "cgload" / "sqlite" / args.accession
        self.sqlite_db = base / "cgload.db"
        self.gffutils_db = args.work / "compare" / "gffutils" / args.accession / "gffutils.db"
        self.tabix_file = next((base / "data").glob("**/annotation.gff3.gz"), None)
        self.tabix = os.environ.get("BENCH_TABIX", "tabix")
        conninfo = os.environ.get("BENCH_PG_CONNINFO", "")
        self.pg_url = None
        if conninfo:
            parts = dict(item.split("=", 1) for item in conninfo.split())
            user = parts.pop("user", "")
            query = "&".join(f"{k}={v}" for k, v in parts.items())
            self.pg_url = (f"postgresql+psycopg://{user}@/cgb_{bc.safe_name(args.accession)}"
                           f"?{query}")

        #: The one assembly in each cgload database, or None: a refused or
        #: failed load leaves none, and its tool is then not timed.
        self.assembly = {
            "cgload-sqlite": self._assembly(f"sqlite:///{self.sqlite_db}")
            if self.sqlite_db.exists() else None,
            "cgload-postgres": self._assembly(self.pg_url) if self.pg_url else None,
        }

    def available(self) -> dict[str, bool]:
        return {
            "cgload-sqlite": self.assembly["cgload-sqlite"] is not None,
            "cgload-postgres": self.assembly["cgload-postgres"] is not None,
            "gffutils": self.gffutils_db.exists(),
            "tabix": self.tabix_file is not None and shutil.which(self.tabix) is not None,
        }

    @staticmethod
    def _assembly(url: str) -> int | None:
        try:
            from sqlalchemy import create_engine, text

            with create_engine(url).connect() as connection:
                ids = connection.execute(text("SELECT assembly_id FROM assembly")).scalars().all()
        except Exception:  # noqa: BLE001 - no database means no timing, which is recorded
            return None
        return ids[0] if len(ids) == 1 else None

    def assembly_arg(self, tool: str) -> str:
        return str(self.assembly[tool])

    def command(self, tool: str, task: str, target, out: Path) -> list[str] | None:
        if tool.startswith("cgload"):
            if task == "gene":
                return ["cgload", "query", "--assembly-id", self.assembly_arg(tool), "--id",
                        target, "--descendants"]
            if task == "region":
                return ["cgload", "query", "--assembly-id", self.assembly_arg(tool),
                        "--region", f"{target[0]}:{target[1]}-{target[2]}"]
            if task == "count":
                return ["cgload", "stats", "--by", "source"]
            return ["cgload", "query", "--assembly-id", self.assembly_arg(tool), "--region",
                    target,
                    "--format", "gff3", "--output", str(out)]
        if tool == "gffutils":
            script = [sys.executable, str(bc.HERE / "gffutils_query.py"), str(self.gffutils_db)]
            if task == "gene":
                return [*script, "gene", target]
            if task == "region":
                return [*script, "region", target[0], str(target[1]), str(target[2])]
            if task == "count":
                return [*script, "count"]
            return [*script, "export", target, str(out)]
        if task == "region":
            return [self.tabix, str(self.tabix_file), f"{target[0]}:{target[1]}-{target[2]}"]
        if task == "sequence":
            return ["sh", "-c", f'"{self.tabix}" "{self.tabix_file}" "{target}" > "{out}"']
        if task == "count":
            return [sys.executable, str(Path(__file__).resolve()), "_scan", str(self.tabix_file)]
        return None


def label_of(target) -> str:
    return target if isinstance(target, str) else f"{target[0]}:{target[1]}-{target[2]}"


def answer_of(tool: str, task: str, stdout: str, out: Path) -> str:
    """What came back, to check that the tools answered the same question: for
    T3 the total count, otherwise how many features (cgload) or records
    (gffutils, tabix) were returned."""
    if task == "sequence":
        if not out.exists():
            return "no output"
        with out.open() as handle:
            return str(sum(1 for line in handle if line.strip() and not line.startswith("#")))
    lines = [line for line in stdout.splitlines() if line.strip()]
    if task == "count":
        if tool.startswith("cgload"):
            return str(_stats_total(stdout))
        return lines[-1] if lines else ""
    data = [line for line in lines if not line.startswith("#")]
    if tool.startswith("cgload"):
        return str(max(len(data) - 1, 0))  # a header row, then one row per feature
    return str(len(data))


def _stats_total(text: str) -> int | str:
    """The total of ``cgload stats --by source``: its 'by source type:' line
    lists each type with its count."""
    for line in text.splitlines():
        if "by source type:" in line:
            items = line.split("by source type:", 1)[1].split(", ")
            return sum(int(item.rsplit(" ", 1)[-1].replace(",", "")) for item in items
                       if item.strip())
    return "no count"


def scan(path: str) -> int:
    """T3 for tabix: count types by reading the whole bgzipped file."""
    counts: dict[str, int] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("#"):
                kind = line.split("\t", 3)[2]
                counts[kind] = counts.get(kind, 0) + 1
    print(sum(counts.values()))
    return 0


def time_commands(args, tools: Tools, targets: dict, table: Path) -> None:
    available = tools.available()
    tasks = [("gene", targets["genes"]), ("region", targets["windows"]), ("count", ["all"]),
             ("sequence", targets["sequences"])]
    out = args.work / "timing" / args.accession / "out.gff3"
    out.parent.mkdir(parents=True, exist_ok=True)
    for tool, present in available.items():
        env = {**os.environ}
        if tool == "cgload-sqlite":
            env["CGLOAD_DB_URL"] = f"sqlite:///{tools.sqlite_db}"
        elif tool == "cgload-postgres":
            env["CGLOAD_DB_URL"] = tools.pg_url or ""
        for task, items in tasks:
            if bc.done(table, accession=args.accession, view="command", tool=tool, task=task):
                continue
            rows = []
            for target in items:
                command = tools.command(tool, task, target, out) if present else None
                label = label_of(target)
                if command is None:
                    rows.append({"accession": args.accession, "view": "command", "tool": tool,
                                 "task": task, "target": label, "run": 0, "first": None,
                                 "caches_dropped": None, "wall_s": None, "exit_code": None,
                                 "answer": "not available" if not present else "not possible"})
                    continue
                for run in range(1, args.repeats + 1):
                    dropped = drop_caches() if run == 1 else False
                    started = time.perf_counter()
                    done = subprocess.run(command, env=env, capture_output=True, text=True,
                                          check=False)
                    wall = time.perf_counter() - started
                    answer = answer_of(tool, task, done.stdout, out)
                    rows.append({"accession": args.accession, "view": "command", "tool": tool,
                                 "task": task, "target": label, "run": run, "first": run == 1,
                                 "caches_dropped": dropped, "wall_s": round(wall, 4),
                                 "exit_code": done.returncode, "answer": answer[:80]})
            for row in rows:
                bc.append_row(table, {name: row.get(name) for name in COLUMNS})
            print(f"{args.accession} command {tool} {task}: {len(rows)} runs")
    out.unlink(missing_ok=True)


def time_library(args, tools: Tools, targets: dict, table: Path) -> None:
    """cgload and gffutils called in one process each, database opened once."""
    for tool in ("cgload-sqlite", "cgload-postgres", "gffutils"):
        if not tools.available()[tool]:
            continue
        if bc.done(table, accession=args.accession, view="library", tool=tool):
            continue
        rows = []
        if tool == "gffutils":
            import gffutils

            db = gffutils.FeatureDB(str(tools.gffutils_db))
            calls = {
                "gene": lambda g, db=db: [db[g], *db.children(db[g])],
                "region": lambda w, db=db: list(db.region(seqid=w[0], start=w[1], end=w[2])),
                "count": lambda _t, db=db: {t: db.count_features_of_type(t)
                                            for t in db.featuretypes()},
            }
            connection = None
        else:
            from cgload import query as cg
            from sqlalchemy import create_engine

            url = f"sqlite:///{tools.sqlite_db}" if tool == "cgload-sqlite" else tools.pg_url
            connection = create_engine(url).connect()
            assembly = tools.assembly[tool]
            calls = {
                "gene": lambda g, c=connection, a=assembly: cg.by_id(
                    c, a, g, descendants=True).features,
                "region": lambda w, c=connection, a=assembly: cg.in_region(
                    c, a, f"{w[0]}:{w[1]}-{w[2]}").features,
                "count": lambda _t, c=connection: _cgload_counts(c),
            }
        for task, items in (("gene", targets["genes"]), ("region", targets["windows"]),
                            ("count", ["all"])):
            for target in items:
                label = label_of(target)
                for run in range(1, args.repeats + 1):
                    started = time.perf_counter()
                    answer = calls[task](target)
                    wall = time.perf_counter() - started
                    rows.append({"accession": args.accession, "view": "library", "tool": tool,
                                 "task": task, "target": label, "run": run, "first": run == 1,
                                 "caches_dropped": False, "wall_s": round(wall, 5),
                                 "exit_code": 0, "answer": len(answer)})
        if connection is not None:
            connection.close()
        for row in rows:
            bc.append_row(table, {name: row.get(name) for name in COLUMNS})
        print(f"{args.accession} library {tool}: {len(rows)} runs")


def _cgload_counts(connection) -> dict:
    from cgload.db.schema import feature
    from sqlalchemy import func, select

    return dict(connection.execute(
        select(feature.c.source_type, func.count()).group_by(feature.c.source_type)).all())


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["_scan"]:
        return scan(argv[1])
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--targets", type=int, default=20)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args(argv)
    os.environ.pop("CGLOAD_DB_URL", None)
    truth = args.results / "truth" / f"{args.accession}.sqlite"
    targets = pick_targets(truth, args.accession, args.targets)
    tools = Tools(args)
    table = args.results / "query_timing.tsv"
    time_commands(args, tools, targets, table)
    time_library(args, tools, targets, table)
    return 0


if __name__ == "__main__":
    sys.exit(main())
