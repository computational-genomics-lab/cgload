#!/usr/bin/env python3
"""The case study: one of the lab's own genomes, with its functional annotation.

    python3 benchmark/case_study.py CONFIG.toml --work DIR --results DIR [--repeats 5]

CONFIG.toml is an ordinary cgload organism config naming the assembly FASTA,
the funannotate GFF3 under [[files.annotation]], and the eggNOG-mapper and
InterProScan outputs under [[files.functional]] (README, "Your own data").

Steps, each recorded:

1. ``cgload load`` into a new SQLite database, and ``cgload verify``.
2. What the functional files added: assignments per analysis, how many
   proteins carry at least one, and how many queries in each file matched no
   feature (cgload stores that rate; D-039).
3. Questions a biologist asks: the features carrying the most frequent Pfam,
   GO, KEGG KO, EC and InterPro accession, each with ``cgload query --term``,
   timed --repeats times (task T5 of the timing protocol).
4. The round trip: ``cgload export`` read by truth.py and compared with the
   input, as for the RefSeq genomes.
5. Where the annotation and FASTA are, for the JBrowse 2 figure
   (benchmark/jbrowse/render.sh makes it from them).

Writes RESULTS/case_study.json and RESULTS/case_study_queries.tsv.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
import tomllib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

ANALYSES = ("Pfam", "GO", "KEGG_ko", "EC", "InterPro")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args(argv)

    work = args.work / "case_study"
    shutil.rmtree(work, ignore_errors=True)
    steps = work / "steps"
    steps.mkdir(parents=True)
    db = work / "cgload.db"
    env = {"CGLOAD_DB_URL": f"sqlite:///{db}"}
    report: dict = {"config": str(args.config.resolve()),
                    "cgload_commit": bc.git_commit(bc.HERE.parent)}

    bc.measured("init", ["cgload", "init"], steps, env=env)
    load = bc.measured("load", ["cgload", "load", "--config", str(args.config),
                                "--data-dir", str(work / "data")], steps, env=env)
    report["load"] = {k: load.get(k) for k in ("exit_code", "wall_s", "cpu_s", "peak_rss_mb")}
    report["load_output"] = bc.text_of(steps, "load")[-4000:]
    if load["exit_code"] != 0:
        report["outcome"] = "refused or failed"
        (args.results / "case_study.json").write_text(json.dumps(report, indent=1) + "\n")
        print(report["load_output"])
        return 1

    verify = bc.measured("verify", ["cgload", "verify", "--assembly-id", "1"], steps, env=env)
    report["verify"] = {"exit_code": verify["exit_code"],
                        "tiers": bc.tiers(bc.text_of(steps, "verify"))}

    with sqlite3.connect(db) as connection:
        report["features_by_type"] = dict(connection.execute(
            "SELECT source_type, count(DISTINCT coalesce(source_id, feature_id)) FROM feature "
            "GROUP BY source_type ORDER BY 2 DESC").fetchall())
        report["functional_runs"] = [dict(zip(("tool", "path", "queries", "matched", "unmatched"),
                                              row, strict=True)) for row in connection.execute(
            "SELECT tool, path, query_count, matched_count, unmatched_count FROM annotation_run")]
        report["assignments_by_analysis"] = dict(connection.execute(
            "SELECT analysis, count(*) FROM protein_annotation GROUP BY analysis "
            "ORDER BY 2 DESC").fetchall())
        report["features_with_any_assignment"] = connection.execute(
            "SELECT count(DISTINCT feature_id) FROM protein_annotation").fetchone()[0]
        top_terms = {}
        for analysis in ANALYSES:
            found = connection.execute(
                "SELECT accession, count(DISTINCT feature_id) FROM protein_annotation "
                "WHERE analysis = ? GROUP BY accession ORDER BY 2 DESC, accession LIMIT 1",
                (analysis,)).fetchone()
            if found:
                top_terms[analysis] = {"accession": found[0], "features": found[1]}
        report["top_terms"] = top_terms

    table = args.results / "case_study_queries.tsv"
    if table.exists():
        table.unlink()
    for analysis, term in top_terms.items():
        command = ["cgload", "query", "--assembly-id", "1", "--term", term["accession"],
                   "--analysis", analysis]
        for run in range(1, args.repeats + 1):
            started = time.perf_counter()
            done = subprocess.run(command, env={**os.environ, **env}, capture_output=True,
                                  text=True, check=False)
            wall = time.perf_counter() - started
            rows = max(len([line for line in done.stdout.splitlines() if line.strip()]) - 1, 0)
            bc.append_row(table, {"analysis": analysis, "term": term["accession"], "run": run,
                                  "wall_s": round(wall, 4), "exit_code": done.returncode,
                                  "rows": rows})

    exported = work / "export.gff3"
    bc.measured("export", ["cgload", "export", "--assembly-id", "1", "--format", "gff3",
                           "--output", str(exported)], steps, env=env)
    with args.config.open("rb") as handle:
        annotation = Path(tomllib.load(handle)["files"]["annotation"][0]["path"])
    if not annotation.is_absolute():
        annotation = args.config.parent / annotation
    truth, stored = work / "truth.sqlite", work / "stored.sqlite"
    bc.measured("truth_input", [sys.executable, str(bc.HERE / "truth.py"), str(annotation),
                                str(truth)], steps)
    bc.measured("truth_export", [sys.executable, str(bc.HERE / "truth.py"), str(exported),
                                 str(stored), "--ignore-attribute", "cgload_linkage_method"], steps)
    fidelity = work / "fidelity.json"
    bc.measured("compare", [sys.executable, str(bc.HERE / "compare.py"), str(truth), str(stored),
                            "--tool", "cgload-sqlite", "--genome", "case_study",
                            "--json", str(fidelity)], steps)
    if fidelity.exists():
        report["round_trip"] = json.loads(fidelity.read_text())["totals"]

    data = work / "data"
    report["for_the_figure"] = {
        "fasta": str(next(data.glob("**/assembly.fa"), "")),
        "annotation": str(next(data.glob("**/annotation.gff3.gz"), "")),
    }
    report["outcome"] = "loaded"
    (args.results / "case_study.json").write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({k: report[k] for k in ("load", "verify", "assignments_by_analysis",
                                             "top_terms")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
