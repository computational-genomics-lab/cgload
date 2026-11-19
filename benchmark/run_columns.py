#!/usr/bin/env python3
"""The round trip again, for GFF3 columns 2, 6 and 8: does cgload keep them?

    python3 benchmark/run_columns.py --inputs DIR --work DIR --results DIR [ACCESSION ...]

The benchmark's round trip (run_cgload.py) compared each feature's sequence,
type, strand, row spans, parents and attributes, but not the program that
wrote each row (column 2), its score (column 6) or its phase (column 8). For
each genome (all 30 of genomes.tsv when none is named) this script repeats
cgload's part of the round trip on SQLite, the same way run_cgload.py does it:

1. ``cgload init`` on a new SQLite file;
2. ``cgload load``, and with ``--repair-envelopes`` if the file is refused;
3. ``cgload export`` to GFF3;
4. ``columns.py`` on the input and the export, row by row.

It appends one row per genome to RESULTS/columns.tsv, keeps each comparison
as RESULTS/columns/ACCESSION.json, and deletes the database, data folder and
export. A genome already in the table is skipped, so a stopped run can be
started again. It changes nothing in the benchmark's other results.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402
import columns  # noqa: E402
from run_cgload import TIMEOUT, write_config  # noqa: E402

COLUMNS = (
    "accession", "cgload_version", "cgload_commit", "outcome", "load_wall_s",
    "rows_input", "rows_output", "rows_compared", "source_differs", "source_encoding_only",
    "score_differs", "score_text_only", "phase_differs", "features_row_count_differs",
    "features_missing", "features_extra", "rows_span_differs",
    "input_score_beyond_six_digits", "input_source_differs_from_first_row",
)


def cgload_version() -> str:
    """The version of the ``cgload`` command that runs, as it reports itself: the
    installed package's metadata can lag behind an editable install."""
    import subprocess

    shown = subprocess.run(["cgload", "--version"], capture_output=True, text=True,
                           check=False).stdout.split()
    return shown[-1] if shown else "unknown"


def run_one(accession: str, args) -> None:
    table = args.results / "columns.tsv"
    if bc.done(table, accession=accession):
        print(f"{accession}: already in {table}, skipped")
        return
    genome = bc.genome(accession)
    inputs = args.inputs / accession
    work = args.work / "columns" / accession
    if work.exists():
        shutil.rmtree(work)
    steps = work / "steps"
    steps.mkdir(parents=True)
    config = work / "organism.toml"
    write_config(config, genome, inputs)
    env = {"CGLOAD_DB_URL": f"sqlite:///{work / 'cgload.db'}"}
    row = dict.fromkeys(COLUMNS)
    row.update({"accession": accession, "cgload_version": cgload_version(),
                "cgload_commit": bc.git_commit(bc.HERE.parent)})

    bc.measured("init", ["cgload", "init"], steps, env=env)
    load_command = ["cgload", "load", "--config", str(config), "--data-dir", str(work / "data")]
    load = bc.measured("load", load_command, steps, timeout=TIMEOUT, env=env)
    row["outcome"] = "loaded"
    if load["exit_code"] != 0:
        load = bc.measured("load_repair", [*load_command, "--repair-envelopes"], steps,
                           timeout=TIMEOUT, env=env)
        row["outcome"] = "loaded_with_repair"
    row["load_wall_s"] = load["wall_s"]
    if load["exit_code"] != 0:
        row["outcome"] = "refused"
    else:
        exported = work / "export.gff3"
        export = bc.measured("export", ["cgload", "export", "--assembly-id", "1",
                                        "--format", "gff3", "--output", str(exported)],
                             steps, timeout=TIMEOUT, env=env)
        if export["exit_code"] != 0:
            row["outcome"] = "export_failed"
        else:
            result = columns.compare(inputs / "genomic.gff", exported,
                                     frozenset({"cgload_linkage_method"}))
            for name in COLUMNS:
                if name in result["counts"]:
                    row[name] = result["counts"][name]
            (args.results / "columns").mkdir(parents=True, exist_ok=True)
            (args.results / "columns" / f"{accession}.json").write_text(
                json.dumps(result, indent=1) + "\n")
    bc.append_row(table, row)
    print(json.dumps(row))
    if row["outcome"] in ("loaded", "loaded_with_repair"):
        shutil.rmtree(work, ignore_errors=True)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accessions", nargs="*")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    args = parser.parse_args(argv)
    accessions = args.accessions or [genome["accession"] for genome in bc.genomes()]
    for accession in accessions:
        run_one(accession, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
