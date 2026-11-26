"""run_blast.py, end to end on a small synthetic genome (PROTOCOL.md, section 14).

run_portal.py loads the genome first, as on the server; then every step of
run_blast.py runs against a real ``cgload serve`` with BLAST+, with the sizes
cut down. The tables it writes are the ones the paper reads.
"""

from __future__ import annotations

import csv
import shutil
import socket
import sys
from pathlib import Path

import pytest

BENCH = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BENCH))
sys.path.insert(0, str(BENCH / "tests"))

pytest.importorskip("flask")
pytest.importorskip("waitress")

import benchcommon as bc  # noqa: E402
import make_synthetic  # noqa: E402
import run_blast  # noqa: E402
import run_portal  # noqa: E402

ACCESSION = run_portal.YEAST


def free_port() -> int:
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.skipif(shutil.which("cgload") is None, reason="cgload is not installed")
@pytest.mark.skipif(shutil.which("blastp") is None or shutil.which("makeblastdb") is None,
                    reason="BLAST+ is not installed")
def test_every_step_runs_and_writes_its_table(tmp_path, monkeypatch):
    genomes = tmp_path / "genomes.tsv"
    genomes.write_text("accession\tclade\torganism\tgenus\tspecies\tstrain\tassembly_name\troles\n"
                       f"{ACCESSION}\ttest\tTestus primus\tTestus\tprimus\treference\tT1\tscale\n")
    monkeypatch.setattr(bc, "GENOMES", genomes)
    make_synthetic.main([str(tmp_path / "inputs" / ACCESSION), "--genes", "12"])
    work = tmp_path / "work"
    assert run_portal.main([
        ACCESSION, "--inputs", str(tmp_path / "inputs"), "--work", str(work),
        "--results", str(tmp_path / "portal"), "--port", str(free_port()), "--threads", "4",
        "--sample", "2", "--translations", "2", "--seconds", "0.2", "--fuzz", "12",
        "--no-demo"]) == 0
    results = tmp_path / "results"
    assert run_blast.main([
        ACCESSION, "--work", str(work), "--results", str(results), "--port", str(free_port()),
        "--proteins", "4", "--pieces", "4", "--own", "2", "--across", "1", "--seconds", "0.5",
        "--fuzz", "40"]) == 0

    build = table(results / "blast_build.tsv")[0]
    assert build["exit_code"] == "0"
    assert build["genome_check"] == build["proteins_check"] == "PASS"
    assert build["sequences_in_database"] == build["sequences_read"] == build["sequences_in_blast"]
    assert build["cds_features"] == build["translated"] == build["proteins_in_blast"]
    assert int(build["genome_letters"]) > 0 and int(build["bytes"]) > 0

    searches = table(results / "blast_searches.tsv")
    by_program = {}
    for row in searches:
        by_program.setdefault((row["program"], row["scope"]), []).append(row)
    assert set(by_program) == {("blastp", "all"), ("blastn", "all"), ("tblastn", "own"),
                               ("blastx", "own"), ("tblastn", "all")}
    assert all(row["status"] == "200" for row in searches if row["kind"] != "skipped_length")
    # A random genome's pieces are unique, so each must be found where it lies.
    assert all(row["found"] == "True" for row in by_program[("blastn", "all")])
    assert all(row["found"] == "True" for row in by_program[("blastp", "all")]
               if row["kind"] != "skipped_length")

    levels = table(results / "blast_concurrency.tsv")
    assert [int(row["clients"]) for row in levels] == [0, 1, 2, 4, 8]
    assert all(int(row["other_errors"]) == 0 and int(row["page_errors"]) == 0 for row in levels)
    assert int(levels[0]["pages"]) > 0 and int(levels[1]["searches_done"]) > 0

    security = table(results / "blast_security.tsv")[0]
    assert int(security["server_errors"]) == 0 and security["failures"] == ""
    assert int(security["missing_security_headers"]) == 0

    readonly = table(results / "blast_readonly.tsv")[0]
    assert readonly["database_unchanged"] == readonly["manifests_unchanged"] == "True"
    assert "found itself" in (results / "blast_summary.txt").read_text()
