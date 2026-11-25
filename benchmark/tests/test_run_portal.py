"""run_portal.py, end to end on a small synthetic genome (PROTOCOL.md, section 13).

Every step runs, against a real ``cgload serve``, with the sizes cut down so
the whole run takes seconds: the tables it writes are the ones the paper
reads, so their columns and the checks that fill them are tested here.
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
import run_portal  # noqa: E402

#: The yeast accession, so that step 8 (the demo) runs too; the genome is synthetic.
ACCESSION = run_portal.YEAST


def free_port() -> int:
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


def table(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.mark.skipif(shutil.which("cgload") is None, reason="cgload is not installed")
def test_every_step_runs_and_writes_its_table(tmp_path, monkeypatch):
    genomes = tmp_path / "genomes.tsv"
    genomes.write_text("accession\tclade\torganism\tgenus\tspecies\tstrain\tassembly_name\troles\n"
                       f"{ACCESSION}\ttest\tTestus primus\tTestus\tprimus\treference\tT1\tscale\n")
    monkeypatch.setattr(bc, "GENOMES", genomes)
    make_synthetic.main([str(tmp_path / "inputs" / ACCESSION), "--genes", "12"])
    results = tmp_path / "results"

    assert run_portal.main([
        ACCESSION, "--inputs", str(tmp_path / "inputs"), "--work", str(tmp_path / "work"),
        "--results", str(results), "--port", str(free_port()), "--threads", "4",
        "--sample", "5", "--translations", "10", "--seconds", "0.5", "--fuzz", "60"]) == 0

    loads = table(results / "portal_loads.tsv")
    assert [row["outcome"] for row in loads] in (["loaded"], ["loaded_with_repair"])
    assert table(results / "portal_index.tsv")[0]["check"] == "PASS"

    requests = table(results / "portal_requests.tsv")
    kinds = {row["kind"] for row in requests}
    assert {"search", "gene_page", "region_page", "protein", "range_annotation",
            "range_assembly"} <= kinds
    assert all(200 <= int(row["status"]) < 400 for row in requests)
    assert all(row["found"] == "True" for row in table(results / "portal_found.tsv"))

    levels = table(results / "portal_concurrency.tsv")
    assert [int(row["clients"]) for row in levels] == [1, 8, 32]
    assert all(int(row["errors"]) == 0 and int(row["requests"]) > 0 for row in levels)

    security = table(results / "portal_security.tsv")[0]
    assert int(security["server_errors"]) == 0
    assert int(security["missing_security_headers"]) == 0
    assert security["failures"] == ""

    translation = table(results / "portal_translation.tsv")[0]
    assert int(translation["checked"]) == 10 and int(translation["failed"]) == 0

    readonly = table(results / "portal_readonly.tsv")[0]
    assert readonly["unchanged"] == "True"
    assert readonly["sha256_before"] == readonly["sha256_after"]

    # The demo records its own copies of the inputs, not this machine's folders.
    import sqlite3
    import tarfile

    with tarfile.open(results / "demo.tar.gz") as archive:
        archive.extractall(tmp_path / "unpacked", filter="data")
    demo = tmp_path / "unpacked" / "demo"
    with sqlite3.connect(demo / "demo.db") as connection:
        paths = [row[0] for row in connection.execute("SELECT path FROM source_file")]
    assert paths and not any(str(tmp_path / "inputs") in path for path in paths)
    assert next(demo.glob("data/*/*/*/annotation.gff3.gz"))

    summary = (results / "portal_summary.txt").read_text()
    assert "database unchanged by serving: True" in summary
    assert "search found the gene: 5 of 5; not found with at most 1,000 matches" in summary
    assert all(int(row["matches"]) >= 1 for row in table(results / "portal_found.tsv"))
