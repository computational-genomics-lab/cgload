"""Tests of the benchmark's own scripts: the reference count and the comparison.

Run with ``pytest benchmark/tests``. Not part of cgload's suite (pyproject's
testpaths), because the benchmark's numbers are only as good as these scripts
and they deserve their own check, not because they are optional.
"""

from __future__ import annotations

import json
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
FIXTURES = BENCH.parent / "tests" / "fixtures"
sys.path.insert(0, str(BENCH))

import compare  # noqa: E402
import gffmodel  # noqa: E402
import truth  # noqa: E402

HARD = (
    "##gff-version 3\n"
    "##sequence-region chr1 1 10000\n"
    "##sequence-region chrM 1 1000\n"
    "chrM\t.\tregion\t1\t1000\t.\t+\t.\tID=chrM:1..1000;Is_circular=true;genome=mitochondrion\n"
    "chr1\t.\tgene\t100\t900\t.\t+\t.\tID=g1;Dbxref=a:1,b:2\n"
    "chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=t1;Parent=g1\n"
    "chr1\t.\tCDS\t100\t200\t.\t+\t0\tID=c1;Parent=t1;Note=first\n"
    "chr1\t.\tCDS\t200\t300\t.\t+\t2\tID=c1;Parent=t1;Note=second\n"
    "chr1\t.\tpseudogene\t2000\t2500\t.\t-\t.\tID=p1\n"
    "chr1\t.\texon\t2000\t2500\t.\t-\t.\tID=e1;Parent=p1\n"
    "chr1\t.\tintron\t201\t299\t.\t+\t.\tParent=t1\n"
    "chrM\t.\tgene\t900\t1100\t.\t?\t.\tID=g2;exception=trans-splicing\n"
    "chrM\t.\tmRNA\t50\t60\t.\t?\t.\tID=t2;Parent=g2\n"
    "chr1\t.\tgene\t9990\t10050\t.\t+\t.\tID=g3\n"
    "chr1\t.\tgene\t500\t400\t.\t+\t.\tID=g4\n"
    "chr1\t.\tmRNA\t3000\t3100\t.\t+\t.\tID=t5;Parent=g9\n"
)


def build(tmp_path: Path, text: str, name: str = "a.gff3") -> tuple[Path, dict]:
    gff = tmp_path / name
    gff.write_text(text)
    out = tmp_path / (name + ".sqlite")
    return out, truth.build(gff, out)


def flags(db: Path) -> dict[str, set[str]]:
    with sqlite3.connect(db) as connection:
        return {key: set(filter(None, value.split(","))) for key, value in
                connection.execute("SELECT key, flags FROM feature")}


def test_the_reference_count_follows_the_specification(tmp_path):
    db, summary = build(tmp_path, HARD)
    assert summary["rows"] == 13
    # c1 is one feature on two rows; the intron, with no ID, is a feature of its own.
    assert summary["features"] == 12
    assert summary["types"]["CDS"] == 1
    found = flags(db)
    assert {"multi_row", "overlapping_rows", "rows_disagree"} <= found["c1"]
    assert "multi_valued" in found["g1"]
    assert "exon_under_pseudogene" in found["e1"]
    assert {"crosses_origin", "trans_spliced", "strand_unknown", "organelle"} <= found["g2"]
    assert "escapes_parent" in found["t2"]
    assert "past_end" in found["g3"]
    assert "inverted" in found["g4"]
    assert "dangling_parent" in found["t5"]
    assert any("no_id" in value for key, value in found.items() if key.startswith("anon:"))


def test_rows_that_disagree_are_compared_row_by_row(tmp_path):
    """A loader that keeps both Notes but gives each to the other row has
    changed the feature, though no value is lost (D-104 was such a case)."""
    db, _ = build(tmp_path, HARD)
    with sqlite3.connect(db) as connection:
        (digest,) = connection.execute("SELECT attrs FROM feature WHERE key = 'c1'").fetchone()
    assert digest == gffmodel.feature_digest([(100, 200, {"Note": ["first"]}),
                                              (200, 300, {"Note": ["second"]})])
    swapped, _ = build(tmp_path, HARD.replace("Note=first", "Note=tmp")
                       .replace("Note=second", "Note=first").replace("Note=tmp", "Note=second"),
                       "swapped.gff3")
    result = compare.compare(db, swapped)
    assert result["field_differences"] == {"attrs": 1}
    assert result["examples"]["changed:attrs"] == ["c1"]


def test_a_file_compared_with_itself_keeps_everything(tmp_path):
    db, _ = build(tmp_path, HARD)
    result = compare.compare(db, db)
    assert result["totals"]["kept"] == result["totals"]["truth"] == 12
    assert result["totals"].get("lost", 0) == 0


def test_the_comparison_sees_each_kind_of_damage(tmp_path):
    original, _ = build(tmp_path, HARD)
    damaged_text = (HARD
                    .replace("chr1\t.\texon\t2000\t2500\t.\t-\t.\tID=e1;Parent=p1\n", "")
                    .replace("ID=g1;Dbxref=a:1,b:2", "ID=g1;Dbxref=a:1")
                    .replace("mRNA\t100\t900\t.\t+", "mRNA\t100\t900\t.\t-"))
    damaged, _ = build(tmp_path, damaged_text, "b.gff3")
    result = compare.compare(original, damaged)
    assert result["totals"]["lost"] == 1
    assert result["field_differences"]["attrs"] == 1
    assert result["field_differences"]["strand"] == 1
    assert result["examples"]["lost"] == ["e1"]


def test_a_reported_change_is_declared_and_a_logged_loss_is_announced(tmp_path):
    original, _ = build(tmp_path, HARD)
    damaged, _ = build(tmp_path, HARD.replace("ID=t1;Parent=g1", "ID=t1;Parent=g1;x=1")
                       .replace("chr1\t.\tgene\t9990\t10050\t.\t+\t.\tID=g3\n", ""), "c.gff3")
    log = tmp_path / "tool.log"
    log.write_text("warning: skipped g3, which runs past the end of chr1\n")
    result = compare.compare(original, damaged, declared={"t1"}, logs=[log])
    assert result["totals"]["changed_declared"] == 1
    assert result["totals"]["lost_announced"] == 1
    assert result["totals"].get("lost_silent", 0) == 0


def test_the_exporters_own_attribute_can_be_left_out(tmp_path):
    plain, _ = build(tmp_path, HARD)
    gff = tmp_path / "exported.gff3"
    gff.write_text(HARD.replace("ID=t1;Parent=g1", "ID=t1;Parent=g1;cgload_linkage_method=x"))
    exported = tmp_path / "exported.sqlite"
    truth.build(gff, exported, frozenset({"cgload_linkage_method"}))
    assert compare.compare(plain, exported)["totals"]["kept"] == 12


def test_measure_records_a_missing_program_as_a_failure(tmp_path):
    record = tmp_path / "m.json"
    code = subprocess.run([sys.executable, str(BENCH / "measure.py"), "--json", str(record),
                           "--", "no-such-program-here"], check=False).returncode
    assert code == 127
    assert json.loads(record.read_text())["exit_code"] == 127


def test_measure_kills_a_command_that_runs_too_long(tmp_path):
    record = tmp_path / "m.json"
    code = subprocess.run([sys.executable, str(BENCH / "measure.py"), "--json", str(record),
                           "--timeout", "1", "--", "sleep", "30"], check=False).returncode
    assert code == 124
    assert json.loads(record.read_text())["timed_out"] is True


@pytest.mark.skipif(shutil.which("cgload") is None, reason="cgload is not installed")
def cgload_round_trip(tmp_path: Path, text: str, name: str = "a.gff3") -> dict:
    """The benchmark's round trip, as run_cgload.py runs it: load, export, and
    the reference count of each file compared."""
    gff = tmp_path / name
    gff.write_text(text)
    lengths = {line.split()[1]: int(line.split()[3]) for line in text.splitlines()
               if line.startswith("##sequence-region")}
    with (tmp_path / "a.fa").open("w") as handle:
        for sequence, length in lengths.items():
            handle.write(f">{sequence}\n" + "ACGT" * (length // 4) + "A" * (length % 4) + "\n")
    (tmp_path / "o.toml").write_text(
        '[organism]\ngenus = "F"\nspecies = "g"\n[assembly]\nname = "a"\nfasta = "a.fa"\n'
        f'[[files.annotation]]\npath = "{gff.name}"\n')
    env = {"CGLOAD_DB_URL": f"sqlite:///{tmp_path / 'g.db'}", "PATH": __import__("os").environ["PATH"]}
    for command in (["cgload", "init"],
                    ["cgload", "load", "--config", str(tmp_path / "o.toml"), "--data-dir",
                     str(tmp_path / "data")],
                    ["cgload", "export", "--assembly-id", "1", "--format", "gff3", "--output",
                     str(tmp_path / "out.gff3")]):
        done = subprocess.run(command, env=env, capture_output=True, text=True, check=False)
        assert done.returncode == 0, done.stdout + done.stderr
    before, after = tmp_path / "t.sqlite", tmp_path / "s.sqlite"
    truth.build(gff, before)
    truth.build(tmp_path / "out.gff3", after, frozenset({"cgload_linkage_method"}))
    return compare.compare(before, after)


def test_cgload_round_trips_a_real_fixture_exactly(tmp_path):
    source = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
    result = cgload_round_trip(tmp_path, source.read_text(), source.name)
    assert result["totals"]["kept"] == result["totals"]["truth"] == 354


def test_cgload_keeps_each_rows_attributes_when_rows_run_highest_first(tmp_path):
    """NCBI writes a minus-strand feature's rows from the highest coordinate
    down. The preflight of 22 October 2026 found the lowest row's own
    attributes replaced by the first row's (D-104)."""
    text = (
        "##gff-version 3\n"
        "##sequence-region chr1 1 10000\n"
        "chr1\tRefSeq\tgene\t6000\t9000\t.\t-\t.\tID=gene-c\n"
        "chr1\tRefSeq\tmRNA\t6000\t9000\t.\t-\t.\tID=rna-c;Parent=gene-c\n"
        "chr1\tRefSeq\texon\t8000\t9000\t.\t-\t.\tID=exon-c-1;Parent=rna-c\n"
        "chr1\tRefSeq\texon\t6000\t6500\t.\t-\t.\tID=exon-c-2;Parent=rna-c\n"
        "chr1\tRefSeq\tCDS\t8000\t8800\t.\t-\t0\tID=cds-c;Parent=rna-c;partial=true;"
        "end_range=8800,.\n"
        "chr1\tRefSeq\tCDS\t6200\t6500\t.\t-\t2\tID=cds-c;Parent=rna-c;partial=true\n"
        "chr1\tRefSeq\tcDNA_match\t8000\t9000\t.\t-\t.\tID=aln-2;Target=NM_3.1 1 1001 +;"
        "Gap=M1001\n"
        "chr1\tRefSeq\tcDNA_match\t6000\t6500\t.\t-\t.\tID=aln-2;Target=NM_3.1 1002 1502 +;"
        "Gap=M501\n"
    )
    result = cgload_round_trip(tmp_path, text)
    assert result["totals"]["kept"] == result["totals"]["truth"] == 6, result["examples"]


TXDB_GFF = (
    "##gff-version 3\n"
    "chr1\t.\tgene\t100\t900\t.\t+\t.\tID=gene-A;Name=A\n"
    "chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=rna-1;Parent=gene-A;transcript_id=NM_1.1\n"
    "chr1\t.\texon\t100\t300\t.\t+\t.\tID=exon-1-1;Parent=rna-1;transcript_id=NM_1.1\n"
    "chr1\t.\texon\t500\t900\t.\t+\t.\tID=exon-1-2;Parent=rna-1;transcript_id=NM_1.1\n"
    "chr1\t.\tCDS\t150\t300\t.\t+\t0\tID=cds-1;Parent=rna-1\n"
    "chr1\t.\tCDS\t500\t800\t.\t+\t0\tID=cds-1;Parent=rna-1\n"
    "chr1\t.\tmRNA\t100\t900\t.\t+\t.\tID=rna-2;Parent=gene-A;transcript_id=NM_2.1\n"
    "chr1\t.\texon\t100\t900\t.\t+\t.\tID=exon-2-1;Parent=rna-2;transcript_id=NM_2.1\n"
    "chr1\t.\tCDS\t150\t800\t.\t+\t0\tID=cds-2;Parent=rna-2\n"
    "chr1\t.\tgene\t2000\t3000\t.\t-\t.\tID=gene-B;Name=B\n"
    "chr1\t.\tmRNA\t2000\t3000\t.\t-\t.\tID=rna-3;Parent=gene-B;transcript_id=NM_3.1\n"
    "chr1\t.\texon\t2000\t3000\t.\t-\t.\tID=exon-3-1;Parent=rna-3;transcript_id=NM_3.1\n"
    "chr1\t.\tpseudogene\t4000\t4500\t.\t+\t.\tID=gene-P;Name=P\n"
    "chr1\t.\texon\t4000\t4500\t.\t+\t.\tID=exon-P-1;Parent=gene-P\n"
    "chr1\t.\tregion\t1\t10000\t.\t+\t.\tID=chr1:1..10000\n"
)


def write_txdb(directory: Path) -> None:
    """What run_txdb.R would write for TXDB_GFF if txdbmaker kept rna-1 and the
    pseudogene, cut rna-2's CDS short and dropped rna-3."""
    directory.mkdir()
    kinds = {"gene": ("gene", "pseudogene"), "transcript": ("mRNA",), "exon": ("exon",),
             "cds": ("CDS",), "stop_codon": ("stop_codon",),
             "imported": ("gene", "pseudogene", "mRNA", "exon", "CDS", "stop_codon")}
    (directory / "types.txt").write_text("".join(
        f"{kind}\ttest\t{value}\n" for kind, values in kinds.items() for value in values))
    (directory / "transcripts.tsv").write_text(
        "tx_id\ttx_name\tseqid\tstrand\tstart\tend\ttype\tgenes\n"
        "1\tNM_1.1\tchr1\t+\t100\t900\tmRNA\tA\n"
        "2\tNM_2.1\tchr1\t+\t100\t900\tmRNA\tA\n"
        "3\t\tchr1\t+\t4000\t4500\tpseudogene\tP\n")
    (directory / "exons.tsv").write_text(
        "tx_id\trank\tseqid\tstrand\tstart\tend\tname\n"
        "1\t1\tchr1\t+\t100\t300\t\n1\t2\tchr1\t+\t500\t900\t\n"
        "2\t1\tchr1\t+\t100\t900\t\n3\t1\tchr1\t+\t4000\t4500\t\n")
    (directory / "cds.tsv").write_text(
        "tx_id\trank\tseqid\tstrand\tstart\tend\tname\n"
        "1\t1\tchr1\t+\t150\t300\t\n1\t2\tchr1\t+\t500\t800\t\n"
        "2\t1\tchr1\t+\t150\t700\t\n")
    (directory / "genes.tsv").write_text("gene_id\nA\nP\n")


def test_the_txdb_comparison_follows_txdbmakers_rules(tmp_path):
    import compare_txdb

    db, _ = build(tmp_path, TXDB_GFF)
    write_txdb(tmp_path / "txdb")
    log = tmp_path / "load.err"
    log.write_text("WARNING: The following transcripts were dropped because ...: rna-3\n")
    result = compare_txdb.compare(db, tmp_path / "txdb", tmp_path / "a.gff3", [log])
    assert result["transcript_names_from"] == "transcript_id"
    assert result["gene_ids_from"] == "Name"
    # rna-1, rna-2, rna-3, and the pseudogene, whose exon hangs directly from it.
    assert result["totals"] == {"truth": 4, "kept": 2, "changed": 1, "changed_silent": 1,
                                "lost": 1, "lost_announced": 1}
    assert result["field_differences"] == {"cds": 1}
    assert result["examples"] == {"changed:cds": ["rna-2"], "lost": ["rna-3"]}
    assert result["genes"] == {"truth": 3, "held": 2, "absent": 1}
    assert result["not_modelled_by_type"] == {"region": 1}
    assert result["by_case"]["exon_under_pseudogene"] == {"truth": 1, "kept": 1}


def test_txdbmaker_names_genes_by_id_when_names_do_not_tell_them_apart(tmp_path):
    import compare_txdb

    gff = tmp_path / "g.gff3"
    gff.write_text("chr1\t.\tgene\t1\t9\t.\t+\t.\tID=gene-X;Name=SHOX\n"
                   "chrY\t.\tgene\t1\t9\t.\t+\t.\tID=gene-X-2;Name=SHOX\n")
    kinds = {kind: set() for kind in compare_txdb.KINDS}
    kinds["gene"] = kinds["imported"] = {"gene"}
    names = compare_txdb.txdb_names(gff, kinds)
    assert names["gene_rule"] == "ID"
    assert names["gene"] == {"gene-X": "X", "gene-X-2": "X-2"}


def test_a_kill_is_judged_against_a_complete_load():
    import kills

    assert kills.state_of(0, 0, (10, 5)) == "empty"
    assert kills.state_of(10, 5, (10, 5)) == "complete"
    assert kills.state_of(10, 2, (10, 5)) == "partial"
    assert kills.state_of("unreadable: x", 0, (10, 5)) == "unreadable"


def test_the_timing_reads_each_tools_answer(tmp_path):
    import query_timing

    stats = "Org / S / a\n  by source type: exon 6,006, gene 2,004, CDS 2,002, region 2\n"
    assert query_timing.answer_of("cgload-sqlite", "count", stats, tmp_path / "x") == "10014"
    rows = "# cgload query: ...\nid\ttype\na\tgene\nb\tmRNA\n"
    assert query_timing.answer_of("cgload-sqlite", "gene", rows, tmp_path / "x") == "2"
    gff_lines = "chr1\t.\tgene\t1\t9\t.\t+\t.\tID=a\nchr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=b\n"
    assert query_timing.answer_of("gffutils", "gene", gff_lines, tmp_path / "x") == "2"
    assert query_timing.answer_of("gffutils", "count", "gene\t1\nmRNA\t1\n2\n",
                                  tmp_path / "x") == "2"
