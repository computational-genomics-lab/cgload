"""``cgload query`` (D-085).

What these assert, stated before they were written:

* Each of the three questions (by ID, by region, by functional term) returns
  whole features, with the right parent and top-level identifiers.
* A region matches rows, not envelopes: a position inside an intron returns the
  gene and the mRNA but not the CDS, whose rows all lie elsewhere.
* Typed accessions reach their stored spelling (K00799 -> ko:K00799), and
  nothing is matched by prefix.
* GFF3 output carries every ancestor, so each Parent resolves within the file.
* The database is unchanged afterwards, and no step depends on the driver's
  bound-parameter limit.
"""

from __future__ import annotations

import hashlib
import random
import sqlite3
from pathlib import Path

import pytest
from cgload.cli.main import cli
from cgload.db.engine import build_engine
from cgload.query import in_region, term_variants
from click.testing import CliRunner
from sqlalchemy import event

FIXTURES = Path(__file__).parent / "fixtures"
FUSARIUM = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
AUGUSTUS = FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3"
EGGNOG = FIXTURES / "eggnog_fgraminearum.emapper.annotations"
INTERPRO = FIXTURES / "interproscan_fgraminearum.tsv"

GENE = "gene-FGSG_11614"  # mRNA rna-XM_011327722.1, six-segment CDS, six exons


def _fasta(path: Path, lengths: dict[str, int]) -> None:
    block = "".join(random.Random(7).choice("ACGT") for _ in range(6000))
    with path.open("w") as handle:
        for name, length in lengths.items():
            sequence = (block * (length // len(block) + 1))[:length]
            handle.write(f">{name}\n")
            for index in range(0, length, 60):
                handle.write(sequence[index : index + 60] + "\n")


def _lengths(gff: Path) -> dict[str, int]:
    """The declared length of each sequence, or 500 bp past its last row."""
    declared: dict[str, int] = {}
    furthest: dict[str, int] = {}
    for line in gff.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            declared[parts[1]] = int(parts[3])
        elif line and not line.startswith("#"):
            columns = line.split("\t")
            furthest[columns[0]] = max(furthest.get(columns[0], 0), int(columns[4]))
    for name, end in furthest.items():
        declared.setdefault(name, end + 500)
    return declared


def _load(directory: Path, gff: Path, functional: list[tuple[Path, str]]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / gff.name).write_text(gff.read_text())
    _fasta(directory / "assembly.fa", _lengths(gff))
    text = (
        '[organism]\ngenus = "Testus"\nspecies = "queryi"\nstrain = "S1"\n\n'
        '[assembly]\nname = "asm1"\nfasta = "assembly.fa"\n\n'
        f'[[files.annotation]]\npath = "{gff.name}"\n'
    )
    for path, tool in functional:
        (directory / path.name).write_text(path.read_text())
        text += f'\n[[files.functional]]\npath = "{path.name}"\ntool = "{tool}"\n'
    (directory / "organism.toml").write_text(text)
    db = directory / "g.db"
    runner = CliRunner()
    assert runner.invoke(cli, ["init", "--db-url", f"sqlite:///{db}"]).exit_code == 0
    loaded = runner.invoke(cli, ["load", "--db-url", f"sqlite:///{db}", "--config",
                                 str(directory / "organism.toml"), "--data-dir",
                                 str(directory / "data")])
    assert loaded.exit_code == 0, loaded.output
    return db


@pytest.fixture(scope="module")
def fusarium(tmp_path_factory):
    return _load(tmp_path_factory.mktemp("fusarium"), FUSARIUM,
                 [(EGGNOG, "eggnog"), (INTERPRO, "interproscan")])


@pytest.fixture(scope="module")
def augustus(tmp_path_factory):
    return _load(tmp_path_factory.mktemp("augustus"), AUGUSTUS, [])


def _query(db: Path, *args: str):
    return CliRunner().invoke(cli, ["query", "--db-url", f"sqlite:///{db}",
                                    "--assembly-id", "1", *args])


def _table(output: str) -> list[dict[str, str]]:
    lines = [line for line in output.splitlines() if line and not line.startswith("#")]
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _gff3_rows(output: str) -> list[list[str]]:
    return [line.split("\t") for line in output.splitlines()
            if line and not line.startswith("#")]


def _attributes(column: str) -> dict[str, str]:
    return dict(pair.split("=", 1) for pair in column.split(";") if pair)


# ---------------------------------------------------------------------------
# by identifier
# ---------------------------------------------------------------------------


def test_an_id_returns_that_feature(fusarium):
    result = _query(fusarium, "--id", "cds-XP_011326024.1")
    assert result.exit_code == 0, result.output
    [row] = _table(result.stdout)
    assert row["type"] == "CDS"
    assert row["segments"] == "6"
    assert (row["start"], row["end"]) == ("3250", "5086")
    assert row["parent"] == "rna-XM_011327722.1"
    assert row["top_level"] == GENE


def test_descendants_return_the_gene_model(fusarium):
    result = _query(fusarium, "--id", GENE, "--descendants")
    assert result.exit_code == 0, result.output
    types = sorted(row["type"] for row in _table(result.stdout))
    assert types == ["CDS", "exon", "exon", "exon", "exon", "exon", "exon", "gene", "mRNA"]
    assert {row["top_level"] for row in _table(result.stdout)} == {GENE}


def test_type_keeps_only_the_named_types(fusarium):
    result = _query(fusarium, "--id", GENE, "--descendants", "--type", "exon")
    assert {row["type"] for row in _table(result.stdout)} == {"exon"}


def test_an_unknown_id_is_an_error(fusarium):
    result = _query(fusarium, "--id", "gene-NOPE")
    assert result.exit_code == 1
    assert "no feature with ID 'gene-NOPE'" in result.output


# ---------------------------------------------------------------------------
# by region
# ---------------------------------------------------------------------------


def test_a_region_matches_rows_not_envelopes(fusarium):
    """3300 lies between the CDS's first two segments. The RefSeq `region` row
    spans the whole chromosome, so it overlaps every position."""

    def found(region: str) -> list[dict[str, str]]:
        rows = _table(_query(fusarium, "--region", region).stdout)
        return [row for row in rows if row["type"] != "region"]

    assert {row["id"] for row in found("NC_026477.1:3300-3300")} == {
        GENE, "rna-XM_011327722.1"}

    by_id = {row["id"]: row for row in found("NC_026477.1:3,260-3,260")}
    assert set(by_id) == {GENE, "rna-XM_011327722.1", "cds-XP_011326024.1",
                          "exon-XM_011327722.1-6"}
    assert by_id["cds-XP_011326024.1"]["segments"] == "6", "the CDS comes back whole"


def test_a_region_with_only_a_start_runs_to_the_end(fusarium):
    to_end = _table(_query(fusarium, "--region", "NC_026477.1:5000000", "--type", "gene").stdout)
    whole = _table(_query(fusarium, "--region", "NC_026477.1", "--type", "gene").stdout)
    assert to_end
    assert {row["id"] for row in to_end} <= {row["id"] for row in whole}
    assert all(int(row["end"]) >= 5_000_000 for row in to_end)


def test_an_end_past_the_sequence_is_cut_and_a_start_past_it_is_refused(fusarium):
    assert _query(fusarium, "--region", "NC_026477.1:8000000-99999999").exit_code == 0
    past = _query(fusarium, "--region", "NC_026477.1:9000000-9000100")
    assert past.exit_code == 2
    assert "8,033,942 bp" in past.output


def test_an_unknown_sequence_is_refused_with_examples(fusarium):
    result = _query(fusarium, "--region", "chrUn:1-100")
    assert result.exit_code == 2
    assert "no sequence named 'chrUn'" in result.output
    assert "NC_026477.1" in result.output


# ---------------------------------------------------------------------------
# by functional term
# ---------------------------------------------------------------------------


def test_term_spellings():
    assert term_variants("K00799") == ["K00799", "ko:K00799"]
    assert term_variants("ko:K00799") == ["K00799", "ko:K00799"]
    assert "ko00010" in term_variants("map00010")
    assert "GO:0005515" in term_variants("go:0005515")
    assert "2.4.1.1" in term_variants("EC:2.4.1.1")
    assert term_variants("PF03198") == ["PF03198"]


def test_a_term_returns_annotated_features_and_their_genes(fusarium):
    result = _query(fusarium, "--term", "K00799")
    assert result.exit_code == 0, result.output
    rows = _table(result.stdout)
    assert rows
    assert {row["accession"] for row in rows} == {"ko:K00799"}
    assert all(row["top_level"].startswith("gene-") for row in rows)
    assert "assignment(s)" in result.stderr


def test_analysis_narrows_a_term(fusarium):
    every = _table(_query(fusarium, "--term", "PF03198").stdout)
    interpro_only = _table(_query(fusarium, "--term", "PF03198", "--analysis", "pfam").stdout)
    assert {row["analysis"] for row in every} == {"Pfam"}
    assert interpro_only
    assert len(interpro_only) <= len(every)


def test_a_term_matches_whole_accessions_only(fusarium):
    result = _query(fusarium, "--term", "PF0319")
    assert result.exit_code == 0
    assert _table(result.stdout) == []
    assert "Stored analyses:" in result.stderr


def test_type_and_term_together_are_refused(fusarium):
    result = _query(fusarium, "--term", "PF03198", "--type", "gene")
    assert result.exit_code == 2
    assert "top_level" in result.output


# ---------------------------------------------------------------------------
# GFF3 output
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("question", [
    ("--id", "cds-XP_011326024.1"),
    ("--region", "NC_026477.1:3260-3260"),
    ("--term", "PF03198"),
])
def test_every_parent_in_gff3_output_resolves(fusarium, question):
    result = _query(fusarium, *question, "--format", "gff3")
    assert result.exit_code == 0, result.output
    rows = _gff3_rows(result.stdout)
    identifiers = {_attributes(row[8]).get("ID") for row in rows}
    parents = {_attributes(row[8])["Parent"] for row in rows if "Parent" in _attributes(row[8])}
    assert parents
    assert parents <= identifiers
    assert "##cgload-export" in result.stdout


TRANS_SPLICED = """\
##gff-version 3
##sequence-region chr1 1 20000
chr1\tRefSeq\tregion\t1\t20000\t.\t+\t.\tID=chr1:1..20000;genome=chromosome
chr1\tRefSeq\tgene\t1000\t1500\t.\t+\t.\tID=gene-ts;exception=trans-splicing;part=1/2
chr1\tRefSeq\tgene\t15000\t15600\t.\t+\t.\tID=gene-ts;exception=trans-splicing;part=2/2
chr1\tRefSeq\tmRNA\t1000\t1500\t.\t+\t.\tID=rna-ts;Parent=gene-ts
chr1\tRefSeq\tmRNA\t15000\t15600\t.\t+\t.\tID=rna-ts;Parent=gene-ts
chr1\tRefSeq\texon\t1000\t1500\t.\t+\t.\tID=exon-ts-1;Parent=rna-ts
chr1\tRefSeq\texon\t15000\t15600\t.\t+\t.\tID=exon-ts-2;Parent=rna-ts
chr1\tRefSeq\tCDS\t1100\t1500\t.\t+\t0\tID=cds-ts;Parent=rna-ts
chr1\tRefSeq\tCDS\t15000\t15500\t.\t+\t2\tID=cds-ts;Parent=rna-ts
"""


def test_descendants_of_a_later_row_include_the_whole_model(tmp_path):
    """A trans-spliced gene and its mRNA are each written on two rows. The
    region overlaps only the second ones; exon-ts-1 hangs from the mRNA's first
    row, far outside the region, and must still come back with --descendants."""
    gff = tmp_path / "ts.gff3"
    gff.write_text(TRANS_SPLICED)
    db = _load(tmp_path / "project", gff, [])
    rows = _table(_query(db, "--region", "chr1:15100-15200", "--descendants").stdout)
    assert {row["id"] for row in rows if row["type"] != "region"} == {
        "gene-ts", "rna-ts", "exon-ts-1", "exon-ts-2", "cds-ts"}


CIRCULAR = """\
##gff-version 3
##sequence-region chrM 1 16569
chrM\tRefSeq\tregion\t1\t16569\t.\t+\t.\tID=chrM:1..16569;Is_circular=true;genome=mitochondrion
chrM\tRefSeq\tgene\t16000\t16700\t.\t+\t.\tID=gene-across;gene_biotype=protein_coding
chrM\tRefSeq\tgene\t5000\t6000\t.\t+\t.\tID=gene-inside;gene_biotype=protein_coding
"""


@pytest.mark.parametrize(("region", "expected"), [
    ("chrM:50-60", {"gene-across"}),          # past the origin: 16,619-16,629 as stored
    ("chrM:16100-16200", {"gene-across"}),
    ("chrM:200-300", set()),                  # after the gene's end at 131
    ("chrM:5500-5500", {"gene-inside"}),
])
def test_a_region_finds_a_feature_across_the_origin(tmp_path, region, expected):
    gff = tmp_path / "mt.gff3"
    gff.write_text(CIRCULAR)
    db = _load(tmp_path / "project", gff, [])
    rows = _table(_query(db, "--region", region).stdout)
    assert {row["id"] for row in rows if row["type"] != "region"} == expected


def test_anonymous_rows_come_with_their_gene(augustus):
    tsv = _table(_query(augustus, "--id", "g1", "--descendants").stdout)
    anonymous = [row for row in tsv if row["id"] == ""]
    assert anonymous, "AUGUSTUS intron and codon rows have no ID"
    assert all(row["top_level"] == "g1" for row in anonymous)
    gff3 = _gff3_rows(_query(augustus, "--id", "g1", "--descendants", "--format", "gff3").stdout)
    assert all("Parent=" in row[8] for row in gff3 if row[2] == "intron")


# ---------------------------------------------------------------------------
# read-only, and independent of the parameter limit
# ---------------------------------------------------------------------------


def test_queries_change_nothing(fusarium):
    before = hashlib.sha256(fusarium.read_bytes()).hexdigest()
    for question in (("--id", GENE, "--descendants"), ("--region", "NC_026477.1"),
                     ("--term", "PF03198", "--format", "gff3")):
        assert _query(fusarium, *question).exit_code == 0
    assert hashlib.sha256(fusarium.read_bytes()).hexdigest() == before


def test_a_whole_sequence_does_not_depend_on_the_parameter_limit(fusarium):
    engine = build_engine(f"sqlite:///{fusarium}")

    @event.listens_for(engine, "connect")
    def _small_limit(dbapi_connection, _record):
        dbapi_connection.setlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER, 10)

    with engine.connect() as connection:
        result = in_region(connection, 1, "NC_026477.1", ancestors=True)
    assert len(result.features) == 354


def test_exactly_one_question(fusarium):
    assert _query(fusarium).exit_code == 2
    both = _query(fusarium, "--id", GENE, "--region", "NC_026477.1")
    assert both.exit_code == 2
    assert "exactly one of" in both.output
