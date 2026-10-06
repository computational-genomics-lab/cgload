"""Audit findings, one test each, written against cgload 2.0.0.dev0.

Every test here asserts the behaviour the README or the GFF3 specification
promises, and every one fails on the code as audited. They are marked
``xfail(strict=True)`` -- the convention this suite already uses for D-059 -- so
the file can be merged without turning CI red, and so that fixing a finding
turns its test into an XPASS failure that forces the marker to be removed.

To see the current failures rather than the expected-failure summary:

    pytest --runxfail tests/test_audit_findings.py

Findings AUDIT-1 to AUDIT-4 are silent corruption: cgload stores a wrong
database or data directory and still reports success. They matter most.
"""

from __future__ import annotations

import random
import sqlite3
from pathlib import Path

import pytest
from cgload.cli.main import cli
from cgload.db.engine import build_engine
from cgload.db.export_source import load_features
from cgload.db.init import create_schema
from cgload.db.schema import feature, sequence_region
from cgload.storage import staged_load
from click.testing import CliRunner
from sqlalchemy import BigInteger, event
from sqlalchemy.dialects import postgresql
from sqlalchemy.schema import CreateTable

FIXTURES = Path(__file__).parent / "fixtures"
EXAMPLE = Path(__file__).parent.parent / "examples" / "minimal"


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _fasta(path: Path, lengths: dict[str, int]) -> Path:
    rnd = random.Random(3)
    with path.open("w") as handle:
        for name, length in lengths.items():
            handle.write(f">{name}\n")
            sequence = "".join(rnd.choice("ACGT") for _ in range(length))
            for index in range(0, length, 60):
                handle.write(sequence[index : index + 60] + "\n")
    return path


def _lengths_of(gff: Path, pad: int = 500) -> dict[str, int]:
    declared: dict[str, int] = {}
    furthest: dict[str, int] = {}
    for line in gff.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            declared[parts[1]] = int(parts[3])
        elif line.startswith("##FASTA"):
            break
        elif line and not line.startswith("#"):
            cols = line.split("\t")
            if len(cols) >= 9:
                furthest[cols[0]] = max(furthest.get(cols[0], 0), int(cols[4]))
    for seqid, end in furthest.items():
        declared.setdefault(seqid, end + pad)
    return declared


def _project(tmp_path: Path, gff_text: str | None = None, *, fixture: Path | None = None,
             lengths: dict[str, int] | None = None) -> Path:
    """Write annotation, FASTA and organism.toml; return the config path."""
    gff = tmp_path / (fixture.name if fixture else "annotation.gff3")
    gff.write_text(fixture.read_text() if fixture else gff_text)
    _fasta(tmp_path / "assembly.fa", lengths or _lengths_of(gff))
    config = tmp_path / "organism.toml"
    config.write_text(
        "[organism]\ngenus = \"Auditus\"\nspecies = \"testii\"\nstrain = \"S1\"\n"
        "ncbi_taxon_id = 1\n\n[assembly]\nname = \"asm1\"\nfasta = \"assembly.fa\"\n\n"
        f"[[files.annotation]]\npath = \"{gff.name}\"\n"
    )
    return config


def _cli(tmp_path: Path, *args: str):
    url = f"sqlite:///{tmp_path / 'g.db'}"
    return CliRunner().invoke(cli, [*args[:1], "--db-url", url, *args[1:]])


def _init_and_load(tmp_path: Path, config: Path, *extra: str):
    assert _cli(tmp_path, "init").exit_code == 0
    return _cli(tmp_path, "load", "--config", str(config),
                "--data-dir", str(tmp_path / "data"), *extra)


def _rows(tmp_path: Path, query: str):
    connection = sqlite3.connect(tmp_path / "g.db")
    try:
        return connection.execute(query).fetchall()
    finally:
        connection.close()


REFSEQ_HEADER = (
    "##gff-version 3\n"
    "##sequence-region NC_000001.1 1 20000\n"
)
NUCLEAR_GENE = (
    "NC_000001.1\tGnomon\tgene\t1000\t3000\t.\t+\t.\tID=gene-A;gene_biotype=protein_coding\n"
    "NC_000001.1\tGnomon\tmRNA\t1000\t3000\t.\t+\t.\tID=rna-A;Parent=gene-A\n"
    "NC_000001.1\tGnomon\texon\t1000\t3000\t.\t+\t.\tID=exon-A-1;Parent=rna-A\n"
    "NC_000001.1\tGnomon\tCDS\t1100\t2900\t.\t+\t0\tID=cds-A;Parent=rna-A\n"
)


# ---------------------------------------------------------------------------
# AUDIT-1  the parent-edge tier is written and tested but never run
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-1: check_parent_edges is never called; "
                   "verify prints three tiers where the README documents four")
def test_verify_reports_the_parent_edge_tier(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "refseq_fgraminearum_NC_026477.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    result = _cli(tmp_path, "verify", "--assembly-id", "1")
    assert "parent edges" in result.output


@pytest.mark.xfail(strict=True, reason="AUDIT-1: a database whose hierarchy is gone "
                   "still verifies PASS (the D-074 failure, still reachable)")
def test_verify_fails_when_every_parent_link_is_lost(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "refseq_fgraminearum_NC_026477.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    connection = sqlite3.connect(tmp_path / "g.db")
    connection.execute("UPDATE feature SET parent_id = NULL")
    connection.commit()
    connection.close()
    result = _cli(tmp_path, "verify", "--assembly-id", "1")
    assert result.exit_code != 0, result.output


# ---------------------------------------------------------------------------
# AUDIT-2  AUGUSTUS rows without an ID lose the parent the file states
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-2: _resolve_parent_ids skips rows with no "
                   "ID; 175 of 322 AUGUSTUS features stored unlinked yet marked "
                   "explicit_parent")
def test_augustus_rows_without_an_id_keep_their_stated_parent(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    unlinked = _rows(tmp_path, """
        SELECT count(*) FROM feature
        WHERE parent_id IS NULL AND linkage_method = 'explicit_parent'""")[0][0]
    assert unlinked == 0, f"{unlinked} rows claim an explicit parent but have none"


# ---------------------------------------------------------------------------
# AUDIT-3  export escapes the commas that separate multiple values
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-3: export writes Dbxref=a%2Cb, which GFF3 "
                   "reads as one value; a round trip loses attribute rows")
def test_export_keeps_multi_valued_attributes_multi_valued(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "refseq_mycmay_GCF_000328475.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    out = tmp_path / "export.gff3"
    assert _cli(tmp_path, "export", "--assembly-id", "1", "--output", str(out)).exit_code == 0
    text = out.read_text()
    assert "Dbxref=GeneID:23561424,GenBank:XM_011387688.1" in text
    assert "%2CGenBank:" not in text


# ---------------------------------------------------------------------------
# AUDIT-4  --force is not crash-safe; verify never checks the data directory
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-4: the superseded directory is moved aside "
                   "when the load starts, so SIGTERM/SIGKILL/OOM mid-load leaves rows "
                   "pointing at nothing")
def test_force_leaves_the_previous_directory_in_place_until_commit(tmp_path):
    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(engine)
    root = tmp_path / "data"
    target = root / "org" / "s1" / "asm1"
    target.mkdir(parents=True)
    (target / "assembly.fa").write_text(">x\nACGT\n")
    with staged_load(engine, root, organism="Org", strain="S1", assembly="asm1", force=True):
        # Any instant here is one at which the process can be killed without
        # Python's exception handling running (SIGTERM, SIGKILL, OOM killer).
        assert (target / "assembly.fa").exists()


@pytest.mark.xfail(strict=True, reason="AUDIT-4: fasta_sha256 is recorded but verify never "
                   "checks the stored assembly still exists")
def test_verify_fails_when_the_stored_assembly_is_missing(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "refseq_mycmay_GCF_000328475.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    for stored in (tmp_path / "data").rglob("assembly.fa"):
        stored.unlink()
    result = _cli(tmp_path, "verify", "--assembly-id", "1")
    assert result.exit_code != 0, result.output


# ---------------------------------------------------------------------------
# AUDIT-5  vertebrate RefSeq mitochondria are refused
# ---------------------------------------------------------------------------

# Verbatim structure from a published RefSeq GFF (chimpanzee GCF_000001515.7,
# ND1): NCBI writes vertebrate mitochondrial CDS directly under the gene.
CHIMP_ND1 = (
    "NC_001643.1\tRefSeq\tregion\t1\t16554\t.\t+\t.\tID=NC_001643.1:1..16554;"
    "Is_circular=true;genome=mitochondrion\n"
    "NC_001643.1\tRefSeq\tgene\t2725\t3681\t.\t+\t.\tID=gene39686;Dbxref=GeneID:807867;"
    "Name=ND1;gbkey=Gene;gene=ND1;gene_biotype=protein_coding;partial=true;"
    "start_range=.,2725\n"
    "NC_001643.1\tRefSeq\tCDS\t2725\t3681\t.\t+\t0\tID=cds80154;Parent=gene39686;"
    "Dbxref=Genbank:NP_008186.1,GeneID:807867;Name=NP_008186.1;gbkey=CDS;gene=ND1;"
    "partial=true;product=NADH dehydrogenase subunit 1;protein_id=NP_008186.1;"
    "start_range=.,2725;transl_table=2\n"
)


@pytest.mark.xfail(strict=True, reason="AUDIT-5: the file-wide transcript-layer rule "
                   "refuses NCBI's standard mitochondrial representation")
def test_vertebrate_refseq_file_with_its_mitochondrion_loads(tmp_path):
    text = REFSEQ_HEADER + "##sequence-region NC_001643.1 1 16554\n" + NUCLEAR_GENE + CHIMP_ND1
    config = _project(tmp_path, text)
    result = _init_and_load(tmp_path, config)
    assert result.exit_code == 0, result.output


# ---------------------------------------------------------------------------
# AUDIT-6  coordinate checks the loader does not make
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-6: features past the end of a linear "
                   "sequence load silently")
def test_feature_past_the_end_of_a_linear_sequence_is_not_silent(tmp_path):
    text = REFSEQ_HEADER + (
        "NC_000001.1\tGnomon\tgene\t19500\t20800\t.\t+\t.\tID=gene-Z;gene_biotype=protein_coding\n"
        "NC_000001.1\tGnomon\tmRNA\t19500\t20800\t.\t+\t.\tID=rna-Z;Parent=gene-Z\n"
        "NC_000001.1\tGnomon\texon\t19500\t20800\t.\t+\t.\tID=exon-Z-1;Parent=rna-Z\n"
        "NC_000001.1\tGnomon\tCDS\t19600\t20700\t.\t+\t0\tID=cds-Z;Parent=rna-Z\n"
    )
    config = _project(tmp_path, text, lengths={"NC_000001.1": 20000})
    result = _init_and_load(tmp_path, config)
    named = any(tag in result.output for tag in ("gene-Z", "rna-Z", "exon-Z", "cds-Z"))
    assert result.exit_code != 0 or named, (
        "a feature ending at 20,800 on a 20,000 bp sequence loaded without comment")


@pytest.mark.xfail(strict=True, reason="AUDIT-6: the envelope check compares numbers, "
                   "not sequences; a child on another sequence is linked anyway")
def test_child_on_a_different_sequence_from_its_parent_is_refused(tmp_path):
    text = REFSEQ_HEADER + "##sequence-region NC_000002.1 1 20000\n" + NUCLEAR_GENE.replace(
        "NC_000001.1\tGnomon\texon", "NC_000002.1\tGnomon\texon")
    config = _project(tmp_path, text)
    result = _init_and_load(tmp_path, config)
    assert result.exit_code != 0, result.output


# ---------------------------------------------------------------------------
# AUDIT-7  a spec-legal ##FASTA section is refused
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-7: ##FASTA ends the feature section per the "
                   "GFF3 spec; the tokenizer reads the sequence lines as rows")
def test_embedded_fasta_section_is_accepted(tmp_path):
    text = REFSEQ_HEADER + NUCLEAR_GENE + "##FASTA\n>NC_000001.1\nACGTACGTACGT\n"
    config = _project(tmp_path, text, lengths={"NC_000001.1": 20000})
    result = _init_and_load(tmp_path, config)
    assert result.exit_code == 0, result.output


# ---------------------------------------------------------------------------
# AUDIT-8  export binds one parameter per row
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-8: load_features puts every row id into one "
                   "IN list; exceeds SQLite's limit (250,000 here, 32,766 or 999 on other "
                   "builds) and PostgreSQL's 65,535")
def test_export_does_not_depend_on_the_bound_parameter_limit(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "refseq_fgraminearum_NC_026477.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")

    @event.listens_for(engine, "connect")
    def _small_limit(dbapi_connection, _record):
        # 400 rows against a limit of 250 stands in for 300,000 rows against
        # 250,000, which is what failed in the audit's scale run.
        dbapi_connection.setlimit(sqlite3.SQLITE_LIMIT_VARIABLE_NUMBER, 250)

    with engine.connect() as connection:
        assert len(load_features(connection, 1)) == 354


# ---------------------------------------------------------------------------
# AUDIT-9  PostgreSQL cannot create the schema; coordinates are 32-bit
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-9: CheckConstraint('end >= start') is raw SQL; "
                   "END is reserved in PostgreSQL, so `cgload init` fails there")
def test_postgresql_ddl_quotes_the_end_column_in_its_check_constraint():
    ddl = str(CreateTable(feature).compile(dialect=postgresql.dialect()))
    assert "CHECK (end >= start)" not in ddl


@pytest.mark.xfail(strict=True, reason="AUDIT-9: INTEGER tops out at 2,147,483,647; "
                   "PostgreSQL rejects a 2.2 Gb sequence with 'integer out of range'")
def test_coordinates_hold_sequences_longer_than_2_1_gb():
    for column in (feature.c.start, feature.c.end, sequence_region.c.length):
        assert isinstance(column.type, BigInteger), column


# ---------------------------------------------------------------------------
# AUDIT-10  README and help-text promises
# ---------------------------------------------------------------------------


@pytest.mark.xfail(strict=True, reason="AUDIT-10: README promises a bgzip'd, tabix-indexed "
                   "annotation in the data directory for JBrowse; nothing writes it")
def test_data_directory_holds_the_indexed_annotation_the_readme_promises(tmp_path):
    for name in ("assembly.fa", "annotation.gff3", "organism.toml"):
        (tmp_path / name).write_text((EXAMPLE / name).read_text())
    assert _init_and_load(tmp_path, tmp_path / "organism.toml").exit_code == 0
    target = next(p for p in (tmp_path / "data").rglob("manifest.json")).parent
    assert (target / "annotation.gff3.gz").exists()
    assert (target / "annotation.gff3.gz.tbi").exists()


@pytest.mark.xfail(strict=True, reason="AUDIT-10: --help says the reference tier is skipped "
                   "without --sequences-in-assembly; the command errors instead")
def test_verify_help_describes_what_happens_without_sequences_in_assembly(tmp_path):
    config = _project(tmp_path, fixture=FIXTURES / "refseq_mycmay_GCF_000328475.gff3")
    assert _init_and_load(tmp_path, config).exit_code == 0
    helptext = " ".join(_cli(tmp_path, "verify", "--help").output.split())
    result = _cli(tmp_path, "verify", "--assembly-id", "1", "--reference-counts",
                  str(FIXTURES / "ncbi_feature_count_GCF_000328475.txt"),
                  "--accession", "GCF_000328475.2")
    says_skipped = "and is skipped" in helptext
    assert not (says_skipped and result.exit_code != 0), (
        "help promises a skip; the command exits with an error")
