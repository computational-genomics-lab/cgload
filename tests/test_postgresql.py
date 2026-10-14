"""The PostgreSQL code path, run against a live server.

The README says MySQL and PostgreSQL work through the same code path as SQLite,
but until this file no test touched a server, and two commands failed on
PostgreSQL 16 the first time one was tried: `init` (AUDIT-9: the CHECK
constraint used the reserved word END) and `stats` (it grouped by the assembly
key alone, which SQLite accepts and PostgreSQL refuses).

Skipped unless CGLOAD_TEST_POSTGRES_URL names a server URL without a password,
for example ``postgresql+psycopg://postgres@127.0.0.1:5432/postgres``. Each
test creates its own database and drops it afterwards. CI runs this file
against a PostgreSQL service container.
"""

from __future__ import annotations

import os
import random
import uuid
from pathlib import Path

import pytest
from cgload.cli.main import cli
from click.testing import CliRunner
from sqlalchemy import create_engine, text
from sqlalchemy.engine import make_url

pytestmark = pytest.mark.requires_server

ADMIN_URL = os.environ.get("CGLOAD_TEST_POSTGRES_URL")
FIXTURES = Path(__file__).parent / "fixtures"
MYCMAY_GFF = FIXTURES / "refseq_mycmay_GCF_000328475.gff3"


@pytest.fixture
def pg_url():
    if not ADMIN_URL:
        pytest.skip("set CGLOAD_TEST_POSTGRES_URL to run the PostgreSQL tests")
    name = f"cgload_test_{uuid.uuid4().hex[:12]}"
    admin = create_engine(ADMIN_URL, isolation_level="AUTOCOMMIT")
    with admin.connect() as connection:
        connection.execute(text(f'CREATE DATABASE "{name}"'))
    try:
        yield make_url(ADMIN_URL).set(database=name).render_as_string(hide_password=False)
    finally:
        with admin.connect() as connection:
            connection.execute(text(f'DROP DATABASE IF EXISTS "{name}" WITH (FORCE)'))
        admin.dispose()


def _project(tmp_path: Path) -> Path:
    gff = tmp_path / MYCMAY_GFF.name
    gff.write_text(MYCMAY_GFF.read_text())
    lengths: dict[str, int] = {}
    for line in gff.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            lengths[parts[1]] = int(parts[3])
        elif line and not line.startswith("#"):
            columns = line.split("\t")
            ends = lengths.get(columns[0], 0)
            lengths[columns[0]] = max(ends, int(columns[4]) + 500)
    rnd = random.Random(5)
    with (tmp_path / "assembly.fa").open("w") as handle:
        for name, length in lengths.items():
            handle.write(f">{name}\n")
            sequence = "".join(rnd.choice("ACGT") for _ in range(length))
            for index in range(0, length, 60):
                handle.write(sequence[index : index + 60] + "\n")
    config = tmp_path / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Pgtest"\nspecies = "exampleii"\nstrain = "S1"\n\n'
        '[assembly]\nname = "asm1"\nfasta = "assembly.fa"\n\n'
        f'[[files.annotation]]\npath = "{gff.name}"\n'
    )
    return config


def _run(*args: str):
    return CliRunner().invoke(cli, list(args))


def test_init_load_verify_export_and_stats(pg_url, tmp_path):
    assert _run("init", "--db-url", pg_url).exit_code == 0
    loaded = _run("load", "--db-url", pg_url, "--config", str(_project(tmp_path)),
                  "--data-dir", str(tmp_path / "data"))
    assert loaded.exit_code == 0, loaded.output

    verified = _run("verify", "--db-url", pg_url, "--assembly-id", "1")
    assert verified.exit_code == 0, verified.output
    assert "source count" in verified.output

    out = tmp_path / "export.gff3"
    exported = _run("export", "--db-url", pg_url, "--assembly-id", "1", "--output", str(out))
    assert exported.exit_code == 0, exported.output
    assert "Dbxref=GeneID:23561424,GenBank:XM_011387688.1" in out.read_text()

    stats = _run("stats", "--db-url", pg_url)
    assert stats.exit_code == 0, stats.output
    assert "Pgtest exampleii" in stats.output


def test_coordinates_past_2_1_gb_are_stored(pg_url):
    """AUDIT-9: INTEGER stops at 2,147,483,647, and PostgreSQL enforces it."""
    assert _run("init", "--db-url", pg_url).exit_code == 0
    engine = create_engine(pg_url)
    try:
        with engine.begin() as connection:
            for statement in (
                "INSERT INTO organism (genus, species) VALUES ('Big', 'genome')",
                "INSERT INTO strain (organism_id, name) VALUES (1, 'S')",
                "INSERT INTO assembly (strain_id, name, vocabulary_version, "
                "vocabulary_hash, cgload_version) VALUES (1, 'a', 1, 'h', 'v')",
                "INSERT INTO sequence_region (assembly_id, source_id, length) "
                "VALUES (1, 'chr1', 2600000000)",
                "INSERT INTO source_file (assembly_id, path, sha256, file_format) "
                "VALUES (1, 'x.gff3', 'h', 'gff3')",
            ):
                connection.execute(text(statement))
            connection.execute(text(
                'INSERT INTO feature (assembly_id, sequence_region_id, source_file_id, '
                'source_id, source_type, feature_type, start, "end", strand) '
                "VALUES (1, 1, 1, 'g', 'gene', 'gene', 2500000000, 2500001000, '+')"
            ))
            stored = connection.execute(text('SELECT "end" FROM feature')).scalar_one()
        assert stored == 2500001000
    finally:
        engine.dispose()
