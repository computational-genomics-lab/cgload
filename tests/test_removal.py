"""Replacing and removing assemblies (D-083).

What these assert, stated before they were written:

* ``load --force`` replaces a loaded assembly. It never did: it deleted only the
  assembly row, and the enforced foreign keys refused that whenever the
  assembly had features, so every replacement stopped with a traceback.
"""

from __future__ import annotations

import random
import sqlite3
from pathlib import Path

import pytest
from cgload.assembly import register_assembly
from cgload.cli.main import cli
from cgload.config import load_config
from cgload.db.engine import build_engine
from click.testing import CliRunner

FIXTURES = Path(__file__).parent / "fixtures"
FUSARIUM = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
MYCMAY = FIXTURES / "refseq_mycmay_GCF_000328475.gff3"
EGGNOG = FIXTURES / "eggnog_fgraminearum.emapper.annotations"

OWNED_TABLES = ("feature", "sequence_region", "source_file", "annotation_run")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _lengths(gff: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line in gff.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            lengths[parts[1]] = max(lengths.get(parts[1], 0), int(parts[3]))
        elif line and not line.startswith("#"):
            columns = line.split("\t")
            lengths[columns[0]] = max(lengths.get(columns[0], 0), int(columns[4]) + 500)
    return lengths


def _fasta(path: Path, lengths: dict[str, int]) -> None:
    block = "".join(random.Random(7).choice("ACGT") for _ in range(6000))
    with path.open("w") as handle:
        for name, length in lengths.items():
            sequence = (block * (length // len(block) + 1))[:length]
            handle.write(f">{name}\n")
            for index in range(0, length, 60):
                handle.write(sequence[index : index + 60] + "\n")


def _project(
    directory: Path,
    gff: Path,
    *,
    genus: str = "Fusarium",
    species: str = "graminearum",
    strain: str = "PH-1",
    assembly: str = "asm1",
    functional: Path | None = None,
) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / gff.name).write_text(gff.read_text())
    _fasta(directory / "assembly.fa", _lengths(gff))
    text = (
        f'[organism]\ngenus = "{genus}"\nspecies = "{species}"\nstrain = "{strain}"\n\n'
        f'[assembly]\nname = "{assembly}"\nfasta = "assembly.fa"\n\n'
        f'[[files.annotation]]\npath = "{gff.name}"\n'
    )
    if functional is not None:
        (directory / functional.name).write_text(functional.read_text())
        text += f'\n[[files.functional]]\npath = "{functional.name}"\n'
    config = directory / "organism.toml"
    config.write_text(text)
    return config


def _cli(db: Path, *args: str, input: str | None = None):
    return CliRunner().invoke(cli, [args[0], "--db-url", f"sqlite:///{db}", *args[1:]],
                              input=input)


def _count(db: Path, sql: str) -> int:
    connection = sqlite3.connect(db)
    try:
        return connection.execute(sql).fetchone()[0]
    finally:
        connection.close()


@pytest.fixture
def fusarium(tmp_path):
    """One assembly with functional annotation, loaded through the CLI."""
    db = tmp_path / "g.db"
    assert _cli(db, "init").exit_code == 0
    config = _project(tmp_path / "fusarium", FUSARIUM, functional=EGGNOG)
    loaded = _cli(db, "load", "--config", str(config), "--data-dir", str(tmp_path / "data"))
    assert loaded.exit_code == 0, loaded.output
    directory = next((tmp_path / "data").rglob("manifest.json")).parent
    return db, config, directory


# ---------------------------------------------------------------------------
# load --force
# ---------------------------------------------------------------------------


def test_force_replaces_a_loaded_assembly(fusarium, tmp_path):
    db, config, _ = fusarium
    features = _count(db, "SELECT count(*) FROM feature")
    annotations = _count(db, "SELECT count(*) FROM protein_annotation")
    assert annotations > 0

    again = _cli(db, "load", "--config", str(config), "--data-dir", str(tmp_path / "data"),
                 "--force")
    assert again.exit_code == 0, again.output
    assert _count(db, "SELECT count(*) FROM assembly") == 1
    assert _count(db, "SELECT count(*) FROM feature") == features
    assert _count(db, "SELECT count(*) FROM protein_annotation") == annotations
    for table in OWNED_TABLES:
        orphans = _count(db, f"SELECT count(*) FROM {table} "
                             "WHERE assembly_id NOT IN (SELECT assembly_id FROM assembly)")
        assert orphans == 0, table
    assembly_id = _count(db, "SELECT assembly_id FROM assembly")
    verified = _cli(db, "verify", "--assembly-id", str(assembly_id))
    assert verified.exit_code == 0, verified.output


def test_register_assembly_with_force_replaces_the_row(tmp_path):
    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    assert _cli(tmp_path / "g.db", "init").exit_code == 0
    config = load_config(_project(tmp_path / "p", MYCMAY, genus="Pgtest", species="exampleii"))
    register_assembly(engine, config, tmp_path / "data")
    register_assembly(engine, config, tmp_path / "data", force=True)
    assert _count(tmp_path / "g.db", "SELECT count(*) FROM assembly") == 1
    assert _count(tmp_path / "g.db", "SELECT count(*) FROM sequence_region") == 2
