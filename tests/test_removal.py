"""Replacing and removing assemblies (D-083, D-084).

What these assert, stated before they were written:

* ``load --force`` replaces a loaded assembly. It never did: it deleted only the
  assembly row, and the enforced foreign keys refused that whenever the
  assembly had features, so every replacement stopped with a traceback.
* ``remove`` deletes every row an assembly owns and its data directory, in one
  transaction, and leaves every other assembly exactly as it was.
* It deletes no file it cannot show belongs to the assembly, and it changes
  nothing when it stops.
* A failure at the last step puts the data directory back.
"""

from __future__ import annotations

import hashlib
import json
import random
import shutil
import sqlite3
from pathlib import Path

import pytest
from cgload.assembly import register_assembly
from cgload.cli.main import cli
from cgload.config import load_config
from cgload.db.engine import build_engine
from cgload.remove import remove_assembly
from click.testing import CliRunner
from sqlalchemy import event

FIXTURES = Path(__file__).parent / "fixtures"
FUSARIUM = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
MYCMAY = FIXTURES / "refseq_mycmay_GCF_000328475.gff3"
EGGNOG = FIXTURES / "eggnog_fgraminearum.emapper.annotations"

OWNED_TABLES = ("feature", "sequence_region", "source_file", "annotation_run")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


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


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


# ---------------------------------------------------------------------------
# remove
# ---------------------------------------------------------------------------


def test_a_dry_run_changes_nothing(fusarium):
    db, _, directory = fusarium
    before = _digest(db)
    result = _cli(db, "remove", "--assembly-id", "1", "--dry-run")
    assert result.exit_code == 0, result.output
    assert "Dry run: nothing was removed." in result.output
    assert "354 features" in result.output
    assert _digest(db) == before
    assert directory.is_dir()


def test_remove_takes_rows_and_files_and_leaves_other_assemblies(fusarium, tmp_path):
    db, _, directory = fusarium
    other = _project(tmp_path / "other", MYCMAY, genus="Pgtest", species="exampleii",
                     strain="S1")
    assert _cli(db, "load", "--config", str(other), "--data-dir",
                str(tmp_path / "data")).exit_code == 0
    other_features = _count(db, "SELECT count(*) FROM feature WHERE assembly_id = 2")

    result = _cli(db, "remove", "--assembly-id", "1", "--yes")
    assert result.exit_code == 0, result.output

    for table in (*OWNED_TABLES, "assembly"):
        assert _count(db, f"SELECT count(*) FROM {table} WHERE assembly_id = 1") == 0, table
    assert _count(db, "SELECT count(*) FROM protein_annotation") == 0
    assert _count(db, "SELECT count(*) FROM feature_attribute WHERE feature_id NOT IN "
                      "(SELECT feature_id FROM feature)") == 0
    assert _count(db, "SELECT count(*) FROM organism WHERE genus = 'Fusarium'") == 0
    assert _count(db, "SELECT count(*) FROM strain WHERE name = 'PH-1'") == 0
    assert not directory.exists()
    assert not directory.parent.parent.exists(), "empty organism directory left behind"
    assert not any((tmp_path / "data" / ".staging").iterdir())

    assert _count(db, "SELECT count(*) FROM feature WHERE assembly_id = 2") == other_features
    verified = _cli(db, "verify", "--assembly-id", "2")
    assert verified.exit_code == 0, verified.output


def test_remove_asks_first(fusarium):
    db, _, directory = fusarium
    before = _digest(db)
    result = _cli(db, "remove", "--assembly-id", "1", input="n\n")
    assert result.exit_code == 1
    assert _digest(db) == before
    assert directory.is_dir()


def test_remove_stops_when_the_files_are_not_where_recorded(fusarium, tmp_path):
    db, _, directory = fusarium
    moved = tmp_path / "moved"
    shutil.move(str(tmp_path / "data"), moved)
    before = _digest(db)

    refused = _cli(db, "remove", "--assembly-id", "1", "--yes")
    assert refused.exit_code == 2
    assert "--data-dir" in refused.output
    assert _digest(db) == before

    done = _cli(db, "remove", "--assembly-id", "1", "--yes", "--data-dir", str(moved))
    assert done.exit_code == 0, done.output
    assert not (moved / directory.relative_to(tmp_path / "data")).exists()
    assert _count(db, "SELECT count(*) FROM assembly") == 0


def test_keep_files_removes_only_the_rows(fusarium):
    db, _, directory = fusarium
    result = _cli(db, "remove", "--assembly-id", "1", "--yes", "--keep-files")
    assert result.exit_code == 0, result.output
    assert _count(db, "SELECT count(*) FROM feature") == 0
    assert (directory / "assembly.fa").is_file()


def test_remove_refuses_a_directory_whose_manifest_names_another_assembly(fusarium):
    db, _, directory = fusarium
    manifest = directory / "manifest.json"
    content = json.loads(manifest.read_text())
    content["assembly"] = "something-else"
    manifest.write_text(json.dumps(content))
    before = _digest(db)

    result = _cli(db, "remove", "--assembly-id", "1", "--yes")
    assert result.exit_code == 2
    assert "does not describe" in result.output
    assert _digest(db) == before
    assert directory.is_dir()


def test_a_failed_commit_puts_the_directory_back(fusarium):
    db, _, directory = fusarium
    engine = build_engine(f"sqlite:///{db}")

    @event.listens_for(engine, "commit")
    def _fail(_connection):
        raise RuntimeError("simulated failure at commit")

    with pytest.raises(RuntimeError, match="simulated"):
        remove_assembly(engine, 1)
    assert (directory / "assembly.fa").is_file()
    assert _count(db, "SELECT count(*) FROM assembly") == 1
    assert _count(db, "SELECT count(*) FROM feature") > 0


def test_remove_keeps_a_strain_that_another_assembly_uses(fusarium, tmp_path):
    db, _, directory = fusarium
    second = _project(tmp_path / "second", MYCMAY, assembly="asm2")
    assert _cli(db, "load", "--config", str(second), "--data-dir",
                str(tmp_path / "data")).exit_code == 0

    assert _cli(db, "remove", "--assembly-id", "1", "--yes").exit_code == 0
    assert _count(db, "SELECT count(*) FROM strain WHERE name = 'PH-1'") == 1
    assert _count(db, "SELECT count(*) FROM organism WHERE genus = 'Fusarium'") == 1
    assert not directory.exists()
    assert (directory.parent / "asm2" / "manifest.json").is_file()


def test_removing_an_unknown_assembly_is_an_error(fusarium):
    db, _, _ = fusarium
    result = _cli(db, "remove", "--assembly-id", "9", "--yes")
    assert result.exit_code == 2
    assert "no assembly with id 9" in result.output
