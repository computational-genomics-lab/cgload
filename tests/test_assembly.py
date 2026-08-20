"""Milestone 2 tests.

What these assert, stated before they were written:

* A good assembly produces a copied FASTA, a .fai, a manifest and one row per
  sequence region.
* A failure anywhere in the load leaves *both* the database and the data
  directory exactly as they were. This is the whole point of D-006 and is the
  test that matters most in this file.
* Malformed assemblies error rather than loading partially: empty file,
  duplicate sequence name, zero-length record.
* Paths stored in the database are relative, so the database stays portable.
* The config loader rejects what it cannot safely guess at.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner
from cgload.assembly import AssemblyExistsError, register_assembly
from cgload.cli.main import cli
from cgload.config import ConfigError, load_config
from cgload.db.engine import build_engine
from cgload.db.init import create_schema
from cgload.db.schema import assembly, organism, sequence_region, source_file, strain
from cgload.fasta import AssemblyReadError, index_and_read
from cgload.storage import DataDirectoryConflictError, slug, staged_load
from sqlalchemy import func, select

REGIONS = [("scaffold_1", 400), ("scaffold_2", 150), ("mito", 60)]


def write_fasta(path: Path, regions=REGIONS) -> Path:
    with path.open("w") as fh:
        for name, n in regions:
            fh.write(f">{name} some description\n")
            seq = "ACGT" * (n // 4 + 1)
            seq = seq[:n]
            for i in range(0, len(seq), 60):
                fh.write(seq[i : i + 60] + "\n")
    return path


def write_config(tmp_path: Path, fasta: Path, **overrides) -> Path:
    body = f"""
[organism]
genus = "Phytophthora"
species = "megasperma"
ncbi_taxon_id = 4788
strain = "{overrides.get('strain', 'CJ26')}"

[assembly]
name = "{overrides.get('name', 'Phyme_CJ26')}"
version = "1.0"
fasta = "{fasta.name}"
"""
    path = tmp_path / "organism.toml"
    path.write_text(body)
    return path


@pytest.fixture
def engine(tmp_path):
    eng = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(eng)
    return eng


@pytest.fixture
def config(tmp_path):
    return load_config(write_config(tmp_path, write_fasta(tmp_path / "tiny.fa")))


# -- the happy path --------------------------------------------------------


def test_assembly_lands_with_index_and_manifest(engine, config, tmp_path):
    result = register_assembly(engine, config, tmp_path / "data")

    assert result.region_count == len(REGIONS)
    assert result.total_length == sum(n for _, n in REGIONS)

    directory = result.data_directory
    assert (directory / "assembly.fa").exists()
    assert (directory / "assembly.fa.fai").exists()

    manifest = json.loads((directory / "manifest.json").read_text())
    assert manifest["regions"] == len(REGIONS)
    assert manifest["files"][0]["mode"] == "copy"
    assert manifest["vocabulary_hash"]


def test_regions_are_stored_with_their_lengths(engine, config, tmp_path):
    register_assembly(engine, config, tmp_path / "data")
    with engine.connect() as conn:
        stored = dict(
            conn.execute(
                select(sequence_region.c.source_id, sequence_region.c.length)
            ).all()
        )
    assert stored == dict(REGIONS)


def test_stored_fasta_path_is_relative(engine, config, tmp_path):
    """An absolute path makes the database non-portable, defeating D-002."""
    register_assembly(engine, config, tmp_path / "data")
    with engine.connect() as conn:
        stored = conn.execute(select(assembly.c.fasta_path)).scalar_one()
    assert not Path(stored).is_absolute()
    assert stored == "phytophthora_megasperma/cj26/phyme_cj26/assembly.fa"


def test_original_file_is_untouched(engine, config, tmp_path):
    before = config.fasta.read_bytes()
    register_assembly(engine, config, tmp_path / "data")
    assert config.fasta.read_bytes() == before


def test_second_assembly_reuses_organism_and_strain(engine, config, tmp_path):
    register_assembly(engine, config, tmp_path / "data")
    second = load_config(
        write_config(tmp_path, config.fasta, name="Phyme_CJ26_v2")
    )
    register_assembly(engine, second, tmp_path / "data")
    with engine.connect() as conn:
        assert conn.execute(select(func.count()).select_from(organism)).scalar_one() == 1
        assert conn.execute(select(func.count()).select_from(strain)).scalar_one() == 1
        assert conn.execute(select(func.count()).select_from(assembly)).scalar_one() == 2


def test_reloading_the_same_assembly_is_refused(engine, config, tmp_path):
    register_assembly(engine, config, tmp_path / "data")
    with pytest.raises((AssemblyExistsError, DataDirectoryConflictError)):
        register_assembly(engine, config, tmp_path / "data")


# -- all or nothing (D-006) ------------------------------------------------


def test_failure_leaves_database_and_data_directory_untouched(engine, config, tmp_path,
                                                              monkeypatch):
    """The test this milestone exists for."""
    data_root = tmp_path / "data"

    import cgload.assembly as assembly_module

    def explode(*_args, **_kwargs):
        raise RuntimeError("simulated failure after files were staged")

    monkeypatch.setattr(assembly_module, "index_and_read", explode)

    with pytest.raises(RuntimeError, match="simulated"):
        register_assembly(engine, config, data_root)

    with engine.connect() as conn:
        for table in (organism, strain, assembly, sequence_region, source_file):
            assert conn.execute(select(func.count()).select_from(table)).scalar_one() == 0

    leftovers = [p for p in data_root.rglob("*") if p.is_file()]
    assert leftovers == []


def test_failed_load_leaves_no_staging_directory(engine, config, tmp_path, monkeypatch):
    data_root = tmp_path / "data"
    import cgload.assembly as assembly_module

    monkeypatch.setattr(
        assembly_module, "index_and_read", lambda *_a, **_k: (_ for _ in ()).throw(OSError("x"))
    )
    with pytest.raises(OSError):
        register_assembly(engine, config, data_root)

    staging = data_root / ".staging"
    assert not staging.exists() or not any(staging.iterdir())


def test_commit_failure_puts_the_directory_back(engine, config, tmp_path):
    """Rename happens inside the transaction; if the commit fails the
    directory must not be left in place with no rows behind it."""
    data_root = tmp_path / "data"
    target = data_root / "phytophthora_megasperma" / "cj26" / "phyme_cj26"

    with pytest.raises(RuntimeError, match="commit failed"), staged_load(
        engine,
        data_root,
        organism="Phytophthora_megasperma",
        strain="CJ26",
        assembly="Phyme_CJ26",
    ) as load:
        (load.staging / "assembly.fa").write_text(">x\nACGT\n")
        raise RuntimeError("commit failed")

    assert not target.exists()


def test_partial_directory_is_not_visible_midway(engine, config, tmp_path):
    """A reader must never see a half-built assembly directory: work happens
    under .staging and appears only at the final rename."""
    data_root = tmp_path / "data"
    seen: list[bool] = []

    with staged_load(
        engine,
        data_root,
        organism="Phytophthora_megasperma",
        strain="CJ26",
        assembly="Phyme_CJ26",
    ) as load:
        (load.staging / "assembly.fa").write_text(">x\nACGT\n")
        seen.append(load.target.exists())

    assert seen == [False]
    assert (data_root / "phytophthora_megasperma" / "cj26" / "phyme_cj26").exists()


# -- malformed assemblies error rather than loading partially --------------


def test_empty_fasta_is_an_error(tmp_path):
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    # pyfaidx rejects it at indexing time; either path must raise ours, not
    # a bare pyfaidx exception leaking through.
    with pytest.raises(AssemblyReadError):
        index_and_read(empty)


def test_duplicate_sequence_name_is_an_error(tmp_path):
    path = write_fasta(tmp_path / "dupe.fa", [("scaffold_1", 40), ("scaffold_1", 40)])
    with pytest.raises(AssemblyReadError):
        index_and_read(path)


def test_zero_length_record_is_an_error(tmp_path):
    path = tmp_path / "zero.fa"
    path.write_text(">scaffold_1\nACGT\n>scaffold_2\n")
    with pytest.raises(AssemblyReadError, match="zero-length"):
        index_and_read(path)


def test_missing_assembly_file_is_an_error(engine, tmp_path):
    config_path = write_config(tmp_path, tmp_path / "absent.fa")
    with pytest.raises(FileNotFoundError):
        register_assembly(engine, load_config(config_path), tmp_path / "data")


# -- directory naming ------------------------------------------------------


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("Phytophthora megasperma", "phytophthora_megasperma"),
        ("CJ26", "cj26"),
        ("strain/with/slashes", "strain_with_slashes"),
        ("  padded  ", "padded"),
    ],
)
def test_slug_is_filesystem_safe(raw, expected):
    assert slug(raw) == expected


def test_slug_rejects_a_name_with_no_safe_form():
    with pytest.raises(ValueError):
        slug("///")


# -- config ----------------------------------------------------------------


def test_strain_defaults_to_reference(tmp_path):
    fasta = write_fasta(tmp_path / "t.fa")
    (tmp_path / "c.toml").write_text(
        f'[organism]\ngenus="G"\nspecies="s"\n[assembly]\nname="A"\nfasta="{fasta.name}"\n'
    )
    assert load_config(tmp_path / "c.toml").strain == "reference"


def test_missing_required_key_names_it(tmp_path):
    (tmp_path / "c.toml").write_text('[organism]\ngenus="G"\n[assembly]\nname="A"\nfasta="x"\n')
    with pytest.raises(ConfigError, match="species"):
        load_config(tmp_path / "c.toml")


def test_non_integer_taxon_id_is_refused(tmp_path):
    (tmp_path / "c.toml").write_text(
        '[organism]\ngenus="G"\nspecies="s"\nncbi_taxon_id="Phytophthora"\n'
        '[assembly]\nname="A"\nfasta="x"\n'
    )
    with pytest.raises(ConfigError, match="integer"):
        load_config(tmp_path / "c.toml")


def test_annotation_files_are_read_as_a_list(tmp_path):
    fasta = write_fasta(tmp_path / "t.fa")
    (tmp_path / "c.toml").write_text(
        f'[organism]\ngenus="G"\nspecies="s"\n[assembly]\nname="A"\nfasta="{fasta.name}"\n'
        '[[files.annotation]]\npath="a.gff3"\n'
        '[[files.annotation]]\npath="b.gff3"\nid_prefix="FUN_"\n'
    )
    cfg = load_config(tmp_path / "c.toml")
    assert len(cfg.annotation) == 2
    assert cfg.annotation[1].id_prefix == "FUN_"


# -- stats -----------------------------------------------------------------


def test_stats_on_an_empty_database(engine, tmp_path):
    result = CliRunner().invoke(cli, ["stats", "--db-url", str(engine.url)])
    assert result.exit_code == 0
    assert "No assemblies loaded" in result.output


def test_stats_reports_the_loaded_assembly(engine, config, tmp_path):
    register_assembly(engine, config, tmp_path / "data")
    result = CliRunner().invoke(cli, ["stats", "--db-url", str(engine.url)])
    assert result.exit_code == 0
    assert "Phytophthora megasperma" in result.output
    assert "3 sequence regions" in result.output
    assert "610 bp" in result.output


def test_stats_is_no_longer_a_stub():
    from cgload.cli.main import PENDING

    assert "stats" not in PENDING
    assert "init" not in PENDING


# -- the wrong-directory failure (regression) ------------------------------


def test_loading_into_an_uninitialised_database_says_run_init(tmp_path):
    """Reported from the field: running the load snippet in a directory where
    `cgload init` had not been run produced a SQLAlchemy missing-table
    traceback. It must name the fix instead."""
    from cgload.db.engine import build_engine as raw_engine
    from cgload.db.init import SchemaMissingError

    blank = raw_engine(f"sqlite:///{tmp_path / 'blank.db'}")
    fasta = write_fasta(tmp_path / "t.fa")
    config = load_config(write_config(tmp_path, fasta))

    with pytest.raises(SchemaMissingError, match="cgload init"):
        register_assembly(blank, config, tmp_path / "data")


def test_missing_sqlite_file_is_refused_rather_than_created(tmp_path):
    """SQLite creates a file on connect, so a mistyped path silently yields an
    empty database. Every command except init must refuse it."""
    from cgload.db.engine import DatabaseNotFoundError
    from cgload.db.engine import build_engine as raw_engine

    absent = tmp_path / "typo.db"
    with pytest.raises(DatabaseNotFoundError, match="cgload init"):
        raw_engine(f"sqlite:///{absent}", must_exist=True)
    assert not absent.exists()


def test_stats_on_a_missing_database_exits_two(tmp_path):
    result = CliRunner().invoke(cli, ["stats", "--db-url", f"sqlite:///{tmp_path / 'nope.db'}"])
    assert result.exit_code == 2
    assert "cgload init" in result.output


# -- accession (keys the optional NCBI cross-check at stage 8) --------------


def test_accession_is_stored_and_manifested(engine, tmp_path):
    fasta = write_fasta(tmp_path / "t.fa")
    (tmp_path / "organism.toml").write_text(
        f'[organism]\ngenus="Phytophthora"\nspecies="megasperma"\nstrain="CJ26"\n'
        f'[assembly]\nname="Phyme_CJ26"\naccession="GCF_000149735.1"\nfasta="{fasta.name}"\n'
    )
    config = load_config(tmp_path / "organism.toml")
    assert config.accession == "GCF_000149735.1"

    result = register_assembly(engine, config, tmp_path / "data")
    with engine.connect() as conn:
        assert conn.execute(select(assembly.c.accession)).scalar_one() == "GCF_000149735.1"
    manifest = json.loads((result.data_directory / "manifest.json").read_text())
    assert manifest["accession"] == "GCF_000149735.1"


def test_accession_is_optional(engine, config, tmp_path):
    register_assembly(engine, config, tmp_path / "data")
    with engine.connect() as conn:
        assert conn.execute(select(assembly.c.accession)).scalar_one() is None
