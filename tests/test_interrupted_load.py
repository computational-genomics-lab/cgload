"""D-095: a load killed between putting its directory in place and committing.

What these assert, stated before they were written:

* A load killed in that window leaves a directory and no rows. The next load
  of the same assembly into the same database replaces the directory, says
  so, and passes `verify`; before D-095 it was refused.
* A directory that another database's load put there, or one whose manifest
  names no database, still stops the load: replacing it could break an
  assembly another database uses.
* Each database gets its own identity at `init`, and every assembly
  directory's manifest records it.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from cgload.cli.main import cli
from cgload.db.engine import build_engine
from cgload.db.schema import assembly, vocabulary_meta
from click.testing import CliRunner
from sqlalchemy import func, select

GFF = (
    "##gff-version 3\n"
    "##sequence-region chr1 1 400\n"
    "chr1\tRefSeq\tgene\t10\t300\t.\t+\t.\tID=gene-a;Name=a\n"
    "chr1\tRefSeq\tmRNA\t10\t300\t.\t+\t.\tID=rna-a;Parent=gene-a\n"
    "chr1\tRefSeq\texon\t10\t300\t.\t+\t.\tID=exon-a-1;Parent=rna-a\n"
)
TARGET = Path("testus_interruptus") / "s1" / "asm1"


def make_inputs(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "a.gff3").write_text(GFF)
    sequence = ("ACGT" * 100)[:400]
    (directory / "a.fa").write_text(
        ">chr1\n" + "\n".join(sequence[i : i + 60] for i in range(0, 400, 60)) + "\n"
    )
    config = directory / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Testus"\nspecies = "interruptus"\nstrain = "S1"\n\n'
        '[assembly]\nname = "asm1"\nfasta = "a.fa"\n\n'
        '[[files.annotation]]\npath = "a.gff3"\n'
    )
    return config


def run(*args: str):
    return CliRunner().invoke(cli, list(args))


def init(db: Path) -> str:
    url = f"sqlite:///{db}"
    assert run("init", "--db-url", url).exit_code == 0
    return url


def load(url: str, config: Path, data: Path):
    return run("load", "--db-url", url, "--config", str(config), "--data-dir", str(data))


def assemblies(url: str) -> int:
    with build_engine(url).connect() as connection:
        return connection.execute(select(func.count()).select_from(assembly)).scalar_one()


KILLED_LOAD = """
import os, sys
from pathlib import Path

target = Path(sys.argv[1])
real_rename = os.rename


def rename(source, destination):
    real_rename(source, destination)
    if Path(destination) == target:
        # Killed with the directory in place and the transaction still open:
        # no exception handler runs, as with SIGKILL or the OOM killer.
        os._exit(9)


os.rename = rename
from cgload.cli.main import cli

cli(sys.argv[2:])
"""


def test_a_load_killed_after_its_directory_landed_does_not_block_the_next(tmp_path):
    config = make_inputs(tmp_path / "in")
    url = init(tmp_path / "g.db")
    data = tmp_path / "data"
    script = tmp_path / "killed_load.py"
    script.write_text(KILLED_LOAD)
    killed = subprocess.run(
        [sys.executable, str(script), str(data.resolve() / TARGET), "load", "--db-url", url,
         "--config", str(config), "--data-dir", str(data)],
        capture_output=True, text=True, check=False,
    )
    assert killed.returncode == 9, killed.stderr
    assert (data / TARGET).is_dir()
    assert assemblies(url) == 0

    again = load(url, config, data)
    assert again.exit_code == 0, again.output
    assert "interrupted before it committed" in again.output
    assert assemblies(url) == 1
    assert run("verify", "--db-url", url, "--assembly-id", "1").exit_code == 0
    # The replaced directory is gone once the new load has committed.
    staging = data / ".staging"
    assert not staging.exists() or not any(staging.iterdir())


def test_a_directory_another_database_loaded_still_blocks(tmp_path):
    config = make_inputs(tmp_path / "in")
    data = tmp_path / "data"
    first = init(tmp_path / "first.db")
    assert load(first, config, data).exit_code == 0
    second = init(tmp_path / "second.db")

    refused = load(second, config, data)
    assert refused.exit_code != 0
    assert "another cgload database" in refused.output
    assert assemblies(second) == 0
    # The first database's assembly is untouched.
    assert run("verify", "--db-url", first, "--assembly-id", "1").exit_code == 0


def test_a_directory_that_names_no_database_still_blocks(tmp_path):
    config = make_inputs(tmp_path / "in")
    url = init(tmp_path / "g.db")
    data = tmp_path / "data"
    (data / TARGET).mkdir(parents=True)
    (data / TARGET / "manifest.json").write_text(json.dumps({"assembly": "asm1"}))

    refused = load(url, config, data)
    assert refused.exit_code != 0
    assert "already exists" in refused.output
    assert (data / TARGET / "manifest.json").exists()


def test_each_database_has_its_own_id_and_manifests_record_it(tmp_path):
    config = make_inputs(tmp_path / "in")
    ids = []
    for name in ("a.db", "b.db"):
        url = init(tmp_path / name)
        with build_engine(url).connect() as connection:
            ids.append(connection.execute(select(vocabulary_meta.c.database_id)).scalar_one())
    assert len(set(ids)) == 2
    assert load(url, config, tmp_path / "data").exit_code == 0
    manifest = json.loads((tmp_path / "data" / TARGET / "manifest.json").read_text())
    assert manifest["database_id"] == ids[1]
