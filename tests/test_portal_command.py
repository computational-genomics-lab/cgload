"""``cgload portal``: genome files to a running portal in one command (D-118).

What these assert, stated before they were written:

* From an organism file and nothing else, one command creates the database,
  loads and checks the genome, indexes it, builds its BLAST databases and is
  ready to serve; run again, it does none of that twice.
* With no database and no genome, it says what to give; a genome that does
  not load stops it before anything is served, in a sentence.
* When the computer has no BLAST+, cgload installs its own, into its own
  folder only: micromamba at a pinned version whose SHA-256 must match, then
  BLAST+ from Bioconda, by an argument list. A download that does not match
  is deleted and never run; a platform Bioconda does not serve, or a failed
  install, is said in words, and the portal goes on without BLAST.
* ``-m network`` installs BLAST+ for real.
"""

from __future__ import annotations

import io
import os
import shutil
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("flask")

from cgload.cli.main import cli  # noqa: E402
from cgload.portal import install  # noqa: E402
from click.testing import CliRunner  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "minimal"
HAVE_BLAST = all(shutil.which(name) for name in ("blastn", "tblastn", "blastp", "blastx",
                                                  "makeblastdb"))


@pytest.fixture
def here(tmp_path, monkeypatch):
    """A folder holding the example genome, as a user's would."""
    for name in ("organism.toml", "assembly.fa", "annotation.gff3"):
        shutil.copy(EXAMPLE / name, tmp_path / name)
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("CGLOAD_DB_URL", raising=False)
    monkeypatch.setenv("CGLOAD_HOME", str(tmp_path / "cgload-home"))
    return tmp_path


def run(*args: str):
    return CliRunner().invoke(cli, ["portal", *args])


def test_one_command_loads_indexes_and_gets_ready(here):
    first = run("organism.toml", "--prepare-only", "--no-blast")
    assert first.exit_code == 0, first.output
    assert "database: created" in first.output
    assert "Demonstratus minimalis / DEMO-1 / Demo_v1" in first.output
    assert "assembly 1: 1 sequence, 12 features; every check passed" in first.output
    assert "search index: assembly 1" in first.output
    assert "Start the portal with: cgload serve" in first.output
    assert (here / "cgload.db").exists() and (here / "data").is_dir()
    again = run("organism.toml", "--prepare-only", "--no-blast")
    assert again.exit_code == 0, again.output
    assert "is loaded already (assembly 1)" in again.output
    assert "search index: in place for 1 genome" in again.output
    assert "database: created" not in again.output
    alone = run("--prepare-only", "--no-blast")  # no genome given: what is there
    assert alone.exit_code == 0 and "in place for 1 genome" in alone.output


def test_without_a_database_or_a_genome_it_says_what_to_give(here):
    result = run("--prepare-only")
    assert result.exit_code == 2
    assert "Give the genomes to load" in result.output and "Traceback" not in result.output
    assert not (here / "cgload.db").exists()


def test_a_genome_that_does_not_load_stops_it_before_serving(here):
    (here / "broken.toml").write_text((here / "organism.toml").read_text()
                                      .replace("assembly.fa", "missing.fa"))
    result = run("broken.toml", "--no-blast")
    assert result.exit_code == 2
    assert "broken.toml" in result.output and "Nothing is served" in result.output
    assert "Traceback" not in result.output and "portal is running" not in result.output


def test_a_failed_blast_install_leaves_the_portal_without_blast(here, monkeypatch):
    from cgload.portal import blast

    monkeypatch.setattr(blast, "find_binaries", lambda directory=None: {})

    def refuse(log=print, **_):
        raise install.InstallError("no network here. Install BLAST+ by hand: ...")

    monkeypatch.setattr(install, "install_blast", refuse)
    result = run("organism.toml", "--prepare-only")
    assert result.exit_code == 0, result.output
    assert "cgload installs it for itself, once" in result.output
    assert "BLAST+: no network here" in result.output
    assert "BLAST databases" not in result.output


@pytest.mark.skipif(not HAVE_BLAST, reason="BLAST+ is not installed")
def test_with_blast_it_builds_the_databases_once(here):
    first = run("organism.toml", "--prepare-only")
    assert first.exit_code == 0, first.output
    assert "BLAST databases: building for 1 genome" in first.output
    assert ": PASS" in first.output
    again = run("--prepare-only")
    assert "BLAST databases: in place for 1 genome" in again.output


# --------------------------------------------------------------------------
# Installing BLAST+
# --------------------------------------------------------------------------


@pytest.mark.parametrize(("system", "machine", "expected"), [
    ("linux", "x86_64", "linux-64"), ("linux", "aarch64", "linux-aarch64"),
    ("darwin", "arm64", "osx-arm64"), ("darwin", "x86_64", "osx-64"),
    ("win32", "AMD64", None), ("linux", "ppc64le", None), ("freebsd14", "amd64", None),
])
def test_platforms_bioconda_serves(monkeypatch, system, machine, expected):
    monkeypatch.setattr(install.sys, "platform", system)
    monkeypatch.setattr(install.platform, "machine", lambda: machine)
    assert install.conda_platform() == expected
    if expected:
        assert expected in install.MICROMAMBA_SHA256 and expected in install.BLAST_VERSIONS


def test_cgloads_folder_is_its_own(monkeypatch, tmp_path):
    monkeypatch.setenv("CGLOAD_HOME", str(tmp_path / "mine"))
    assert install.home() == tmp_path / "mine"
    assert install.blast_directory() == tmp_path / "mine" / "blast" / "bin"
    monkeypatch.delenv("CGLOAD_HOME")
    monkeypatch.setattr(install.sys, "prefix", str(tmp_path / "env"))
    monkeypatch.setattr(install.sys, "base_prefix", "/usr")
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "data"))
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    (tmp_path / "env").mkdir()
    in_env = tmp_path / "env" / "share" / "cgload"
    assert install.home() == in_env
    # A service with a read-only file system cannot write to the environment;
    # it must still find what an install put there.
    monkeypatch.setattr(install, "_writable", lambda place: False)
    assert in_env / "blast" / "bin" in install.blast_directories()
    assert install.blast_directories()[0] == in_env / "blast" / "bin"


class _Response(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def test_a_download_is_kept_only_if_its_checksum_matches(monkeypatch, tmp_path):
    import hashlib

    body = b"#!/bin/sh\necho micromamba\n"
    monkeypatch.setattr(install.urllib.request, "urlopen",
                        lambda url, timeout=0: _Response(body))
    target = tmp_path / "bin" / "micromamba"
    install.download("https://example.org/mm", target, hashlib.sha256(body).hexdigest())
    assert target.read_bytes() == body and os.access(target, os.X_OK)
    other = tmp_path / "bin" / "tampered"
    with pytest.raises(install.InstallError, match="does not have the SHA-256"):
        install.download("https://example.org/mm", other, "0" * 64)
    assert not other.exists() and list((tmp_path / "bin").iterdir()) == [target]
    with pytest.raises(install.InstallError, match="HTTPS"):
        install.download("http://example.org/mm", other, "0" * 64)


def test_blast_is_installed_by_an_argument_list_into_cgloads_folder(monkeypatch, tmp_path):
    monkeypatch.setenv("CGLOAD_HOME", str(tmp_path))
    monkeypatch.setattr(install, "conda_platform", lambda: "linux-64")
    fetched, ran = [], []

    def fake_download(url, target, sha256, **_):
        fetched.append((url, sha256))
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("micromamba")

    def fake_run(args, **kwargs):
        ran.append((args, kwargs))
        (tmp_path / "blast" / "bin").mkdir(parents=True)
        (tmp_path / "blast" / "bin" / "makeblastdb").write_text("")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(install, "download", fake_download)
    monkeypatch.setattr(install, "sha256_of", lambda path: "x")
    monkeypatch.setattr(install.subprocess, "run", fake_run)
    assert install.install_blast(log=lambda line: None) == tmp_path / "blast" / "bin"
    assert fetched == [(install.MICROMAMBA_URL.format(version=install.MICROMAMBA_VERSION,
                                                      platform="linux-64"),
                        install.MICROMAMBA_SHA256["linux-64"])]
    args, kwargs = ran[0]
    assert args[1:] == ["create", "--yes", "--quiet", "--prefix", str(tmp_path / "blast"),
                        "--override-channels", "--channel", "conda-forge", "--channel",
                        "bioconda", "blast=2.17.0"]
    assert not kwargs.get("shell") and kwargs["env"]["MAMBA_ROOT_PREFIX"] == str(tmp_path /
                                                                                 "mamba")
    assert install.install_blast(log=lambda line: None) == tmp_path / "blast" / "bin"
    assert len(ran) == 1  # once installed, never again


def test_a_failed_or_impossible_install_says_what_to_do(monkeypatch, tmp_path):
    monkeypatch.setenv("CGLOAD_HOME", str(tmp_path))
    monkeypatch.setattr(install, "conda_platform", lambda: None)
    with pytest.raises(install.InstallError, match="by hand"):
        install.install_blast(log=lambda line: None)
    monkeypatch.setattr(install, "conda_platform", lambda: "osx-arm64")
    monkeypatch.setattr(install, "download", lambda url, target, sha256, **_: (
        target.parent.mkdir(parents=True, exist_ok=True), target.write_text("mm")))
    monkeypatch.setattr(install, "sha256_of", lambda path: "x")

    def fail(args, **kwargs):
        (tmp_path / "blast").mkdir(exist_ok=True)
        return SimpleNamespace(returncode=1, stdout="", stderr="could not solve\nno blast")

    monkeypatch.setattr(install.subprocess, "run", fail)
    with pytest.raises(install.InstallError, match="no blast.*by hand"):
        install.install_blast(log=lambda line: None)
    assert not (tmp_path / "blast").exists()  # no half-made environment is left


@pytest.mark.network
@pytest.mark.skipif(install.conda_platform() is None, reason="Bioconda has no BLAST+ here")
def test_blast_installs_for_real(tmp_path, monkeypatch):
    monkeypatch.setenv("CGLOAD_HOME", str(tmp_path))
    directory = install.install_blast(log=lambda line: None)
    import subprocess

    shown = subprocess.run([str(directory / "blastp"), "-version"], capture_output=True,
                           text=True, check=True)
    assert install.BLAST_VERSIONS[install.conda_platform()] in shown.stdout
    from cgload.portal.blast import find_binaries

    monkeypatch.setenv("PATH", os.defpath)
    monkeypatch.setattr(sys, "executable", str(tmp_path / "python"))
    assert Path(find_binaries()["makeblastdb"]).parent == directory.resolve()
