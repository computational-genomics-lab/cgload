"""BLAST+ for the portal, installed by cgload when the computer has none (D-118).

BLAST+ is not a Python package, so pip cannot install it with cgload. When
``cgload portal`` finds no BLAST+ on the computer, it installs one, once, into
cgload's own folder (:func:`home`), and nothing anywhere else:

1. micromamba, one program from mamba-org's releases on GitHub, at a pinned
   version whose SHA-256 for each platform is written below. A download whose
   checksum differs is deleted before it can run.
2. micromamba makes an environment holding BLAST+ from Bioconda, at a pinned
   version, and checks every package it fetches against the channel's own
   checksums.

No shell profile, PATH or system package is changed; the portal finds BLAST+
in that folder by itself, as do ``cgload serve`` and ``cgload blastdb``.
Removing the folder removes it. Where Bioconda has no BLAST+ (Windows, and
processors other than x86-64 and ARM), :class:`InstallError` says how to
install BLAST+ by hand.
"""

from __future__ import annotations

import hashlib
import os
import platform
import shutil
import subprocess  # micromamba only, with an argument list  # nosec B404
import sys
import tempfile
import urllib.request
from pathlib import Path

MICROMAMBA_VERSION = "2.10.0-0"
MICROMAMBA_URL = ("https://github.com/mamba-org/micromamba-releases/releases/download/"
                  "{version}/micromamba-{platform}")
#: SHA-256 of each micromamba 2.10.0-0 program, as its release publishes them.
MICROMAMBA_SHA256 = {
    "linux-64": "701f39c06f8e6a670a0e186c8479b9e5c0502a2f0cda6cef65d83bdd2875e8bd",
    "linux-aarch64": "3bb56786e21b9b24024f599fe9185b4a59dc3f9f59937a42311ee583adbbb23c",
    "osx-64": "8421462a1ef1df9f7acd114f9e2a24c1332ce4f4e29f44465b4acb44e2755172",
    "osx-arm64": "c547639f5caf742236641c362746abb4a3505774753ecd08fe341643670c4c01",
}
#: The BLAST+ that Bioconda builds for each platform; 2.17.0 is the one the
#: paper's benchmark ran. Bioconda's ARM Linux build is 2.16.0.
BLAST_VERSIONS = {"linux-64": "2.17.0", "osx-64": "2.17.0", "osx-arm64": "2.17.0",
                  "linux-aarch64": "2.16.0"}
CHANNELS = ("conda-forge", "bioconda")
BY_HAND = ("Install BLAST+ by hand: conda install -c conda-forge -c bioconda blast, "
           "or sudo apt install ncbi-blast+, or NCBI's own programs from "
           "https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ on your PATH")


class InstallError(RuntimeError):
    """BLAST+ could not be installed; the message says why and what to do."""


def conda_platform() -> str | None:
    """This computer's platform as conda names it, or None if Bioconda has no
    BLAST+ for it."""
    machine = platform.machine().lower()
    if sys.platform.startswith("linux"):
        return {"x86_64": "linux-64", "amd64": "linux-64", "aarch64": "linux-aarch64",
                "arm64": "linux-aarch64"}.get(machine)
    if sys.platform == "darwin":
        return {"x86_64": "osx-64", "arm64": "osx-arm64"}.get(machine)
    return None


def homes() -> list[Path]:
    """Where cgload keeps the programs it installs, in the order it looks:
    ``CGLOAD_HOME`` alone if set; else ``<environment>/share/cgload`` in a
    virtual or conda environment, which goes when the environment goes, then
    the user's data folder."""
    if os.environ.get("CGLOAD_HOME"):
        return [Path(os.environ["CGLOAD_HOME"]).expanduser()]
    places = []
    prefix = Path(sys.prefix)
    if sys.prefix != sys.base_prefix or (prefix / "conda-meta").is_dir():
        places.append(prefix / "share" / "cgload")
    if sys.platform == "darwin":
        places.append(Path.home() / "Library" / "Application Support" / "cgload")
    elif os.name == "nt":
        places.append(Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
                      / "cgload")
    else:
        places.append(Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
                      / "cgload")
    return places


def _writable(place: Path) -> bool:
    """Whether cgload could make or change ``place``: its nearest existing
    folder can be written."""
    while not place.exists() and place != place.parent:
        place = place.parent
    return os.access(place, os.W_OK)


def home() -> Path:
    """The folder an install goes into: the first of :func:`homes` cgload may
    write to. Looking does not depend on it: a service whose file system is
    read-only still finds what was installed (:func:`blast_directories`)."""
    places = homes()
    return next((place for place in places if _writable(place)), places[-1])


def blast_directory() -> Path:
    """Where an install puts cgload's own BLAST+ programs."""
    return home() / "blast" / "bin"


def blast_directories() -> list[Path]:
    """Every folder where cgload's own BLAST+ programs may be, in order."""
    return [place / "blast" / "bin" for place in homes()]


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, target: Path, sha256: str, *, timeout: float = 120) -> None:
    """Fetch ``url`` over HTTPS into ``target``, keeping it only if its SHA-256
    is ``sha256``."""
    if not url.startswith("https://"):
        raise InstallError(f"refusing to download over anything but HTTPS: {url}")
    target.parent.mkdir(parents=True, exist_ok=True)
    handle, name = tempfile.mkstemp(dir=target.parent, prefix=".download-")
    partial = Path(name)
    try:
        with os.fdopen(handle, "wb") as out, \
                urllib.request.urlopen(url, timeout=timeout) as response:  # noqa: S310 - https only, checked above  # nosec B310
            shutil.copyfileobj(response, out, 1 << 20)
        found = sha256_of(partial)
        if found != sha256:
            raise InstallError(f"the download of {url} does not have the SHA-256 cgload "
                               f"expects ({found[:16]}…, not {sha256[:16]}…); it was "
                               f"deleted and not run")
        partial.chmod(0o755)
        partial.replace(target)
    except OSError as exc:
        raise InstallError(f"could not download {url}: {exc}") from None
    finally:
        partial.unlink(missing_ok=True)


def micromamba_command(micromamba: Path, prefix: Path, blast_version: str) -> list[str]:
    """The one command that installs BLAST+: an argument list, no shell."""
    command = [str(micromamba), "create", "--yes", "--quiet", "--prefix", str(prefix),
               "--override-channels"]
    for channel in CHANNELS:
        command += ["--channel", channel]
    return command + [f"blast={blast_version}"]


def install_blast(log=print, *, timeout: float = 1800) -> Path:
    """Install BLAST+ into cgload's folder; the folder of its programs.

    Takes a few minutes and a few hundred megabytes, once. Raises
    InstallError, with what to do instead, when it cannot."""
    target = conda_platform()
    if target is None:
        raise InstallError(f"cgload cannot install BLAST+ on this computer "
                           f"({sys.platform}, {platform.machine()}). {BY_HAND}.")
    root = home()
    prefix = root / "blast"
    if (blast_directory() / "makeblastdb").exists():
        return blast_directory()
    micromamba = root / "micromamba" / MICROMAMBA_VERSION / "micromamba"
    if not micromamba.exists() or sha256_of(micromamba) != MICROMAMBA_SHA256[target]:
        log(f"  downloading micromamba {MICROMAMBA_VERSION}, which installs BLAST+ from "
            f"Bioconda")
        download(MICROMAMBA_URL.format(version=MICROMAMBA_VERSION, platform=target),
                 micromamba, MICROMAMBA_SHA256[target])
    version = BLAST_VERSIONS[target]
    log(f"  installing BLAST+ {version} into {prefix} (once; a few minutes)")
    shutil.rmtree(prefix, ignore_errors=True)  # a half-made environment from before
    environment = dict(os.environ, MAMBA_ROOT_PREFIX=str(root / "mamba"))
    try:
        done = subprocess.run(  # noqa: S603 - a fixed argument list; no shell  # nosec B603
            micromamba_command(micromamba, prefix, version), env=environment,
            capture_output=True, text=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        shutil.rmtree(prefix, ignore_errors=True)
        raise InstallError(f"micromamba could not install BLAST+ ({exc}). {BY_HAND}.") from None
    if done.returncode != 0 or not (blast_directory() / "makeblastdb").exists():
        shutil.rmtree(prefix, ignore_errors=True)
        said = " ".join((done.stderr or done.stdout).strip().splitlines()[-3:])[:400]
        raise InstallError(f"micromamba could not install BLAST+: {said}. {BY_HAND}.")
    # The packages are unpacked into the environment; the downloads are not needed.
    shutil.rmtree(root / "mamba" / "pkgs", ignore_errors=True)
    return blast_directory()
