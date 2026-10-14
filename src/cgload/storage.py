"""The data directory, and the transaction that spans it and the database.

Sequence lives on disk (D-004), so "succeeds completely or changes nothing"
is only true if the filesystem takes part in the transaction. This module is
that mechanism.

Layout, relative to the data directory root::

    <root>/
      .staging/<uuid>/            work in progress; never visible to a reader
      <organism>/<strain>/<assembly>/
        assembly.fa               the copied assembly
        assembly.fa.fai           the index (D-004)
        manifest.json             what was loaded, from where, under which
                                  vocabulary and cgload version
        annotation.gff3.gz        the normalised annotation, bgzipped (D-082)
        annotation.gff3.gz.tbi    its tabix index, for JBrowse 2, IGV, tabix
        annotation/               declared proteins, where the dialect has them

**Ordering (amends D-006).** The decision as written stages files, commits the
database, then renames. That leaves a window: a crash between commit and
rename leaves rows pointing at a directory that is not there. Here the rename
happens *inside* the open database transaction, immediately before commit, and
is reversed if the commit fails. The remaining window is one rename plus one
commit, and the failure it leaves -- a directory in place with no rows -- is
detectable and harmless, where the reverse is not.

**Under --force (AUDIT-4, amends D-073)** the directory being replaced stays
where it is for the whole load and is moved aside only in that final step,
next to the rename of the new one. Moving it aside when the load *started*
meant a process killed mid-load (SIGTERM, SIGKILL, the OOM killer -- none of
which run Python's exception handling) left the old rows pointing at a
directory that was gone. The window that remains is the two renames plus the
commit; `verify` checks the stored assembly's checksum, so a crash inside it is
detected rather than silent, and the old directory is still under ``.staging``.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

from sqlalchemy import Connection, Engine

CHUNK = 1024 * 1024


class DataDirectoryConflictError(RuntimeError):
    """The target directory for this assembly already exists."""


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while chunk := fh.read(CHUNK):
            h.update(chunk)
    return h.hexdigest()


def slug(text: str) -> str:
    """Filesystem-safe directory component.

    Deliberately strict: a species name with a slash or a space becomes a path
    bug on someone else's machine, and case-insensitive filesystems make
    'CJ26' and 'cj26' the same directory.
    """
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip()).strip("._-")
    if not cleaned:
        raise ValueError(f"{text!r} has no filesystem-safe form")
    return cleaned.lower()


@dataclass
class StagedLoad:
    """A directory being built, plus the database transaction it belongs to."""

    staging: Path
    target: Path
    connection: Connection
    manifest: dict = field(default_factory=dict)
    #: The data directory root, absolute. Recorded on the assembly row so that
    #: `verify` can find the stored files (D-079).
    data_root: Path | None = None

    def place(self, source: Path, name: str, *, copy: bool = True) -> tuple[Path, str]:
        """Copy (default) or hard-link a file into the staging directory.

        Returns the path relative to the data root, and the SHA-256. Copying
        is the default because the database and its data directory must stay
        self-contained: a referenced file that the user later edits or deletes
        makes every stored coordinate silently wrong.
        """
        destination = self.staging / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if copy:
            shutil.copy2(source, destination)
        else:
            os.link(source, destination)
        return destination, sha256_of(destination)

    def relative(self, path: Path) -> str:
        """Path as stored in the database: relative to the data root, POSIX
        separators. An absolute path makes the database non-portable."""
        root = self.target.parent.parent.parent
        return (self.target.relative_to(root) / path.relative_to(self.staging)).as_posix()


@contextmanager
def staged_load(
    engine: Engine,
    data_root: Path,
    *,
    organism: str,
    strain: str,
    assembly: str,
    force: bool = False,
) -> Iterator[StagedLoad]:
    """Run a load. Either the rows and the files both land, or neither does."""
    data_root = Path(data_root).resolve()
    target = data_root / slug(organism) / slug(strain) / slug(assembly)

    # Under --force the existing directory is replaced, never deleted up front
    # and never moved before the end of the load.
    #
    # Deleting it at the start broke the atomicity guarantee outright (D-073):
    # the rmtree ran before the transaction opened and no exception path undid
    # it. Moving it aside at the start was the next version, and it still left a
    # killed process with the old rows pointing at a directory that was gone
    # (AUDIT-4). It is now moved aside only immediately before the new
    # directory is renamed into place, inside the open transaction.
    superseded: Path | None = None
    if target.exists() and not force:
        raise DataDirectoryConflictError(
            f"{target} already exists. Use a different assembly name, or "
            f"--force to replace it."
        )

    staging_root = data_root / ".staging"
    staging_root.mkdir(parents=True, exist_ok=True)
    # Same filesystem as the target, so the final rename is atomic.
    staging = Path(tempfile.mkdtemp(dir=staging_root))

    connection = engine.connect()
    transaction = connection.begin()
    load = StagedLoad(staging=staging, target=target, connection=connection, data_root=data_root)
    renamed = False
    try:
        yield load

        (staging / "manifest.json").write_text(json.dumps(load.manifest, indent=2, default=str))
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            # Same filesystem as the target, so this and the restore are atomic.
            superseded = Path(tempfile.mkdtemp(dir=staging_root)) / "superseded"
            os.rename(target, superseded)
        os.rename(staging, target)
        renamed = True
        transaction.commit()
        # Only now is the superseded directory unrecoverable, and only after
        # the database has agreed to the replacement.
        if superseded is not None:
            shutil.rmtree(superseded.parent, ignore_errors=True)
    except BaseException:
        if renamed:
            # Commit failed after the directory landed. Put it back so the
            # data directory matches the database, which is now unchanged.
            os.rename(target, staging)
        if superseded is not None and not target.exists():
            os.rename(superseded, target)
            shutil.rmtree(superseded.parent, ignore_errors=True)
        transaction.rollback()
        shutil.rmtree(staging, ignore_errors=True)
        raise
    finally:
        connection.close()
