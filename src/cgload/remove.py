"""``cgload remove``: take one assembly out of the database and the data directory.

The load is transactional across the database and the data directory
(``storage.staged_load``), and so is removal (``storage.staged_removal``):
either the rows and the files both go, or neither does (D-084).

Removal never guesses which files to delete. The assembly row records where its
data directory was at load time (D-079) and where the FASTA sits inside it. The
directory is deleted only if its ``manifest.json`` names this assembly and
records the same FASTA checksum as the database. A directory that cannot be
found, or that holds something else, stops the removal before anything changes;
``--keep-files`` removes the rows alone, on purpose.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

from sqlalchemy import Engine, func, select

from cgload.db.removal import Footprint, delete_assembly_rows, footprint
from cgload.db.schema import (
    annotation_run,
    assembly,
    feature,
    organism,
    sequence_region,
    source_file,
    strain,
)
from cgload.storage import staged_removal


class RemovalRefusedError(RuntimeError):
    """The removal would delete files cgload cannot show belong to the assembly."""


@dataclass(frozen=True)
class RemovalPlan:
    assembly_id: int
    label: str
    footprint: Footprint
    #: The data directory that will be deleted, or None when files are kept.
    directory: Path | None
    data_root: Path | None
    #: Where the files are, when they are deliberately kept.
    kept_directory: Path | None = None


def _label(row) -> str:
    version = f" {row.version}" if row.version else ""
    return f"{row.genus} {row.species} / {row.strain} / {row.name}{version}"


def plan_removal(
    engine: Engine,
    assembly_id: int,
    *,
    data_dir: Path | None = None,
    keep_files: bool = False,
) -> RemovalPlan:
    """Work out what removing ``assembly_id`` would delete, changing nothing."""
    with engine.connect() as connection:
        row = connection.execute(
            select(
                assembly.c.name,
                assembly.c.version,
                assembly.c.fasta_path,
                assembly.c.fasta_sha256,
                assembly.c.data_root,
                strain.c.name.label("strain"),
                organism.c.genus,
                organism.c.species,
            )
            .select_from(assembly.join(strain).join(organism))
            .where(assembly.c.assembly_id == assembly_id)
        ).one_or_none()
        if row is None:
            raise LookupError(f"no assembly with id {assembly_id}")
        counted = footprint(connection, assembly_id)

    label = _label(row)
    root = Path(data_dir) if data_dir is not None else (
        Path(row.data_root) if row.data_root else None
    )
    directory = None
    if row.fasta_path and root is not None:
        directory = root / PurePosixPath(row.fasta_path).parent

    if keep_files:
        return RemovalPlan(assembly_id, label, counted, None, root, kept_directory=directory)

    if directory is None:
        raise RemovalRefusedError(
            f"the database does not record where assembly {assembly_id}'s data "
            f"directory is. Pass --data-dir to remove its files too, or "
            f"--keep-files to remove only the database rows."
        )
    _check_ownership(directory, root, row, label)
    return RemovalPlan(assembly_id, label, counted, directory, root)


def _check_ownership(directory: Path, root: Path, row, label: str) -> None:
    """Refuse unless ``directory`` is this assembly's, by its own manifest."""
    resolved, resolved_root = directory.resolve(), root.resolve()
    if resolved == resolved_root or not resolved.is_relative_to(resolved_root):
        raise RemovalRefusedError(
            f"{directory} is not inside the data directory {root}; nothing was removed."
        )
    if not directory.is_dir():
        raise RemovalRefusedError(
            f"the files of {label} are not at {directory}. If the data directory "
            f"moved, pass --data-dir. To remove only the database rows, pass "
            f"--keep-files. Nothing was removed."
        )
    manifest_path = directory / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, ValueError):
        raise RemovalRefusedError(
            f"{directory} has no readable manifest.json, so cgload cannot show it "
            f"holds {label}. Nothing was removed."
        ) from None

    recorded = {
        (entry.get("role"), entry.get("sha256"))
        for entry in manifest.get("files", [])
        if isinstance(entry, dict)
    }
    same_assembly = (
        manifest.get("organism") == f"{row.genus} {row.species}"
        and manifest.get("strain") == row.strain
        and manifest.get("assembly") == row.name
    )
    same_fasta = row.fasta_sha256 is None or ("assembly", row.fasta_sha256) in recorded
    if not (same_assembly and same_fasta):
        raise RemovalRefusedError(
            f"the manifest in {directory} does not describe {label} as loaded in "
            f"this database (names or assembly checksum differ). Nothing was removed."
        )


def remove_assembly(
    engine: Engine,
    assembly_id: int,
    *,
    data_dir: Path | None = None,
    keep_files: bool = False,
) -> RemovalPlan:
    """Remove the assembly's rows and, unless kept, its data directory.

    Returns the plan that was carried out. The rows are counted again inside the
    transaction after the delete, and anything left over stops the removal.
    """
    plan = plan_removal(engine, assembly_id, data_dir=data_dir, keep_files=keep_files)
    with staged_removal(engine, plan.directory, data_root=plan.data_root) as connection:
        delete_assembly_rows(connection, assembly_id, prune=True)
        left = {
            table.name: connection.execute(
                select(func.count()).select_from(table).where(table.c.assembly_id == assembly_id)
            ).scalar_one()
            for table in (assembly, feature, sequence_region, source_file, annotation_run)
        }
        if any(left.values()):
            raise RuntimeError(
                f"rows of assembly {assembly_id} remained after the delete: "
                + ", ".join(f"{name} {count}" for name, count in left.items() if count)
                + ". Nothing was removed. This is a bug in cgload; please report it."
            )
    return plan
