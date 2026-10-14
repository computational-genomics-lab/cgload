"""``cgload remove`` -- take an assembly out, rows and files together (D-084)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import click

from cgload.db.engine import (
    CredentialInArgumentError,
    DatabaseNotFoundError,
    build_engine,
)
from cgload.db.init import (
    SchemaMissingError,
    VocabularyMismatchError,
    assert_vocabulary_current,
)
from cgload.remove import RemovalPlan, RemovalRefusedError, plan_removal, remove_assembly

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"


def _describe(plan: RemovalPlan) -> str:
    counted = plan.footprint
    lines = [
        f"assembly {plan.assembly_id}: {plan.label}",
        f"  {counted.features:,} features in {counted.rows:,} rows, "
        f"{counted.attributes:,} attributes",
        f"  {counted.functional_assignments:,} functional assignments from "
        f"{counted.annotation_runs:,} file(s)",
        f"  {counted.sequence_regions:,} sequence regions, {counted.source_files:,} "
        f"source file record(s)",
    ]
    if plan.directory is not None:
        lines.append(f"  data directory {plan.directory}")
    elif plan.kept_directory is not None:
        lines.append(f"  data directory kept: {plan.kept_directory}")
    else:
        lines.append("  data directory kept (location not recorded)")
    return "\n".join(lines)


@click.command("remove")
@click.option("--assembly-id", type=int, required=True, help="Which assembly to remove.")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option(
    "--data-dir",
    type=click.Path(file_okay=False, path_type=Path),
    default=None,
    help="Data directory holding the assembly's files. Defaults to where it was "
    "when the assembly was loaded; pass this if the directory has moved since.",
)
@click.option(
    "--keep-files",
    is_flag=True,
    help="Remove only the database rows and leave the data directory alone.",
)
@click.option("--dry-run", is_flag=True, help="Show what would be removed, and change nothing.")
@click.option("--yes", "-y", is_flag=True, help="Do not ask for confirmation.")
def remove_command(
    assembly_id: int,
    db_url: str | None,
    data_dir: Path | None,
    keep_files: bool,
    dry_run: bool,
    yes: bool,
) -> None:
    """Remove an assembly: its rows and its data directory, or nothing.

    Every feature, attribute, functional assignment, sequence region and file
    record of the assembly is deleted in one transaction, and the data
    directory goes with it. The strain and organism go too when no other
    assembly uses them.
    """
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL
    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        assert_vocabulary_current(engine)
        plan = plan_removal(engine, assembly_id, data_dir=data_dir, keep_files=keep_files)
    except (
        CredentialInArgumentError,
        DatabaseNotFoundError,
        SchemaMissingError,
        VocabularyMismatchError,
        LookupError,
        RemovalRefusedError,
    ) as exc:
        click.echo(f"cgload remove: {exc}", err=True)
        sys.exit(2)

    click.echo(_describe(plan))
    if dry_run:
        click.echo("Dry run: nothing was removed.")
        return
    if not yes:
        click.confirm("Remove it?", abort=True)

    try:
        remove_assembly(engine, assembly_id, data_dir=data_dir, keep_files=keep_files)
    except (LookupError, RemovalRefusedError) as exc:
        click.echo(f"cgload remove: {exc}", err=True)
        sys.exit(2)
    click.echo(f"Removed assembly {assembly_id}.")
