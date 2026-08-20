"""``cgload load`` -- the real thing (milestone 4)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import click

from cgload.assembly import AssemblyExistsError
from cgload.config import ConfigError, load_config
from cgload.db.engine import (
    CredentialInArgumentError,
    DatabaseNotFoundError,
    build_engine,
)
from cgload.db.init import SchemaMissingError, VocabularyMismatchError
from cgload.fasta import AssemblyReadError
from cgload.functional import FunctionalLoadError
from cgload.loader import LoadError, load_organism
from cgload.parsers.functional import FunctionalParseError
from cgload.parsers.tokenizer import TokenizeError
from cgload.storage import DataDirectoryConflictError

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"

#: Everything a bad input can raise. Each of these already names the file and,
#: where it can, the line -- so the handler prints the message and exits 2
#: rather than showing a traceback.
INPUT_ERRORS = (
    AssemblyExistsError,
    AssemblyReadError,
    ConfigError,
    CredentialInArgumentError,
    DataDirectoryConflictError,
    DatabaseNotFoundError,
    FileNotFoundError,
    FunctionalLoadError,
    FunctionalParseError,
    LoadError,
    SchemaMissingError,
    TokenizeError,
    VocabularyMismatchError,
)


@click.command("load")
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Organism TOML: the assembly, the annotation files, and the taxon id.",
)
@click.option(
    "--data-dir",
    default="data",
    type=click.Path(file_okay=False, path_type=Path),
    help="Where the copied assembly, its index and the annotation live.",
)
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option(
    "--link",
    is_flag=True,
    help="Hard-link the assembly instead of copying it. Saves disk, but a file the "
    "user later edits or deletes makes every stored coordinate wrong (D-018).",
)
@click.option(
    "--repair-envelopes",
    is_flag=True,
    help="Recompute a parent's extent from the features it contains, where the "
    "file's own extent excludes them. Needed for NCBI's trans-spliced organellar "
    "genes; every change is printed and written into the manifest (D-055).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Replace an assembly already loaded under this name. Destroys its features.",
)
def load_command(
    config_path: Path,
    data_dir: Path,
    db_url: str | None,
    link: bool,
    repair_envelopes: bool,
    force: bool,
) -> None:
    """Load one organism: assembly and annotation, together or not at all."""
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL

    try:
        config = load_config(config_path)
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        report = load_organism(
            engine,
            config,
            data_dir,
            copy=not link,
            force=force,
            repair_envelopes=repair_envelopes,
        )
    except INPUT_ERRORS as exc:
        click.echo(f"cgload load: {exc}", err=True)
        sys.exit(2)

    click.echo(
        f"Loaded {config.genus} {config.species} / {config.strain} / "
        f"{config.assembly_name} (assembly_id={report.assembly_id})"
    )
    click.echo(f"  {report.region_count:,} sequence regions, {report.total_length:,} bp")

    for entry in report.files:
        click.echo(
            f"  {entry.path.name}: {entry.dialect} ({entry.confidence}) -- "
            f"{entry.feature_count:,} features"
        )
        shown = ", ".join(
            f"{name} {count:,}" for name, count in sorted(entry.counts_by_type.items())
        )
        click.echo(f"    {shown}")
        if entry.protein_count:
            click.echo(
                f"    {entry.protein_count:,} declared proteins checked against the CDS "
                f"spans and written to annotation/proteins.faa"
            )

    for item in report.envelope_repairs:
        click.echo(
            f"  repaired envelope: {item.feature_type} {item.feature_id} "
            f"{item.old_start:,}-{item.old_end:,} -> {item.new_start:,}-{item.new_end:,}",
            err=True,
        )

    for entry in report.functional:
        version = f" {entry.tool_version}" if entry.tool_version else ""
        click.echo(
            f"  {entry.path.name}: {entry.tool}{version} -- {entry.hit_count:,} "
            f"annotations on {entry.matched_count:,}/{entry.query_count:,} proteins"
        )
        shown = ", ".join(
            f"{name} {count:,}" for name, count in sorted(entry.analyses.items())
        )
        click.echo(f"    {shown}")

    for warning in report.warnings:
        click.echo(f"  warning: {warning}", err=True)

    click.echo(f"  data directory: {report.data_directory}")
