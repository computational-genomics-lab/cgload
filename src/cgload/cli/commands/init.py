"""``cgload init`` -- create the schema (milestone 1)."""

from __future__ import annotations

import os
import sys

import click

from cgload.db.engine import CredentialInArgumentError, build_engine
from cgload.db.init import SchemaAlreadyPresentError, create_schema

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"


@click.command("init")
@click.option(
    "--db-url",
    default=None,
    help=f"SQLAlchemy URL. Defaults to ${ENV_DB_URL}, then {DEFAULT_DB_URL}. "
    "A URL containing a password is refused here (D-003); use the environment.",
)
@click.option("--force", is_flag=True, help="Drop and recreate. Destroys existing data.")
def init_command(db_url: str | None, force: bool) -> None:
    """Create the schema and seed the controlled vocabulary."""
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL

    try:
        engine = build_engine(resolved, from_environment=from_environment)
    except CredentialInArgumentError as exc:
        click.echo(f"cgload init: {exc}", err=True)
        sys.exit(2)

    try:
        result = create_schema(engine, force=force)
    except SchemaAlreadyPresentError as exc:
        click.echo(f"cgload init: {exc}", err=True)
        sys.exit(2)

    safe = engine.url.render_as_string(hide_password=True)
    click.echo(f"Created {len(result.tables)} tables in {safe}")
    click.echo(
        f"Seeded vocabulary version {result.vocabulary_version} "
        f"({result.vocabulary_terms} terms, {result.vocabulary_hash[:12]})"
    )
