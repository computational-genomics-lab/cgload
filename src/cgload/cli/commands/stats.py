"""``cgload stats`` -- what is in this database (milestone 2)."""

from __future__ import annotations

import os
import sys

import click
from sqlalchemy import func, select

from cgload.db.counts import counts_by
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
from cgload.db.schema import assembly, feature, organism, sequence_region, strain

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"


@click.command("stats")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option(
    "--by",
    type=click.Choice(["both", "source", "canonical"]),
    default="both",
    help="Which feature-type column to count by. 'source' is the term the file "
    "used and agrees with `verify`; 'canonical' is the normalised term and "
    "agrees with `query`. Default shows both, because they can differ.",
)
def stats_command(db_url: str | None, by: str) -> None:
    """Counts per organism, strain, assembly and feature type."""
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL
    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        assert_vocabulary_current(engine)
    except (
        CredentialInArgumentError,
        DatabaseNotFoundError,
        SchemaMissingError,
        VocabularyMismatchError,
    ) as exc:
        click.echo(f"cgload stats: {exc}", err=True)
        sys.exit(2)

    with engine.connect() as conn:
        rows = conn.execute(
            select(
                organism.c.genus,
                organism.c.species,
                strain.c.name.label("strain"),
                assembly.c.assembly_id,
                assembly.c.name.label("assembly"),
                assembly.c.version,
                func.count(sequence_region.c.sequence_region_id).label("regions"),
                func.coalesce(func.sum(sequence_region.c.length), 0).label("bases"),
            )
            .select_from(
                organism.join(strain).join(assembly).outerjoin(
                    sequence_region, sequence_region.c.assembly_id == assembly.c.assembly_id
                )
            )
            .group_by(assembly.c.assembly_id)
            .order_by(organism.c.genus, organism.c.species, strain.c.name)
        ).all()

        if not rows:
            click.echo("No assemblies loaded.")
            return

        for row in rows:
            version = f" {row.version}" if row.version else ""
            click.echo(
                f"{row.genus} {row.species} / {row.strain} / {row.assembly}{version}\n"
                f"  {row.regions:,} sequence regions, {row.bases:,} bp"
            )
            _feature_counts(conn, row.assembly_id, by)


def _feature_counts(conn, assembly_id: int, by: str) -> None:
    columns = {"source": [("source", feature.c.source_type)],
               "canonical": [("canonical", feature.c.feature_type)]}
    columns["both"] = columns["source"] + columns["canonical"]

    for label, column in columns[by]:
        counts = conn.execute(counts_by(column, assembly_id)).all()
        if not counts:
            click.echo("  no features loaded")
            return
        shown = ", ".join(f"{name} {n:,}" for name, n in counts)
        click.echo(f"  by {label} type: {shown}")
