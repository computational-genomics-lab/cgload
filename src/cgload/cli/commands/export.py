"""``cgload export`` -- write the database back out (milestone 9)."""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

import click
from sqlalchemy import select

from cgload.db.engine import (
    CredentialInArgumentError,
    DatabaseNotFoundError,
    build_engine,
)
from cgload.db.export_source import load_context, load_features
from cgload.db.init import SchemaMissingError, VocabularyMismatchError
from cgload.db.schema import annotation_run, protein_annotation
from cgload.export import export_gff3

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"

#: Types no source dialect in the registry emits for every assembly, recorded in
#: the export header so absence stays distinguishable from loss.
KNOWN_OPTIONAL = ("exon", "intron", "start_codon", "stop_codon")


@click.command("export")
@click.option("--assembly-id", type=int, required=True, help="Which assembly to write out.")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option(
    "--format",
    "output_format",
    type=click.Choice(["gff3", "tsv"]),
    default="gff3",
    help="gff3 writes the cgload dialect (D-052). tsv writes functional "
    "annotation, one assignment per row.",
)
@click.option(
    "--output",
    type=click.Path(dir_okay=False, path_type=Path),
    default=None,
    help="Where to write. Defaults to standard output.",
)
def export_command(
    assembly_id: int, db_url: str | None, output_format: str, output: Path | None
) -> None:
    """Write a loaded assembly back out."""
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL

    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        with engine.connect() as connection:
            if output_format == "gff3":
                text, summary = _gff3(connection, assembly_id)
                report = (
                    f"{summary.features:,} features in {summary.rows:,} rows"
                    + (
                        f", {summary.anonymous_features:,} without an identifier"
                        if summary.anonymous_features
                        else ""
                    )
                )
            else:
                text, count = _tsv(connection, assembly_id)
                report = f"{count:,} functional assignments"
    except LookupError as exc:
        click.echo(f"cgload export: {exc}", err=True)
        sys.exit(2)
    except (
        CredentialInArgumentError,
        DatabaseNotFoundError,
        SchemaMissingError,
        VocabularyMismatchError,
    ) as exc:
        click.echo(f"cgload export: {exc}", err=True)
        sys.exit(2)

    if output is None:
        click.echo(text, nl=False)
    else:
        output.write_text(text)
        click.echo(f"Wrote {report} to {output}", err=True)


def _gff3(connection, assembly_id: int) -> tuple[str, object]:
    context = load_context(connection, assembly_id)
    features = load_features(connection, assembly_id)
    if not features:
        raise LookupError(
            f"assembly {assembly_id} has no features stored, so there is nothing "
            f"to export."
        )
    present = {item.source_type for item in features}
    return export_gff3(
        features,
        version=context.cgload_version,
        vocabulary_hash=context.vocabulary_hash,
        declared_absent=[name for name in KNOWN_OPTIONAL if name not in present],
        sequence_regions=context.sequence_regions,
    )


def _tsv(connection, assembly_id: int) -> tuple[str, int]:
    """Functional annotation as TSV.

    A header comment records that GO counts include inherited ancestors (D-045);
    anyone summing this column needs to know that before they do.
    """
    from io import StringIO

    from cgload.db.schema import feature

    rows = connection.execute(
        select(
            feature.c.source_id,
            feature.c.source_type,
            protein_annotation.c.query_id,
            annotation_run.c.tool,
            protein_annotation.c.analysis,
            protein_annotation.c.accession,
            protein_annotation.c.description,
            protein_annotation.c.protein_start,
            protein_annotation.c.protein_end,
            protein_annotation.c.evalue,
            protein_annotation.c.score,
        )
        .select_from(
            protein_annotation.join(
                feature, feature.c.feature_id == protein_annotation.c.feature_id
            ).join(
                annotation_run,
                annotation_run.c.annotation_run_id
                == protein_annotation.c.annotation_run_id,
            )
        )
        .where(feature.c.assembly_id == assembly_id)
        .order_by(feature.c.source_id, protein_annotation.c.analysis,
                  protein_annotation.c.accession)
    ).all()

    buffer = StringIO()
    buffer.write(
        "# cgload functional annotation export\n"
        "# GO rows are ontology terms including inherited ancestors, not "
        "independent functions\n"
    )
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "feature_id", "feature_type", "query_id", "tool", "analysis",
            "accession", "description", "protein_start", "protein_end",
            "evalue", "score",
        ]
    )
    for row in rows:
        writer.writerow(["" if value is None else value for value in row])
    return buffer.getvalue(), len(rows)
