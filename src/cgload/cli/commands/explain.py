"""``cgload explain`` -- why is this feature the way it is (milestone 8).

The companion to ``verify``. Where ``verify`` answers *did we lose anything*,
``explain`` answers *why is this particular thing here* — which file it came
from, which dialect that file was recognised as, how its parent was determined,
and whether that determination was stated by the file or computed by cgload.

It recomputes nothing. Every value it prints was recorded at load time, which is
what makes it an audit rather than a second opinion: a recomputed answer could
differ from the stored one and neither would be checkable.
"""

from __future__ import annotations

import os
import sys

import click
from sqlalchemy import select

from cgload.db.engine import (
    CredentialInArgumentError,
    DatabaseNotFoundError,
    build_engine,
)
from cgload.db.init import SchemaMissingError, VocabularyMismatchError
from cgload.db.schema import (
    annotation_run,
    assembly,
    feature,
    feature_attribute,
    organism,
    protein_annotation,
    sequence_region,
    source_file,
    strain,
    vocabulary_meta,
)
from cgload.db.vocabulary import COMPUTED_LINKAGE

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"

#: How each linkage method was arrived at, in words. `explain` exists to make
#: this distinction visible: the difference between a parent the file stated and
#: one cgload worked out is the whole substance of D-007 and of the paper.
LINKAGE_PROSE: dict[str, str] = {
    "explicit_parent": "the file stated it (GFF3 Parent attribute)",
    "locus_tag": "the file stated it (parent and child share a locus tag)",
    "qualifier": "the file stated it (/transcript_id on both features)",
    "note": "cgload read it from the /note text naming the parent transcript",
    "structural": "cgload computed it: the only transcript containing this feature",
    "positional": "cgload computed it from order of appearance",
    "synthesised": "cgload created the parent, which the file omitted",
}


@click.command("explain")
@click.argument("identifier")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option(
    "--assembly-id",
    type=int,
    default=None,
    help="Restrict the search. Needed when one identifier occurs in several "
    "assemblies, which cgload refuses to disambiguate for you.",
)
def explain_command(identifier: str, db_url: str | None, assembly_id: int | None) -> None:
    """Show where a feature came from and how its parent was determined."""
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL

    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        with engine.connect() as connection:
            output = _explain(connection, identifier, assembly_id)
    except (
        CredentialInArgumentError,
        DatabaseNotFoundError,
        SchemaMissingError,
        VocabularyMismatchError,
    ) as exc:
        click.echo(f"cgload explain: {exc}", err=True)
        sys.exit(2)

    if output is None:
        click.echo(
            f"cgload explain: no feature with source id {identifier!r}"
            + (f" in assembly {assembly_id}" if assembly_id else ""),
            err=True,
        )
        sys.exit(2)

    click.echo(output)


def _explain(connection, identifier: str, assembly_id: int | None) -> str | None:
    condition = [feature.c.source_id == identifier]
    if assembly_id is not None:
        condition.append(feature.c.assembly_id == assembly_id)

    segments = connection.execute(
        select(feature).where(*condition).order_by(feature.c.segment_index)
    ).all()
    if not segments:
        return None

    assemblies = {row.assembly_id for row in segments}
    if len(assemblies) > 1:
        return (
            f"{identifier!r} occurs in {len(assemblies)} assemblies "
            f"({', '.join(str(a) for a in sorted(assemblies))}). Pass --assembly-id; "
            f"cgload will not choose one for you."
        )

    first = segments[0]
    lines: list[str] = []

    where = connection.execute(
        select(
            organism.c.genus,
            organism.c.species,
            strain.c.name.label("strain"),
            assembly.c.name.label("assembly"),
            assembly.c.vocabulary_version,
            assembly.c.cgload_version,
        )
        .select_from(assembly.join(strain).join(organism))
        .where(assembly.c.assembly_id == first.assembly_id)
    ).one()
    region = connection.execute(
        select(sequence_region.c.source_id).where(
            sequence_region.c.sequence_region_id == first.sequence_region_id
        )
    ).scalar_one()
    origin = connection.execute(
        select(
            source_file.c.path,
            source_file.c.detected_dialect,
            source_file.c.detection_confidence,
            source_file.c.id_prefix,
            source_file.c.sha256,
        ).where(source_file.c.source_file_id == first.source_file_id)
    ).one()

    lines.append(f"{first.source_type} {identifier}")
    lines.append(
        f"  organism    {where.genus} {where.species} / {where.strain} / {where.assembly}"
    )
    lines.append(
        f"  location    {region}:{min(s.start for s in segments):,}-"
        f"{max(s.end for s in segments):,} ({first.strand})"
    )
    if len(segments) > 1:
        spans = ", ".join(f"{s.start:,}-{s.end:,}" for s in segments)
        lines.append(f"  segments    {len(segments)} in coordinate order: {spans}")

    lines.append(f"  type        {first.source_type} in the file, stored canonically "
                 f"as {first.feature_type}")
    lines.append(f"  from file   {origin.path}")
    lines.append(
        f"              recognised as {origin.detected_dialect} "
        f"({origin.detection_confidence} confidence)"
        + (f", id_prefix {origin.id_prefix!r}" if origin.id_prefix else "")
    )
    lines.append(f"              sha256 {origin.sha256[:16]}...")

    # The point of the command.
    if first.parent_id is None:
        lines.append("  parent      none (this is a top-level feature)")
    else:
        parent = connection.execute(
            select(feature.c.source_id, feature.c.source_type).where(
                feature.c.feature_id == first.parent_id
            )
        ).one()
        method = first.linkage_method or "unrecorded"
        prose = LINKAGE_PROSE.get(method, "method not in the recorded vocabulary")
        computed = method in COMPUTED_LINKAGE
        lines.append(
            f"  parent      {parent.source_type} {parent.source_id}"
        )
        lines.append(f"  linked by   {method} — {prose}")
        lines.append(
            "              this link was computed by cgload, not read from the file"
            if computed
            else "              this link was stated by the file"
        )

    attributes = connection.execute(
        select(feature_attribute.c.key, feature_attribute.c.value)
        .where(feature_attribute.c.feature_id == first.feature_id)
        .order_by(feature_attribute.c.key)
    ).all()
    if attributes:
        lines.append("  attributes")
        for key, value in attributes[:12]:
            shown = value if len(value) <= 70 else value[:67] + "..."
            lines.append(f"              {key} = {shown}")
        if len(attributes) > 12:
            lines.append(f"              ... and {len(attributes) - 12} more")

    functional = connection.execute(
        select(
            protein_annotation.c.analysis,
            protein_annotation.c.accession,
            annotation_run.c.tool,
        )
        .select_from(
            protein_annotation.join(
                annotation_run,
                annotation_run.c.annotation_run_id
                == protein_annotation.c.annotation_run_id,
            )
        )
        .where(protein_annotation.c.feature_id == first.feature_id)
    ).all()
    if functional:
        by_analysis: dict[str, int] = {}
        for row in functional:
            by_analysis[row.analysis] = by_analysis.get(row.analysis, 0) + 1
        shown = ", ".join(f"{name} {count}" for name, count in sorted(by_analysis.items()))
        tools = ", ".join(sorted({row.tool for row in functional}))
        lines.append(f"  functional  {len(functional)} assignments from {tools}: {shown}")
        if "GO" in by_analysis:
            lines.append(
                "              GO counts are ontology terms including inherited "
                "ancestors, not independent functions (D-045)"
            )

    version = connection.execute(select(vocabulary_meta.c.version)).scalar_one()
    lines.append(
        f"  recorded by cgload {where.cgload_version} under vocabulary version {version}"
    )
    return "\n".join(lines)
