"""``cgload query`` -- read-only questions about one assembly (D-085)."""

from __future__ import annotations

import csv
import os
import sys
from io import StringIO
from pathlib import Path

import click

from cgload.db.engine import (
    CredentialInArgumentError,
    DatabaseNotFoundError,
    build_engine,
)
from cgload.db.export_source import load_context
from cgload.db.init import (
    SchemaMissingError,
    VocabularyMismatchError,
    assert_vocabulary_current,
)
from cgload.export import export_gff3
from cgload.query import QueryError, QueryResult, by_id, in_region, stored_analyses, with_term

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"

FEATURE_COLUMNS = [
    "id", "type", "feature_type", "seqid", "start", "end", "strand", "segments",
    "parent", "top_level",
]
TERM_COLUMNS = ["analysis", "accession", "description", "protein_start", "protein_end", "evalue"]


@click.command("query")
@click.option("--assembly-id", type=int, required=True, help="Which assembly to ask about.")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option("--id", "feature_id", default=None,
              help="A feature's ID, exactly as written in its file.")
@click.option("--region", default=None,
              help="chr, chr:start-end or chr:start (to the end of chr); 1-based, inclusive.")
@click.option("--term", default=None,
              help="An accession that eggNOG-mapper or InterProScan assigned: GO:0005524, "
              "K00799, PF00069, 2.7.11.1, IPR000719, ...")
@click.option("--analysis", default=None,
              help="With --term, only this analysis: GO, KEGG_ko, Pfam, EC, InterPro, ...")
@click.option("--descendants", is_flag=True,
              help="Also every feature below each match: a gene's transcripts, exons and CDS.")
@click.option("--type", "types", multiple=True,
              help="Keep only features of this canonical type, e.g. gene. Repeatable. "
              "Not with --term.")
@click.option("--format", "output_format", type=click.Choice(["tsv", "gff3"]), default="tsv",
              help="tsv: one row per feature (per assignment, with --term). gff3: the "
              "features plus their ancestors, so every Parent resolves within the file.")
@click.option("--output", type=click.Path(dir_okay=False, path_type=Path), default=None,
              help="Where to write. Defaults to standard output.")
def query_command(
    assembly_id: int,
    db_url: str | None,
    feature_id: str | None,
    region: str | None,
    term: str | None,
    analysis: str | None,
    descendants: bool,
    types: tuple[str, ...],
    output_format: str,
    output: Path | None,
) -> None:
    """Find features by ID, by region, or by functional term. Changes nothing.

    Answers are whole features: a feature is returned when any of its rows
    overlaps the region, and with all its rows. A summary goes to standard
    error, so standard output holds only the answer.
    """
    asked = [name for name, value in (("--id", feature_id), ("--region", region),
                                       ("--term", term)) if value is not None]
    if len(asked) != 1:
        _fail("give exactly one of --id, --region or --term")
    if analysis and term is None:
        _fail("--analysis applies only to --term")
    if types and term is not None:
        _fail("--type does not apply to --term. The annotated feature is the CDS or "
              "mRNA the functional file named; the top_level column names its gene.")

    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL
    gff3 = output_format == "gff3"
    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        assert_vocabulary_current(engine)
        with engine.connect() as connection:
            if feature_id is not None:
                result = by_id(connection, assembly_id, feature_id, descendants=descendants,
                               types=types, ancestors=gff3)
            elif region is not None:
                result = in_region(connection, assembly_id, region, descendants=descendants,
                                   types=types, ancestors=gff3)
            else:
                result = with_term(connection, assembly_id, term, analysis=analysis,
                                   descendants=descendants, ancestors=gff3)
            if gff3:
                text, summary = _gff3(connection, assembly_id, result)
            else:
                text, summary = _tsv(result, term=term is not None)
            if term is not None and not result.hits:
                present = stored_analyses(connection, assembly_id)
                summary += (". Stored analyses: " + ", ".join(present) if present
                            else ". No functional annotation is loaded for this assembly")
    except LookupError as exc:
        click.echo(f"cgload query: {exc}", err=True)
        sys.exit(1)
    except (
        CredentialInArgumentError,
        DatabaseNotFoundError,
        SchemaMissingError,
        VocabularyMismatchError,
        QueryError,
    ) as exc:
        _fail(str(exc))

    if output is None:
        click.echo(text, nl=False)
    else:
        output.write_text(text)
    click.echo(f"{summary}", err=True)


def _fail(message: str) -> None:
    click.echo(f"cgload query: {message}", err=True)
    sys.exit(2)


def _feature_row(item, result: QueryResult) -> list:
    return [
        item.source_id or "",
        item.source_type,
        item.feature_type,
        item.seqid,
        item.start,
        item.end,
        item.strand,
        len(item.segments),
        item.parent_source_id or "",
        result.top_level.get(item.row_id) or "",
    ]


def _tsv(result: QueryResult, *, term: bool) -> tuple[str, str]:
    buffer = StringIO()
    buffer.write(f"# cgload query: {result.description}\n")
    if any(hit.analysis == "GO" for hit in result.hits):
        buffer.write("# GO matches include terms eggNOG-mapper inherited from more specific "
                     "ones (D-045)\n")
    writer = csv.writer(buffer, delimiter="\t", lineterminator="\n")
    writer.writerow(FEATURE_COLUMNS + (TERM_COLUMNS if term else []))

    hits_by_row: dict[int, list] = {}
    for hit in result.hits:
        hits_by_row.setdefault(hit.row_id, []).append(hit)

    for item in result.features:
        base = _feature_row(item, result)
        if not term:
            writer.writerow(base)
            continue
        for hit in hits_by_row.get(item.row_id, [None]):
            extra = ["", "", "", "", "", ""] if hit is None else [
                hit.analysis, hit.accession, hit.description or "",
                "" if hit.protein_start is None else hit.protein_start,
                "" if hit.protein_end is None else hit.protein_end,
                "" if hit.evalue is None else f"{hit.evalue:g}",
            ]
            writer.writerow(base + extra)

    if term:
        annotated = len({hit.row_id for hit in result.hits})
        summary = f"{len(result.hits):,} assignment(s) on {annotated:,} feature(s)"
        if len(result.features) > annotated:
            summary += f", and {len(result.features) - annotated:,} feature(s) below them"
    else:
        summary = f"{len(result.features):,} feature(s)"
    return buffer.getvalue(), summary


def _gff3(connection, assembly_id: int, result: QueryResult) -> tuple[str, str]:
    context = load_context(connection, assembly_id)
    present = {item.seqid for item in result.features}
    text, written = export_gff3(
        result.features,
        version=context.cgload_version,
        vocabulary_hash=context.vocabulary_hash,
        sequence_regions={
            name: length for name, length in context.sequence_regions.items() if name in present
        },
        comments=[f"cgload query: {result.description}, with their ancestors"],
    )
    return text, f"{written.features:,} feature(s) in {written.rows:,} row(s)"
