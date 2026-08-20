"""``cgload verify`` -- prove the load lost nothing (milestone 8)."""

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
from cgload.db.init import SchemaMissingError, VocabularyMismatchError
from cgload.db.verification import (
    annotation_files,
    compare_with_database,
    count_database_for_file,
    sequences_loaded,
)
from cgload.verify import (
    VerifyError,
    compare_across_formats,
    compare_with_reference,
    count_source_naively,
    exit_code,
    parse_ncbi_feature_counts,
    render,
)

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"


@click.command("verify")
@click.option("--assembly-id", type=int, required=True, help="Which assembly to check.")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option(
    "--against",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="A second file describing the same assembly in another format, for the "
    "cross-parser check. GenBank against GFF3 is the strongest available.",
)
@click.option(
    "--reference-counts",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="NCBI's published feature_count file, for the reference check.",
)
@click.option(
    "--accession",
    default=None,
    help="Assembly accession to select within the reference-counts file. Required "
    "with --reference-counts: the file lists more than one, differing by a digit.",
)
@click.option(
    "--assembly-unit",
    default=None,
    help="Assembly unit to compare against within the reference-counts file, by "
    "accession or name (e.g. GCF_902167144.1, or 'Primary Assembly'). A file may "
    "carry an 'all' rollup plus one row-set per unit; the rollup is used by "
    "default. Name the unit when only part of an assembly was loaded, or the "
    "comparison is against totals that include sequences you did not load.",
)
@click.option(
    "--sequences-in-assembly",
    type=int,
    default=None,
    help="How many sequences the published assembly has. Without it the reference "
    "check cannot tell a partial load from a real discrepancy, and is skipped.",
)
def verify_command(
    assembly_id: int,
    db_url: str | None,
    against: Path | None,
    reference_counts: Path | None,
    accession: str | None,
    assembly_unit: str | None,
    sequences_in_assembly: int | None,
) -> None:
    """Count what was stored and compare it against what the files said."""
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL

    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        with engine.connect() as connection:
            results = run_verification(
                connection,
                assembly_id,
                against=against,
                reference_counts=reference_counts,
                accession=accession,
                assembly_unit=assembly_unit,
                sequences_in_assembly=sequences_in_assembly,
            )
    except (
        CredentialInArgumentError,
        DatabaseNotFoundError,
        SchemaMissingError,
        VerifyError,
        VocabularyMismatchError,
    ) as exc:
        click.echo(f"cgload verify: {exc}", err=True)
        sys.exit(2)

    click.echo(render(results))
    sys.exit(exit_code(results))


def run_verification(
    connection,
    assembly_id: int,
    *,
    against: Path | None = None,
    reference_counts: Path | None = None,
    accession: str | None = None,
    assembly_unit: str | None = None,
    sequences_in_assembly: int | None = None,
) -> list:
    """Run every applicable check and return one result per tier.

    Every tier is returned, including skipped ones. A check that vanishes when it
    cannot run is indistinguishable from one that passed (D-033).
    """
    files = annotation_files(connection, assembly_id)
    if not files:
        raise VerifyError(
            f"assembly {assembly_id} has no annotation files recorded, so there is "
            f"nothing to check it against."
        )

    results = []
    combined = None
    for source_file_id, path in files:
        source_path = Path(path)
        if not source_path.exists():
            results.append(
                _skipped(
                    f"source count [{source_path.name}]",
                    f"the input file is no longer at {path}; the check needs the "
                    f"original to count independently",
                )
            )
            continue
        source_counts = count_source_naively(source_path)
        stored = count_database_for_file(connection, assembly_id, source_file_id)
        label = (
            "source count"
            if len(files) == 1
            else f"source count [{source_path.name}]"
        )
        results.append(compare_with_database(source_counts, stored, label=label))
        combined = source_counts if combined is None else combined + source_counts

    if against is not None:
        primary = Path(files[0][1])
        if primary.exists():
            results.append(compare_across_formats(primary, against))
        else:
            results.append(
                _skipped("cross-parser", "the originally loaded file is no longer present")
            )
    else:
        results.append(
            _skipped(
                "cross-parser",
                "no second format supplied; pass --against with a GenBank or GFF3 "
                "file describing the same assembly",
            )
        )

    if reference_counts is not None and combined is not None:
        if accession is None:
            raise VerifyError(
                "--reference-counts needs --accession: the file lists several "
                "assemblies whose accessions differ by one digit, and selecting the "
                "wrong one would report every feature as missing."
            )
        published = parse_ncbi_feature_counts(
            reference_counts, accession, assembly_unit=assembly_unit
        )
        loaded = sequences_loaded(connection, assembly_id)
        if sequences_in_assembly is None:
            # Never default this to `loaded`. The precondition becomes
            # `loaded == loaded`, always true, so the tier can no longer skip --
            # a partial load is then compared against whole-assembly totals and
            # either passes by luck or fails for a reason cgload already
            # classifies as innocent. A precondition that cannot fail is D-033
            # inverted (D-067).
            raise VerifyError(
                "--reference-counts needs --sequences-in-assembly: without it the "
                "partial-load precondition cannot be evaluated, and the tier would "
                f"compare the {loaded} sequence(s) you loaded against totals that may "
                f"describe the whole assembly. Pass the published sequence count for "
                f"the unit you are comparing against, or pass {loaded} to state "
                f"deliberately that the counts describe exactly what you loaded."
            )
        results.append(
            compare_with_reference(
                combined,
                published,
                sequences_loaded=loaded,
                sequences_in_assembly=sequences_in_assembly,
            )
        )
    else:
        results.append(
            _skipped(
                "reference (NCBI)",
                "no published counts supplied; pass --reference-counts and --accession",
            )
        )

    return results


def _skipped(name: str, detail: str):
    from cgload.verify import Outcome, TierResult

    return TierResult(name=name, outcome=Outcome.SKIPPED, detail=detail)
