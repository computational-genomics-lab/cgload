"""``cgload portal``: from genome files to a running portal, in one command (D-118).

It runs, in order, what the separate commands do, and skips what is done:

1. the database: created if there is none and genomes are given;
2. each genome given (an organism TOML, as for ``cgload load``): loaded, with
   every check of ``cgload load``, unless it is already in the database;
3. the search index of every genome that has none (``cgload index``);
4. BLAST+: found on the computer, or installed into cgload's own folder
   (``cgload.portal.install``), unless ``--no-blast``;
5. the BLAST databases of every genome that has none, or older ones
   (``cgload blastdb``);
6. the portal (``cgload serve``), read-only, with BLAST.

A step that fails says why in a sentence. A genome that does not load stops
the command before anything is served; a missing index or BLAST leaves the
portal running without them, and says so. Run it again at any time: it does
only what is not done yet. ``--prepare-only`` stops before step 6, for a
portal that a service starts with ``cgload serve``.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import click

from cgload.cli.commands.load import INPUT_ERRORS
from cgload.cli.commands.serve import (
    DEFAULT_DB_URL,
    ENV_DB_URL,
    blast_jobs_for,
    given,
    run_server,
    serving_options,
)
from cgload.db.engine import CredentialInArgumentError, DatabaseNotFoundError, build_engine
from cgload.db.init import (
    SchemaMissingError,
    VocabularyMismatchError,
    assert_vocabulary_current,
    create_schema,
)


def _say(text: str) -> None:
    click.echo(text, err=True)


def _plural(number: int, word: str) -> str:
    return f"{number:,} {word}{'' if number == 1 else 's'}"


def loaded_assembly(engine, config) -> int | None:
    """The assembly id under which this organism file's assembly is already
    loaded, or None: matched as ``cgload load`` matches it, by organism,
    strain, assembly name and version."""
    from sqlalchemy import select

    from cgload.db.schema import assembly, organism, strain

    version = assembly.c.version
    with engine.connect() as connection:
        return connection.execute(
            select(assembly.c.assembly_id)
            .join(strain, strain.c.strain_id == assembly.c.strain_id)
            .join(organism, organism.c.organism_id == strain.c.organism_id)
            .where(organism.c.genus == config.genus, organism.c.species == config.species,
                   strain.c.name == config.strain, assembly.c.name == config.assembly_name,
                   version.is_(None) if config.assembly_version is None
                   else version == config.assembly_version)
        ).scalar_one_or_none()


def _database(resolved: str, from_environment: bool, *, create: bool):
    """The database's engine; created first if it is missing and ``create``."""
    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        assert_vocabulary_current(engine)
        return engine
    except (DatabaseNotFoundError, SchemaMissingError) as exc:
        if not create:
            _say(f"cgload portal: {exc}. Give the genomes to load, as organism files: "
                 f"cgload portal organism.toml")
            sys.exit(2)
    engine = build_engine(resolved, from_environment=from_environment)
    result = create_schema(engine)
    _say(f"  database: created, {len(result.tables)} tables, at "
         f"{engine.url.render_as_string(hide_password=True)}")
    return engine


def _load(engine, paths: tuple[Path, ...], data_dir: Path) -> None:
    from cgload.config import load_config
    from cgload.loader import load_organism

    for path in paths:
        try:
            config = load_config(path)
            name = f"{config.genus} {config.species} / {config.strain} / {config.assembly_name}"
            existing = loaded_assembly(engine, config)
            if existing is not None:
                _say(f"  load: {name} is loaded already (assembly {existing})")
                continue
            _say(f"  load: {name} ...")
            report = load_organism(engine, config, data_dir)
        except INPUT_ERRORS as exc:
            _say(f"cgload portal: {path}: {exc}")
            _say("  Nothing is served until every genome given loads; the genomes loaded "
                 "before this one stay loaded.")
            sys.exit(2)
        checks = ", ".join(sorted({tier.name for tier in report.verification})) or "counts"
        _say(f"    assembly {report.assembly_id}: {_plural(report.region_count, 'sequence')}, "
             f"{_plural(report.feature_count, 'feature')}; every check passed ({checks})")
        for warning in report.warnings:
            _say(f"    warning: {warning}")


def _index(engine) -> None:
    from sqlalchemy import select

    from cgload.db.schema import assembly
    from cgload.portal.search import build_index, indexed_assemblies

    with engine.connect() as connection:
        every = connection.execute(select(assembly.c.assembly_id)).scalars().all()
        done = indexed_assemblies(connection)
    todo = [assembly_id for assembly_id in every if assembly_id not in done]
    if not todo:
        if every:
            _say(f"  search index: in place for {_plural(len(every), 'genome')}")
        return
    for report in build_index(engine, todo):
        if report.ok:
            _say(f"  search index: assembly {report.assembly_id}, "
                 f"{_plural(report.terms, 'term')}; every count matches")
        else:
            detail = "; ".join(f"{count.key} {count.in_attributes:,} against "
                               f"{count.in_index:,}" for count in report.failures())
            _say(f"  search index: assembly {report.assembly_id} failed its check ({detail}); "
                 f"names there are not searchable")


def _blast(resolved: str, from_environment: bool) -> str | None:
    """Make BLAST ready; None if it is, else why not, in a few words."""
    from cgload.portal.blast import build, databases_of, find_binaries
    from cgload.portal.data import PortalData, read_only_engine
    from cgload.portal.install import InstallError, install_blast

    binaries = find_binaries()
    if "makeblastdb" not in binaries or len(binaries) < 5:
        _say("  BLAST+: not on this computer; cgload installs it for itself, once "
             "(about 900 MB)")
        try:
            install_blast(log=_say)
        except InstallError as exc:
            _say(f"  BLAST+: {exc}")
            return "BLAST+ could not be installed (above)."
        binaries = find_binaries()
        if "makeblastdb" not in binaries:
            return "BLAST+ was installed, but its programs are not where cgload looks."
    portal = PortalData(read_only_engine(resolved, from_environment=from_environment))
    try:
        genomes = portal.genomes()
        current = databases_of(genomes)
        todo = [genome for genome in genomes if genome.assembly_id not in current]
        if not todo:
            if genomes:
                _say(f"  BLAST databases: in place for {_plural(len(genomes), 'genome')}")
            return None
        _say(f"  BLAST databases: building for {_plural(len(todo), 'genome')}")
        results = build(portal, todo, binaries["makeblastdb"],
                        log_line=lambda line: _say(f"    {line}"))
    finally:
        portal.close()
    if all(result.check != "PASS" for result in results) and not current:
        return "no genome's BLAST databases passed their check (above)."
    return None


@click.command("portal")
@click.argument("organisms", nargs=-1,
                type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--db-url", default=None,
              help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}, then {DEFAULT_DB_URL}.")
@click.option("--data-dir", default="data", show_default=True,
              type=click.Path(file_okay=False, path_type=Path),
              help="Where genomes being loaded go, as for cgload load.")
@click.option("--no-blast", is_flag=True,
              help="No BLAST: neither BLAST+ nor its databases are installed or built.")
@click.option("--prepare-only", is_flag=True,
              help="Do everything but start the portal: for a portal that a service "
              "starts with cgload serve.")
@serving_options
@click.pass_context
def portal_command(ctx: click.Context, organisms: tuple[Path, ...], db_url: str | None,
                   data_dir: Path, no_blast: bool, prepare_only: bool, host: str, port: int,
                   threads: int, title: str, url_prefix: str, blast_jobs: int,
                   blast_threads: int, blast_timeout: int) -> None:
    """Load genomes, index them, set up BLAST, and serve the portal: one command.

    ORGANISMS are organism files (TOML), as cgload load reads them; none to
    serve what the database holds. Each step is skipped when it is done
    already, so the same command starts the portal again later. BLAST+ is
    installed into cgload's own folder if the computer has none.
    """
    from_environment = db_url is None and ENV_DB_URL in os.environ
    resolved = db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL
    jobs_given = given(ctx, "blast_jobs")
    try:
        jobs = blast_jobs_for(threads, blast_jobs, given=jobs_given)
    except click.UsageError as exc:
        _say(f"cgload portal: {exc.message}")
        sys.exit(2)
    steps = [f"{_plural(len(organisms), 'genome')} to load" if organisms else None,
             "the search index", None if no_blast or not jobs else "BLAST",
             None if prepare_only else "the portal"]
    _say("cgload portal: " + ", then ".join(step for step in steps if step) + ".")
    try:
        engine = _database(resolved, from_environment, create=bool(organisms))
        _load(engine, organisms, data_dir)
        _index(engine)
    except (CredentialInArgumentError, VocabularyMismatchError) as exc:
        _say(f"cgload portal: {exc}")
        sys.exit(2)
    note = "--no-blast." if no_blast else None
    if not (no_blast or not jobs):
        note = _blast(resolved, from_environment)
    if prepare_only:
        _say("cgload portal: ready. Start the portal with: cgload serve")
        return
    run_server("cgload portal", db_url, None, host=host, port=port,
               port_given=given(ctx, "port"), threads=threads, title=title,
               url_prefix=url_prefix,
               blast_jobs=0 if no_blast else blast_jobs, jobs_given=jobs_given or no_blast,
               blast_threads=blast_threads, blast_timeout=blast_timeout, blast_note=note)
