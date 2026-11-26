"""``cgload serve`` and ``cgload index``: the portal (D-111, D-112)."""

from __future__ import annotations

import ipaddress
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

DEFAULT_DB_URL = "sqlite:///cgload.db"
ENV_DB_URL = "CGLOAD_DB_URL"


def _reach(host: str) -> str:
    """"here" for this computer only, "everywhere" for every interface, or "named"."""
    if host == "localhost":
        return "here"
    try:
        address = ipaddress.ip_address(host)
    except ValueError:
        return "named"
    return "here" if address.is_loopback else "everywhere" if address.is_unspecified else "named"


def where_it_listens(host: str, port: int, prefix: str = "") -> list[str]:
    """What ``cgload serve`` tells its user about who can open the portal."""
    path = prefix or "/"
    def url(name: str) -> str:
        return f"http://{f'[{name}]' if ':' in name else name}:{port}{path}"
    elsewhere = url("<this computer's address>")
    beyond = ("  Listening beyond this computer: put a web server with HTTPS in front "
              "(docs/portal.md).")
    reach = _reach(host)
    if reach == "here":
        return [f"on {url(host)}, for this computer only.",
                f"  For the other computers on your network: cgload serve --host 0.0.0.0, "
                f"then open {elsewhere} on them."]
    if reach == "everywhere":
        here = url("::1" if ":" in host else "127.0.0.1")
        return [f"on port {port} of every network interface: {elsewhere} from other "
                f"computers, {here} from this one.", beyond]
    return [f"on {url(host)}", beyond]


def _resolve(db_url: str | None) -> tuple[str, bool]:
    from_environment = db_url is None and ENV_DB_URL in os.environ
    return db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL, from_environment


@click.command("serve")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option("--data-dir", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="The data directory the assemblies were loaded into, if it has moved "
              "since the load. The portal reads assemblies and annotation files from it.")
@click.option("--host", default="127.0.0.1", show_default=True,
              help="Address to listen on. 127.0.0.1 serves this computer only; 0.0.0.0 "
              "serves every network interface, best behind a web server that adds HTTPS.")
@click.option("--port", type=click.IntRange(1, 65535), default=8000, show_default=True)
@click.option("--threads", type=click.IntRange(1, 64), default=8, show_default=True,
              help="Requests served at once.")
@click.option("--title", default="cgload portal", show_default=True,
              help="The portal's name, shown on every page.")
@click.option("--url-prefix", default="",
              help="The path the portal is served under, behind a web server: /genomes.")
def serve_command(db_url: str | None, data_dir: Path | None, host: str, port: int,
                  threads: int, title: str, url_prefix: str) -> None:
    """Serve the database as a read-only web portal.

    People search genomes by name, open genes with their transcripts drawn to
    scale, browse regions in a genome browser and download sequences; programs
    get the same answers as JSON under /api/v1/. The database is opened
    read-only and nothing is ever written.
    """
    try:
        from waitress import serve

        from cgload.portal.data import PortalData, read_only_engine
        from cgload.portal.web import create_app
    except ImportError as exc:
        click.echo(f"cgload serve: the portal needs its extra packages ({exc.name} is "
                   f"missing). Install them with: pip install \"cgload[portal]\"", err=True)
        sys.exit(2)

    resolved, from_environment = _resolve(db_url)
    try:
        assert_vocabulary_current(build_engine(resolved, from_environment=from_environment,
                                               must_exist=True))
        engine = read_only_engine(resolved, from_environment=from_environment)
        portal = PortalData(engine, data_dir=data_dir)
    except (CredentialInArgumentError, DatabaseNotFoundError, SchemaMissingError,
            VocabularyMismatchError, FileNotFoundError, ValueError) as exc:
        click.echo(f"cgload serve: {exc}", err=True)
        sys.exit(2)

    genomes = portal.genomes()
    unindexed = [genome for genome in genomes if genome.indexed is None]
    missing = [genome for genome in genomes if genome.fasta is None]
    where = where_it_listens(host, port, url_prefix or "")
    click.echo(f"cgload serve: {len(genomes)} assembl{'y' if len(genomes) == 1 else 'ies'}, "
               f"read-only, {where[0]}", err=True)
    for line in where[1:]:
        click.echo(line, err=True)
    if unindexed:
        click.echo(f"  {len(unindexed)} without a search index: names are not searchable "
                   f"there until `cgload index` is run.", err=True)
    if missing:
        click.echo(f"  {len(missing)} whose assembly file is not found: give --data-dir "
                   f"if the data directory has moved.", err=True)
    app = create_app(portal, title=title)
    serve(app, host=host, port=port, threads=threads, ident="cgload",
          url_prefix=url_prefix or "", clear_untrusted_proxy_headers=True)


@click.command("index")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option("--assembly-id", "assembly_ids", type=int, multiple=True,
              help="Index this assembly only. Repeatable; every assembly by default.")
def index_command(db_url: str | None, assembly_ids: tuple[int, ...]) -> None:
    """Build the portal's search index of names, and check it by counting.

    For each assembly, every ID, name, synonym, locus tag, protein and
    transcript ID, cross-reference and word of a product goes into one
    indexed table. The index is then checked as a load is: for each attribute,
    the features that carry it are counted in the database and in the index,
    and an assembly whose counts differ is left without an index.
    """
    from cgload.portal.search import build_index

    resolved, from_environment = _resolve(db_url)
    try:
        engine = build_engine(resolved, from_environment=from_environment, must_exist=True)
        assert_vocabulary_current(engine)
        reports = build_index(engine, assembly_ids or None)
    except (CredentialInArgumentError, DatabaseNotFoundError, SchemaMissingError,
            VocabularyMismatchError) as exc:
        click.echo(f"cgload index: {exc}", err=True)
        sys.exit(2)
    except LookupError as exc:
        click.echo(f"cgload index: {exc}", err=True)
        sys.exit(1)
    failed = 0
    for report in reports:
        if report.ok:
            click.echo(f"assembly {report.assembly_id}: {report.features:,} features, "
                       f"{report.terms:,} terms; every key's count matches: PASS")
            continue
        failed += 1
        detail = "; ".join(f"{count.key}: {count.in_attributes:,} features carry it, "
                           f"{count.in_index:,} are indexed" for count in report.failures())
        click.echo(f"assembly {report.assembly_id}: FAIL, no index kept. {detail}", err=True)
    sys.exit(1 if failed else 0)
