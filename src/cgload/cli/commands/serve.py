"""``cgload serve`` and ``cgload index``: the portal (D-111, D-112)."""

from __future__ import annotations

import errno
import ipaddress
import logging
import os
import socket
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


TRIES = 10  # ports tried after the default one, when it is taken


class PortalStartError(Exception):
    """A reason the portal could not start, in words for its user."""


def _reach(host: str) -> str:
    """"here" for this computer only, "everywhere" for every interface, or "named"."""
    if host == "localhost":
        return "here"
    try:
        address = ipaddress.ip_address(host)
    except ValueError:
        return "named"
    return "here" if address.is_loopback else "everywhere" if address.is_unspecified else "named"


def network_addresses() -> list[str]:
    """This computer's IPv4 addresses that other computers may use, best first.

    Read without sending anything: a UDP socket "connected" to a private and a
    public address reports the local address the system would use, and the
    host name's own addresses are added. Loopback addresses are left out.
    """
    found: list[str] = []
    for target in ("10.255.255.255", "192.0.2.1"):  # private; TEST-NET-1 (never routed)
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_DGRAM) as probe:
                probe.connect((target, 9))
                found.append(probe.getsockname()[0])
        except OSError:
            continue
    try:
        for *_, address in socket.getaddrinfo(socket.gethostname(), None, socket.AF_INET):
            found.append(address[0])
    except OSError:
        pass
    usable = []
    for address in found:
        try:
            parsed = ipaddress.ip_address(address)
        except ValueError:
            continue
        if not (parsed.is_loopback or parsed.is_unspecified or parsed.is_link_local) \
                and address not in usable:
            usable.append(address)
    return usable[:3]


def where_it_listens(host: str, port: int, prefix: str = "",
                     addresses: list[str] | None = None) -> list[str]:
    """What ``cgload serve`` tells its user about where the portal opens."""
    path = prefix or "/"

    def url(name: str) -> str:
        return f"http://{f'[{name}]' if ':' in name else name}:{port}{path}"

    reach = _reach(host)
    if reach == "here":
        return [f"  Open it on this computer:   {url(host)}",
                "  Only this computer can open it. For the other computers on your "
                "network, leave out --host.",
                "  Stop it with Ctrl-C."]
    if reach == "named":
        return [f"  Open it at:   {url(host)}", "  Stop it with Ctrl-C."]
    others = addresses or ["<this computer's address>"]
    lines = [f"  Open it on this computer:   {url('::1' if ':' in host else '127.0.0.1')}"]
    lines += [f"  or from another computer:   {url(name)}" for name in others]
    lines += ["  Anyone who can reach this computer can open it. To keep it to this "
              "computer: --host 127.0.0.1;",
              "  to serve it beyond your network, put a web server with HTTPS in front "
              "(docs/portal.md).",
              "  Stop it with Ctrl-C."]
    return lines


def why_not(exc: OSError, host: str, port: int) -> str:
    """One sentence for an error from opening the port."""
    code = exc.errno
    if code == errno.EADDRINUSE or getattr(exc, "winerror", None) == 10048:
        return (f"port {port} is already in use on this computer, perhaps by another "
                f"cgload serve. Stop that program, or choose another port: "
                f"cgload serve --port {port + 1 if port < 65535 else 8000}")
    if code == errno.EACCES or getattr(exc, "winerror", None) == 10013:
        return (f"port {port} needs administrator rights on this computer. Use a port "
                f"above 1023 (the default is 8000), or put a web server in front "
                f"(docs/portal.md).")
    if code == errno.EADDRNOTAVAIL or getattr(exc, "winerror", None) == 10049:
        return (f"this computer has no address {host}. Leave out --host to listen on "
                f"all of its addresses.")
    if isinstance(exc, socket.gaierror):
        return f"{host} is not an address or a name this computer can look up."
    return f"could not listen on {host}:{port}: {exc.strerror or exc}."


def open_server(app, host: str, port: int, *, port_given: bool, threads: int,
                url_prefix: str):
    """Open the portal's port, or say why not. When the default port is taken,
    the next free one of the following ten is used; a port the user chose is
    used or refused, never changed."""
    from waitress import create_server

    candidates = [port] if port_given else range(port, min(port + TRIES, 65535) + 1)
    for candidate in candidates:
        try:
            server = create_server(app, host=host, port=candidate, threads=threads,
                                   ident="cgload", url_prefix=url_prefix,
                                   clear_untrusted_proxy_headers=True)
        except OSError as exc:
            in_use = exc.errno == errno.EADDRINUSE or getattr(exc, "winerror", None) == 10048
            if in_use and not port_given and candidate != candidates[-1]:
                continue
            if in_use and not port_given:
                raise PortalStartError(
                    f"ports {port} to {candidate} are all in use on this computer. "
                    f"Choose a free one: cgload serve --port 9000") from None
            raise PortalStartError(why_not(exc, host, candidate)) from None
        return server, candidate
    raise PortalStartError(f"no port to try from {port}")  # pragma: no cover


def _resolve(db_url: str | None) -> tuple[str, bool]:
    from_environment = db_url is None and ENV_DB_URL in os.environ
    return db_url or os.environ.get(ENV_DB_URL) or DEFAULT_DB_URL, from_environment


@click.command("serve")
@click.option("--db-url", default=None, help=f"SQLAlchemy URL; defaults to ${ENV_DB_URL}.")
@click.option("--data-dir", type=click.Path(file_okay=False, path_type=Path), default=None,
              help="The data directory the assemblies were loaded into, if it has moved "
              "since the load. The portal reads assemblies and annotation files from it.")
# Every interface by default, as a lab server needs (D-115); the start-up lines say who
# can then open the portal, and --host 127.0.0.1 keeps it to this computer.
@click.option("--host", default="0.0.0.0", show_default=True,  # noqa: S104  # nosec B104
              help="Address to listen on. The default lets every computer that can reach "
              "this one open the portal; 127.0.0.1 keeps it to this computer, as behind "
              "a web server.")
@click.option("--port", type=click.IntRange(1, 65535), default=8000, show_default=True,
              help="If the default is taken, the next free port is used.")
@click.option("--threads", type=click.IntRange(1, 64), default=8, show_default=True,
              help="Requests served at once.")
@click.option("--title", default="cgload portal", show_default=True,
              help="The portal's name, shown on every page.")
@click.option("--url-prefix", default="",
              help="The path the portal is served under, behind a web server: /genomes.")
@click.pass_context
def serve_command(ctx: click.Context, db_url: str | None, data_dir: Path | None, host: str,
                  port: int, threads: int, title: str, url_prefix: str) -> None:
    """Serve the database as a read-only web portal.

    People search genomes by name, open genes with their transcripts drawn to
    scale, browse regions in a genome browser and download sequences; programs
    get the same answers as JSON under /api/v1/. The database is opened
    read-only and nothing is ever written.
    """
    try:
        import waitress  # noqa: F401

        from cgload.portal.data import PortalData, read_only_engine
        from cgload.portal.web import create_app
    except ImportError as exc:
        click.echo(f"cgload serve: the portal needs Flask and waitress, which install with "
                   f"cgload, but {exc.name} is missing here. Install them with: "
                   f"pip install \"flask>=3.0\" \"waitress>=3.0\"", err=True)
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
    app = create_app(portal, title=title)
    port_given = ctx.get_parameter_source("port") is not click.core.ParameterSource.DEFAULT
    try:
        server, port_used = open_server(app, host, port, port_given=port_given,
                                        threads=threads, url_prefix=url_prefix or "")
    except PortalStartError as exc:
        click.echo(f"cgload serve: {exc}", err=True)
        portal.close()
        sys.exit(1)

    count = f"{len(genomes)} genome{'' if len(genomes) == 1 else 's'}"
    click.echo(f"cgload serve: the portal is running, read-only, with {count}.", err=True)
    if port_used != port:
        click.echo(f"  Port {port} was in use, so it is on port {port_used}.", err=True)
    addresses = network_addresses() if _reach(host) == "everywhere" else []
    for line in where_it_listens(host, port_used, url_prefix or "", addresses):
        click.echo(line, err=True)
    unindexed = [genome for genome in genomes if genome.indexed is None]
    missing = [genome for genome in genomes if genome.fasta is None]
    if not genomes:
        click.echo("  The database has no genomes yet: load one with `cgload load`, then "
                   "run `cgload index`.", err=True)
    if unindexed:
        click.echo(f"  {len(unindexed)} without a search index: names are not searchable "
                   f"there until `cgload index` is run.", err=True)
    if missing:
        click.echo(f"  {len(missing)} whose assembly file is not found: give --data-dir "
                   f"if the data directory has moved.", err=True)
    logging.basicConfig()  # waitress reports a long queue of requests as a warning
    server.run()  # returns on Ctrl-C, after the requests in hand are answered
    portal.close()
    click.echo("cgload serve: stopped.", err=True)


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
