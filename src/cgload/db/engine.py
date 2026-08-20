"""Engine construction and the credential rule of D-003.

Credentials never appear in a CLI argument or a tracked file. ``--db-url`` is
accepted for SQLite paths and for server URLs that carry no password; a URL
with an inline password is refused, because it lands in shell history and in
the process table. Server credentials come from ``CGLOAD_DB_URL`` or from a
gitignored config file.
"""

from __future__ import annotations

from pathlib import Path

from sqlalchemy import Engine, create_engine, event
from sqlalchemy.engine import URL, make_url


class CredentialInArgumentError(ValueError):
    """A password was passed where D-003 forbids one."""


class DatabaseNotFoundError(FileNotFoundError):
    """The SQLite file does not exist.

    SQLite creates a database file on first connect, so a mistyped path
    produces an empty database rather than an error, and the failure surfaces
    later as a confusing missing-table traceback. Every command except `init`
    therefore requires the file to exist already.
    """


def validated_url(db_url: str, *, from_environment: bool = False) -> URL:
    """Apply D-003 to a URL. Separate from engine construction so the rule can
    be tested without the server driver being installed."""
    url = make_url(db_url)
    if url.password and not from_environment:
        raise CredentialInArgumentError(
            "the database URL contains a password. Set CGLOAD_DB_URL in the "
            "environment instead; a password in an argument is recorded in "
            "shell history and visible in the process table (D-003)."
        )
    return url


def build_engine(
    db_url: str, *, from_environment: bool = False, must_exist: bool = False
) -> Engine:
    url = validated_url(db_url, from_environment=from_environment)
    if must_exist and url.get_backend_name() == "sqlite" and url.database not in (None, ":memory:"):
        path = Path(url.database)
        if not path.exists():
            raise DatabaseNotFoundError(
                f"no database at {path.resolve()}. Run `cgload init` first, or point "
                f"--db-url at the right file. (SQLite would otherwise create an empty "
                f"one here and fail later with a confusing missing-table error.)"
            )
    engine = create_engine(url, future=True)
    if engine.dialect.name == "sqlite":
        # Off by default in SQLite, so every ForeignKey in the schema would be
        # decorative rather than enforced -- and the loader's correctness
        # claims rest on them.
        @event.listens_for(engine, "connect")
        def _fk_on(dbapi_connection, _record):  # pragma: no cover - trivial
            cur = dbapi_connection.cursor()
            cur.execute("PRAGMA foreign_keys=ON")
            cur.close()

    return engine
