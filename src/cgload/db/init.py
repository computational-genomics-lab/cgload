"""Create the schema and seed the vocabulary."""

from __future__ import annotations

from dataclasses import dataclass

from sqlalchemy import Engine, insert, inspect, select

from cgload import __version__
from cgload.db import vocabulary as vocab
from cgload.db.schema import metadata, vocabulary, vocabulary_meta


class SchemaAlreadyPresentError(RuntimeError):
    """init would overwrite an existing database."""


class VocabularyMismatchError(RuntimeError):
    """The stored vocabulary is not the one this cgload would apply."""


class SchemaMissingError(RuntimeError):
    """The database has no cgload schema.

    Raised in place of the driver's missing-table error: a user who ran a
    command in the wrong directory needs to be told to run `init`, not handed
    a SQLAlchemy traceback.
    """


@dataclass(frozen=True)
class InitResult:
    tables: tuple[str, ...]
    vocabulary_terms: int
    vocabulary_version: int
    vocabulary_hash: str


def existing_tables(engine: Engine) -> set[str]:
    return set(inspect(engine).get_table_names())


def create_schema(engine: Engine, *, force: bool = False) -> InitResult:
    """Create every table and seed the vocabulary.

    Refuses a database that already holds cgload tables unless ``force``.
    P1: overwriting a populated database on a mistyped path is exactly the
    silent destructive behaviour this project rejects.
    """
    present = existing_tables(engine) & set(metadata.tables)
    if present and not force:
        raise SchemaAlreadyPresentError(
            f"{len(present)} cgload tables already exist. Point --db-url at a "
            f"new database, or pass --force to recreate (this drops data)."
        )
    if present and force:
        metadata.drop_all(engine)

    metadata.create_all(engine)

    rows = vocab.seed_rows()
    with engine.begin() as conn:
        conn.execute(insert(vocabulary), rows)
        conn.execute(
            insert(vocabulary_meta).values(
                vocabulary_meta_id=1,
                version=vocab.VOCABULARY_VERSION,
                content_hash=vocab.content_hash(),
            )
        )

    return InitResult(
        tables=tuple(sorted(metadata.tables)),
        vocabulary_terms=len(rows),
        vocabulary_version=vocab.VOCABULARY_VERSION,
        vocabulary_hash=vocab.content_hash(),
    )


def assert_vocabulary_current(engine: Engine) -> None:
    """Guard every later command. Called by load, verify, query and export.

    A vocabulary row insert changes what ``feature_type`` means. Rows already
    stored are not re-mapped, so continuing against a changed vocabulary would
    silently mix two linkage semantics in one database.
    """
    if "vocabulary_meta" not in existing_tables(engine):
        where = engine.url.render_as_string(hide_password=True)
        raise SchemaMissingError(
            f"{where} has no cgload schema. Run `cgload init` first -- and check "
            f"you are in the directory where you ran it, since the default database "
            f"is ./cgload.db."
        )

    with engine.connect() as conn:
        row = conn.execute(
            select(vocabulary_meta.c.version, vocabulary_meta.c.content_hash)
        ).first()
    if row is None:
        raise SchemaMissingError(
            "the vocabulary table is empty; this database was not created by "
            "`cgload init`. Recreate it with `cgload init --force`."
        )
    if row.content_hash != vocab.content_hash():
        raise VocabularyMismatchError(
            f"this database was created under vocabulary version {row.version} "
            f"({row.content_hash[:12]}), and cgload {__version__} would apply "
            f"version {vocab.VOCABULARY_VERSION} ({vocab.content_hash()[:12]}). "
            f"Run `cgload remap` to rewrite feature_type across affected rows, "
            f"or use a cgload whose vocabulary matches."
        )
