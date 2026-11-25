"""The portal on PostgreSQL (D-111, D-112).

Two things differ from SQLite and are what these check: PostgreSQL compares
text by the database's collation, which would make a prefix search written as a
range miss names unless the term column is byte-ordered (``COLLATE "C"``); and
the portal's read-only engine relies on ``default_transaction_read_only``.

Skipped unless CGLOAD_TEST_POSTGRES_URL names a server URL without a password.
CI runs it against a PostgreSQL service container.
"""

from __future__ import annotations

import os
import uuid

import pytest

pytest.importorskip("flask")

from cgload.db.engine import build_engine  # noqa: E402
from cgload.db.schema import feature, feature_attribute  # noqa: E402
from cgload.portal.data import PortalData, read_only_engine  # noqa: E402
from cgload.portal.search import NAME_KEYS, PART_TYPES, check_index, search  # noqa: E402
from cgload.portal.web import create_app  # noqa: E402
from sqlalchemy import create_engine, select, text  # noqa: E402
from sqlalchemy.engine import make_url  # noqa: E402
from sqlalchemy.exc import DBAPIError  # noqa: E402

from portal_db import build  # noqa: E402

pytestmark = pytest.mark.requires_server
ADMIN_URL = os.environ.get("CGLOAD_TEST_POSTGRES_URL")


@pytest.fixture(scope="module")
def pg(tmp_path_factory):
    if not ADMIN_URL:
        pytest.skip("set CGLOAD_TEST_POSTGRES_URL to run the PostgreSQL tests")
    name = f"cgload_portal_{uuid.uuid4().hex[:12]}"
    admin = create_engine(ADMIN_URL, isolation_level="AUTOCOMMIT")
    with admin.connect() as connection:
        connection.execute(text(f'CREATE DATABASE "{name}"'))
    url = make_url(ADMIN_URL).set(database=name).render_as_string(hide_password=False)
    try:
        build(tmp_path_factory.mktemp("portal-pg"), url)
        yield url
    finally:
        with admin.connect() as connection:
            connection.execute(text(f'DROP DATABASE IF EXISTS "{name}" WITH (FORCE)'))
        admin.dispose()


def test_index_passes_its_count_check(pg):
    engine = build_engine(pg)
    with engine.connect() as connection:
        for assembly_id in (1, 2):
            assert check_index(connection, assembly_id).ok


def test_every_name_finds_its_feature(pg):
    engine = build_engine(pg)
    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.assembly_id, feature.c.feature_id, feature_attribute.c.value)
            .join(feature, feature.c.feature_id == feature_attribute.c.feature_id)
            .where(feature.c.assembly_id.in_([1, 2]), feature.c.segment_index == 0,
                   feature.c.feature_type.not_in(PART_TYPES),
                   feature_attribute.c.key.in_(NAME_KEYS))).all()
        missed = [row.value for row in rows if row.feature_id not in {
            match.feature_id for match in search(connection, row.assembly_id, row.value.upper(),
                                                 limit=5000).matches}]
    assert rows and missed == []


def test_a_prefix_finds_every_name_that_starts_with_it(pg):
    """Under an ICU or libc collation, '_' and '-' sort as if absent, and a
    range over such a column misses names; the term column is COLLATE "C"."""
    engine = build_engine(pg)
    with engine.connect() as connection:
        names = connection.execute(
            select(feature_attribute.c.feature_id, feature_attribute.c.value)
            .join(feature, feature.c.feature_id == feature_attribute.c.feature_id)
            .where(feature.c.assembly_id == 2, feature.c.segment_index == 0,
                   feature.c.feature_type.not_in(PART_TYPES),
                   feature_attribute.c.key.in_(NAME_KEYS))).all()
        for prefix in ("fgsg_1", "xp_0113", "xm_0113"):
            expected = {row.feature_id for row in names
                        if row.value.lower().startswith(prefix) and row.value.lower() != prefix}
            found = {match.feature_id for match in search(connection, 2, prefix,
                                                          limit=5000).matches}
            assert expected and expected <= found, prefix


def test_the_portal_engine_cannot_write(pg):
    engine = read_only_engine(pg)
    with engine.connect() as connection, pytest.raises(DBAPIError):
        connection.execute(text("DELETE FROM feature"))


def test_pages_on_postgresql(pg, tmp_path):
    portal = PortalData(read_only_engine(pg))
    client = create_app(portal, count_in_background=False).test_client()
    assert client.get("/").status_code == 200
    response = client.get("/g/2/search?q=FGSG_11614")
    assert response.status_code == 200 and "FGSG_11614" in response.get_data(as_text=True)
    row = portal.row_of_id(portal.genome(2), "gene-FGSG_11614")
    assert client.get(f"/g/2/f/{row}").status_code == 200
    portal.close()
