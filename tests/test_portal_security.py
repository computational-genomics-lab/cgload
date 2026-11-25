"""The portal against hostile input (D-111).

A genome portal has two kinds of input an attacker controls: what a visitor
sends, and what the annotation file said, which reaches every page that shows
the feature. These tests send both and assert:

* nothing a visitor sends makes the server fail (no 5xx) or reach SQL as
  text: classic injection strings come back as ordinary searches with no
  matches;
* no file path outside the four an assembly serves can be read;
* text from the file is escaped wherever it is shown, in pages, in the gene
  model drawing and in links, and a cross-reference becomes a link only when
  its accession has the expected form;
* every response carries the security headers, and the policy allows no
  inline script and no script from another site;
* the portal answers only GET and HEAD.

Plus a seeded fuzz: thousands of random requests over every route, none of
which may produce a server error.
"""

from __future__ import annotations

import random
import string

import pytest

flask = pytest.importorskip("flask")

from cgload.db.schema import feature, source_file  # noqa: E402
from cgload.portal.data import PortalData, file_name, read_only_engine  # noqa: E402
from cgload.portal.web import CONTENT_SECURITY_POLICY, create_app, xref_link  # noqa: E402
from sqlalchemy import select  # noqa: E402

from portal_db import HOSTILE_ID, HOSTILE_NAME, HOSTILE_PRODUCT, build  # noqa: E402

INJECTIONS = [
    "' OR '1'='1",
    "' OR 1=1 --",
    "\" OR \"\"=\"",
    "1; DROP TABLE feature; --",
    "x' UNION SELECT source_id, 1 FROM feature --",
    "%' AND 1=1 AND '%'='",
    "\\'; SELECT pg_sleep(5); --",
    "<script>alert(1)</script>",
    "\"><svg onload=alert(1)>",
    "javascript:alert(1)",
    "{{7*7}}",
    "${7*7}",
    "%00",
    "\x00",
    "../../../../etc/passwd",
    "chr1:1-99999999999999999999999",
    "NT_033778.4:-5--1",
    "NT_033778.4:1-2:3-4",
    "A" * 5000,
]


@pytest.fixture(scope="module")
def app(tmp_path_factory):
    db, _ = build(tmp_path_factory.mktemp("portal-security"))
    portal = PortalData(read_only_engine(f"sqlite:///{db}"))
    portal.count_features()
    application = create_app(portal, count_in_background=False)
    application.config["PROPAGATE_EXCEPTIONS"] = False
    yield application
    portal.close()


@pytest.fixture(scope="module")
def client(app):
    return app.test_client()


@pytest.fixture(scope="module")
def hostile_rows(app):
    portal = app.extensions["cgload_portal"]
    with portal.engine.connect() as connection:
        return dict(connection.execute(
            select(feature.c.source_id, feature.c.feature_id)
            .where(feature.c.assembly_id == 3, feature.c.segment_index == 0)).all())


def routes(value: str) -> list[str]:
    from urllib.parse import quote

    q = quote(value, safe="")
    return [
        f"/search?q={q}", f"/g/1/search?q={q}", f"/g/2/search?q={q}",
        f"/g/1/region?r={q}", f"/g/1/region.fa?r={q}", f"/api/v1/genomes/2/search?q={q}",
        f"/api/v1/genomes/1/region?r={q}", f"/g/1/locate?q={q}", f"/g/2/search?q=a&page={q}",
        f"/g/1/id/{q}", f"/g/1/files/{q}",
    ]


@pytest.mark.parametrize("value", INJECTIONS)
def test_hostile_requests_never_fail_the_server(client, value):
    for url in routes(value):
        response = client.get(url)
        assert response.status_code < 500, (url, response.status_code)
        body = response.get_data(as_text=True)
        assert "<script>alert(1)</script>" not in body
        assert "<svg onload" not in body
        assert "Traceback" not in body


def test_sql_injection_is_only_ever_a_search_for_that_text(app):
    portal = app.extensions["cgload_portal"]
    genome = portal.genome(2)
    for value in INJECTIONS[:7]:
        result, rows, total = portal.search(genome, value)
        assert total == 0, value


def test_no_file_outside_the_four_can_be_read(client):
    for name in ("../manifest.json", "..%2Fmanifest.json", "%2e%2e%2fportal.db",
                 "/etc/passwd", "..\\..\\portal.db", "assembly.fa%00.fai", "manifest.json",
                 ".", "..", "assembly.fa/"):
        response = client.get(f"/g/1/files/{name}")
        assert response.status_code in (404, 308), name
        if response.status_code == 308:  # Werkzeug normalising a trailing slash
            assert client.get(response.headers["Location"]).status_code == 404


def test_text_from_the_file_is_escaped_everywhere(client, hostile_rows):
    gene = hostile_rows[HOSTILE_ID]
    rna = hostile_rows["rna-1"]
    for url in (f"/g/3/f/{gene}", f"/g/3/f/{rna}", "/g/3/search?q=rna-1",
                "/g/3/region?r=chrX:1-5000"):
        body = client.get(url, follow_redirects=True).get_data(as_text=True)
        assert HOSTILE_NAME not in body
        assert HOSTILE_PRODUCT not in body
        assert "<img src=x" not in body
        assert "&lt;script&gt;" in body or "&#34;&gt;&lt;img" in body or url.endswith("5000")


def test_the_drawing_escapes_its_labels(client, hostile_rows):
    body = client.get(f"/g/3/f/{hostile_rows[HOSTILE_ID]}").get_data(as_text=True)
    drawing = body[body.index("<svg"): body.index("</svg>")]
    assert "<script" not in drawing and "&lt;script" in drawing


def test_a_cross_reference_links_only_in_its_proper_form():
    assert xref_link("GeneID:672") == "https://www.ncbi.nlm.nih.gov/gene/672"
    assert xref_link("GeneID:javascript:alert(1)") is None
    assert xref_link("HGNC:HGNC:1100").endswith("hgnc_id/HGNC:1100")
    assert xref_link("Unknown:1") is None


def test_pages_name_input_files_without_the_loading_machines_paths(app, client):
    """The database keeps where each input file was on the machine that loaded
    it; a public page shows its name only, not that machine's directories."""
    portal = app.extensions["cgload_portal"]
    with portal.engine.connect() as connection:
        paths = connection.execute(select(source_file.c.assembly_id, source_file.c.path)).all()
    assert paths
    for assembly_id, path in paths:
        body = client.get(f"/g/{assembly_id}/").get_data(as_text=True)
        assert path.replace("\\", "/").rsplit("/", 1)[0] not in body
        assert file_name(path) in body
    assert file_name("C:\\genomes\\fusarium\\genomic.gff") == "genomic.gff"


def test_hostile_ids_route_safely(client, hostile_rows):
    from urllib.parse import quote

    response = client.get("/g/3/id/" + quote(HOSTILE_ID, safe=""))
    assert response.status_code == 302
    assert response.headers["Location"].endswith(f"/g/3/f/{hostile_rows[HOSTILE_ID]}")


def test_every_response_carries_the_security_headers(client):
    for url in ("/", "/g/1/", "/g/1/f/1", "/api/v1/genomes", "/g/1/files/assembly.fa.fai",
                "/nowhere", "/g/9/", "/static/style.css"):
        response = client.get(url)
        assert response.headers["Content-Security-Policy"] == CONTENT_SECURITY_POLICY
        assert response.headers["X-Content-Type-Options"] == "nosniff"
        assert response.headers["X-Frame-Options"] == "DENY"
    policy = dict(part.split(" ", 1) for part in CONTENT_SECURITY_POLICY.split("; "))
    assert policy["script-src"] == "'self'"
    assert policy["object-src"] == "'none'" and policy["frame-ancestors"] == "'none'"


def test_pages_carry_no_inline_script(client):
    for url in ("/", "/g/1/", "/g/1/f/7", "/g/1/region?r=NT_033778.4:1-9000", "/about"):
        body = client.get(url).get_data(as_text=True)
        for tag in body.split("<script")[1:]:
            assert tag.lstrip().startswith("src="), url
        assert " on" + "load=" not in body and "javascript:" not in body


def test_only_get_and_head(client):
    for method in ("post", "put", "delete", "patch"):
        for url in ("/", "/g/1/", "/api/v1/genomes", "/g/1/files/assembly.fa"):
            assert getattr(client, method)(url).status_code == 405


def test_limits_on_what_a_request_may_say(client):
    assert client.get("/g/1/search?q=" + "a" * 201).status_code == 400
    assert client.get("/g/1/region?r=" + "a" * 301).status_code == 400
    assert client.get("/g/1/search?q=a&page=0").status_code == 400
    assert client.get("/g/1/search?q=a&page=1001").status_code == 400
    assert client.get("/g/1/region.fa?r=NT_033778.4").status_code == 200


def test_seeded_fuzz_finds_no_server_error(client):
    """3,000 requests of random text, bytes and numbers over every route."""
    rng = random.Random(20261125)
    alphabet = string.printable + "äé漢字‮​<>\"'%;&#\\"
    failures = []
    for _ in range(200):
        length = rng.choice([0, 1, 2, 5, 20, 120, 199, 200, 201, 299])
        value = "".join(rng.choice(alphabet) for _ in range(length))
        number = str(rng.choice([0, 1, 2, 7, 99, 10**9, 2**63, -1]))
        urls = routes(value) + [f"/g/{number}/", f"/g/1/f/{number}",
                                f"/api/v1/genomes/{number}/features/{number}",
                                f"/g/1/f/{number}/protein.fa"]
        assert len(urls) == 15  # 200 rounds of 15: 3,000 requests
        for url in urls:
            status = client.get(url).status_code
            if status >= 500:
                failures.append((url, status))
    assert failures == []
