"""The portal (D-111 to D-114): what it shows, and that what it shows is right.

What these assert, stated before they were written:

* The search index holds, for each attribute it indexes, every feature that
  carries the attribute: ``cgload index`` counts both sides and refuses an
  index that falls short. Every name in the fixtures finds its feature.
* A search finds a feature by its exact ID, by a name in any case, by the
  start of a name, by every word of its product, by a cross-reference and by
  a functional accession, best matches first.
* A protein the portal shows is the protein NCBI wrote in the GenBank record:
  the three Dscam1 isoforms, 21 to 23 exons on the minus strand, translate to
  their ``/translation`` exactly.
* A page shows a feature with its whole family; a download is the family's
  rows as GFF3, every Parent resolving inside it.
* The portal changes nothing: its engine cannot write, and the database file
  is the same after every page has been read.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil

import pytest

flask = pytest.importorskip("flask")

from cgload.cli.main import cli  # noqa: E402
from cgload.db.engine import build_engine  # noqa: E402
from cgload.db.schema import feature, feature_attribute  # noqa: E402
from cgload.export import round_trip_counts  # noqa: E402
from cgload.portal.data import PortalData, read_only_engine  # noqa: E402
from cgload.portal.search import (  # noqa: E402
    EXACT_ID,
    EXACT_NAME,
    FUNCTIONAL,
    NAME_KEYS,
    PART_TYPES,
    check_index,
    search,
)
from cgload.portal.web import create_app  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from sqlalchemy import select, text  # noqa: E402
from sqlalchemy.exc import OperationalError  # noqa: E402

from portal_db import DSCAM, FUSARIUM, build  # noqa: E402

FLY, FUSARIUM_ID, HOSTILE = 1, 2, 3


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build(tmp_path_factory.mktemp("portal"))


@pytest.fixture(scope="module")
def portal(built):
    db, _ = built
    data = PortalData(read_only_engine(f"sqlite:///{db}"))
    data.count_features()
    yield data
    data.close()


@pytest.fixture(scope="module")
def client(portal):
    return create_app(portal, count_in_background=False).test_client()


def page(client, url: str, status: int = 200) -> str:
    response = client.get(url)
    assert response.status_code == status, (url, response.status_code)
    return response.get_data(as_text=True)


# --------------------------------------------------------------------------
# The index, and searching it
# --------------------------------------------------------------------------


def test_the_index_is_checked_by_counting(built):
    db, _ = built
    engine = build_engine(f"sqlite:///{db}")
    with engine.connect() as connection:
        for assembly_id in (FLY, FUSARIUM_ID):
            report = check_index(connection, assembly_id)
            assert report.ok, report.failures()
            assert report.features > 0 and report.terms > report.features
            name = next(count for count in report.counts if count.key == "Name")
            assert name.in_attributes == name.in_index


def test_an_index_that_falls_short_is_refused(built, tmp_path):
    """A term deleted from the index makes the counts differ."""
    db, _ = built
    copy = tmp_path / "copy.db"
    shutil.copy(db, copy)
    engine = build_engine(f"sqlite:///{copy}")
    with engine.begin() as connection:
        victim = connection.execute(text(
            "SELECT feature_id FROM search_term WHERE assembly_id = 2 AND field = 'Name' "
            "LIMIT 1")).scalar_one()
        connection.execute(text("DELETE FROM search_term WHERE feature_id = :v AND field = "
                                "'Name'"), {"v": victim})
        report = check_index(connection, FUSARIUM_ID)
    assert not report.ok
    assert [count.key for count in report.failures()] == ["Name"]


def test_index_command_reports_pass_and_ignores_unknown_assemblies(built):
    db, _ = built
    runner = CliRunner()
    result = runner.invoke(cli, ["index", "--db-url", f"sqlite:///{db}", "--assembly-id", "2"])
    assert result.exit_code == 0, result.output
    assert "PASS" in result.output
    missing = runner.invoke(cli, ["index", "--db-url", f"sqlite:///{db}",
                                  "--assembly-id", "99"])
    assert missing.exit_code == 1
    assert "no assembly with id 99" in missing.output


def test_every_name_in_the_files_finds_its_feature(built):
    """End to end: each value of each name key, searched as typed, returns the
    feature that carries it."""
    db, _ = built
    engine = build_engine(f"sqlite:///{db}")
    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.assembly_id, feature.c.feature_id, feature_attribute.c.key,
                   feature_attribute.c.value)
            .join(feature, feature.c.feature_id == feature_attribute.c.feature_id)
            .where(feature.c.assembly_id.in_([FLY, FUSARIUM_ID]),
                   feature.c.segment_index == 0,
                   feature.c.feature_type.not_in(PART_TYPES),
                   feature_attribute.c.key.in_(NAME_KEYS))).all()
        assert len(rows) > 300
        missed = []
        for row in rows:
            found = search(connection, row.assembly_id, row.value.upper(), limit=5000)
            if row.feature_id not in {match.feature_id for match in found.matches}:
                missed.append((row.key, row.value))
    assert missed == []


def test_search_ranks_an_exact_id_first(portal):
    genome = portal.genome(FUSARIUM_ID)
    result, rows, total = portal.search(genome, "gene-FGSG_11614")
    assert rows[0][0].rank == EXACT_ID
    assert rows[0][1].source_id == "gene-FGSG_11614"


def test_search_by_name_is_case_blind_and_puts_genes_first(portal):
    genome = portal.genome(FUSARIUM_ID)
    _, rows, total = portal.search(genome, "fgsg_11614")
    assert rows[0][0].rank == EXACT_NAME
    assert rows[0][1].feature_type == "gene"
    assert {item.feature_type for _, item, _ in rows} >= {"gene", "mRNA", "CDS"}


def test_search_by_the_start_of_a_name(portal):
    genome = portal.genome(FUSARIUM_ID)
    _, rows, total = portal.search(genome, "FGSG_116")
    assert total >= 3
    assert all(any(value.lower().startswith("fgsg_116") for values in item.attributes.values()
                   for value in values) or (item.source_id or "").lower().startswith(
                       "fgsg_116") for _, item, _ in rows)


def test_search_by_words_of_a_product(portal):
    genome = portal.genome(FLY)
    _, rows, _ = portal.search(genome, "cell adhesion molecule")
    assert rows
    assert all("cell adhesion molecule" in " ".join(item.attributes.get("product", [])).lower()
               for _, item, _ in rows)
    # Each match names the attribute its words were found in.
    assert {match.field for match, _, _ in rows} == {"product"}


def test_search_by_a_cross_reference_and_its_number(portal):
    genome = portal.genome(FLY)
    _, by_whole, _ = portal.search(genome, "GeneID:35652")
    _, by_number, _ = portal.search(genome, "35652")
    assert by_whole and {item.row_id for _, item, _ in by_whole} == {
        item.row_id for _, item, _ in by_number}


def test_search_by_a_functional_accession_without_its_prefix(portal):
    genome = portal.genome(FUSARIUM_ID)
    result, rows, _ = portal.search(genome, "K00799")
    assert rows and all(match.rank == FUNCTIONAL for match, _, _ in rows)


def test_without_an_index_names_are_not_searched_but_ids_are(portal):
    genome = portal.genome(HOSTILE)
    assert genome.indexed is None
    result, rows, _ = portal.search(genome, "rna-1")
    assert [item.source_id for _, item, _ in rows] == ["rna-1"]
    result, rows, _ = portal.search(genome, "RNA-1")
    assert rows == [] and not result.indexed


def test_search_pages_do_not_overlap(portal):
    genome = portal.genome(FUSARIUM_ID)
    _, first, total = portal.search(genome, "fgsg", per_page=10, page=1)
    _, second, _ = portal.search(genome, "fgsg", per_page=10, page=2)
    assert total > 10
    assert not {item.row_id for _, item, _ in first} & {item.row_id for _, item, _ in second}


# --------------------------------------------------------------------------
# Sequences
# --------------------------------------------------------------------------


def test_proteins_are_the_proteins_ncbi_wrote(portal):
    """The GenBank record carries each CDS's translation; the portal's own
    translation of the stored rows on the record's sequence must equal it. Two
    CDS are cut by the edge of the region (``<1..96``); their proteins end with
    what the portal translates."""
    genome = portal.genome(FLY)
    with portal.engine.connect() as connection:
        rows = connection.execute(select(feature.c.feature_id).where(
            feature.c.assembly_id == FLY, feature.c.feature_type == "CDS",
            feature.c.segment_index == 0)).scalars().all()
    assert len(rows) == 5
    complete = 0
    for row in rows:
        family = portal.family(genome, row)
        item = family.selected
        declared = "".join(item.attributes["translation"][0].split())
        protein = portal.sequences_of(family, item)["protein"].rstrip("*")
        if item.start == 1:
            assert declared.endswith(protein)
        else:
            assert protein == declared
            assert len(item.segments) >= 21 and item.strand == "-"
            complete += 1
    assert complete == 3


def test_spliced_transcript_is_its_exons_in_order(portal):
    genome = portal.genome(FLY)
    row = portal.row_of_id(genome, "NM_001273835.1")
    family = portal.family(genome, row)
    item = family.selected
    spliced = portal.sequences_of(family, item)["spliced"]
    assert len(spliced) == sum(segment.end - segment.start + 1 for segment in item.segments)
    genomic = portal.sequences_of(family, item)["genomic"]
    assert len(genomic) == item.end - item.start + 1


# --------------------------------------------------------------------------
# Pages and downloads
# --------------------------------------------------------------------------


def test_home_lists_every_genome_with_its_counts(client):
    body = page(client, "/")
    for name in ("Drosophila melanogaster", "Fusarium graminearum", "Hostilia testis"):
        assert name in body
    assert "counting" not in body


def test_genome_page_shows_the_load_receipt(client):
    body = page(client, "/g/2/")
    assert "source count" in body
    assert "pass" in body
    assert FUSARIUM.name in body
    assert "eggnog" in body and "interproscan" in body


def test_feature_page_shows_the_family_drawing_and_annotation(client, portal):
    genome = portal.genome(FUSARIUM_ID)
    row = portal.row_of_id(genome, "rna-XM_011327722.1")
    body = page(client, f"/g/2/f/{row}")
    assert "<svg" in body and 'class="model"' in body
    assert "FGSG_11614" in body
    for member in ("XP_011326024.1", "exon-XM_011327722.1-1"):
        assert member in body
    assert "PF03198" in body and "https://www.ebi.ac.uk/interpro/entry/pfam/PF03198" in body
    assert 'id="igv"' in body


def test_a_later_row_of_a_feature_opens_the_feature(client, portal):
    genome = portal.genome(FUSARIUM_ID)
    with portal.engine.connect() as connection:
        later = connection.execute(select(feature.c.feature_id).where(
            feature.c.assembly_id == FUSARIUM_ID, feature.c.source_id == "cds-XP_011326024.1",
            feature.c.segment_index == 3)).scalar_one()
    family = portal.family(genome, later)
    assert family.selected.source_id == "cds-XP_011326024.1"
    assert len(family.selected.segments) == 6


def test_id_route_redirects_to_the_feature(client):
    response = client.get("/g/2/id/gene-FGSG_11614")
    assert response.status_code == 302
    assert re.search(r"/g/2/f/\d+$", response.headers["Location"])
    assert client.get("/g/2/id/no-such-id").status_code == 404


def test_family_gff3_download_is_whole_and_self_contained(client, portal):
    genome = portal.genome(FUSARIUM_ID)
    row = portal.row_of_id(genome, "gene-FGSG_11614")
    response = client.get(f"/g/2/f/{row}/family.gff3")
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    ids = set(re.findall(r"\tID=([^;\n]+)", body))
    parents = set(re.findall(r";Parent=([^;\n]+)", body))
    assert parents <= ids
    counts = round_trip_counts(body)
    assert counts["gene"] == 1 and counts["CDS"] == 1 and counts["exon"] == 6


def test_fasta_downloads(client, portal):
    genome = portal.genome(FLY)
    row = portal.row_of_id(genome, "NP_001260764.1")
    protein = page(client, f"/g/1/f/{row}/protein.fa")
    assert protein.startswith(">NP_001260764.1 protein")
    assert len("".join(protein.splitlines()[1:]).rstrip("*")) == 1947
    page(client, f"/g/1/f/{row}/spliced.fa", status=404)  # a CDS is not spliced twice
    page(client, f"/g/1/f/{row}/everything.fa", status=404)


def test_region_fasta_is_the_assembly_slice(client, built):
    _, data = built
    fasta = next(path for path in data.rglob("assembly.fa") if "drosophila" in str(path))
    bases = "".join(fasta.read_text().splitlines()[1:])
    body = page(client, "/g/1/region.fa?r=NT_033778.4:101-160")
    assert "".join(body.splitlines()[1:]) == bases[100:160]


def test_region_page_lists_the_genes_there(client):
    body = page(client, "/g/1/region?r=NT_033778.4:5000-6000")
    assert "Dscam1" in body
    assert "CG12164" not in body.split("top-level feature")[1]
    page(client, "/g/1/region?r=NT_033778.4:9-1", status=400)
    page(client, "/g/1/region?r=nowhere:1-10", status=400)


def test_files_are_served_with_ranges_and_only_those_four(client):
    response = client.get("/g/1/files/assembly.fa", headers={"Range": "bytes=0-9"})
    assert response.status_code == 206
    assert response.data.startswith(b">NT_033778")
    assert client.get("/g/1/files/annotation.gff3.gz.tbi").status_code == 200
    for name in ("manifest.json", "organism.toml", "portal.db", "..%2F..%2Fportal.db"):
        assert client.get(f"/g/1/files/{name}").status_code == 404


def test_json_api(client):
    genomes = json.loads(page(client, "/api/v1/genomes"))["genomes"]
    assert [genome["organism"] for genome in genomes][:2] == [
        "Drosophila melanogaster", "Fusarium graminearum"]
    assert genomes[0]["files"]["assembly.fa"].endswith("/files/assembly.fa")
    found = json.loads(page(client, "/api/v1/genomes/2/search?q=FGSG_11614"))
    assert found["total"] >= 3 and found["results"][0]["feature"]["type"] == "gene"
    row = found["results"][0]["feature"]["row"]
    detail = json.loads(page(client, f"/api/v1/genomes/2/features/{row}"))
    assert len(detail["family"]) == 9
    assert any(hits for hits in detail["annotations"].values())
    region = json.loads(page(client, "/api/v1/genomes/1/region?r=NT_033778.4:1-80000"))
    assert {item["name"] for item in region["features"]} >= {"Dscam1", "CG12164"}
    located = json.loads(page(client, "/g/1/locate?q=DSCAM1"))
    assert located[0]["chromosome"] == "NT_033778.4"
    assert json.loads(page(client, "/api/v1/genomes/9/search?q=x", status=404))["status"] == 404


# --------------------------------------------------------------------------
# Read-only
# --------------------------------------------------------------------------


def test_the_portal_cannot_write(portal):
    with portal.engine.connect() as connection, pytest.raises(OperationalError):
        connection.execute(text("DELETE FROM feature"))


def test_reading_every_page_changes_nothing(built, client, portal):
    db, _ = built
    before = hashlib.sha256(db.read_bytes()).hexdigest()
    urls = ["/", "/about", "/search?q=dscam", "/g/1/", "/g/2/", "/g/3/",
            "/g/2/search?q=fgsg", "/g/1/region?r=NT_033778.4", "/api/v1/genomes"]
    with portal.engine.connect() as connection:
        rows = connection.execute(select(feature.c.feature_id, feature.c.assembly_id).where(
            feature.c.segment_index == 0)).all()
    urls += [f"/g/{row.assembly_id}/f/{row.feature_id}" for row in rows]
    for url in urls:
        assert client.get(url).status_code == 200, url
    assert hashlib.sha256(db.read_bytes()).hexdigest() == before


def test_removing_an_assembly_removes_its_index(built, tmp_path):
    db, data = built
    copy_db = tmp_path / "copy.db"
    shutil.copy(db, copy_db)
    result = CliRunner().invoke(cli, ["remove", "--db-url", f"sqlite:///{copy_db}",
                                      "--assembly-id", "2", "--keep-files", "--yes"])
    assert result.exit_code == 0, result.output
    engine = build_engine(f"sqlite:///{copy_db}")
    with engine.connect() as connection:
        left = connection.execute(text(
            "SELECT COUNT(*) FROM search_term WHERE assembly_id = 2")).scalar_one()
        meta = connection.execute(text(
            "SELECT COUNT(*) FROM search_index WHERE assembly_id = 2")).scalar_one()
    assert (left, meta) == (0, 0)


def test_fixture_paths_exist():
    assert DSCAM.exists() and FUSARIUM.exists()


def test_waitress_serves_many_requests_at_once(portal):
    """The production server, not the test client: 8 clients, 200 requests,
    pages, JSON and byte ranges, every one answered."""
    waitress = pytest.importorskip("waitress")
    import threading
    import urllib.request
    from concurrent.futures import ThreadPoolExecutor

    app = create_app(portal, count_in_background=False)
    server = waitress.create_server(app, host="127.0.0.1", port=0, threads=4)
    port = server.effective_port
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    genome = portal.genome(FUSARIUM_ID)
    row = portal.row_of_id(genome, "gene-FGSG_11614")
    paths = ["/", "/g/2/", f"/g/2/f/{row}", "/g/2/search?q=fgsg", "/api/v1/genomes",
             "/g/1/region?r=NT_033778.4:1-9000", "/g/1/locate?q=dscam1"]

    def fetch(index: int) -> int:
        path = paths[index % len(paths)]
        request = urllib.request.Request(f"http://127.0.0.1:{port}{path}")
        if index % 10 == 0:
            request = urllib.request.Request(
                f"http://127.0.0.1:{port}/g/1/files/assembly.fa",
                headers={"Range": "bytes=100-199"})
        with urllib.request.urlopen(request, timeout=30) as response:  # noqa: S310
            body = response.read()
            if index % 10 == 0:
                assert response.status == 206 and len(body) == 100
            return response.status

    try:
        with ThreadPoolExecutor(max_workers=8) as pool:
            statuses = list(pool.map(fetch, range(200)))
    finally:
        server.close()
    assert set(statuses) <= {200, 206} and len(statuses) == 200
