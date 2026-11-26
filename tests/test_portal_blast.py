"""BLAST in the portal (D-116).

What these assert, stated before they were written:

* A query is one sequence of its alphabet's letters, 10 to 10,000 long; the
  program, task, E-value, number of hits and genomes each come from a fixed
  list. Anything else is refused in a sentence, with status 400.
* No text a visitor sends reaches BLAST+ as an argument: the argument list is
  made of fixed choices and of paths the portal found itself, and the query
  goes to BLAST's standard input under the header ``>query``.
* Searches beyond those running and waiting are turned away at once, with 503
  and Retry-After.
* BLAST's output becomes hits whose links open the CDS or the region, and
  text in a hit's title is escaped like all text from the files.
* POST is answered on the two BLAST URLs only, 405 elsewhere as before; a body
  over 64 KiB is refused with 413.
* A build that fails leaves no database, not even an older one, and a load
  newer than its databases hides them.
* With BLAST+ installed: ``cgload blastdb`` builds databases whose counts match
  the database (every sequence; a protein for every CDS); each Dscam1
  isoform's protein finds itself first, at 100% identity over its length;
  tblastn finds it on the genome inside the CDS's span; blastn finds a piece of
  the genome where it lies, on either strand; blastx finds the CDS's protein.
"""

from __future__ import annotations

import json
import shutil
import threading
import time
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace

import pytest

flask = pytest.importorskip("flask")

from cgload.cli.main import cli  # noqa: E402
from cgload.db.schema import feature  # noqa: E402
from cgload.portal import blast as blasting  # noqa: E402
from cgload.portal.blast import (  # noqa: E402
    Blast,
    BlastBusy,
    Databases,
    Hit,
    Hsp,
    QueryError,
    clean_query,
    databases_of,
    make_request,
    ncbi_url,
    parse,
)
from cgload.portal.data import PortalData, read_only_engine  # noqa: E402
from cgload.portal.web import CONTENT_SECURITY_POLICY, create_app  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from sqlalchemy import select  # noqa: E402

from portal_db import build  # noqa: E402

PROTEIN = "MKVLAAGIVGLLLAQSAWAEDPKTSRELIEKAVRDF"
DNA = "ATGAAAGTTCTGGCAGCAGGTATTGTTGGTCTGCTGCTG"
HOSTILE = ['<script>alert(1)</script>', '"><img src=x onerror=alert(1)>', "' OR '1'='1",
           "1; DROP TABLE feature; --", "{{7*7}}", "../../etc/passwd", "\x00", "A" * 40_000]
HAVE_BLAST = all(shutil.which(name) for name in ("blastn", "tblastn", "blastp", "blastx",
                                                  "makeblastdb"))
needs_blast = pytest.mark.skipif(not HAVE_BLAST, reason="BLAST+ is not installed")


# --------------------------------------------------------------------------
# What a visitor may ask
# --------------------------------------------------------------------------


def test_a_query_may_be_fasta_or_bare_sequence_with_genbank_numbering():
    sequence, title = clean_query(">my query here\nACGT ACGT\n  61 acgtu\r\n", "nucleotide")
    assert sequence == "ACGTACGTACGTT" and title == "my"
    sequence, title = clean_query("mkvl aagi vgll laqs\n*", "protein")
    assert sequence == "MKVLAAGIVGLLLAQS" and title == "query"
    _, title = clean_query('>"><script>x\n' + PROTEIN, "protein")
    assert set(title) <= set("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._:|-")


@pytest.mark.parametrize(("text", "alphabet", "words"), [
    (">a\nACGTACGTACGT\n>b\nACGTACGTACGT", "nucleotide", "One sequence at a time"),
    (">only a header\n", "nucleotide", "no sequence"),
    ("", "protein", "no sequence"),
    ("A" * 10_001, "nucleotide", "at most 10,000"),
    ("ACGTACGT", "nucleotide", "at least 10"),
    (PROTEIN, "nucleotide", "It looks like a protein"),
    ("MKVL<script>AAGIVG", "protein", "this one also has"),
    ("ACGT;rm -rf /;ACGTACGT", "nucleotide", "this one also has"),
])
def test_what_is_not_one_sequence_of_its_alphabet_is_refused_in_a_sentence(text, alphabet, words):
    with pytest.raises(QueryError) as caught:
        clean_query(text, alphabet)
    assert words in str(caught.value)
    assert str(caught.value)[0].isupper() and str(caught.value).endswith(".")


def test_every_choice_comes_from_a_fixed_list():
    def available(program):
        return [1, 2]

    base = {"query": PROTEIN, "program": "blastp"}
    request = make_request(base, available)
    assert request.assembly_ids == (1, 2) and request.task is None
    assert request.evalue == blasting.DEFAULT_EVALUE and request.max_hits == 25
    # A form sends a task whatever the program; only blastn has one.
    assert make_request({**base, "task": "anything"}, available).task is None
    request = make_request({"query": DNA, "program": "blastn", "genome": "2"}, available)
    assert request.task == "megablast" and request.assembly_ids == (2,)
    # What a form sends with flat=False, and what JSON sends.
    request = make_request({"query": [PROTEIN], "program": ["blastp"], "genome": ["1", "2"],
                            "max_hits": ["100"], "evalue": ["1e-10"]}, available)
    assert request.assembly_ids == (1, 2) and request.max_hits == 100
    assert make_request({**base, "genome": [2], "max_hits": 10}, available).assembly_ids == (2,)
    for bad in ({"program": "psiblast"}, {"program": "blastp -remote"},
                {"program": "blastn", "query": DNA, "task": "dc-megablast"},
                {"evalue": "1e-3"}, {"evalue": "0.001; rm -rf /"}, {"evalue": "-1"},
                {"max_hits": "500"}, {"max_hits": "-1"}, {"max_hits": "25 -out /tmp/x"},
                {"genome": "3"}, {"genome": "1 OR 1=1"}, {"genome": ["1", "x"]},
                {"genome": None}, {"query": "M" * 30_001}):
        with pytest.raises(QueryError):
            make_request({**base, **bad}, available)
    with pytest.raises(QueryError, match="No genome here"):
        make_request(base, lambda program: [])


# --------------------------------------------------------------------------
# What BLAST is given, and what comes back
# --------------------------------------------------------------------------


class _NoGenomes:
    def genomes(self):
        return []


def _blast(tmp_path: Path, jobs: int = 1, wait: float = 0.05) -> Blast:
    names = ("blastn", "tblastn", "blastp", "blastx")
    blast = Blast(_NoGenomes(), binaries={name: str(tmp_path / name) for name in names},
                  jobs=jobs, threads=2, timeout=5, wait=wait)
    blast.databases = {1: Databases(1, tmp_path / "a" / "blast", True, True, {}),
                       2: Databases(2, tmp_path / "b" / "blast", True, True, {})}
    return blast


def test_the_query_is_never_an_argument(tmp_path, monkeypatch):
    blast = _blast(tmp_path)
    request = make_request({"query": ">x';$(reboot)\n" + PROTEIN, "program": "blastp",
                            "evalue": "1e-5", "max_hits": "10"}, blast.searchable)
    seen = {}

    def run(args, **kwargs):
        seen.update(args=args, **kwargs)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(blasting.subprocess, "run", run)
    assert blast.search(request).hits == []
    args = seen["args"]
    assert args[0] == str(tmp_path / "blastp")
    assert seen["input"] == f">query\n{PROTEIN}\n"
    assert not seen.get("shell") and seen["cwd"] == tmp_path
    assert args[args.index("-db") + 1] == "a/blast/proteins b/blast/proteins"
    assert "BLASTDB" not in seen["env"]
    fixed = {"-db", "-outfmt", blasting.OUTFMT, "-evalue", "1e-5", "-max_target_seqs", "10",
             "-max_hsps", str(blasting.MAX_HSPS), "-num_threads", "2"}
    assert set(args[1:]) - fixed == {args[args.index("-db") + 1]}
    assert not any(PROTEIN in arg or "reboot" in arg for arg in args)


def test_databases_are_named_from_their_common_folder(tmp_path):
    paths = [tmp_path / "a" / "x" / "blast" / "genome", tmp_path / "b" / "blast" / "genome"]
    assert blasting._groups(paths) == [(tmp_path, ["a/x/blast/genome", "b/blast/genome"])]
    # One database keeps a folder in its name, so BLASTDB is never searched for it.
    assert blasting._groups(paths[:1]) == [(tmp_path / "a" / "x", ["blast/genome"])]
    spaced = [tmp_path / "a b" / "blast" / "genome", tmp_path / "c" / "blast" / "genome"]
    assert blasting._groups(spaced) == [(path.parent.parent, ["blast/genome"])
                                        for path in spaced]


def test_searches_beyond_those_running_and_waiting_are_turned_away_at_once(tmp_path,
                                                                            monkeypatch):
    blast = _blast(tmp_path, jobs=1, wait=5)
    entered, gate = threading.Event(), threading.Event()

    def run(request, group):
        entered.set()
        gate.wait(10)
        return ""

    monkeypatch.setattr(blast, "_run", run)
    request = make_request({"query": PROTEIN}, blast.searchable)
    outcomes = []

    def search():
        try:
            blast.search(request)
            outcomes.append("done")
        except BlastBusy:
            outcomes.append("busy")

    first = threading.Thread(target=search)
    first.start()
    assert entered.wait(5)
    second = threading.Thread(target=search)  # waits for the one slot
    second.start()
    time.sleep(0.2)
    started = time.monotonic()
    with pytest.raises(BlastBusy):  # a third may not even wait
        blast.search(request)
    assert time.monotonic() - started < 1
    gate.set()
    first.join(10)
    second.join(10)
    assert outcomes == ["done", "done"]

    blast = _blast(tmp_path, jobs=1, wait=0.1)
    gate.clear()
    entered.clear()
    monkeypatch.setattr(blast, "_run", run)
    first = threading.Thread(target=search)
    first.start()
    assert entered.wait(5)
    with pytest.raises(BlastBusy):  # a second waits its time, then is turned away
        blast.search(request)
    gate.set()
    first.join(10)


def _line(title: str, qstart=1, qend=36, sstart=1, send=36, evalue="1e-20", bits="70.1",
          strand="N/A", slen=400) -> str:
    return "\t".join(map(str, ["query", title, "100.000", 36, 0, 0, qstart, qend, sstart, send,
                               evalue, bits, 36, slen, strand, PROTEIN, PROTEIN]))


def test_output_becomes_hits_with_their_hsps_in_order():
    text = "\n".join([
        _line("p1.77 NP_001260764.1 down syndrome cell adhesion molecule 1"),
        _line("p1.77 NP_001260764.1 down syndrome cell adhesion molecule 1", evalue="1e-5"),
        _line("p2.5 FGSG_1 hypothetical protein", evalue="1e-3"),
        _line("p2.6"),
        "not\ta\tline",
        _line("gNOT A TITLE"),
        _line("p1.x broken row"),
        "\t".join(["query", "p1.9 bad numbers"] + ["x"] * 15),
    ])
    hits = parse(text, "proteins")
    assert [(hit.assembly_id, hit.subject, len(hit.hsps)) for hit in hits] == [
        (1, "77", 2), (2, "5", 1), (2, "6", 1)]
    assert hits[0].name == "NP_001260764.1"
    assert hits[0].product == "down syndrome cell adhesion molecule 1"
    assert hits[2].name == "row 6" and hits[2].product == ""
    assert hits[0].best.strand is None and hits[0].best.evalue == 1e-20
    genome = parse(_line("g3 NT_033778.4", strand="minus", sstart=900, send=800), "genome")
    assert genome[0].name == "NT_033778.4" and genome[0].best.strand == "minus"


def test_alignment_rows_count_positions_on_each_side():
    gapped = Hsp(pident=90.0, length=6, mismatches=0, gap_opens=1, qstart=1, qend=5,
                 sstart=10, send=15, evalue=1e-5, bitscore=20.0, strand=None,
                 qseq="MKV-LL", sseq="MKVALL")
    assert gapped.match_line() == "|||-||".replace("-", " ")
    assert gapped.blocks("blastp", width=3) == [(1, "MKV", "|||", "MKV", 10, 3, 12),
                                                (4, "-LL", " ||", "ALL", 13, 5, 15)]
    tblastn = Hsp(pident=100.0, length=10, mismatches=0, gap_opens=0, qstart=1, qend=10,
                  sstart=300, send=271, evalue=1e-5, bitscore=20.0, strand="minus",
                  qseq="MKVLAAGIVG", sseq="MKVLAAGIVG")
    assert tblastn.blocks("tblastn", width=5) == [(1, "MKVLA", "|||||", "MKVLA", 300, 5, 286),
                                                  (6, "AGIVG", "|||||", "AGIVG", 285, 10, 271)]
    blastx = Hsp(pident=100.0, length=4, mismatches=0, gap_opens=0, qstart=12, qend=1,
                 sstart=5, send=8, evalue=1e-5, bitscore=20.0, strand=None,
                 qseq="MKVL", sseq="MKVL")
    assert blastx.blocks("blastx") == [(12, "MKVL", "||||", "MKVL", 5, 1, 8)]


def test_a_genome_hit_spans_its_hsps_on_one_strand_nearby():
    def hsp(start, end, strand="plus"):
        return Hsp(100.0, 10, 0, 0, 1, 10, start, end, 1e-9, 30.0, strand, "A", "A")

    hit = Hit("genome", 1, "chr1", "", 10**7,
              [hsp(5000, 5100), hsp(7000, 7200), hsp(9000, 9100, "minus"), hsp(900_000, 900_100)])
    assert blasting.span_of(hit, within=100_000) == (5000, 7200)


def test_the_ncbi_link_carries_an_accession_or_a_short_sequence():
    link = ncbi_url("blastp", PROTEIN, "NP_001260764.1")
    assert link.startswith("https://blast.ncbi.nlm.nih.gov/Blast.cgi?")
    assert "QUERY=NP_001260764.1" in link and "PROGRAM=blastp" in link
    link = ncbi_url("blastp", PROTEIN, '"><script>')
    assert f"QUERY={PROTEIN}" in link and "<" not in link
    assert ncbi_url("blastn", "A" * 4001) is None
    assert ncbi_url("psiblast", PROTEIN) is None


# --------------------------------------------------------------------------
# Which databases are used
# --------------------------------------------------------------------------


def _genome(directory: Path, assembly_id: int = 1, loaded: str = "2026-11-26T10:00:00"):
    return SimpleNamespace(assembly_id=assembly_id, directory=directory, label="G",
                           loaded_at=datetime.fromisoformat(loaded))


def _manifest(directory: Path, **changes) -> None:
    blast_dir = directory / "blast"
    blast_dir.mkdir(parents=True, exist_ok=True)
    for name in ("genome.nsq", "proteins.psq"):
        (blast_dir / name).write_text("x")
    manifest = {"assembly_id": 1, "loaded_at": "2026-11-26T10:00:00",
                "genome": {"check": "PASS"}, "proteins": {"check": "PASS"}}
    manifest.update(changes)
    (blast_dir / "manifest.json").write_text(json.dumps(manifest))


def test_only_databases_that_passed_and_match_the_load_are_used(tmp_path):
    _manifest(tmp_path)
    found = databases_of([_genome(tmp_path)])
    assert found[1].genome and found[1].proteins
    assert databases_of([_genome(tmp_path, loaded="2026-11-27T09:00:00")]) == {}
    assert databases_of([_genome(tmp_path, assembly_id=2)]) == {}
    _manifest(tmp_path, proteins={"check": "FAIL"})
    found = databases_of([_genome(tmp_path)])
    assert found[1].genome and not found[1].proteins
    (tmp_path / "blast" / "genome.nsq").unlink()
    _manifest(tmp_path, proteins={"check": "FAIL"})
    (tmp_path / "blast" / "genome.nsq").unlink()
    assert databases_of([_genome(tmp_path)]) == {}
    (tmp_path / "blast" / "manifest.json").write_text("{not json")
    assert databases_of([_genome(tmp_path)]) == {}
    assert databases_of([SimpleNamespace(directory=None)]) == {}


# --------------------------------------------------------------------------
# Through the portal, without BLAST+: BLAST's answer is made up here
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build(tmp_path_factory.mktemp("portal-blast"))


@pytest.fixture(scope="module")
def portal(built):
    db, _ = built
    data = PortalData(read_only_engine(f"sqlite:///{db}"))
    yield data
    data.close()


@pytest.fixture(scope="module")
def cds_rows(portal):
    with portal.engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.assembly_id, feature.c.source_id, feature.c.feature_id)
            .where(feature.c.feature_type == "CDS", feature.c.segment_index == 0)).all()
    return {(assembly_id, source_id): row for assembly_id, source_id, row in rows}


@pytest.fixture
def made_up(portal, cds_rows, tmp_path, monkeypatch):
    """An app whose BLAST answers with made-up output, hostile titles included."""
    names = ("blastn", "tblastn", "blastp", "blastx")
    blast = Blast(portal, binaries={name: str(tmp_path / name) for name in names}, jobs=1,
                  wait=0.05)
    blast.databases = {genome.assembly_id: Databases(genome.assembly_id, tmp_path, True, True,
                                                     {}) for genome in portal.genomes()}
    dscam = cds_rows[(1, "NP_001260764.1")] if (1, "NP_001260764.1") in cds_rows \
        else next(row for (aid, _), row in cds_rows.items() if aid == 1)
    hostile = next(row for (aid, _), row in cds_rows.items() if aid == 3)
    answers = {
        "proteins": "\n".join([
            _line(f"p1.{dscam} NP_001260764.1 down syndrome cell adhesion molecule 1"),
            _line(f"p3.{hostile} {HOSTILE[0]} {HOSTILE[1]}", evalue="1e-3")]),
        "genome": _line("g1 NT_033778.4", sstart=19280, send=18054, strand="minus",
                        slen=73976),
    }
    monkeypatch.setattr(blast, "_run", lambda request, group: answers[
        blasting.PROGRAMS[request.program][1]])
    app = create_app(portal, count_in_background=False, blast=blast)
    app.config["PROPAGATE_EXCEPTIONS"] = False
    return SimpleNamespace(app=app, client=app.test_client(), blast=blast, dscam=dscam,
                           hostile=hostile)


def test_the_blast_page_offers_the_programs_and_the_genomes(made_up):
    body = made_up.client.get("/blast").get_data(as_text=True)
    assert 'action="/blast"' in body and 'method="post"' in body
    for program in ("blastn", "tblastn", "blastp", "blastx"):
        assert f'value="{program}"' in body
    assert "every genome here (3)" in body
    assert '<a href="/blast">BLAST</a>' in made_up.client.get("/").get_data(as_text=True)


def test_hits_link_to_their_feature_or_region_and_their_titles_are_escaped(made_up):
    response = made_up.client.post("/blast", data={"query": PROTEIN, "program": "blastp"})
    assert response.status_code == 200
    body = response.get_data(as_text=True)
    assert f'href="/g/1/f/{made_up.dscam}"' in body and f'href="/g/3/f/{made_up.hostile}"' in body
    assert "down syndrome cell adhesion molecule 1" in body
    assert HOSTILE[0] not in body and "<img src=x" not in body and "&lt;script&gt;" in body
    assert response.headers["Content-Security-Policy"] == CONTENT_SECURITY_POLICY
    for tag in body.split("<script")[1:]:
        assert tag.lstrip().startswith("src=")
    body = made_up.client.post("/blast", data={"query": PROTEIN, "program": "tblastn",
                                               "genome": "1"}).get_data(as_text=True)
    assert 'href="/g/1/region?r=NT_033778.4:17554-19780"' in body
    assert "minus strand" in body


def test_the_api_answers_in_json(made_up):
    response = made_up.client.post("/api/v1/blast",
                                   json={"query": PROTEIN, "program": "blastp", "genome": [1, 3]})
    assert response.status_code == 200
    answer = response.get_json()
    assert answer["program"] == "blastp" and answer["genomes"] == [1, 3]
    first = answer["hits"][0]
    assert first["url"] == f"/g/1/f/{made_up.dscam}" and first["name"] == "NP_001260764.1"
    assert first["hsps"][0]["identity"] == 100.0
    refused = made_up.client.post("/api/v1/blast", json={"query": PROTEIN, "program": "psiblast"})
    assert refused.status_code == 400 and "program" in refused.get_json()["error"]
    assert made_up.client.post("/api/v1/blast", json=[1, 2]).status_code == 400
    assert made_up.client.post("/api/v1/blast", data="{not json",
                               content_type="application/json").status_code == 400


def test_a_full_portal_turns_searches_away_with_503(made_up):
    blast = made_up.blast
    for _ in range(2 * blast.jobs):
        blast._admitted.acquire()
    try:
        page = made_up.client.post("/blast", data={"query": PROTEIN})
        assert page.status_code == 503 and page.headers["Retry-After"] == "30"
        assert "as many as it allows" in page.get_data(as_text=True)
        api = made_up.client.post("/api/v1/blast", json={"query": PROTEIN})
        assert api.status_code == 503 and api.headers["Retry-After"] == "30"
        assert made_up.client.get("/g/1/").status_code == 200  # the pages still answer
    finally:
        for _ in range(2 * blast.jobs):
            blast._admitted.release()
    assert made_up.client.post("/blast", data={"query": PROTEIN}).status_code == 200


def test_post_is_answered_on_the_blast_urls_only_and_bodies_are_bounded(made_up):
    client = made_up.client
    for url in ("/", "/g/1/", "/api/v1/genomes", "/g/1/files/assembly.fa", "/search"):
        assert client.post(url).status_code == 405, url
    for method in ("put", "delete", "patch"):
        assert getattr(client, method)("/blast").status_code == 405
        assert getattr(client, method)("/api/v1/blast").status_code == 405
    assert client.get("/api/v1/blast").status_code == 405
    assert client.post("/blast", data={"query": "M" * 70_000}).status_code == 413


@pytest.mark.parametrize("value", HOSTILE)
def test_hostile_values_are_refused_without_a_server_error(made_up, value):
    for field in ("query", "program", "genome", "evalue", "max_hits", "task"):
        for response in (made_up.client.post("/blast", data={"query": PROTEIN, field: value}),
                         made_up.client.post("/api/v1/blast", json={"query": PROTEIN,
                                                                     field: value})):
            assert response.status_code in (200, 400, 413), (field, response.status_code)
            body = response.get_data(as_text=True)
            assert "<script>alert(1)</script>" not in body and "<img src=x" not in body
            assert "Traceback" not in body


def test_seeded_fuzz_of_blast_requests_finds_no_server_error(made_up):
    """600 searches of random text in random fields, as a form and as JSON."""
    import random
    import string

    rng = random.Random(20261126)
    alphabet = string.printable + "äé漢字‮​<>\"'%;&#\\" + "ACGTMKV" * 20
    fields = ("query", "program", "genome", "evalue", "max_hits", "task", "other")
    failures = []
    for round_ in range(300):
        values = {"query": PROTEIN if round_ % 3 else ""}
        for name in rng.sample(fields, rng.randint(1, len(fields))):
            length = rng.choice([0, 1, 5, 20, 120, 1000])
            values[name] = "".join(rng.choice(alphabet) for _ in range(length))
        for response in (made_up.client.post("/blast", data=values),
                         made_up.client.post("/api/v1/blast", json=values)):
            if response.status_code >= 500:
                failures.append((values, response.status_code))
    assert failures == []


def test_gene_pages_offer_blast_here_and_at_ncbi(made_up):
    body = made_up.client.get(f"/g/1/f/{made_up.dscam}").get_data(as_text=True)
    assert f'href="/blast?from=1:{made_up.dscam}:protein"' in body
    assert "blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&amp;PROGRAM=blastp" in body
    prefilled = made_up.client.get(f"/blast?from=1:{made_up.dscam}:protein")
    assert prefilled.status_code == 200
    assert "&gt;NP_001260764.1\nMNMPNERLKWLMLFAAVALIACGSQ" in prefilled.get_data(as_text=True)
    assert made_up.client.get("/blast?from=1:x:protein").status_code == 400
    assert made_up.client.get("/blast?from=1:999999:protein").status_code == 404
    assert made_up.client.get("/blast?from=9:1:protein").status_code == 404


def test_without_blast_the_page_says_so_and_searches_are_refused(portal):
    client = create_app(portal, count_in_background=False).test_client()
    body = client.get("/blast").get_data(as_text=True)
    assert "BLAST is not set up on this portal" in body and "<form method=\"post\"" not in body
    assert '<a href="/blast">BLAST</a>' not in client.get("/").get_data(as_text=True)
    response = client.post("/api/v1/blast", json={"query": PROTEIN})
    assert response.status_code == 503
    assert client.post("/blast", data={"query": PROTEIN}).status_code == 503


def test_a_failed_build_leaves_no_database(built, portal, tmp_path):
    genome = portal.genome(3)
    old = genome.directory / "blast"
    old.mkdir(exist_ok=True)
    (old / "genome.nsq").write_text("from an earlier load")
    results = blasting.build(portal, [genome], str(tmp_path / "no-makeblastdb"),
                             log_line=lambda line: None)
    assert results[0].check == "FAIL" and results[0].problems
    assert not old.exists() and not (genome.directory / "blast.new").exists()


def test_the_blast_options_must_leave_threads_for_the_pages(built):
    from cgload.cli.commands.serve import blast_jobs_for

    db, data = built
    result = CliRunner().invoke(cli, ["serve", "--db-url", f"sqlite:///{db}", "--threads", "4",
                                      "--blast-jobs", "2"])
    assert result.exit_code == 2
    assert "--blast-jobs is at most 1" in result.output and "Traceback" not in result.output
    # The default is lowered to fit, never refused: --threads alone still starts.
    assert blast_jobs_for(8, 2, given=False) == 2
    assert blast_jobs_for(4, 2, given=False) == 1
    assert blast_jobs_for(2, 2, given=False) == 0
    assert blast_jobs_for(8, 3, given=True) == 3
    assert blast_jobs_for(8, 0, given=True) == 0


# --------------------------------------------------------------------------
# With BLAST+
# --------------------------------------------------------------------------


@pytest.fixture(scope="module")
def blasted(tmp_path_factory):
    if not HAVE_BLAST:
        pytest.skip("BLAST+ is not installed")
    db, data = build(tmp_path_factory.mktemp("portal-blast-real"))
    result = CliRunner().invoke(cli, ["blastdb", "--db-url", f"sqlite:///{db}"])
    assert result.exit_code == 0, result.output
    assert result.output.count(": PASS") == 3
    portal = PortalData(read_only_engine(f"sqlite:///{db}"))
    yield SimpleNamespace(db=db, data=data, portal=portal,
                          blast=Blast(portal, jobs=2, threads=2, timeout=60))
    portal.close()


def _cds(portal, assembly_id: int):
    """Each CDS feature's row, span and protein, as the portal shows them."""
    with portal.engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.feature_id).where(feature.c.assembly_id == assembly_id,
                                               feature.c.feature_type == "CDS",
                                               feature.c.segment_index == 0)).scalars().all()
    genome = portal.genome(assembly_id)
    found = []
    for row in rows:
        family = portal.family(genome, row)
        sequences = portal.sequences_of(family, family.selected)
        found.append((row, family.selected, sequences["protein"].rstrip("*"), sequences["cds"]))
    return found


@needs_blast
def test_blastdb_counts_every_sequence_and_a_protein_for_every_cds(blasted):
    for genome in blasted.portal.genomes():
        manifest = json.loads((genome.directory / "blast" / "manifest.json").read_text())
        assert manifest["genome"]["sequences"] == manifest["genome"]["in_database"] \
            == manifest["genome"]["in_blast"] == genome.sequences
        assert manifest["proteins"]["translated"] == manifest["proteins"]["cds_features"] \
            == manifest["proteins"]["in_blast"] == len(_cds(blasted.portal, genome.assembly_id))
        assert manifest["genome"]["check"] == manifest["proteins"]["check"] == "PASS"
        assert manifest["loaded_at"] == genome.loaded_at.isoformat()
    assert sorted(blasted.blast.databases) == [1, 2, 3] and blasted.blast.ready


@needs_blast
def test_each_dscam1_isoform_finds_itself_first(blasted):
    blast = blasted.blast
    isoforms = [entry for entry in _cds(blasted.portal, 1) if len(entry[2]) > 1000]
    assert len(isoforms) == 3
    for row, _, protein, _ in isoforms:
        result = blast.search(make_request({"query": protein, "program": "blastp"},
                                           blast.searchable))
        best = result.hits[0]
        assert (best.assembly_id, best.subject) == (1, str(row))
        assert best.best.pident == 100.0 and best.best.length == len(protein)
        assert best.best.qstart == best.best.sstart == 1
        assert {hit.assembly_id for hit in result.hits} == {1}


@needs_blast
def test_tblastn_finds_the_protein_on_its_genome_inside_the_cds(blasted):
    blast = blasted.blast
    row, item, protein, _ = max(_cds(blasted.portal, 1), key=lambda entry: len(entry[2]))
    result = blast.search(make_request({"query": protein, "program": "tblastn", "genome": "1"},
                                       blast.searchable))
    hit = result.hits[0]
    assert hit.subject == item.seqid and hit.best.strand == "minus"
    for hsp in hit.hsps:
        assert item.start <= min(hsp.sstart, hsp.send) <= max(hsp.sstart, hsp.send) <= item.end


@needs_blast
def test_blastn_finds_a_piece_of_the_genome_where_it_lies(blasted):
    from cgload.portal.sequence import reverse_complement

    blast, portal = blasted.blast, blasted.portal
    genome = portal.genome(1)
    seqid = "NT_033778.4"
    piece = portal.sequence(genome).fetch(seqid, 20_001, 20_300)
    for query, strand, start, end in ((piece, "plus", 20_001, 20_300),
                                      (reverse_complement(piece), "minus", 20_300, 20_001)):
        for task in ("megablast", "blastn"):
            result = blast.search(make_request({"query": query, "program": "blastn",
                                                "task": task}, blast.searchable))
            best = result.hits[0]
            assert (best.assembly_id, best.subject, best.best.strand) == (1, seqid, strand)
            assert (best.best.sstart, best.best.send) == (start, end)
            assert best.best.pident == 100.0 and best.best.length == 300


@needs_blast
def test_blastx_finds_the_protein_of_a_coding_sequence(blasted):
    """On the fly: the Fusarium slice's made-up sequence translates to runs of
    one amino acid, which blastx's low-complexity filter rightly masks."""
    blast = blasted.blast
    row, _, protein, cds = max(_cds(blasted.portal, 1), key=lambda entry: len(entry[2]))
    assert len(cds) <= blasting.MAX_QUERY
    result = blast.search(make_request({"query": cds, "program": "blastx"}, blast.searchable))
    best = result.hits[0]
    assert (best.assembly_id, best.subject) == (1, str(row))
    assert best.best.pident == 100.0 and best.best.length == len(protein)


@needs_blast
def test_a_search_through_the_portal_links_to_the_hit(blasted):
    blast = blasted.blast
    row, _, protein, _ = max(_cds(blasted.portal, 1), key=lambda entry: len(entry[2]))
    client = create_app(blasted.portal, count_in_background=False, blast=blast).test_client()
    body = client.post("/blast", data={"query": protein, "program": "blastp"}).get_data(
        as_text=True)
    assert f'href="/g/1/f/{row}"' in body and "100.0%" in body
    answer = client.post("/api/v1/blast", json={"query": protein}).get_json()
    assert answer["hits"][0]["url"] == f"/g/1/f/{row}"


@needs_blast
def test_building_again_replaces_the_databases_and_leaves_nothing_behind(blasted):
    result = CliRunner().invoke(cli, ["blastdb", "--db-url", f"sqlite:///{blasted.db}",
                                      "--assembly-id", "1"])
    assert result.exit_code == 0 and "PASS" in result.output
    directory = blasted.portal.genome(1).directory
    assert (directory / "blast" / "manifest.json").exists()
    assert not (directory / "blast.new").exists() and not (directory / "blast.old").exists()
    unknown = CliRunner().invoke(cli, ["blastdb", "--db-url", f"sqlite:///{blasted.db}",
                                       "--assembly-id", "99"])
    assert unknown.exit_code == 1 and "no assembly 99" in unknown.output
