#!/usr/bin/env python3
"""BLAST in the portal, measured on the benchmark's genomes (PROTOCOL.md, section 14).

    python3 benchmark/run_blast.py --work DIR --results DIR [ACCESSION ...]

Added with cgload 0.4.0, after ``run_portal.py``, whose database it reads:
``--work`` is that run's work directory, holding ``portal.db`` and ``data/``
with the genomes loaded. Nothing there changes but each assembly's ``blast/``
directory, which ``cgload blastdb`` writes. Every search goes through
``cgload serve`` with its default BLAST settings, over HTTP on this machine:

1. **build**: ``cgload blastdb`` for each genome in turn, timed; the counts and
   checks its manifest records, and the databases' size; ``blast_build.tsv``.
2. **searches**, one at a time, each asking whether a query finds itself:

   * blastp: ``--proteins`` CDS per genome, chosen with the accession as seed,
     their proteins as the portal serves them, against every genome;
   * blastn (megablast): ``--pieces`` stretches of 500 bp per genome, at seeded
     places with no N, against every genome;
   * tblastn: ``--own`` of those proteins against their own genome, and
     ``--across`` proteins against every genome;
   * blastx: ``--own`` of those CDS against their own genome's proteins.

   The rank of the query's own feature or place among the hits, its identity
   and length, the hits that tie with it, the status and the time;
   ``blast_searches.tsv``.
3. **load**: blastp searches from 1, 2, 4 and 8 clients at once for
   ``--seconds`` each, while one more client opens gene pages, and the pages
   alone first; searches done and turned away (503), and the times of both;
   ``blast_concurrency.tsv``.
4. **security probe**: hostile values in each field of both BLAST URLs, as a
   form and as JSON, ``--fuzz`` seeded random searches, and POST, PUT and DELETE
   where the portal answers only GET; ``blast_security.tsv``. The database
   file's SHA-256, and every BLAST manifest's, before step 2 and after step 4;
   ``blast_readonly.tsv``.

A stopped run started again skips the genomes already built; the other steps
run again.
"""

from __future__ import annotations

import argparse
import contextlib
import csv
import hashlib
import json
import os
import platform
import random
import shutil
import statistics
import string
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402
from run_cgload import TIMEOUT  # noqa: E402
from run_portal import (  # noqa: E402
    HUMAN,
    INJECTIONS,
    SECURITY_HEADERS,
    YEAST,
    Server,
    cgload_version,
    fingerprint,
    http,
    loaded_accessions,
    log,
    percentile,
)

FLY = "GCF_000001215.4"
PIECE = 500
#: Values a visitor could send to reach a shell or BLAST's own options.
COMMAND_LIKE = ["-remote", "blastp -remote", "; rm -rf /", "$(reboot)", "`id`",
                "| cat /etc/passwd", "-out /tmp/x", "1e-3 -outfmt 5", "-db /etc/passwd",
                "all; DROP TABLE feature", "../../../../etc/passwd", "\n>second\nMKV"]
FIELDS = ("query", "program", "genome", "evalue", "max_hits", "task")


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def post(base: str, path: str, payload: dict, *, as_json: bool = True, timeout: float = 900
         ) -> tuple[int, int, float, dict, bytes]:
    """(status, bytes, milliseconds, headers, body) of one POST."""
    if as_json:
        data, kind = json.dumps(payload).encode(), "application/json"
    else:
        data = urllib.parse.urlencode(payload, doseq=True).encode()
        kind = "application/x-www-form-urlencoded"
    request = urllib.request.Request(base + path, data=data,  # noqa: S310 - http only
                                     headers={"Content-Type": kind}, method="POST")
    started = time.perf_counter()
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
            body = response.read()
            status, found = response.status, dict(response.headers)
    except urllib.error.HTTPError as error:
        body = error.read()
        status, found = error.code, dict(error.headers)
    except Exception:  # noqa: BLE001 - a refused or reset connection is counted, not raised
        return -1, 0, (time.perf_counter() - started) * 1000, {}, b""
    return status, len(body), (time.perf_counter() - started) * 1000, found, body


def blast(base: str, payload: dict) -> tuple[int, float, dict]:
    status, _, ms, _, body = post(base, "/api/v1/blast", payload)
    try:
        answer = json.loads(body) if body else {}
    except ValueError:
        answer = {}
    return status, ms, answer


def fasta_body(body: bytes) -> str:
    return "".join(body.decode(errors="replace").splitlines()[1:]).strip()


def blast_version() -> str:
    binary = shutil.which("blastp")
    if binary is None:
        return "not installed"
    shown = subprocess.run([binary, "-version"], capture_output=True, text=True, check=False)
    return shown.stdout.splitlines()[0].split(":", 1)[-1].strip() if shown.stdout else "unknown"


def portal_of(url: str, data: Path):
    from cgload.portal.data import PortalData, read_only_engine

    return PortalData(read_only_engine(url), data_dir=data)


def manifests_digest(portal) -> str:
    digest = hashlib.sha256()
    for genome in sorted(portal.genomes(), key=lambda item: item.assembly_id):
        path = genome.directory / "blast" / "manifest.json" if genome.directory else None
        if path is not None and path.exists():
            digest.update(path.read_bytes())
    return digest.hexdigest()


# --------------------------------------------------------------------------
# 1. build
# --------------------------------------------------------------------------


def build_all(ids: dict[str, int], url: str, data: Path, args, table: Path) -> None:
    steps = args.work / "steps-blast"
    portal = portal_of(url, data)
    for accession, assembly_id in ids.items():
        if bc.done(table, accession=accession):
            continue
        log(f"blastdb {accession}")
        record = bc.measured(f"blastdb-{accession}",
                             ["cgload", "blastdb", "--db-url", url, "--data-dir", str(data),
                              "--assembly-id", str(assembly_id)], steps, timeout=TIMEOUT)
        directory = portal.genome(assembly_id).directory / "blast"
        manifest = {}
        with contextlib.suppress(OSError, ValueError):
            manifest = json.loads((directory / "manifest.json").read_text())
        genome, proteins = manifest.get("genome", {}), manifest.get("proteins", {})
        bc.append_row(table, {
            "accession": accession, "assembly_id": assembly_id,
            "exit_code": record["exit_code"], "wall_s": record["wall_s"],
            "peak_rss_mb": record["peak_rss_mb"],
            "sequences_in_database": genome.get("in_database"),
            "sequences_read": genome.get("sequences"), "sequences_in_blast": genome.get("in_blast"),
            "genome_letters": genome.get("letters"), "genome_check": genome.get("check", "none"),
            "cds_features": proteins.get("cds_features"), "translated": proteins.get("translated"),
            "proteins_in_blast": proteins.get("in_blast"),
            "with_internal_stop": proteins.get("with_internal_stop"),
            "protein_letters": proteins.get("letters"),
            "proteins_check": proteins.get("check", "none"),
            "bytes": bc.tree_bytes(directory) if directory.exists() else 0,
            "problems": bc.first_lines(bc.text_of(steps, f"blastdb-{accession}"), 4)
            if record["exit_code"] != 0 else ""})
    portal.close()


# --------------------------------------------------------------------------
# 2. searches
# --------------------------------------------------------------------------


def sample_cds(url: str, assembly_id: int, accession: str, size: int) -> list[int]:
    from cgload.db.schema import feature
    from sqlalchemy import create_engine, select

    engine = create_engine(url)
    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.feature_id)
            .where(feature.c.assembly_id == assembly_id, feature.c.feature_type == "CDS",
                   feature.c.segment_index == 0)
            .order_by(feature.c.feature_id)).scalars().all()
    engine.dispose()
    return random.Random(f"{accession} blastp").sample(rows, min(size, len(rows)))


def sequences_of(url: str, assembly_id: int) -> list[tuple[str, int]]:
    from cgload.db.schema import sequence_region
    from sqlalchemy import create_engine, select

    engine = create_engine(url)
    with engine.connect() as connection:
        rows = connection.execute(
            select(sequence_region.c.source_id, sequence_region.c.length)
            .where(sequence_region.c.assembly_id == assembly_id)
            .order_by(sequence_region.c.sequence_region_id)).all()
    engine.dispose()
    return [(name, length) for name, length in rows]


def _row(program: str, accession: str, query: str, length: int, status: int, ms: float,
         answer: dict, scope: str, own_rank, own: dict | None, ties: int, kind: str,
         *, at_top: bool = True) -> dict:
    """One search. It found itself when its own feature or place is among the
    hits at 100% identity, and, where ``at_top``, with the top score (ties
    included: an identical copy elsewhere scores the same)."""
    hits = answer.get("hits", [])
    top = hits[0] if hits else None
    return {
        "program": program, "accession": accession, "scope": scope, "query": query,
        "kind": kind, "length": length, "status": status, "ms": round(ms, 1),
        "server_s": answer.get("seconds"), "hits": len(hits),
        "top_genome": top["genome"] if top else None,
        "top_subject": top["subject"] if top else None,
        "top_bits": top["hsps"][0]["bitscore"] if top else None,
        "own_rank": own_rank, "own_identity": own["identity"] if own else None,
        "own_length": own["length"] if own else None, "ties_at_top": ties,
        "found": own_rank is not None and own is not None and own["identity"] == 100.0
        and top is not None and (not at_top or own["bitscore"] == top["hsps"][0]["bitscore"]),
        "error": answer.get("error", "")[:200] if status != 200 else "",
    }


def _ties(hits: list[dict]) -> int:
    if not hits:
        return 0
    best = hits[0]["hsps"][0]["bitscore"]
    return sum(1 for hit in hits for hsp in hit["hsps"] if hsp["bitscore"] == best)


def search_proteins(base: str, ids: dict[str, int], url: str, args, table: Path
                    ) -> dict[str, list[dict]]:
    """blastp, and the proteins and CDS that tblastn and blastx then use."""
    chosen: dict[str, list[dict]] = {}
    for accession, assembly_id in ids.items():
        log(f"blastp {accession}")
        chosen[accession] = []
        for row in sample_cds(url, assembly_id, accession, args.proteins):
            status, _, _, _, body = http(base, f"/g/{assembly_id}/f/{row}/protein.fa")
            protein = fasta_body(body).rstrip("*") if status == 200 else ""
            status, _, _, _, body = http(base, f"/api/v1/genomes/{assembly_id}/features/{row}")
            item = {}
            if status == 200:
                with contextlib.suppress(ValueError, KeyError):
                    item = json.loads(body)["feature"]
            kind = ("stop_inside" if "*" in protein else
                    "ambiguous" if "X" in protein else "plain")
            if not 10 <= len(protein) <= 10_000:
                bc.append_row(table, _row("blastp", accession, str(row), len(protein), 0, 0,
                                          {}, "all", None, None, 0, "skipped_length"))
                continue
            status, ms, answer = blast(base, {"query": protein, "program": "blastp",
                                              "genome": "all", "max_hits": 100})
            hits = answer.get("hits", [])
            rank = next((index for index, hit in enumerate(hits, start=1)
                         if hit["genome"] == assembly_id and hit["subject"] == str(row)), None)
            own = hits[rank - 1]["hsps"][0] if rank else None
            bc.append_row(table, _row("blastp", accession, str(row), len(protein), status, ms,
                                      answer, "all", rank, own, _ties(hits), kind))
            chosen[accession].append({"row": row, "protein": protein, "kind": kind,
                                      "seqid": item.get("seqid"), "start": item.get("start"),
                                      "end": item.get("end")})
    return chosen


def search_pieces(base: str, ids: dict[str, int], url: str, args, table: Path) -> None:
    for accession, assembly_id in ids.items():
        log(f"blastn {accession}")
        sequences = [(name, length) for name, length in sequences_of(url, assembly_id)
                     if length >= 4 * PIECE]
        if not sequences:
            continue
        rng = random.Random(f"{accession} blastn")
        weights = [length for _, length in sequences]
        for _ in range(args.pieces):
            for _attempt in range(20):  # a stretch with no N
                seqid, length = rng.choices(sequences, weights=weights)[0]
                start = rng.randint(1, length - PIECE + 1)
                region = urllib.parse.quote(f"{seqid}:{start}-{start + PIECE - 1}", safe="")
                status, _, _, _, body = http(base, f"/g/{assembly_id}/region.fa?r={region}")
                bases = fasta_body(body).upper() if status == 200 else ""
                if len(bases) == PIECE and "N" not in bases:
                    break
            else:
                continue
            status, ms, answer = blast(base, {"query": bases, "program": "blastn",
                                              "task": "megablast", "genome": "all",
                                              "max_hits": 100})
            hits = answer.get("hits", [])
            rank, own = None, None
            for index, hit in enumerate(hits, start=1):
                if hit["genome"] != assembly_id or hit["subject"] != seqid:
                    continue
                own = next((hsp for hsp in hit["hsps"] if hsp["subject_start"] == start
                            and hsp["subject_end"] == start + PIECE - 1), None)
                if own is not None:
                    rank = index
                    break
            if own is not None and own["length"] != PIECE:
                own = dict(own, identity=-own["identity"])  # its place, but not all of it
            bc.append_row(table, _row("blastn", accession, f"{seqid}:{start}", PIECE, status, ms,
                                      answer, "all", rank, own, _ties(hits), "piece"))


def _inside(hit: dict, query: dict) -> dict | None:
    """The best HSP of a genome hit that lies in the CDS's span, at 100%."""
    if hit["subject"] != query["seqid"]:
        return None
    for hsp in hit["hsps"]:
        low, high = sorted((hsp["subject_start"], hsp["subject_end"]))
        if query["start"] <= low and high <= query["end"] and hsp["identity"] == 100.0:
            return hsp
    return None


def search_translated(base: str, ids: dict[str, int], chosen: dict[str, list[dict]], args,
                      table: Path) -> None:
    for accession, assembly_id in ids.items():
        # Proteins with no stop or X inside first; the others only if too few.
        queries = sorted((query for query in chosen.get(accession, []) if query["seqid"]),
                         key=lambda query: query["kind"] != "plain")[:args.own]
        log(f"tblastn and blastx {accession}")
        for query in queries:
            status, ms, answer = blast(base, {"query": query["protein"], "program": "tblastn",
                                              "genome": str(assembly_id), "max_hits": 25})
            hits = answer.get("hits", [])
            rank, own = None, None
            for index, hit in enumerate(hits, start=1):
                own = _inside(hit, query)
                if own is not None:
                    rank = index
                    break
            bc.append_row(table, _row("tblastn", accession, str(query["row"]),
                                      len(query["protein"]), status, ms, answer, "own", rank,
                                      own, _ties(hits), query["kind"], at_top=False))
            status, _, _, _, body = http(base, f"/g/{assembly_id}/f/{query['row']}/cds.fa")
            cds = fasta_body(body) if status == 200 else ""
            if not 10 <= len(cds) <= 10_000:
                continue
            status, ms, answer = blast(base, {"query": cds, "program": "blastx",
                                              "genome": str(assembly_id), "max_hits": 25})
            hits = answer.get("hits", [])
            rank = next((index for index, hit in enumerate(hits, start=1)
                         if hit["subject"] == str(query["row"])), None)
            own = hits[rank - 1]["hsps"][0] if rank else None
            bc.append_row(table, _row("blastx", accession, str(query["row"]), len(cds), status,
                                      ms, answer, "own", rank, own, _ties(hits), query["kind"]))
    # tblastn against every genome: the slowest search a visitor can ask for.
    preferred = [accession for accession in (HUMAN, FLY, YEAST) if accession in chosen]
    preferred += [accession for accession in chosen if accession not in preferred]
    across = []
    for accession in preferred:
        query = next(iter(sorted((query for query in chosen[accession] if query["seqid"]),
                                 key=lambda query: query["kind"] != "plain")), None)
        if query is not None:
            across.append((accession, query))
        if len(across) == args.across:
            break
    for accession, query in across:
        log(f"tblastn {accession} against every genome")
        status, ms, answer = blast(base, {"query": query["protein"], "program": "tblastn",
                                          "genome": "all", "max_hits": 25})
        hits = answer.get("hits", [])
        rank, own = None, None
        for index, hit in enumerate(hits, start=1):
            own = _inside(hit, query) if hit["genome"] == ids[accession] else None
            if own is not None:
                rank = index
                break
        bc.append_row(table, _row("tblastn", accession, str(query["row"]), len(query["protein"]),
                                  status, ms, answer, "all", rank, own, _ties(hits), query["kind"],
                                  at_top=False))


# --------------------------------------------------------------------------
# 3. load
# --------------------------------------------------------------------------


def load_test(base: str, chosen: dict[str, list[dict]], ids: dict[str, int], levels: list[int],
              seconds: float, table: Path) -> None:
    proteins = [query["protein"] for queries in chosen.values() for query in queries
                if len(query["protein"]) <= 2000]
    pages = [f"/g/{ids[accession]}/f/{query['row']}" for accession, queries in chosen.items()
             for query in queries]
    if not proteins or not pages:
        return

    def run(level: int) -> dict:
        stop = threading.Event()
        lock = threading.Lock()
        searches: list[tuple[int, float]] = []
        page_times: list[tuple[int, float]] = []

        def searcher(seed: int) -> None:
            rng = random.Random(seed)
            while not stop.is_set():
                status, ms, _ = blast(base, {"query": rng.choice(proteins), "program": "blastp",
                                             "genome": "all"})
                with lock:
                    searches.append((status, ms))
                if status == 503:  # turned away: a person would try again shortly
                    stop.wait(1.0)

        def reader() -> None:
            rng = random.Random(level)
            while not stop.is_set():
                status, _, ms, _, _ = http(base, rng.choice(pages))
                with lock:
                    page_times.append((status, ms))

        workers = [threading.Thread(target=searcher, args=(level * 100 + index,))
                   for index in range(level)] + [threading.Thread(target=reader)]
        for worker in workers:
            worker.start()
        time.sleep(seconds)
        stop.set()
        for worker in workers:
            worker.join()
        done = [ms for status, ms in searches if status == 200]
        served = [ms for status, ms in page_times if status == 200]
        return {
            "clients": level, "seconds": seconds, "searches": len(searches),
            "searches_done": len(done),
            "turned_away_503": sum(1 for status, _ in searches if status == 503),
            "other_errors": sum(1 for status, _ in searches if status not in (200, 503)),
            "search_p50_ms": statistics.median(done) if done else None,
            "search_p95_ms": percentile(done, 0.95),
            "pages": len(page_times), "page_errors": len(page_times) - len(served),
            "page_p50_ms": statistics.median(served) if served else None,
            "page_p95_ms": percentile(served, 0.95), "page_max_ms": max(served, default=None),
        }

    for level in levels:
        log(f"load test: {level} searching clients")
        bc.append_row(table, run(level))


# --------------------------------------------------------------------------
# 4. security probe
# --------------------------------------------------------------------------


def security(base: str, assembly_id: int, protein: str, fuzz: int, table: Path) -> None:
    log("security probe")
    statuses: dict[int, int] = {}
    failures: list[str] = []
    missing = 0
    sent = 0
    base_payload = {"query": protein, "program": "blastp", "genome": str(assembly_id)}

    def check(label: str, status: int, headers: dict, body: bytes) -> None:
        nonlocal missing, sent
        sent += 1
        statuses[status] = statuses.get(status, 0) + 1
        if status >= 0 and not all(name in headers for name in SECURITY_HEADERS):
            missing += 1
        busy = status == 503 and "Retry-After" in headers
        if status < 0 or (status >= 500 and not busy) or b"<script>alert(1)</script>" in body \
                or b"<svg onload" in body:
            failures.append(f"{status} {label[:120]}")

    for value in INJECTIONS + COMMAND_LIKE:
        for field in FIELDS:
            payload = {**base_payload, field: value}
            for path, as_json in (("/blast", False), ("/api/v1/blast", True)):
                status, _, _, headers, body = post(base, path, payload, as_json=as_json)
                check(f"{path} {field}={value!r}", status, headers, body)
    rng = random.Random(20261126)
    alphabet = string.printable + "äé漢字‮<>\"'%;&#\\" + "ACDEFGHIKLMNPQRSTVWY" * 4
    for _ in range(fuzz // 2):
        payload = dict(base_payload)
        for field in rng.sample(FIELDS, rng.randint(1, len(FIELDS))):
            payload[field] = "".join(rng.choice(alphabet)
                                     for _ in range(rng.choice([0, 1, 8, 64, 600])))
        for path, as_json in (("/blast", False), ("/api/v1/blast", True)):
            status, _, _, headers, body = post(base, path, payload, as_json=as_json)
            check(f"{path} fuzz", status, headers, body)
    for method, path in (("POST", "/"), ("PUT", "/"), ("DELETE", "/"),
                         ("POST", f"/g/{assembly_id}/"), ("POST", "/api/v1/genomes"),
                         ("PUT", "/blast"), ("DELETE", "/api/v1/blast"),
                         ("GET", "/api/v1/blast")):
        request = urllib.request.Request(base + path, method=method)  # noqa: S310
        try:
            with urllib.request.urlopen(request, timeout=30):  # noqa: S310 - http only
                status = 200
        except urllib.error.HTTPError as error:
            status = error.code
        sent += 1
        statuses[status] = statuses.get(status, 0) + 1
        if status != 405:
            failures.append(f"{method} {path} answered {status}")
    bc.append_row(table, {
        "requests": sent, "statuses": json.dumps(dict(sorted(statuses.items()))),
        "server_errors": sum(count for status, count in statuses.items()
                             if status >= 500 and status != 503) + statuses.get(-1, 0),
        "busy_503": statuses.get(503, 0), "missing_security_headers": missing,
        "failures": "; ".join(failures[:20])})


# --------------------------------------------------------------------------


def environment(results: Path, url: str, args) -> None:
    info = {
        "started": time.strftime("%Y-%m-%d"), "cgload": cgload_version(),
        "commit": bc.git_commit(bc.HERE.parent), "blast": blast_version(),
        "python": platform.python_version(), "platform": platform.platform(),
        "cpus": os.cpu_count(), "database": url,
        "serve": {"threads": args.threads, "blast_jobs": 2, "blast_threads": 4,
                  "blast_timeout_s": 120, "note": "cgload serve's defaults"},
    }
    with contextlib.suppress(OSError, StopIteration):
        info["memory_kb"] = next(int(line.split()[1]) for line in
                                 Path("/proc/meminfo").read_text().splitlines()
                                 if line.startswith("MemTotal"))
    (results / "blast_environment.json").write_text(json.dumps(info, indent=1) + "\n")


def summary(results: Path) -> str:
    def read(name: str) -> list[dict]:
        path = results / name
        if not path.exists():
            return []
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    lines = []
    build = read("blast_build.tsv")
    def number(value: str) -> float:
        return float(value) if value not in ("", "NA", None) else 0.0

    passed = sum(1 for row in build if row["genome_check"] == row["proteins_check"] == "PASS")
    proteins = sum(number(row["proteins_in_blast"]) for row in build)
    seconds = sum(number(row["wall_s"]) for row in build)
    gigabytes = sum(number(row["bytes"]) for row in build) / 1e9
    lines.append(f"build: {len(build)} genomes, both checks PASS {passed}, {proteins:,.0f} "
                 f"proteins, {seconds:,.0f} s, {gigabytes:.2f} GB")
    searches = read("blast_searches.tsv")
    for program, scope in sorted({(row["program"], row["scope"]) for row in searches}):
        rows = [row for row in searches if row["program"] == program and row["scope"] == scope
                and row["kind"] != "skipped_length"]
        answered = [row for row in rows if row["status"] == "200"]
        times = [float(row["ms"]) / 1000 for row in answered]
        found = sum(1 for row in answered if row["found"] == "True")
        plain = [row for row in answered if row["kind"] in ("plain", "piece")]
        plain_found = sum(1 for row in plain if row["found"] == "True")
        unique_missed = sum(1 for row in plain if row["found"] != "True"
                            and int(row["ties_at_top"]) <= 1)
        unanswered = sorted({row["status"] for row in rows if row["status"] != "200"})
        lines.append(
            f"{program:<8} {scope:<4} n={len(rows):>4}  answered {len(answered):>4}  "
            f"found itself {found}  (plain {plain_found} of {len(plain)}; missed with no tie "
            f"{unique_missed})  "
            + (f"median {statistics.median(times):.1f} s, p95 {percentile(times, 0.95):.1f} s, "
               f"max {max(times):.1f} s" if times else "no times")
            + f"  not answered: {', '.join(unanswered) or 'none'}")
    for row in read("blast_concurrency.tsv"):
        lines.append(f"load {row['clients']:>2} clients: {row['searches_done']} searches done, "
                     f"{row['turned_away_503']} turned away, {row['other_errors']} errors; "
                     f"search p50 {row['search_p50_ms']} ms; pages {row['pages']}, "
                     f"p50 {row['page_p50_ms']} ms, p95 {row['page_p95_ms']} ms")
    for row in read("blast_security.tsv"):
        lines.append(f"security probe: {row['requests']} requests, statuses {row['statuses']}, "
                     f"server errors {row['server_errors']}, missing headers "
                     f"{row['missing_security_headers']}, failures: {row['failures'] or 'none'}")
    for row in read("blast_readonly.tsv"):
        lines.append(f"database unchanged: {row['database_unchanged']}; BLAST manifests "
                     f"unchanged: {row['manifests_unchanged']}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accessions", nargs="*")
    parser.add_argument("--work", type=Path, required=True,
                        help="run_portal.py's work directory: portal.db and data/")
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--port", type=int, default=8791)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--proteins", type=int, default=20, help="blastp queries per genome")
    parser.add_argument("--pieces", type=int, default=10, help="blastn queries per genome")
    parser.add_argument("--own", type=int, default=3,
                        help="tblastn and blastx queries per genome, against that genome")
    parser.add_argument("--across", type=int, default=3,
                        help="tblastn queries against every genome")
    parser.add_argument("--seconds", type=float, default=60, help="per load-test level")
    parser.add_argument("--fuzz", type=int, default=600, help="random searches in step 4")
    args = parser.parse_args(argv)
    if shutil.which("blastp") is None or shutil.which("makeblastdb") is None:
        raise SystemExit("BLAST+ is not on PATH")
    database = args.work / "portal.db"
    if not database.exists():
        raise SystemExit(f"no {database}: run run_portal.py first")
    args.results.mkdir(parents=True, exist_ok=True)
    url = f"sqlite:///{database}"
    data = args.work / "data"
    ids = loaded_accessions(url)
    accessions = args.accessions or [genome["accession"] for genome in bc.genomes()]
    ids = {accession: ids[accession] for accession in accessions if accession in ids}
    environment(args.results, url, args)

    build_all(ids, url, data, args, args.results / "blast_build.tsv")
    for name in ("blast_searches.tsv", "blast_concurrency.tsv", "blast_security.tsv",
                 "blast_readonly.tsv"):
        (args.results / name).unlink(missing_ok=True)
    portal = portal_of(url, data)
    log("fingerprints")
    before = fingerprint(database), manifests_digest(portal)
    with Server(url, data, args.port, args.threads, args.results / "serve-blast.log") as server:
        chosen = search_proteins(server.base, ids, url, args, args.results / "blast_searches.tsv")
        search_pieces(server.base, ids, url, args, args.results / "blast_searches.tsv")
        search_translated(server.base, ids, chosen, args, args.results / "blast_searches.tsv")
        load_test(server.base, chosen, ids, [0, 1, 2, 4, 8], args.seconds,
                  args.results / "blast_concurrency.tsv")
        smallest = min(ids.values(), key=lambda aid: sum(length for _, length
                                                         in sequences_of(url, aid)))
        protein = next((query["protein"] for queries in chosen.values() for query in queries),
                       "MKVLAAGIVGLLLAQSAWAEDPKTSRELIEKAVRDF")
        security(server.base, smallest, protein[:300], args.fuzz,
                 args.results / "blast_security.tsv")
    after = fingerprint(database), manifests_digest(portal)
    portal.close()
    bc.append_row(args.results / "blast_readonly.tsv", {
        "database": database.name, "sha256_before": before[0], "sha256_after": after[0],
        "database_unchanged": before[0] == after[0], "manifests_before": before[1],
        "manifests_after": after[1], "manifests_unchanged": before[1] == after[1]})
    text = summary(args.results)
    (args.results / "blast_summary.txt").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
