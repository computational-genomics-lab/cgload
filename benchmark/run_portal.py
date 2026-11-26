#!/usr/bin/env python3
"""The portal, measured as its users meet it (PROTOCOL.md, section 13).

    python3 benchmark/run_portal.py --inputs DIR --work DIR --results DIR [ACCESSION ...]

Added after the benchmark's runs, with the portal (cgload 0.3.0); it changes
nothing in their results. One database holds all 30 genomes, as a lab's portal
would, and every step reads it through ``cgload serve``, the production
server, over HTTP on this machine:

1. **load**: ``cgload init``, then ``cgload load`` of each genome into the one
   SQLite database (with ``--repair-envelopes`` where the file is refused, as
   in the benchmark); ``portal_loads.tsv``.
2. **index**: ``cgload index``'s build and count check, per genome, timed;
   ``portal_index.tsv``.
3. **pages**: for a fixed random sample of genes per genome (seeded by the
   accession), one client asks what a user asks: search by the gene's name,
   open its page, its region, its protein, and two byte ranges of the files
   the genome browser reads; every response's status, size and time go to
   ``portal_requests.tsv``. Whether the search found the gene, and on which
   result page, goes to ``portal_found.tsv``.
4. **load test**: the same requests from 1, 8 and 32 clients at once for 60 s
   each; ``portal_concurrency.tsv``.
5. **translations**: up to 2,000 CDS per genome, chosen at random, translated
   as the portal shows them, and classified: no stop codon inside; a stop
   inside where the file declares a translation exception, a partial CDS or a
   pseudogene; or a stop inside with none of these; ``portal_translation.tsv``.
6. **security probe**: injection strings and 3,000 seeded random requests sent
   to the running server; every status counted, and every response checked
   for the security headers; ``portal_security.tsv``. The database file's
   SHA-256 is taken after indexing and again after steps 3 to 6 and the
   translations, which read it; ``portal_readonly.tsv``.
7. **PostgreSQL** (with ``--postgres-conninfo``): steps 1 to 3 for the genomes named
   by ``--postgres-genomes``, on PostgreSQL; ``portal_requests_postgres.tsv``.
8. **demo**: a small database of *S. cerevisiae* and of human chromosome 17,
   with their files, packed as ``demo.tar.gz``, for the paper's screenshots.

Every step appends its rows as it goes; a stopped run started again skips the
genomes already loaded and indexed and repeats the rest, which is quick.
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
import tarfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402
from run_cgload import TIMEOUT, write_config  # noqa: E402

SECURITY_HEADERS = ("Content-Security-Policy", "X-Content-Type-Options", "X-Frame-Options")
INJECTIONS = [
    "' OR '1'='1", "' OR 1=1 --", "1; DROP TABLE feature; --",
    "x' UNION SELECT source_id, 1 FROM feature --", "<script>alert(1)</script>",
    "\"><svg onload=alert(1)>", "../../../../etc/passwd", "{{7*7}}", "%00", "A" * 5000,
]
YEAST = "GCF_000146045.2"
HUMAN = "GCF_000001405.40"
HUMAN_CHR17 = "NC_000017.11"


# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------


def log(message: str) -> None:
    print(time.strftime("%H:%M:%S"), message, flush=True)


def cgload_version() -> str:
    shown = subprocess.run(["cgload", "--version"], capture_output=True, text=True,
                           check=False).stdout.split()
    return shown[-1] if shown else "unknown"


def http(base: str, path: str, *, headers: dict | None = None, timeout: float = 120
         ) -> tuple[int, int, float, dict, bytes]:
    """(status, bytes, milliseconds, headers, body) of one GET; redirects not followed."""
    request = urllib.request.Request(base + path, headers=headers or {})  # noqa: S310
    opener = urllib.request.build_opener(_NoRedirect)
    started = time.perf_counter()
    try:
        with opener.open(request, timeout=timeout) as response:  # noqa: S310 - http only
            body = response.read()
            status, found = response.status, dict(response.headers)
    except urllib.error.HTTPError as error:
        body = error.read()
        status, found = error.code, dict(error.headers)
    except Exception:  # noqa: BLE001 - a refused or reset connection is counted, not raised
        return -1, 0, (time.perf_counter() - started) * 1000, {}, b""
    return status, len(body), (time.perf_counter() - started) * 1000, found, body


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):  # noqa: ARG002
        return None


def quoted(text: str) -> str:
    return urllib.parse.quote(text, safe="")


def fingerprint(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def percentile(values: list[float], share: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(round(share * (len(ordered) - 1))))]


# --------------------------------------------------------------------------
# 1. load, 2. index
# --------------------------------------------------------------------------


def loaded_accessions(url: str) -> dict[str, int]:
    from cgload.db.schema import assembly
    from sqlalchemy import create_engine, select

    engine = create_engine(url)
    try:
        with engine.connect() as connection:
            return {row.accession: row.assembly_id for row in connection.execute(
                select(assembly.c.accession, assembly.c.assembly_id))}
    except Exception:  # noqa: BLE001 - no schema yet
        return {}
    finally:
        engine.dispose()


def load_all(accessions: list[str], args, url: str, data: Path, table: Path) -> dict[str, int]:
    steps = args.work / "steps"
    if not loaded_accessions(url):
        bc.measured("init", ["cgload", "init"], steps, env={"CGLOAD_DB_URL": url})
    for accession in accessions:
        if accession in loaded_accessions(url):
            continue
        genome = bc.genome(accession)
        inputs = args.inputs / accession
        config = args.work / "configs" / f"{accession}.toml"
        config.parent.mkdir(parents=True, exist_ok=True)
        write_config(config, genome, inputs)
        command = ["cgload", "load", "--config", str(config), "--data-dir", str(data)]
        log(f"load {accession}")
        record = bc.measured(f"load-{accession}", command, steps, timeout=TIMEOUT,
                             env={"CGLOAD_DB_URL": url})
        outcome = "loaded"
        if record["exit_code"] != 0:
            record = bc.measured(f"load-{accession}-repair", [*command, "--repair-envelopes"],
                                 steps, timeout=TIMEOUT, env={"CGLOAD_DB_URL": url})
            outcome = "loaded_with_repair" if record["exit_code"] == 0 else "refused"
        bc.append_row(table, {"accession": accession, "backend": url.split(":", 1)[0],
                              "outcome": outcome, "wall_s": record["wall_s"],
                              "peak_rss_mb": record["peak_rss_mb"]})
    return loaded_accessions(url)


def index_all(ids: dict[str, int], url: str, table: Path) -> None:
    from cgload.portal.search import build_index, indexed_assemblies
    from sqlalchemy import create_engine

    engine = create_engine(url)
    with engine.connect() as connection:
        done = set(indexed_assemblies(connection))
    for accession, assembly_id in ids.items():
        if assembly_id in done:
            continue
        log(f"index {accession}")
        started = time.perf_counter()
        report = build_index(engine, [assembly_id])[0]
        bc.append_row(table, {
            "accession": accession, "backend": url.split(":", 1)[0],
            "assembly_id": assembly_id, "features": report.features, "terms": report.terms,
            "keys_checked": len(report.counts), "check": "PASS" if report.ok else "FAIL",
            "failures": "; ".join(f"{c.key}:{c.in_attributes}/{c.in_index}"
                                  for c in report.failures()),
            "seconds": round(time.perf_counter() - started, 2)})
    engine.dispose()


# --------------------------------------------------------------------------
# 3. pages
# --------------------------------------------------------------------------


class Server:
    """``cgload serve`` in a child process, for the length of a ``with``."""

    def __init__(self, url: str, data: Path, port: int, threads: int, log_path: Path) -> None:
        self.command = ["cgload", "serve", "--host", "127.0.0.1", "--db-url", url,
                        "--data-dir", str(data),
                        "--port", str(port), "--threads", str(threads)]
        self.base = f"http://127.0.0.1:{port}"
        self.log_path = log_path
        self.process: subprocess.Popen | None = None

    def __enter__(self) -> Server:
        handle = self.log_path.open("a")
        env = dict(os.environ)
        env.pop("CGLOAD_DB_URL", None)
        self.process = subprocess.Popen(self.command, stdout=handle, stderr=handle, env=env)
        for _ in range(600):
            if self.process.poll() is not None:
                raise SystemExit(f"cgload serve stopped; see {self.log_path}")
            if http(self.base, "/healthz", timeout=5)[0] == 200:
                return self
            time.sleep(0.5)
        raise SystemExit(f"cgload serve did not answer; see {self.log_path}")

    def __exit__(self, *exc) -> None:
        if self.process is not None:
            self.process.terminate()
            try:
                self.process.wait(timeout=30)
            except subprocess.TimeoutExpired:
                self.process.kill()


def sample_genes(base: str, assembly_id: int, accession: str, size: int, url: str) -> list[dict]:
    """Up to ``size`` named genes of one genome, chosen with the accession as seed,
    each with its page row, name, location and first CDS (if any)."""
    from cgload.db.schema import feature, feature_attribute, sequence_region
    from sqlalchemy import create_engine, select

    engine = create_engine(url)
    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.feature_id, feature_attribute.c.value, sequence_region.c.source_id,
                   feature.c.start, feature.c.end)
            .join(feature_attribute, feature_attribute.c.feature_id == feature.c.feature_id)
            .join(sequence_region,
                  sequence_region.c.sequence_region_id == feature.c.sequence_region_id)
            .where(feature.c.assembly_id == assembly_id, feature.c.segment_index == 0,
                   feature.c.feature_type == "gene", feature_attribute.c.key == "Name")
            .order_by(feature.c.feature_id)
        ).all()
    engine.dispose()
    # The same genes on every run and backend: the accession seeds the choice.
    chosen = random.Random(accession).sample(rows, min(size, len(rows)))
    genes = []
    for row in chosen:
        status, _, _, _, body = http(base, f"/api/v1/genomes/{assembly_id}/features/"
                                           f"{row.feature_id}")
        cds = None
        if status == 200:
            family = json.loads(body)["family"]
            cds = next((item["row"] for item in family if item["type"] == "CDS"), None)
        genes.append({"row": row.feature_id, "name": row.value, "seqid": row.source_id,
                      "start": row.start, "end": row.end, "cds": cds})
    return genes


def requests_for(assembly_id: int, gene: dict, sizes: dict[str, int]) -> list[tuple[str, str,
                                                                                    dict]]:
    """The requests a user makes about one gene: (kind, path, headers)."""
    flank = 5000
    start = max(gene["start"] - flank, 1)
    end = min(gene["end"] + flank, start + 100_000)
    region = quoted(f"{gene['seqid']}:{start}-{end}")
    out = [
        ("search", f"/g/{assembly_id}/search?q={quoted(gene['name'])}", {}),
        ("gene_page", f"/g/{assembly_id}/f/{gene['row']}", {}),
        ("region_page", f"/g/{assembly_id}/region?r={region}", {}),
    ]
    if gene["cds"] is not None:
        out.append(("protein", f"/g/{assembly_id}/f/{gene['cds']}/protein.fa", {}))
    rng = random.Random(gene["row"])
    for name in ("annotation.gff3.gz", "assembly.fa"):
        size = sizes.get(name)
        if size:
            offset = rng.randrange(0, max(size - 65536, 1))
            out.append((f"range_{name.split('.')[0]}", f"/g/{assembly_id}/files/{name}",
                        {"Range": f"bytes={offset}-{offset + 65535}"}))
    return out


def file_sizes(base: str, assembly_id: int) -> dict[str, int]:
    sizes = {}
    for name in ("annotation.gff3.gz", "assembly.fa"):
        status, _, _, headers, _ = http(base, f"/g/{assembly_id}/files/{name}",
                                        headers={"Range": "bytes=0-0"})
        if status == 206 and "Content-Range" in headers:
            sizes[name] = int(headers["Content-Range"].rsplit("/", 1)[1])
    return sizes


def found_on_page(base: str, assembly_id: int, gene: dict, pages: int = 20
                  ) -> tuple[int | None, int | None]:
    """The result page on which a search by the gene's name lists it (None if
    none of the first ``pages`` does), and how many features the search found."""
    total = None
    for page in range(1, pages + 1):
        status, _, _, _, body = http(base, f"/api/v1/genomes/{assembly_id}/search?q="
                                           f"{quoted(gene['name'])}&page={page}")
        if status != 200:
            return None, total
        answer = json.loads(body)
        total = answer["total"] if total is None else total
        if any(item["feature"]["row"] == gene["row"] for item in answer["results"]):
            return page, total
        if len(answer["results"]) < 50:
            return None, total
    return None, total


def pages(base: str, ids: dict[str, int], args, url: str, backend: str, out_requests: Path,
          out_found: Path | None) -> list[tuple[str, str, dict]]:
    workload: list[tuple[str, str, dict]] = []
    for accession, assembly_id in ids.items():
        log(f"pages {backend} {accession}")
        genes = sample_genes(base, assembly_id, accession, args.sample, url)
        sizes = file_sizes(base, assembly_id)
        rows = []
        for gene in genes:
            for kind, path, headers in requests_for(assembly_id, gene, sizes):
                status, size, ms, _, _ = http(base, path, headers=headers)
                rows.append({"accession": accession, "backend": backend, "kind": kind,
                             "status": status, "bytes": size, "ms": round(ms, 2)})
                workload.append((kind, path, headers))
            if out_found is not None:
                page, matches = found_on_page(base, assembly_id, gene)
                bc.append_row(out_found, {"accession": accession, "row": gene["row"],
                                          "name": gene["name"], "found": page is not None,
                                          "page": page, "matches": matches})
        for row in rows:
            bc.append_row(out_requests, row)
    return workload


# --------------------------------------------------------------------------
# 4. load test
# --------------------------------------------------------------------------


def load_test(base: str, workload: list[tuple[str, str, dict]], levels: list[int],
              seconds: float, table: Path) -> None:
    for clients in levels:
        log(f"load test, {clients} clients, {seconds:.0f} s")
        stop = time.perf_counter() + seconds
        lock = threading.Lock()
        times: list[float] = []
        statuses: dict[int, int] = {}

        def client(seed: int, stop=stop, lock=lock, times=times, statuses=statuses) -> None:
            rng = random.Random(seed)
            while time.perf_counter() < stop:
                _, path, headers = rng.choice(workload)
                status, _, ms, _, _ = http(base, path, headers=headers)
                with lock:
                    times.append(ms)
                    statuses[status] = statuses.get(status, 0) + 1

        started = time.perf_counter()
        with ThreadPoolExecutor(max_workers=clients) as pool:
            list(pool.map(client, range(clients)))
        elapsed = time.perf_counter() - started
        errors = sum(count for status, count in statuses.items() if status >= 500 or status < 0)
        bc.append_row(table, {
            "clients": clients, "seconds": round(elapsed, 1), "requests": len(times),
            "per_second": round(len(times) / elapsed, 1), "errors": errors,
            "statuses": json.dumps(dict(sorted(statuses.items()))),
            "p50_ms": round(percentile(times, 0.5) or 0, 1),
            "p95_ms": round(percentile(times, 0.95) or 0, 1),
            "p99_ms": round(percentile(times, 0.99) or 0, 1)})


# --------------------------------------------------------------------------
# 5. translations
# --------------------------------------------------------------------------


def translations(ids: dict[str, int], url: str, data: Path, per_genome: int, table: Path
                 ) -> None:
    from cgload.db.schema import feature
    from cgload.portal.data import PortalData, exceptions, read_only_engine
    from sqlalchemy import select

    portal = PortalData(read_only_engine(url), data_dir=data)
    for accession, assembly_id in ids.items():
        log(f"translations {accession}")
        genome = portal.genome(assembly_id)
        with portal.engine.connect() as connection:
            rows = sorted(connection.execute(select(feature.c.feature_id).where(
                feature.c.assembly_id == assembly_id, feature.c.feature_type == "CDS",
                feature.c.segment_index == 0)).scalars())
        chosen = random.Random(accession).sample(rows, min(per_genome, len(rows)))
        counts = {"checked": 0, "no_internal_stop": 0, "stop_declared_exception": 0,
                  "stop_partial": 0, "stop_pseudo": 0, "stop_unexplained": 0, "failed": 0}
        examples = []
        for row in chosen:
            try:
                family = portal.family(genome, row)
                item = family.selected
                protein = portal.sequences_of(family, item)["protein"]
            except Exception as exc:  # noqa: BLE001 - counted, with the first few kept
                counts["failed"] += 1
                if len(examples) < 5:
                    examples.append(f"{row}: {exc}")
                continue
            counts["checked"] += 1
            if "*" not in protein.rstrip("*"):
                counts["no_internal_stop"] += 1
            elif exceptions(item):
                counts["stop_declared_exception"] += 1
            elif any(key in item.attributes for key in ("start_range", "end_range", "partial")):
                counts["stop_partial"] += 1
            elif "pseudo" in item.attributes or family.top.attributes.get("pseudo"):
                counts["stop_pseudo"] += 1
            else:
                counts["stop_unexplained"] += 1
                if len(examples) < 5:
                    examples.append(item.source_id or str(row))
        bc.append_row(table, {"accession": accession, "cds_in_genome": len(rows), **counts,
                              "examples": "; ".join(examples)})
    portal.close()


# --------------------------------------------------------------------------
# 6. security probe
# --------------------------------------------------------------------------


def security(base: str, assembly_id: int, table: Path, fuzz: int = 3000) -> None:
    log("security probe")
    paths = []
    for value in INJECTIONS:
        q = quoted(value)
        paths += [f"/search?q={q}", f"/g/{assembly_id}/search?q={q}",
                  f"/g/{assembly_id}/region?r={q}", f"/api/v1/genomes/{assembly_id}/search?q={q}",
                  f"/g/{assembly_id}/id/{q}", f"/g/{assembly_id}/files/{q}"]
    rng = random.Random(20261125)
    alphabet = string.printable + "äé漢字‮<>\"'%;&#\\"
    for _ in range(fuzz // 6):
        value = quoted("".join(rng.choice(alphabet) for _ in range(rng.choice([1, 8, 64, 199]))))
        number = str(rng.choice([0, 1, 7, 10**9, 2**63]))
        paths += [f"/g/{assembly_id}/search?q={value}", f"/g/{assembly_id}/region?r={value}",
                  f"/g/{number}/f/{number}", f"/api/v1/genomes/{number}/features/{number}",
                  f"/g/{assembly_id}/locate?q={value}", f"/g/{assembly_id}/files/{value}"]
    statuses: dict[int, int] = {}
    missing_headers = 0
    failures = []
    for path in paths:
        status, _, _, headers, body = http(base, path)
        statuses[status] = statuses.get(status, 0) + 1
        if status >= 0 and not all(name in headers for name in SECURITY_HEADERS):
            missing_headers += 1
        if status >= 500 or status < 0 or b"<script>alert(1)</script>" in body:
            failures.append(f"{status} {path[:120]}")
    for method in ("POST", "PUT", "DELETE"):
        request = urllib.request.Request(base + "/", method=method)  # noqa: S310
        try:
            urllib.request.urlopen(request, timeout=30)  # noqa: S310 - http only
            status = 200
        except urllib.error.HTTPError as error:
            status = error.code
        statuses[status] = statuses.get(status, 0) + 1
        if status != 405:
            failures.append(f"{method} / answered {status}")
    bc.append_row(table, {"requests": len(paths) + 3,
                          "statuses": json.dumps(dict(sorted(statuses.items()))),
                          "server_errors": sum(c for s, c in statuses.items() if s >= 500 or s < 0),
                          "missing_security_headers": missing_headers,
                          "failures": "; ".join(failures[:20])})


# --------------------------------------------------------------------------
# 8. demo
# --------------------------------------------------------------------------


def demo(args, results: Path) -> None:
    """Yeast whole, and human chromosome 17, each in a small database of its own."""
    from pyfaidx import Fasta

    target = args.work / "demo"
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    url = f"sqlite:///{target / 'demo.db'}"
    data = target / "data"
    subprocess.run(["cgload", "init", "--db-url", url], check=True)
    for accession, seqid in ((YEAST, None), (HUMAN, HUMAN_CHR17)):
        inputs = args.inputs / accession
        if not (inputs / "genomic.gff").exists():
            continue
        genome = bc.genome(accession)
        part = target / accession
        part.mkdir()
        # The inputs are copied beside the demo first, so the paths its database
        # and receipts record are the demo's own, not this machine's input folders.
        fasta, gff = part / "assembly.fa", part / "genomic.gff"
        if seqid is None:
            shutil.copyfile(inputs / "assembly.fa", fasta)
            shutil.copyfile(inputs / "genomic.gff", gff)
        else:
            with Fasta(str(inputs / "assembly.fa"), as_raw=True) as handle, fasta.open("w") as out:
                bases = str(handle[seqid][:])
                out.write(f">{seqid}\n")
                for index in range(0, len(bases), 60):
                    out.write(bases[index : index + 60] + "\n")
            with (inputs / "genomic.gff").open() as source, gff.open("w") as out:
                for line in source:
                    if line.startswith("##sequence-region") and seqid not in line:
                        continue
                    if line.startswith("#") or line.startswith(seqid + "\t"):
                        out.write(line)
        config = part / "organism.toml"
        write_config(config, genome, part)
        if seqid is not None:
            text = config.read_text().replace(genome["assembly_name"],
                                              f"{genome['assembly_name']} chromosome 17")
            config.write_text(text)
        command = ["cgload", "load", "--config", str(config), "--db-url", url, "--data-dir",
                   str(data)]
        if subprocess.run(command, check=False).returncode != 0:
            subprocess.run([*command, "--repair-envelopes"], check=False)
    subprocess.run(["cgload", "index", "--db-url", url], check=False)
    for path in target.iterdir():
        if path.is_dir() and path.name.startswith("GCF_"):
            shutil.rmtree(path)
    with tarfile.open(results / "demo.tar.gz", "w:gz") as archive:
        archive.add(target, arcname="demo")


# --------------------------------------------------------------------------


def environment(results: Path, url: str) -> None:
    import importlib.metadata as metadata

    info = {
        "started": time.strftime("%Y-%m-%d"),
        "cgload": cgload_version(), "commit": bc.git_commit(bc.HERE.parent),
        "python": platform.python_version(), "platform": platform.platform(),
        "processor": platform.processor(), "cpus": os.cpu_count(), "database": url,
        "packages": {name: _version(metadata, name) for name in
                     ("flask", "werkzeug", "waitress", "jinja2", "sqlalchemy", "pyfaidx",
                      "psycopg")},
    }
    with contextlib.suppress(OSError, StopIteration):
        info["memory_kb"] = next(int(line.split()[1]) for line in
                                 Path("/proc/meminfo").read_text().splitlines()
                                 if line.startswith("MemTotal"))
    (results / "portal_environment.json").write_text(json.dumps(info, indent=1) + "\n")


def _version(metadata, name: str) -> str | None:
    try:
        return metadata.version(name)
    except metadata.PackageNotFoundError:
        return None


def summary(results: Path) -> str:
    lines = []

    def read(name: str) -> list[dict]:
        path = results / name
        if not path.exists():
            return []
        return list(csv.DictReader(path.open(), delimiter="\t"))

    loads = read("portal_loads.tsv")
    lines.append(f"loads: {len(loads)}, " + ", ".join(
        f"{outcome} {sum(1 for row in loads if row['outcome'] == outcome)}"
        for outcome in sorted({row["outcome"] for row in loads})))
    index = read("portal_index.tsv")
    passed = sum(row["check"] == "PASS" for row in index)
    lines.append(f"index: {len(index)} genomes, check PASS {passed}"
                 f", {sum(int(r['features']) for r in index):,} features, "
                 f"{sum(int(r['terms']) for r in index):,} terms, "
                 f"{sum(float(r['seconds']) for r in index):,.0f} s")
    requests = read("portal_requests.tsv")
    for backend in sorted({row["backend"] for row in requests}):
        lines.append(f"pages ({backend}):")
        for kind in sorted({row["kind"] for row in requests if row["backend"] == backend}):
            times = [float(row["ms"]) for row in requests
                     if row["backend"] == backend and row["kind"] == kind]
            bad = sum(1 for row in requests if row["backend"] == backend and row["kind"] == kind
                      and not 200 <= int(row["status"]) < 400)
            lines.append(f"  {kind:<20} n={len(times):>6}  median {statistics.median(times):8.1f}"
                         f" ms  p95 {percentile(times, 0.95):8.1f} ms  max {max(times):9.1f} ms"
                         f"  not 2xx/3xx {bad}")
    found = read("portal_found.tsv")
    if found:
        hits = sum(row["found"] == "True" for row in found)
        missed = [row for row in found if row["found"] != "True"]
        few = sum(1 for row in missed if row["matches"] and int(row["matches"]) <= 1000)
        lines.append(f"search found the gene: {hits} of {len(found)}; not found with at most "
                     f"1,000 matches (a failure, section 13): {few}")
    for row in read("portal_concurrency.tsv"):
        lines.append(f"load test {row['clients']:>2} clients: {row['requests']} requests, "
                     f"{float(row['per_second']):.1f}/s, errors {row['errors']}, "
                     f"p50 {float(row['p50_ms']):.1f} ms, p95 {float(row['p95_ms']):.1f} ms, "
                     f"p99 {float(row['p99_ms']):.1f} ms")
    translation = read("portal_translation.tsv")
    if translation:
        total = {key: sum(int(row[key]) for row in translation) for key in
                 ("checked", "no_internal_stop", "stop_declared_exception", "stop_partial",
                  "stop_pseudo", "stop_unexplained", "failed")}
        lines.append("translations: " + ", ".join(f"{key} {value:,}"
                                                  for key, value in total.items()))
    for row in read("portal_security.tsv"):
        lines.append(f"security probe: {row['requests']} requests, statuses {row['statuses']},"
                     f" server errors {row['server_errors']}, missing headers "
                     f"{row['missing_security_headers']}, failures: {row['failures'] or 'none'}")
    for row in read("portal_readonly.tsv"):
        lines.append(f"database unchanged by serving: {row['unchanged']}")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accessions", nargs="*")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--port", type=int, default=8790)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--sample", type=int, default=200, help="genes per genome")
    parser.add_argument("--translations", type=int, default=2000, help="CDS per genome")
    parser.add_argument("--seconds", type=float, default=60, help="per load-test level")
    parser.add_argument("--fuzz", type=int, default=3000, help="random requests in step 6")
    parser.add_argument("--postgres-conninfo", default=None,
                        help="libpq connection string of a server for step 7, as the "
                        "benchmark's BENCH_PG_CONNINFO; a new database is made there")
    parser.add_argument("--postgres-genomes", nargs="*", default=[HUMAN])
    parser.add_argument("--no-demo", action="store_true")
    args = parser.parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=True)
    args.results.mkdir(parents=True, exist_ok=True)
    accessions = args.accessions or [genome["accession"] for genome in bc.genomes()]

    url = f"sqlite:///{args.work / 'portal.db'}"
    data = args.work / "data"
    environment(args.results, url)
    ids = load_all(accessions, args, url, data, args.results / "portal_loads.tsv")
    ids = {accession: ids[accession] for accession in accessions if accession in ids}
    index_all(ids, url, args.results / "portal_index.tsv")
    for name in ("portal_requests.tsv", "portal_found.tsv", "portal_concurrency.tsv",
                 "portal_translation.tsv", "portal_security.tsv", "portal_readonly.tsv"):
        (args.results / name).unlink(missing_ok=True)
    database = args.work / "portal.db"
    log("database fingerprint")
    before = fingerprint(database)
    with Server(url, data, args.port, args.threads, args.results / "serve.log") as server:
        workload = pages(server.base, ids, args, url, "sqlite",
                         args.results / "portal_requests.tsv", args.results / "portal_found.tsv")
        load_test(server.base, workload, [1, 8, 32], args.seconds,
                  args.results / "portal_concurrency.tsv")
        security(server.base, next(iter(ids.values())), args.results / "portal_security.tsv",
                 fuzz=args.fuzz)
    translations(ids, url, data, args.translations, args.results / "portal_translation.tsv")
    after = fingerprint(database)
    bc.append_row(args.results / "portal_readonly.tsv", {
        "database": database.name, "bytes": database.stat().st_size, "sha256_before": before,
        "sha256_after": after, "unchanged": before == after})

    if args.postgres_conninfo:
        from run_cgload import postgres

        wanted = [accession for accession in args.postgres_genomes if accession in accessions]
        pg_url = postgres(args.postgres_conninfo, "cgload_portal")
        pg_data = args.work / "data-postgres"
        shutil.rmtree(pg_data, ignore_errors=True)
        pg_ids = load_all(wanted, args, pg_url, pg_data, args.results / "portal_loads.tsv")
        pg_ids = {accession: pg_ids[accession] for accession in wanted if accession in pg_ids}
        index_all(pg_ids, pg_url, args.results / "portal_index.tsv")
        with Server(pg_url, pg_data, args.port + 1, args.threads,
                    args.results / "serve-postgres.log") as server:
            pages(server.base, pg_ids, args, pg_url, "postgresql",
                  args.results / "portal_requests.tsv", None)

    if not args.no_demo:
        log("demo")
        demo(args, args.results)
    text = summary(args.results)
    (args.results / "portal_summary.txt").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
