#!/usr/bin/env python3
"""Make the paper's two Zenodo drafts, each with a reserved DOI, and upload their files.

    export ZENODO_TOKEN=...          # scopes deposit:write and deposit:actions
    python3 paper/deposit_zenodo.py --repo REPO --results RESULTS_TARBALL

1. Software: the source archive of the release tag (git archive), MIT.
2. Dataset: the benchmark results tarball as benchmark/run.sh packed it,
   superseded results and the PostgreSQL server log included, and the census
   results (survey/results at the same tag), under CC0, which GigaScience
   asks for where a repository offers it.

Nothing is published. A draft can be edited or discarded at zenodo.org, and
its DOI is only reserved; publishing makes it permanent, and is done by hand
when the paper is submitted.

This makes NEW records. The records made for v0.2.0 on 10 November 2026 exist
already: for a later tag, open the software record at zenodo.org and choose
"New version" (or, while it is still a draft, replace its file), with the
archive that `git archive --format=tar.gz --prefix=cgload-X.Y.Z/ vX.Y.Z` makes.
The results did not change after v0.2.0, so the dataset record needs no new
file. The token is read from the environment and never
printed or written anywhere. Python standard library only.
"""

from __future__ import annotations

import argparse
import io
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request
from pathlib import Path

CREATORS = [
    {"name": "Upadhyay, Aditya",
     "affiliation": "CSIR-Indian Institute of Chemical Biology, Kolkata; AcSIR, Ghaziabad"},
    {"name": "Panda, Arijit", "affiliation": "Mayo Clinic, Rochester, MN"},
    {"name": "Tripathy, Sucheta",
     "affiliation": "CSIR-Indian Institute of Chemical Biology, Kolkata; AcSIR, Ghaziabad"},
]
REPOSITORY = "https://github.com/computational-genomics-lab/cgload"


class Zenodo:
    def __init__(self, base: str, token: str):
        self.base = base.rstrip("/")
        self.token = token

    def call(self, method: str, url: str, body: bytes | None = None,
             content_type: str = "application/json") -> dict:
        if not url.startswith("http"):
            url = self.base + url
        # The token goes only to Zenodo over https (or to a local test server).
        if not url.startswith(("https://", "http://127.0.0.1:", "http://localhost:")):
            raise SystemExit(f"refusing to send the token to {url.split('?')[0]}")
        request = urllib.request.Request(url, data=body, method=method, headers={  # noqa: S310
            "Authorization": f"Bearer {self.token}", "Content-Type": content_type})
        try:
            with urllib.request.urlopen(request, timeout=600) as response:  # noqa: S310
                text = response.read().decode()
        except urllib.error.HTTPError as error:
            detail = error.read().decode(errors="replace")[:2000]
            raise SystemExit(f"Zenodo said {error.code} to {method} {url.split('?')[0]}:\n"
                             f"{detail}") from None
        return json.loads(text) if text else {}

    def draft(self, metadata: dict, files: list[Path]) -> dict:
        created = self.call("POST", "/deposit/depositions",
                            json.dumps({"metadata": metadata}).encode())
        bucket = created["links"]["bucket"]
        for path in files:
            with path.open("rb") as handle:
                self.call("PUT", f"{bucket}/{path.name}", handle.read(),
                          "application/octet-stream")
            print(f"  uploaded {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
        return created


def reserved_doi(deposition: dict) -> str:
    reserved = (deposition.get("metadata") or {}).get("prereserve_doi") or {}
    return reserved.get("doi") or f"10.5281/zenodo.{deposition['id']} (expected)"


def source_archive(repo: Path, tag: str, out: Path) -> Path:
    version = tag.removeprefix("v")
    target = out / f"cgload-{version}.tar.gz"
    subprocess.run(["git", "-C", str(repo), "archive", "--format=tar.gz",  # noqa: S603, S607
                    f"--prefix=cgload-{version}/", "-o", str(target), tag], check=True)
    return target


def census_archive(repo: Path, tag: str, out: Path) -> Path:
    """survey/results as the tag has it, so that the deposit matches the release."""
    listing = subprocess.run(["git", "-C", str(repo), "ls-tree", "-r", "--name-only",  # noqa: S603, S607
                              tag, "survey/results"], check=True, capture_output=True,
                             text=True).stdout.split()
    if not listing:
        raise SystemExit(f"{tag} has no survey/results")
    target = out / "cgload-census-results.tar.gz"
    with tarfile.open(target, "w:gz") as archive:
        for name in listing:
            blob = subprocess.run(["git", "-C", str(repo), "show", f"{tag}:{name}"],  # noqa: S603, S607
                                  check=True, capture_output=True).stdout
            info = tarfile.TarInfo(name)
            info.size = len(blob)
            info.mode = 0o644
            archive.addfile(info, io.BytesIO(blob))
    return target


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", type=Path, required=True, help="the cgload clone")
    parser.add_argument("--results", type=Path, required=True,
                        help="cgload-benchmark-results.tar.gz from the last run")
    parser.add_argument("--tag", default="v0.2.1")
    parser.add_argument("--zenodo", default="https://zenodo.org/api",
                        help="use https://sandbox.zenodo.org/api to rehearse")
    args = parser.parse_args(argv)
    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        raise SystemExit("set ZENODO_TOKEN first (scopes deposit:write, deposit:actions)")
    if not args.results.is_file():
        raise SystemExit(f"no results tarball at {args.results}")
    head = subprocess.run(["git", "-C", str(args.repo), "rev-parse", "--verify", "-q",  # noqa: S603, S607
                           f"{args.tag}^{{commit}}"], capture_output=True, text=True)
    if head.returncode:
        raise SystemExit(f"{args.repo} has no tag {args.tag}: create and push it first")
    version = args.tag.removeprefix("v")
    zenodo = Zenodo(args.zenodo, token)

    with tempfile.TemporaryDirectory() as scratch:
        scratch = Path(scratch)
        print(f"software draft: cgload {version}")
        software = zenodo.draft({
            "upload_type": "software",
            "title": f"cgload {version}",
            "creators": CREATORS,
            "description": (
                "<p>cgload loads a genome assembly and its annotation (GFF3, GenBank, "
                "eggNOG-mapper and InterProScan output) into SQLite, PostgreSQL or MySQL "
                "in one transaction, and checks the stored features and parent links "
                "against an independent count of the input files. A load that falls "
                "short is rolled back.</p>"
                f"<p>Source archive of tag {args.tag} of {REPOSITORY}.</p>"),
            "version": version,
            "license": "mit",
            "access_right": "open",
            "keywords": ["genome annotation", "GFF3", "database", "data integrity"],
            "related_identifiers": [{"identifier": f"{REPOSITORY}/tree/{args.tag}",
                                     "relation": "isIdenticalTo"}],
            "prereserve_doi": True,
        }, [source_archive(args.repo, args.tag, scratch)])
        software_doi = reserved_doi(software)

        results = scratch / "cgload-benchmark-results.tar.gz"
        results.symlink_to(args.results.resolve())
        print("dataset draft: benchmark and census results")
        dataset = zenodo.draft({
            "upload_type": "dataset",
            "title": f"Benchmark of cgload {version} and three other GFF3 loaders "
                     "on 30 RefSeq genomes: results",
            "creators": CREATORS,
            "description": (
                "<p>Results of the benchmark protocol in benchmark/PROTOCOL.md of "
                f"{REPOSITORY}: cgload on SQLite and PostgreSQL, gffutils, txdbmaker "
                "and python-chado on 30 RefSeq genomes, compared feature by feature with "
                "a reference reading of each file; planted faults; interrupted loads; "
                "query timing; and the census of the 30 files.</p>"
                "<p>cgload-benchmark-results.tar.gz is the results folder as run.sh "
                "packed it, with the results that later repeats replaced "
                "(results/superseded) and the PostgreSQL server log; the reference "
                "databases are left out, since benchmark/truth.py rebuilds them from "
                "the input files. cgload-census-results.tar.gz is survey/results at "
                f"tag {args.tag}. PROTOCOL.md section 12 says which tag each result "
                "comes from.</p>"),
            "version": "benchmark-run-2.2",
            "license": "cc0-1.0",
            "access_right": "open",
            "keywords": ["genome annotation", "GFF3", "benchmark", "RefSeq"],
            "related_identifiers": [
                {"identifier": software_doi.split()[0], "relation": "isSupplementTo"},
                {"identifier": REPOSITORY, "relation": "isSupplementTo"}],
            "prereserve_doi": True,
        }, [results, census_archive(args.repo, args.tag, scratch)])

    print()
    print(f"software  DOI {software_doi}  draft {software['links']['html']}")
    print(f"dataset   DOI {reserved_doi(dataset)}  draft {dataset['links']['html']}")
    print("Both are drafts: check them at zenodo.org, and publish when the paper is submitted.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
