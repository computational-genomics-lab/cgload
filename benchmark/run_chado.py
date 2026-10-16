#!/usr/bin/env python3
"""Chado in the loader comparison: python-chado's GFF3 loader, then read it back.

    python3 benchmark/run_chado.py setup  --dsn DSN --schema default_schema.sql --so so.obo
    python3 benchmark/run_chado.py load   --dsn DSN GFF --genus G --species S
    python3 benchmark/run_chado.py extract --dsn DSN STORED.sqlite

Runs in its own environment: python-chado 2.3.9 needs SQLAlchemy 1.4 and
Biopython 1.81, and cgload needs SQLAlchemy 2. DSN is a libpq connection string
for one empty database per genome, e.g. "dbname=chado_x host=/tmp/sock port=5432".

* ``setup`` loads the Chado 1.31 schema with psql (the GMOD repository's
  ``default_schema.sql``) and then the Sequence Ontology from an OBO file into
  the ``sequence`` vocabulary, which is where python-chado looks up every
  feature type and the ``part_of`` relationship. GMOD's own Docker image loads
  the schema the same way and gets the same nine errors on PostgreSQL 14 and
  later: helper functions for interval overlap use the ``box @ box`` operator,
  which PostgreSQL no longer has. No table the loader uses is affected; the
  errors are counted and reported.
* ``load`` is the timed step: ``load_gff`` with ``landmark_type="region"``
  (sequences are created as features of that type, as they have no other
  record), ``no_seq_compute=True`` (no sequence is loaded, so none can be
  computed), and ``add_only=True`` (the database is new; python-chado documents
  this as faster).
* ``extract`` reads the database with SQL: a feature's rows are its
  ``featureloc`` rows, its parents its ``feature_relationship`` objects. The
  landmark features ``load`` created for the sequences are left out.
  Attributes are not compared: python-chado spreads them over several tables.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gffmodel  # noqa: E402

STRAND = {1: "+", -1: "-", 0: ".", None: "."}


def _connect(dsn: str):
    import psycopg2

    return psycopg2.connect(dsn)


def parse_obo(path: Path) -> list[dict]:
    stanzas, current = [], None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("["):
                current = {"kind": line.strip("[]"), "synonym": []}
                stanzas.append(current)
                continue
            if current is None or not line or line.startswith("!"):
                continue
            key, _, value = line.partition(": ")
            if key == "synonym":
                found = re.match(r'"((?:[^"\\]|\\.)*)"\s+(\w+)', value)
                if found:
                    current["synonym"].append(
                        (found.group(1).replace('\\"', '"'), found.group(2).lower()))
            elif key in ("id", "name", "def", "is_obsolete") and key not in current:
                current[key] = value
    return stanzas


def load_ontology(dsn: str, obo: Path) -> int:
    """The Sequence Ontology into cv 'sequence', with its synonyms."""
    import psycopg2

    connection = _connect(dsn)
    cursor = connection.cursor()

    def one(sql, args=()):
        cursor.execute(sql, args)
        return cursor.fetchone()[0]

    db_so = one("INSERT INTO db (name) VALUES ('SO') RETURNING db_id")
    db_internal = one("INSERT INTO db (name) VALUES ('internal') ON CONFLICT (name) "
                      "DO UPDATE SET name = EXCLUDED.name RETURNING db_id")
    cv_sequence = one("INSERT INTO cv (name) VALUES ('sequence') RETURNING cv_id")
    cv_synonym = one("INSERT INTO cv (name) VALUES ('synonym_type') RETURNING cv_id")
    synonym_type = {}
    for name in ("exact", "narrow", "broad", "related"):
        xref = one("INSERT INTO dbxref (db_id, accession) VALUES (%s, %s) RETURNING dbxref_id",
                   (db_internal, "synonym_type:" + name))
        synonym_type[name] = one(
            "INSERT INTO cvterm (cv_id, name, dbxref_id) VALUES (%s, %s, %s) RETURNING cvterm_id",
            (cv_synonym, name, xref))
    loaded = 0
    for stanza in parse_obo(obo):
        if stanza["kind"] not in ("Term", "Typedef") or "id" not in stanza or "name" not in stanza:
            continue
        accession = stanza["id"]
        if accession.startswith("SO:"):
            accession = accession.split(":", 1)[1]
        xref = one("INSERT INTO dbxref (db_id, accession) VALUES (%s, %s) RETURNING dbxref_id",
                   (db_so, accession))
        cursor.execute("SAVEPOINT term")
        try:
            term = one(
                "INSERT INTO cvterm (cv_id, name, definition, dbxref_id, is_obsolete, "
                "is_relationshiptype) VALUES (%s, %s, %s, %s, %s, %s) RETURNING cvterm_id",
                (cv_sequence, stanza["name"], stanza.get("def", "")[:10000], xref,
                 1 if stanza.get("is_obsolete") == "true" else 0,
                 1 if stanza["kind"] == "Typedef" else 0))
        except psycopg2.Error:
            # An obsolete term may share its name with a current one.
            cursor.execute("ROLLBACK TO SAVEPOINT term")
            continue
        loaded += 1
        seen = set()
        for text, scope in stanza["synonym"]:
            if text in seen:
                continue
            seen.add(text)
            cursor.execute(
                "INSERT INTO cvtermsynonym (cvterm_id, synonym, type_id) VALUES (%s, %s, %s)",
                (term, text[:1024], synonym_type.get(scope, synonym_type["exact"])))
    connection.commit()
    connection.close()
    return loaded


def setup(args) -> int:
    done = subprocess.run(
        [args.psql, "-q", "-v", "ON_ERROR_STOP=0", "-d", args.dsn, "-f", str(args.schema)],
        capture_output=True, text=True, check=False)
    errors = [line for line in done.stderr.splitlines() if "ERROR" in line]
    terms = load_ontology(args.dsn, args.so)
    print(json.dumps({"schema_errors": len(errors), "schema_error_kinds": sorted(
        {re.sub(r"^.*ERROR:\s*", "", line) for line in errors}), "ontology_terms": terms}))
    return 0


def load(args) -> int:
    from chado import ChadoInstance

    url = args.url or "postgresql+psycopg2:///?" + args.dsn.replace(" ", "&")
    instance = ChadoInstance(dburl=url)
    organism = instance.organism.add_organism(genus=args.genus, species=args.species,
                                              common=args.genus, abbr=args.genus[:8])
    analysis = instance.analysis.add_analysis(name="benchmark", program="truth",
                                              programversion="1", sourcename=str(args.gff))
    started = time.monotonic()
    result = instance.feature.load_gff(str(args.gff), analysis["analysis_id"],
                                       organism["organism_id"], landmark_type="region",
                                       no_seq_compute=True, add_only=True, quiet=True)
    import chado

    print(json.dumps({"tool": "python-chado", "version": getattr(chado, "__version__", "?"),
                      "seconds": round(time.monotonic() - started, 3), "result": result}))
    return 0


def extract(args) -> dict:
    connection = _connect(args.dsn)
    cursor = connection.cursor()
    cursor.execute(
        "SELECT f.feature_id, f.uniquename, t.name, src.uniquename, fl.fmin, fl.fmax, "
        "fl.strand FROM feature f JOIN cvterm t ON t.cvterm_id = f.type_id "
        "LEFT JOIN featureloc fl ON fl.feature_id = f.feature_id "
        "LEFT JOIN feature src ON src.feature_id = fl.srcfeature_id")
    features: dict[int, dict] = {}
    for feature_id, uniquename, ftype, seqid, fmin, fmax, strand in cursor.fetchall():
        entry = features.setdefault(feature_id, {"key": uniquename, "type": ftype,
                                                 "seqid": seqid, "segs": [], "strands": set()})
        if fmin is not None:
            entry["segs"].append((fmin + 1, fmax))
            entry["strands"].add(STRAND.get(strand, "."))
            entry["seqid"] = entry["seqid"] or seqid
    cursor.execute("SELECT subject_id, object_id FROM feature_relationship")
    parents: dict[int, set[int]] = defaultdict(set)
    for child, parent in cursor.fetchall():
        parents[child].add(parent)
    connection.close()

    landmarks = {fid for fid, entry in features.items() if not entry["segs"]}
    anonymous_seen: dict[str, int] = {}

    def rows():
        for feature_id, entry in features.items():
            if feature_id in landmarks:
                continue
            names = ",".join(sorted(features[p]["key"] for p in parents.get(feature_id, ())
                                    if p in features and p not in landmarks))
            key = entry["key"]
            if key is None:
                key = gffmodel.unique_anon(gffmodel.anon_key(
                    entry["seqid"], entry["type"], entry["segs"][0][0], entry["segs"][0][1],
                    ",".join(sorted(entry["strands"])), names), anonymous_seen)
            low = min(a for a, _ in entry["segs"])
            high = max(b for _, b in entry["segs"])
            yield (key, entry["seqid"], entry["type"], ",".join(sorted(entry["strands"])),
                   low, high, len(entry["segs"]), gffmodel.segs_text(entry["segs"]), names,
                   "", 1, "")

    out = gffmodel.create(args.out)
    written = gffmodel.write_features(out, rows())
    gffmodel.set_meta(out, tool="python-chado", features=written, landmarks=len(landmarks))
    out.commit()
    out.close()
    return {"features": written, "landmarks": len(landmarks)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("setup")
    one.add_argument("--dsn", required=True)
    one.add_argument("--schema", type=Path, required=True)
    one.add_argument("--so", type=Path, required=True)
    one.add_argument("--psql", default="psql")
    two = commands.add_parser("load")
    two.add_argument("--dsn", required=True)
    two.add_argument("--url", help="SQLAlchemy URL, if the DSN cannot be translated")
    two.add_argument("gff", type=Path)
    two.add_argument("--genus", required=True)
    two.add_argument("--species", required=True)
    three = commands.add_parser("extract")
    three.add_argument("--dsn", required=True)
    three.add_argument("out", type=Path)
    args = parser.parse_args(argv)
    if args.command == "setup":
        return setup(args)
    if args.command == "load":
        return load(args)
    print(json.dumps(extract(args)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
