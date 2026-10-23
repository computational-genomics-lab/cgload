#!/usr/bin/env python3
"""gffutils in the loader comparison: build its database, then read it back.

    python3 benchmark/run_gffutils.py load GFF DB.sqlite --strategy create_unique
    python3 benchmark/run_gffutils.py extract DB.sqlite STORED.sqlite

``load`` calls ``gffutils.create_db`` the way its documentation describes, with
the merge strategy named; it is the step that is timed. ``extract`` reads the
finished database with plain SQL and writes the comparison table described in
gffmodel.py. It does not use gffutils to read, so it reports what the database
holds rather than what gffutils' API chooses to show.

How a gffutils database is read back:

* each row of its ``features`` table is one record; the record's own ``ID``
  attribute (kept even when gffutils renames the record, as ``create_unique``
  does for a second row with the same ID) says which feature of the file it
  holds, so records sharing one are one feature held as several records;
* a record's parents are the level-1 rows of ``relations``, translated back to
  the features those parent records hold;
* a record with no ``ID`` attribute is keyed exactly as truth.py keys a row
  without one.
"""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gffmodel  # noqa: E402


def load(gff: Path, db: Path, strategy: str) -> int:
    import gffutils  # imported here so `extract` runs without it

    if db.exists():
        db.unlink()
    started = time.monotonic()
    gffutils.create_db(str(gff), str(db), merge_strategy=strategy, force=True,
                       keep_order=False, verbose=False)
    print(json.dumps({"tool": "gffutils", "version": gffutils.__version__,
                      "merge_strategy": strategy,
                      "seconds": round(time.monotonic() - started, 3)}))
    return 0


def extract(db: Path, out: Path) -> dict:
    source = sqlite3.connect(db)
    records: dict[str, tuple] = {}
    key_of: dict[str, str] = {}
    anonymous_seen: dict[str, int] = {}
    pending_anonymous = []
    for record_id, seqid, ftype, start, end, strand, attributes in source.execute(
        "SELECT id, seqid, featuretype, start, end, strand, attributes FROM features "
        "ORDER BY seqid, start, end"
    ):
        parsed = json.loads(attributes) if attributes else {}
        parsed = {key: [unquote(v) for v in values] for key, values in parsed.items()}
        records[record_id] = (seqid, ftype, int(start), int(end), strand, parsed)
        original = parsed.get("ID", [None])[0]
        if original is None:
            pending_anonymous.append(record_id)
        else:
            key_of[record_id] = original
    parents_of: dict[str, set[str]] = defaultdict(set)
    for parent, child in source.execute("SELECT parent, child FROM relations WHERE level = 1"):
        parents_of[child].add(parent)
    source.close()

    def parent_keys(record_id):
        return {key_of.get(p, p) for p in parents_of.get(record_id, ())}

    for record_id in pending_anonymous:
        seqid, ftype, start, end, strand, parsed = records[record_id]
        parents = ",".join(sorted(parent_keys(record_id)))
        key_of[record_id] = gffmodel.unique_anon(
            gffmodel.anon_key(seqid, ftype, start, end, strand, parents), anonymous_seen)

    grouped: dict[str, list[str]] = defaultdict(list)
    for record_id, key in key_of.items():
        grouped[key].append(record_id)

    def rows():
        for key, members in grouped.items():
            first = records[members[0]]
            strands, segs, parents, rows_of = set(), [], set(), []
            for record_id in members:
                seqid, ftype, start, end, strand, parsed = records[record_id]
                strands.add(strand)
                segs.append((start, end))
                parents |= parent_keys(record_id)
                rows_of.append((start, end, parsed))
            low = min(min(a, b) for a, b in segs)
            high = max(max(a, b) for a, b in segs)
            yield (key, first[0], first[1], ",".join(sorted(strands)), low, high, len(segs),
                   gffmodel.segs_text(segs), ",".join(sorted(parents)),
                   gffmodel.feature_digest(rows_of), len(members), "")

    connection = gffmodel.create(out)
    written = gffmodel.write_features(connection, rows())
    gffmodel.set_meta(connection, tool="gffutils", records=len(records), features=written)
    connection.commit()
    connection.close()
    return {"records": len(records), "features": written}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("load")
    one.add_argument("gff", type=Path)
    one.add_argument("db", type=Path)
    one.add_argument("--strategy", default="create_unique",
                     choices=["error", "warning", "merge", "create_unique", "replace"])
    two = commands.add_parser("extract")
    two.add_argument("db", type=Path)
    two.add_argument("out", type=Path)
    args = parser.parse_args(argv)
    if args.command == "load":
        return load(args.gff, args.db, args.strategy)
    print(json.dumps(extract(args.db, args.out)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
