#!/usr/bin/env python3
"""The reference count: what a GFF3 file holds, read without any loader's code.

    python3 benchmark/truth.py INPUT.gff3[.gz] OUTPUT.sqlite [--summary OUT.json]

Reads the file and writes one row per feature into OUTPUT.sqlite (the table
described in gffmodel.py), plus the sequences and a summary. Every loader's
result is later compared with this.

Deliberately simple and standard-library only. It shares no code with cgload,
gffutils, txdbmaker or Chado, so it cannot repeat their mistakes. It follows the
GFF3 specification and nothing else: rows sharing an ID form one feature, a row
without an ID is a feature of its own, `Parent` names the parents, and
`##sequence-region` and the `region` rows give each sequence's length, whether
it is circular, and which genome it belongs to.

The same program reads the GFF3 that a loader writes back out (cgload's
export), so the input and the output of a round trip are measured one way.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gffmodel  # noqa: E402

TRANS_SPLICING = ("trans-splicing", "trans_splicing")


def md5_of_content(path: Path) -> str:
    """MD5 of the file's bytes, decompressed if gzipped: the survey's checksum."""
    digest = hashlib.md5(usedforsecurity=False)
    opener = gzip.open if _is_gzip(path) else open
    with opener(path, "rb") as raw:
        for chunk in iter(lambda: raw.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _is_gzip(path: Path) -> bool:
    with open(path, "rb") as probe:
        return probe.read(2) == b"\x1f\x8b"


class Feature:
    __slots__ = ("seqid", "ftype", "strands", "segs", "parents", "digest", "disagree",
                 "flags", "line", "lo", "hi")

    def __init__(self, seqid, ftype, parents, digest, line):
        self.seqid = seqid
        self.ftype = ftype
        #: Every strand its rows give. One, normally; a trans-spliced feature may
        #: give two, and the comparison needs to see whether a loader kept both.
        self.strands: set[str] = set()
        self.segs: list[tuple[int, int]] = []
        self.parents = parents
        self.digest = digest
        self.disagree = False
        self.flags: set[str] = set()
        self.line = line
        self.lo = None
        self.hi = None


def read(path: Path, ignore: frozenset[str] = frozenset()):
    """One pass: features, sequences and per-row facts.

    ``ignore`` names attributes to leave out, for an attribute a loader adds on
    export that the input never had (cgload writes ``cgload_linkage_method``)."""
    features: dict[str, Feature] = {}
    anonymous_seen: dict[str, int] = {}
    sequences: dict[str, dict] = {}
    rows = 0
    id_conflicts = 0
    non_numeric = 0
    with gffmodel.open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4 and parts[3].isdigit():
                    sequences.setdefault(parts[1], {"length": None, "circular": 0, "genome": ""})
                    sequences[parts[1]]["length"] = int(parts[3])
                continue
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                continue
            rows += 1
            seqid, _source, ftype, start, end, _score, strand, _phase, column9 = columns
            try:
                start_n, end_n = int(start), int(end)
            except ValueError:
                non_numeric += 1
                continue
            attributes = gffmodel.parse_attributes(column9)
            for name in ignore:
                attributes.pop(name, None)
            if ftype == "region":
                entry = sequences.setdefault(seqid, {"length": None, "circular": 0, "genome": ""})
                if attributes.get("Is_circular", [""])[0].lower() == "true":
                    entry["circular"] = 1
                genome = attributes.get("genome", [""])[0]
                if genome:
                    entry["genome"] = genome
            parents = ",".join(sorted(set(attributes.get("Parent", []))))
            digest = gffmodel.attrs_digest(attributes)
            identifier = attributes.get("ID", [None])[0]
            if identifier is None:
                key = gffmodel.unique_anon(
                    gffmodel.anon_key(seqid, ftype, start_n, end_n, strand, parents),
                    anonymous_seen,
                )
            else:
                key = identifier
            feature = features.get(key)
            if feature is None:
                feature = Feature(seqid, ftype, parents, digest, number)
                features[key] = feature
                if identifier is None:
                    feature.flags.add("no_id")
            else:
                if feature.seqid != seqid or feature.ftype != ftype:
                    id_conflicts += 1
                if digest != feature.digest:
                    feature.disagree = True
            feature.strands.add(strand)
            feature.segs.append((start_n, end_n))
            low, high = min(start_n, end_n), max(start_n, end_n)
            feature.lo = low if feature.lo is None else min(feature.lo, low)
            feature.hi = high if feature.hi is None else max(feature.hi, high)
            if start_n > end_n:
                feature.flags.add("inverted")
            if strand == "?":
                feature.flags.add("strand_unknown")
            exceptions = " ".join(attributes.get("exception", [])).lower()
            if any(word in exceptions for word in TRANS_SPLICING):
                feature.flags.add("trans_spliced")
            if any(len(values) > 1 for key2, values in attributes.items()
                   if key2 not in ("ID", "Parent")):
                feature.flags.add("multi_valued")
    return features, sequences, rows, id_conflicts, non_numeric


def digest_disagreeing(path: Path, features: dict[str, Feature],
                       ignore: frozenset[str] = frozenset()) -> dict[str, str]:
    """Second pass, only for features whose rows carry different attributes:
    their digest is taken row by row (``gffmodel.feature_digest``), so a loader
    that keeps every value but gives it to the wrong row is caught."""
    wanted = {key for key, feature in features.items() if feature.disagree}
    if not wanted:
        return {}
    rows: dict[str, list] = {key: [] for key in wanted}
    for _number, columns in gffmodel.iter_rows(path):
        attributes = gffmodel.parse_attributes(columns[8])
        for name in ignore:
            attributes.pop(name, None)
        identifier = attributes.get("ID", [None])[0]
        if identifier in rows:
            rows[identifier].append((columns[3], columns[4], attributes))
    return {key: gffmodel.feature_digest(values) for key, values in rows.items()}


def flag(features: dict[str, Feature], sequences: dict[str, dict]) -> None:
    for feature in features.values():
        sequence = sequences.get(feature.seqid, {})
        length = sequence.get("length")
        circular = sequence.get("circular", 0)
        if len(feature.segs) > 1:
            feature.flags.add("multi_row")
            ordered = sorted((min(a, b), max(a, b)) for a, b in feature.segs)
            if any(ordered[i][0] <= ordered[i - 1][1] for i in range(1, len(ordered))):
                feature.flags.add("overlapping_rows")
        if length is not None:
            if circular:
                if feature.lo > length:
                    feature.flags.add("past_end")
                elif feature.hi > length:
                    feature.flags.add("crosses_origin")
            elif feature.hi > length:
                feature.flags.add("past_end")
        if sequence.get("genome", "") in gffmodel.ORGANELLE_GENOMES:
            feature.flags.add("organelle")
        if feature.disagree:
            feature.flags.add("rows_disagree")
        if not feature.parents:
            continue
        names = feature.parents.split(",")
        if len(names) > 1:
            feature.flags.add("multi_parent")
        for name in names:
            parent = features.get(name)
            if parent is None:
                feature.flags.add("dangling_parent")
                continue
            if feature.lo < parent.lo or feature.hi > parent.hi:
                feature.flags.add("escapes_parent")
            if feature.ftype == "exon" and parent.ftype == "pseudogene":
                feature.flags.add("exon_under_pseudogene")


def build(path: Path, out: Path, ignore: frozenset[str] = frozenset()) -> dict:
    features, sequences, rows, id_conflicts, non_numeric = read(path, ignore)
    merged = digest_disagreeing(path, features, ignore)
    flag(features, sequences)

    connection = gffmodel.create(out)

    def records():
        for key, f in features.items():
            yield (
                key, f.seqid, f.ftype, ",".join(sorted(f.strands)), f.lo, f.hi, len(f.segs),
                gffmodel.segs_text(f.segs), f.parents, merged.get(key, f.digest), 1,
                ",".join(sorted(f.flags)),
            )

    gffmodel.write_features(connection, records())
    connection.executemany(
        "INSERT INTO sequence VALUES (?, ?, ?, ?)",
        [(name, s["length"], s["circular"], s["genome"]) for name, s in sequences.items()],
    )
    types = Counter(f.ftype for f in features.values())
    with_parent = Counter(f.ftype for f in features.values() if f.parents)
    hard = Counter(name for f in features.values() for name in f.flags)
    summary = {
        "input": str(path),
        "md5": md5_of_content(path),
        "rows": rows,
        "features": len(features),
        "parent_links": sum(len(f.parents.split(",")) for f in features.values() if f.parents),
        "features_with_parent": sum(with_parent.values()),
        "sequences": len(sequences),
        "circular_sequences": sum(1 for s in sequences.values() if s["circular"]),
        "organelle_sequences": sum(
            1 for s in sequences.values() if s["genome"] in gffmodel.ORGANELLE_GENOMES
        ),
        "id_conflicts": id_conflicts,
        "non_numeric_rows": non_numeric,
        "types": dict(sorted(types.items())),
        "with_parent_by_type": dict(sorted(with_parent.items())),
        "hard_cases": {name: hard.get(name, 0) for name in gffmodel.HARD_CASES},
    }
    gffmodel.set_meta(connection, **{k: json.dumps(v) for k, v in summary.items()})
    connection.commit()
    connection.close()
    return summary


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("gff", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--summary", type=Path, help="also write the summary here as JSON")
    parser.add_argument("--ignore-attribute", action="append", default=[], metavar="KEY",
                        help="leave this attribute out (repeatable); for attributes a loader "
                             "adds when it writes the file back out")
    args = parser.parse_args(argv)
    summary = build(args.gff, args.out, frozenset(args.ignore_attribute))
    text = json.dumps(summary, indent=1)
    if args.summary:
        args.summary.write_text(text + "\n")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
