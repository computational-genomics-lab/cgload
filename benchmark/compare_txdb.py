#!/usr/bin/env python3
"""Compare a TxDb (txdbmaker) with the reference count of the same file.

    python3 benchmark/compare_txdb.py TRUTH.sqlite TXDB_DIR --gff GFF \\
        [--log LOG ...] --genome ACC [--json OUT] [--cases-tsv TSV]

TXDB_DIR is what run_txdb.R wrote. A TxDb models transcripts, their exons and
CDS, and which gene each transcript belongs to; it has no place for other
feature types or for attributes. So the comparison is made per transcript, on
the fields a TxDb holds:

* the transcript itself: sequence, strand, start and end;
* its gene: the gene the TxDb links it to, against the gene the file names as
  its Parent;
* its exons: the exon rows under it in the file, against the TxDb's exons
  (only when the file gives it exon rows: a TxDb makes exons for a transcript
  that has none, from its CDS or its span, and those are not compared);
* its CDS: the CDS rows under it in the file, against the TxDb's CDS pieces.

Which features become transcripts follows txdbmaker's own rules, read from
the type lists run_txdb.R took from the package (types.txt): a feature whose
type is a transcript type, and a gene whose exons or CDS hang directly from it.

A TxDb does not keep the file's IDs. It names a transcript by its
``transcript_id`` attribute when the file uses that attribute, otherwise by its
``Name``, and a gene by its ``gene_id``, ``Name`` or ID, by the rules of
makeTxDbFromGRanges, which ``txdb_names`` repeats. Each transcript of the file
is paired with one of the TxDb's in stages: same name and span; same name,
when only one of each has it; same span and type; same span, when only one of
each has it. A transcript left unpaired is lost; a TxDb transcript left
unpaired is extra. Exon and CDS features under no transcript, which txdbmaker
drops, are counted apart.

A lost or changed transcript is announced when txdbmaker's own output (--log)
names it or a feature below it; txdbmaker lists the transcripts it drops by
ID. A transcript belongs to a hard case when it, or a feature below it, has
that flag in the reference count. Genes are counted separately: a gene is held
when the TxDb has its gene ID. Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import gffmodel  # noqa: E402

#: GRanges has no '?': rtracklayer reads it as '*', which means "no strand".
TXDB_STRAND = {"*": ".", "+": "+", "-": "-"}
KINDS = ("gene", "transcript", "exon", "cds", "stop_codon", "imported")


def read_types(path: Path) -> tuple[dict[str, set[str]], set[str]]:
    kinds: dict[str, set[str]] = {kind: set() for kind in KINDS}
    sources = set()
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        kind, source, value = line.split("\t")
        kinds[kind].add(value)
        sources.add(source)
    return kinds, sources


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _drop_prefix(ids: list[str | None], level: str) -> list[str | None]:
    """makeTxDbFromGRanges's .drop_ID_prefix: strip 'level:' or 'level-' when
    every ID that has a value starts with it."""
    for separator in (":", "-"):
        prefix = level + separator
        if all(i.startswith(prefix) for i in ids if i is not None):
            return [i[len(prefix):] if i is not None else None for i in ids]
    return ids


def _pattern(values: list) -> list[int]:
    """R's match(x, x): where each value first occurs."""
    first: dict = {}
    return [first.setdefault(value, index) for index, value in enumerate(values)]


def txdb_names(gff: Path, kinds: dict[str, set[str]]) -> dict:
    """The names a TxDb gives the file's transcripts and genes.

    makeTxDbFromGFF imports only the rows whose type it models, and then:
    a feature's Name is its ``Name``, or, for a type none of whose rows has a
    Name, its ID with a 'type:' or 'type-' prefix dropped; a transcript's name
    is its ``transcript_id`` if any imported row has that attribute, otherwise
    its Name; a gene's ID in the TxDb is its ``gene_id``, else its Name (each
    only when every gene has one and they tell genes apart as the IDs do),
    else its ID with a 'gene:' or 'gene-' prefix dropped.
    """
    rows = []
    has_transcript_id = has_gene_id = False
    for _number, columns in gffmodel.iter_rows(gff):
        ftype = columns[2]
        if ftype not in kinds["imported"]:
            continue
        attributes = gffmodel.parse_attributes(columns[8])
        has_transcript_id |= "transcript_id" in attributes
        has_gene_id |= "gene_id" in attributes

        def value(key, attributes=attributes):
            found = ",".join(attributes.get(key, []))
            return found or None

        rows.append((ftype, value("ID"), value("Name"), value("transcript_id"),
                     value("gene_id")))

    names = [name for _t, _i, name, _x, _g in rows]
    by_type: dict[str, list[int]] = defaultdict(list)
    for index, (ftype, *_rest) in enumerate(rows):
        by_type[ftype].append(index)
    for ftype, indices in by_type.items():
        if all(names[i] is None for i in indices):
            for i, name in zip(indices, _drop_prefix([rows[i][1] for i in indices], ftype),
                               strict=True):
                names[i] = name

    tx_name: dict[str, str | None] = {}
    for index, (ftype, row_id, _name, transcript_id, _g) in enumerate(rows):
        if row_id is None or row_id in tx_name:
            continue
        if ftype in kinds["transcript"] or ftype in kinds["gene"]:
            tx_name[row_id] = transcript_id if has_transcript_id else names[index]

    gene_rows = [index for index, row in enumerate(rows) if row[0] in kinds["gene"]]
    gene_ids = [rows[i][1] for i in gene_rows]
    candidate = None
    rule = "ID"
    if has_gene_id and all(rows[i][4] is not None for i in gene_rows):
        candidate, rule = [rows[i][4] for i in gene_rows], "gene_id"
    elif all(names[i] is not None for i in gene_rows):
        candidate, rule = [names[i] for i in gene_rows], "Name"
    if candidate is None or _pattern(gene_ids) != _pattern(candidate):
        candidate, rule = _drop_prefix(gene_ids, "gene"), "ID"
    gene_name: dict[str, str] = {}
    for gene_id, name in zip(gene_ids, candidate, strict=True):
        if gene_id is not None and name is not None:
            gene_name.setdefault(gene_id, name)
    return {"transcript": tx_name, "gene": gene_name,
            "transcript_rule": "transcript_id" if has_transcript_id else "Name",
            "gene_rule": rule}


def _log_text(paths) -> str:
    return "\n".join(Path(p).read_text(errors="replace") for p in paths if Path(p).exists())


def compare(truth: Path, txdb: Path, gff: Path, logs=(), examples: int = 5) -> dict:
    kinds, sources = read_types(txdb / "types.txt")
    names = txdb_names(gff, kinds)
    log_text = _log_text(logs)

    connection = sqlite3.connect(truth)
    features = {}
    children: dict[str, list[str]] = defaultdict(list)
    for key, seqid, ftype, strand, lo, hi, segs, parents, flags in connection.execute(
        'SELECT key, seqid, type, strand, start, "end", segs, parents, flags FROM feature'
    ):
        owners = parents.split(",") if parents else []
        features[key] = {"seqid": seqid, "type": ftype, "strand": strand, "start": lo,
                         "end": hi, "segs": gffmodel.parse_segs(segs), "parents": owners,
                         "flags": {name for name in (flags or "").split(",") if name}}
        for owner in owners:
            children[owner].append(key)
    connection.close()

    def kids(key: str, wanted: set[str]) -> list[str]:
        return [child for child in children.get(key, ()) if features[child]["type"] in wanted]

    genes = {key for key, f in features.items() if f["type"] in kinds["gene"]}
    as_transcript = {key for key in genes if kids(key, kinds["exon"] | kinds["cds"])}
    transcripts = sorted(key for key, f in features.items()
                         if f["type"] in kinds["transcript"] or key in as_transcript)
    modelled = kinds["gene"] | kinds["transcript"] | kinds["exon"] | kinds["cds"]
    not_modelled = Counter(f["type"] for f in features.values() if f["type"] not in modelled)
    in_transcripts = set(transcripts)
    orphan_parts = Counter(
        f["type"] for f in features.values()
        if f["type"] in kinds["exon"] | kinds["cds"]
        and not any(p in in_transcripts for p in f["parents"]))

    expected = {}
    for key in transcripts:
        feature = features[key]
        below = kids(key, kinds["exon"] | kinds["cds"])
        exons = sorted(seg for child in below if features[child]["type"] in kinds["exon"]
                       for seg in features[child]["segs"])
        cds = sorted(seg for child in below if features[child]["type"] in kinds["cds"]
                     for seg in features[child]["segs"])
        owners = [key] if key in as_transcript else [p for p in feature["parents"] if p in genes]
        expected[key] = {
            **feature, "name": names["transcript"].get(key), "exons": exons, "cds": cds,
            "genes": {names["gene"][p] for p in owners if p in names["gene"]},
            "below": below,
            "cases": set(feature["flags"]).union(*(features[c]["flags"] for c in below)),
        }

    held = read_tsv(txdb / "transcripts.tsv")
    parts = {}
    for name in ("exons", "cds"):
        grouped = defaultdict(list)
        for row in read_tsv(txdb / f"{name}.tsv"):
            grouped[row["tx_id"]].append((int(row["start"]), int(row["end"])))
        parts[name] = {tx: sorted(ranges) for tx, ranges in grouped.items()}
    for row in held:
        row["start"], row["end"] = int(row["start"]), int(row["end"])

    pairs, stages = _pair(expected, held)
    totals = Counter()
    by_case: dict[str, Counter] = defaultdict(Counter)
    field_diff = Counter()
    shown: dict[str, list[str]] = defaultdict(list)
    for key, want in expected.items():
        totals["truth"] += 1
        for case in want["cases"]:
            by_case[case]["truth"] += 1
        row = pairs.get(key)
        diffs = []
        if row is None:
            state = "lost"
        else:
            if row["seqid"] != want["seqid"]:
                diffs.append("seqid")
            if TXDB_STRAND.get(row["strand"], row["strand"]) != want["strand"]:
                diffs.append("strand")
            if (row["start"], row["end"]) != (want["start"], want["end"]):
                diffs.append("span")
            if {g for g in row["genes"].split(",") if g} != want["genes"]:
                diffs.append("gene")
            if want["exons"] and parts["exons"].get(row["tx_id"], []) != want["exons"]:
                diffs.append("exons")
            if parts["cds"].get(row["tx_id"], []) != want["cds"]:
                diffs.append("cds")
            state = "changed" if diffs else "kept"
        totals[state] += 1
        loud = state != "kept" and any(name in log_text for name in (key, *want["below"]))
        if state != "kept":
            totals[f"{state}_{'announced' if loud else 'silent'}"] += 1
        for name in diffs:
            field_diff[name] += 1
        for case in want["cases"]:
            by_case[case][state] += 1
            if state != "kept":
                by_case[case][f"{state}_{'announced' if loud else 'silent'}"] += 1
        label = state if not diffs else f"{state}:{'+'.join(diffs)}"
        if state != "kept" and len(shown[label]) < examples:
            shown[label].append(key)

    held_genes = {row["gene_id"] for row in read_tsv(txdb / "genes.tsv")}
    gene_totals = Counter()
    for key in genes:
        gene_totals["truth"] += 1
        gene_totals["held" if names["gene"].get(key) in held_genes else "absent"] += 1
    paired_rows = {id(row) for row in pairs.values()}
    return {
        "tool": "txdbmaker",
        "types_from": sorted(sources),
        "transcript_names_from": names["transcript_rule"],
        "gene_ids_from": names["gene_rule"],
        "totals": dict(totals),
        "extra_transcripts": sum(1 for row in held if id(row) not in paired_rows),
        "paired_by": dict(stages),
        "field_differences": dict(field_diff),
        "by_case": {name: dict(by_case[name]) for name in gffmodel.HARD_CASES if name in by_case},
        "genes": dict(gene_totals),
        "parts_under_no_transcript": dict(orphan_parts),
        "not_modelled_by_type": dict(not_modelled),
        "examples": dict(shown),
    }


def _pair(expected: dict, held: list[dict]) -> tuple[dict[str, dict], Counter]:
    """Pair each transcript of the file with at most one of the TxDb's."""
    pairs: dict[str, dict] = {}
    used: set[int] = set()
    stages = Counter()

    def run(stage: str, key_of_file, key_of_txdb, unique: bool) -> None:
        pool: dict = defaultdict(list)
        for index, row in enumerate(held):
            if index not in used:
                pool[key_of_txdb(row)].append(index)
        waiting: dict = defaultdict(list)
        for key, want in expected.items():
            if key not in pairs:
                waiting[key_of_file(want)].append(key)
        for signature, keys in waiting.items():
            if None in signature:
                continue
            candidates = pool.get(signature, [])
            if not candidates or (unique and (len(keys) > 1 or len(candidates) > 1)):
                continue
            for key, index in zip(sorted(keys), candidates, strict=False):
                pairs[key] = held[index]
                used.add(index)
                stages[stage] += 1

    def name(row):
        return row.get("tx_name") or None

    typed = all(row.get("type") for row in held)
    run("name and span", lambda w: (w["name"], w["seqid"], w["start"], w["end"]),
        lambda r: (name(r), r["seqid"], r["start"], r["end"]), unique=False)
    run("name", lambda w: (w["name"],), lambda r: (name(r),), unique=True)
    if typed:
        run("span and type", lambda w: (w["seqid"], w["start"], w["end"], w["type"]),
            lambda r: (r["seqid"], r["start"], r["end"], r["type"]), unique=False)
    run("span", lambda w: (w["seqid"], w["start"], w["end"]),
        lambda r: (r["seqid"], r["start"], r["end"]), unique=True)
    return pairs, stages


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("truth", type=Path)
    parser.add_argument("txdb", type=Path)
    parser.add_argument("--gff", type=Path, required=True, help="the file txdbmaker read")
    parser.add_argument("--log", type=Path, action="append", default=[],
                        help="txdbmaker's output, to tell announced from silent (repeatable)")
    parser.add_argument("--genome", default="")
    parser.add_argument("--json", type=Path)
    parser.add_argument("--cases-tsv", type=Path, help="append one row per hard case here")
    args = parser.parse_args(argv)
    result = compare(args.truth, args.txdb, args.gff, args.log)
    result["genome"] = args.genome
    if args.json:
        args.json.write_text(json.dumps(result, indent=1) + "\n")
    if args.cases_tsv:
        gffmodel.locked_append(
            args.cases_tsv,
            "genome\ttool\tcase\ttruth\tkept\tchanged\tchanged_declared\tlost"
            "\tchanged_silent\tlost_silent\tfields_changed",
            ("\t".join(map(str, (args.genome, "txdbmaker", case, counts.get("truth", 0),
                                 counts.get("kept", 0), counts.get("changed", 0), 0,
                                 counts.get("lost", 0), counts.get("changed_silent", 0),
                                 counts.get("lost_silent", 0), "")))
             for case, counts in result["by_case"].items()))
    totals = result["totals"]
    print(f"txdbmaker {args.genome}: {totals.get('truth', 0):,} transcripts in the file; "
          f"kept {totals.get('kept', 0):,}, changed {totals.get('changed', 0):,}, "
          f"lost {totals.get('lost', 0):,} ({totals.get('lost_silent', 0):,} silently); "
          f"extra {result['extra_transcripts']:,}; genes held "
          f"{result['genes'].get('held', 0):,} of {result['genes'].get('truth', 0):,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
