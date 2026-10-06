#!/usr/bin/env python3
"""Per-genome structural-defect tally for RefSeq GFF3, one row per genome.

Extends cgload's preflight_gff3.py in four ways that the survey needs and the
single-file version does not have:

  1. Per-genome tallies emitted as a machine-readable row, not a printout.
  2. Molecule class per defect (nuclear / mitochondrion / chloroplast /
     plasmid / unplaced), read from the `genome=` attribute on RefSeq `region`
     rows.  No assembly_report needed.
  3. Denominators and applicability.  A check that could not run is reported
     as NA, never as zero.  This is the difference between a rate and a
     collection of hits.
  4. Attribution.  Rows carrying `exception=`, `part=`, `partial=` or
     `pseudo=true` are counted separately, because RefSeq documents several of
     these as intentional conventions.  A defect that a reviewer can explain
     away must be separable from one they cannot.

Usage:
    python cgsurvey.py --gff genomic.gff.gz --accession GCF_000001405.40 \
        --species "Homo sapiens" --clade vertebrate --out-dir results/

Emits results/<accession>.json (full detail) and appends a row to
results/survey.tsv.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import re
import sys
from collections import defaultdict

ATTR = {
    "ID": re.compile(r"(?:^|;)ID=([^;]*)"),
    "Parent": re.compile(r"(?:^|;)Parent=([^;]*)"),
    "Is_circular": re.compile(r"(?:^|;)Is_circular=([^;]*)"),
    "genome": re.compile(r"(?:^|;)genome=([^;]*)"),
    "exception": re.compile(r"(?:^|;)exception=([^;]*)"),
    "part": re.compile(r"(?:^|;)part=([^;]*)"),
    "partial": re.compile(r"(?:^|;)partial=([^;]*)"),
    "pseudo": re.compile(r"(?:^|;)pseudo=([^;]*)"),
    "gene_biotype": re.compile(r"(?:^|;)gene_biotype=([^;]*)"),
}

GENE_LEVEL = {"gene", "pseudogene"}
ORGANELLE = {"mitochondrion", "chloroplast", "plastid", "apicoplast",
             "kinetoplast", "plasmid"}


def get(attrs: str, key: str) -> str | None:
    m = ATTR[key].search(attrs)
    return m.group(1) if m else None


def compact_flags(attrs: str) -> tuple:
    """(excuse tuple, is_pseudo). Keeps memory flat on multi-million-row genomes."""
    excuses = []
    for key in ("exception", "part", "partial"):
        v = get(attrs, key)
        if v:
            excuses.append(f"{key}={v}")
    pseudo = ((get(attrs, "pseudo") or "").lower() == "true"
              or (get(attrs, "gene_biotype") or "") == "pseudogene")
    return (tuple(excuses), pseudo)


def fasta_lengths(path: str) -> dict[str, int]:
    """Sequence lengths from <fasta>.fai if present, else by scanning the FASTA."""
    fai = path + ".fai"
    lengths: dict[str, int] = {}
    if os.path.exists(fai):
        with open(fai) as fh:
            for line in fh:
                p = line.split("\t")
                lengths[p[0]] = int(p[1])
        return lengths
    name, n = None, 0
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def circular_seqids(refseq_gff: str) -> set[str]:
    """Seqids whose RefSeq region row carries Is_circular=true."""
    out: set[str] = set()
    with opener(refseq_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\r\n").split("\t", 8)
            if len(c) == 9 and c[2] == "region" and \
                    (get(c[8], "Is_circular") or "").lower() == "true":
                out.add(c[0])
    return out


def molecule_map(refseq_gff: str) -> dict[str, str]:
    """seqid -> molecule class, borrowed from a RefSeq GFF of the same assembly."""
    out: dict[str, str] = {}
    with opener(refseq_gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\r\n").split("\t", 8)
            if len(c) == 9 and c[2] == "region":
                out[c[0]] = (get(c[8], "genome") or "nuclear").lower()
    return out


def opener(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, errors="replace")


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def survey(path: str, lengths: dict[str, int] | None = None,
           molecules: dict[str, str] | None = None,
           circular: set[str] | None = None) -> dict:
    # FASTA lengths take precedence over ##sequence-region, so a pipeline GFF
    # with no headers is checked exactly as hard as a RefSeq one.
    regions: dict[str, int] = dict(lengths or {})
    molecule: dict[str, str] = dict(molecules or {})
    lengths_source = "fasta" if lengths else "sequence-region"
    circ: set[str] = set(circular or ())
    header_provenance: list[str] = []
    col2_sources: dict[str, int] = defaultdict(int)
    features: dict[str, tuple] = {}       # ID -> (line, seqid, type, s, e, strand, attrs)
    children: dict[str, list] = defaultdict(list)
    duplicate_ids = 0

    n_rows = 0
    n_with_id = 0
    n_with_parent = 0
    n_edges = 0
    types = defaultdict(int)
    seqids_seen: set[str] = set()

    problems: list[dict] = []

    def note(kind, line, seqid, detail, **extra):
        rec = {"kind": kind, "line": line, "seqid": seqid,
               "molecule": molecule.get(seqid, "unknown"), "detail": detail}
        rec.update(extra)
        problems.append(rec)

    with opener(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4 and parts[1] not in (lengths or {}):
                    try:
                        regions[parts[1]] = int(parts[3])
                    except ValueError:
                        pass
                continue
            if line.startswith("#!annotation-source") or line.startswith("#!annotation-date") \
                    or line.startswith("#!genome-build"):
                header_provenance.append(line[2:].strip())
                continue
            if line.startswith("#"):
                continue
            columns = line.rstrip("\r\n").split("\t")
            if len(columns) < 9:
                continue
            seqid, src_col, ftype, raw_start, raw_end, _score, strand = columns[:7]
            attrs = columns[8]
            n_rows += 1
            types[ftype] += 1
            col2_sources[src_col] += 1
            seqids_seen.add(seqid)

            # RefSeq emits one `region` row per sequence carrying genome=...
            if ftype == "region":
                if seqid not in (molecules or {}):
                    g = get(attrs, "genome")
                    molecule[seqid] = (g or "nuclear").lower()
                if (get(attrs, "Is_circular") or "").lower() == "true":
                    circ.add(seqid)

            try:
                start, end = int(raw_start), int(raw_end)
            except ValueError:
                note("non-numeric", number, seqid, f"{raw_start}..{raw_end}")
                continue

            if start > end:
                note("inverted", number, seqid,
                     f"{ftype} start {start:,} exceeds end {end:,}", ftype=ftype)

            if strand == "?":
                note("strand-?", number, seqid, f"{ftype} strand '?'", ftype=ftype)
            elif strand not in ("+", "-", "."):
                note("strand-invalid", number, seqid,
                     f"{strand!r} is not + - . or ?", ftype=ftype)

            length = regions.get(seqid)
            if length is not None and (start > length or end > length):
                # GFF3 spec and NCBI both extend `end` into virtual space for a
                # row that crosses the origin of a circular sequence. That is a
                # convention, and counting it as a defect would inflate exactly
                # the organellar category the survey reports on.
                if seqid in circ and start <= length and end <= 2 * length:
                    note("origin-spanning", number, seqid,
                         f"{ftype} {start:,}..{end:,} crosses origin of circular "
                         f"{seqid} (length {length:,})", ftype=ftype, source=src_col)
                else:
                    note("out-of-range", number, seqid,
                         f"{ftype} {start:,}..{end:,} exceeds {seqid} length {length:,}",
                         ftype=ftype, source=src_col)

            flags = compact_flags(attrs)
            fid = get(attrs, "ID")
            if fid:
                n_with_id += 1
                if fid in features:
                    duplicate_ids += 1
                features[fid] = (number, seqid, ftype, start, end, strand, flags)

            parent = get(attrs, "Parent")
            if parent:
                n_with_parent += 1
                for one in parent.split(","):
                    one = one.strip()
                    n_edges += 1
                    children[one].append((number, ftype, start, end, seqid, strand, flags,
                                          sys.intern(src_col)))

    # ---- escapes-parent, with attribution ----------------------------------
    dangling = 0
    for parent_id, kids in children.items():
        if parent_id not in features:
            dangling += len(kids)
            continue
        pline, pseq, ptype, pstart, pend, pstrand, pattrs = features[parent_id]
        for kline, ktype, kstart, kend, kseq, kstrand, kattrs, ksrc in kids:
            cross_seq = kseq != pseq
            outside = kstart < pstart or kend > pend
            if not (outside or cross_seq):
                continue
            excused = [f"child:{e}" for e in kattrs[0]] + [f"parent:{e}" for e in pattrs[0]]
            note("escapes-parent" if not cross_seq else "escapes-parent-seqid",
                 kline, kseq,
                 f"{ktype} {kstart:,}..{kend:,} outside {ptype} {parent_id} "
                 f"{pstart:,}..{pend:,} (line {pline})",
                 ftype=ktype, parent_type=ptype, parent_line=pline, source=ksrc,
                 excused_by=excused or None)

    # ---- gene-level exons on pseudogenes (RefSeq CONVENTION, not a defect) --
    pseudo_gene_exons = 0
    gene_level_exons = 0
    for parent_id, kids in children.items():
        if parent_id not in features:
            continue
        _pl, pseq, ptype, _ps, _pe, _pst, pattrs = features[parent_id]
        if ptype not in GENE_LEVEL:
            continue
        for kline, ktype, _ks, _ke, kseq, _kst, _ka, _ksrc in kids:
            if ktype != "exon":
                continue
            gene_level_exons += 1
            is_pseudo = ptype == "pseudogene" or pattrs[1]
            if is_pseudo:
                pseudo_gene_exons += 1
                note("pseudogene-gene-exon", kline, kseq,
                     f"exon parented directly by {ptype} {parent_id}",
                     ftype="exon", parent_type=ptype)

    # ---- roll up -----------------------------------------------------------
    by_kind = defaultdict(int)
    by_kind_organellar = defaultdict(int)
    by_kind_excused = defaultdict(int)
    for p in problems:
        by_kind[p["kind"]] += 1
        if p.get("molecule") in ORGANELLE:
            by_kind_organellar[p["kind"]] += 1
        if p.get("excused_by"):
            by_kind_excused[p["kind"]] += 1

    seq_region_coverage = (
        len(seqids_seen & set(regions)) / len(seqids_seen) if seqids_seen else 0.0
    )

    by_kind_source: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for p in problems:
        if p.get("source"):
            by_kind_source[p["kind"]][p["source"]] += 1

    return {
        "provenance": {
            "header": header_provenance,
            "column2_sources": dict(sorted(col2_sources.items(), key=lambda x: -x[1])),
            "circular_sequences": sorted(circ),
        },
        "counts_by_source": {k: dict(v) for k, v in by_kind_source.items()},
        "file_md5": md5(path),
        "file_bytes": os.path.getsize(path),
        "denominators": {
            "rows": n_rows,
            "features_with_ID": n_with_id,
            "rows_with_Parent": n_with_parent,
            "parent_child_edges": n_edges,
            "resolvable_edges": n_edges - dangling,
            "dangling_edges": dangling,
            "duplicate_IDs": duplicate_ids,
            "sequences": len(seqids_seen),
            "sequences_with_sequence_region": len(seqids_seen & set(regions)),
            "organellar_sequences": sum(1 for v in molecule.values() if v in ORGANELLE),
            "feature_types": len(types),
        },
        "applicability": {
            # out-of-range can only run where a ##sequence-region exists
            "out-of-range": "full" if seq_region_coverage == 1.0
            else ("partial" if seq_region_coverage > 0 else "NA"),
            "out-of-range_coverage": round(seq_region_coverage, 4),
            "lengths_source": lengths_source,
            "molecule_class_coverage": (round(
                len(seqids_seen & set(molecule)) / len(seqids_seen), 4)
                if seqids_seen else 0.0),
            "escapes-parent": "full" if n_edges else "NA",
        },
        "counts": dict(by_kind),
        "counts_organellar": dict(by_kind_organellar),
        "counts_with_documented_exception": dict(by_kind_excused),
        "convention_gene_level_exons": {
            "exons_parented_by_gene_level": gene_level_exons,
            "of_which_pseudogene": pseudo_gene_exons,
        },
        "molecule_classes": dict(sorted(
            (k, sum(1 for v in molecule.values() if v == k))
            for k in set(molecule.values()))),
        "problems": problems,
    }


TSV_COLUMNS = [
    "arm", "dialect", "clade", "annotation_source", "top_col2_source", "species", "accession", "assembly_name", "annotation_release",
    "downloaded", "file_md5", "rows", "parent_child_edges", "sequences",
    "organellar_sequences",
    "escapes_parent", "escapes_parent_organellar", "escapes_parent_excused",
    "escapes_parent_seqid",
    "inverted", "out_of_range", "origin_spanning_convention", "out_of_range_applicability",
    "molecule_class_coverage", "strand_q", "strand_invalid", "non_numeric",
    "pseudogene_gene_exons", "gene_level_exons",
    "dangling_edges", "duplicate_ids", "any_defect",
]


def row_for(meta: dict, r: dict) -> list[str]:
    c, o, x = r["counts"], r["counts_organellar"], r["counts_with_documented_exception"]
    d = r["denominators"]
    hard = sum(c.get(k, 0) for k in
               ("escapes-parent", "escapes-parent-seqid", "inverted",
                "out-of-range", "strand-invalid", "non-numeric"))
    return [
        meta.get("arm", ""), meta.get("dialect", ""), meta.get("clade", ""),
        next((h for h in r["provenance"]["header"] if h.startswith("annotation-source")),
             "").replace("annotation-source", "").strip() or "NA",
        next(iter(r["provenance"]["column2_sources"]), "NA"), meta.get("species", ""), meta.get("accession", ""),
        meta.get("assembly_name", ""), meta.get("annotation_release", ""),
        meta.get("downloaded", ""), r["file_md5"],
        d["rows"], d["parent_child_edges"], d["sequences"], d["organellar_sequences"],
        c.get("escapes-parent", 0), o.get("escapes-parent", 0), x.get("escapes-parent", 0),
        c.get("escapes-parent-seqid", 0),
        c.get("inverted", 0), c.get("out-of-range", 0), c.get("origin-spanning", 0),
        r["applicability"]["out-of-range"],
        r["applicability"]["molecule_class_coverage"],
        c.get("strand-?", 0), c.get("strand-invalid", 0), c.get("non-numeric", 0),
        r["convention_gene_level_exons"]["of_which_pseudogene"],
        r["convention_gene_level_exons"]["exons_parented_by_gene_level"],
        d["dangling_edges"], d["duplicate_IDs"],
        "yes" if hard else "no",
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff")
    ap.add_argument("--fasta", help="assembly FASTA; supplies lengths when the GFF lacks ##sequence-region")
    ap.add_argument("--molecules-from", help="RefSeq GFF of the same assembly; supplies organellar labels")
    ap.add_argument("--arm", default="refseq", choices=["refseq", "pipeline"])
    ap.add_argument("--dialect", default="refseq", help="refseq / funannotate / braker")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--accession", default="")
    ap.add_argument("--species", default="")
    ap.add_argument("--clade", default="")
    ap.add_argument("--assembly-name", default="")
    ap.add_argument("--annotation-release", default="")
    ap.add_argument("--downloaded", default="")
    ap.add_argument("--out-dir", default="results")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if not args.gff:
        ap.error("--gff is required")
    os.makedirs(args.out_dir, exist_ok=True)
    lengths = fasta_lengths(args.fasta) if args.fasta else None
    molecules = molecule_map(args.molecules_from) if args.molecules_from else None
    circular = circular_seqids(args.molecules_from) if args.molecules_from else None
    result = survey(args.gff, lengths, molecules, circular)
    if result["applicability"]["out-of-range"] != "full":
        print(f"WARNING: out-of-range coverage "
              f"{result['applicability']['out-of-range_coverage']:.0%}; pass --fasta",
              file=sys.stderr)
    meta = {"arm": args.arm, "dialect": args.dialect,
        "clade": args.clade, "species": args.species, "accession": args.accession,
        "assembly_name": args.assembly_name,
        "annotation_release": args.annotation_release, "downloaded": args.downloaded,
    }
    stem = (args.accession or os.path.basename(args.gff).split(".")[0]) + f".{args.arm}.{args.dialect}"
    with open(os.path.join(args.out_dir, f"{stem}.json"), "w") as fh:
        json.dump({"meta": meta, **result}, fh, indent=1)

    tsv = os.path.join(args.out_dir, "survey.tsv")
    new = not os.path.exists(tsv)
    with open(tsv, "a") as fh:
        if new:
            fh.write("\t".join(TSV_COLUMNS) + "\n")
        fh.write("\t".join(str(v) for v in row_for(meta, result)) + "\n")

    print(f"{args.species or args.gff}: "
          f"{sum(result['counts'].values()):,} flagged rows "
          f"across {len(result['counts'])} kind(s); "
          f"{result['denominators']['rows']:,} rows, "
          f"{result['denominators']['parent_child_edges']:,} edges")
    for k, v in sorted(result["counts"].items(), key=lambda x: -x[1]):
        org = result["counts_organellar"].get(k, 0)
        exc = result["counts_with_documented_exception"].get(k, 0)
        print(f"  {k:26s} {v:>8,}   organellar {org:>6,}   documented-exception {exc:>6,}")
    return 0


SELFTEST_GFF = (
    "##gff-version 3\n##sequence-region NC_1.1 1 1000\n##sequence-region NC_M.1 1 500\n"
    "NC_1.1\tRefSeq\tregion\t1\t1000\t.\t+\t.\tID=r1;genome=chromosome\n"
    "NC_M.1\tRefSeq\tregion\t1\t500\t.\t+\t.\tID=r2;genome=mitochondrion;Is_circular=true\n"
    "NC_1.1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=g0\n"
    "NC_1.1\tRefSeq\tmRNA\t100\t200\t.\t+\t.\tID=t0;Parent=g0\n"
    "NC_1.1\tRefSeq\texon\t160\t250\t.\t+\t.\tID=e1;Parent=t0\n"
    "NC_1.1\tRefSeq\tpseudogene\t300\t400\t.\t+\t.\tID=g1;pseudo=true\n"
    "NC_1.1\tRefSeq\texon\t300\t350\t.\t+\t.\tID=e2;Parent=g1\n"
    "NC_1.1\tRefSeq\tgene\t500\t450\t.\t?\t.\tID=g2\n"
    "NC_1.1\tRefSeq\tgene\t900\t1200\t.\t+\t.\tID=g3\n"
    "NC_M.1\tRefSeq\tgene\t10\t490\t.\t+\t.\tID=mg;exception=trans-splicing\n"
    "NC_M.1\tRefSeq\tmRNA\t10\t495\t.\t+\t.\tID=mt;Parent=mg\n"
    "NC_M.1\tRefSeq\tgene\t480\t520\t.\t+\t.\tID=wrap\n"
)
EXPECTED = {"escapes-parent": 2, "inverted": 1, "strand-?": 1,
            "out-of-range": 1, "origin-spanning": 1, "pseudogene-gene-exon": 1}


def selftest() -> int:
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".gff", delete=False) as fh:
        fh.write(SELFTEST_GFF)
    r = survey(fh.name)
    os.unlink(fh.name)
    ok = (r["counts"] == EXPECTED
          and r["counts_organellar"].get("escapes-parent") == 1
          and r["counts_with_documented_exception"].get("escapes-parent") == 1)
    print("SELFTEST PASS" if ok else f"SELFTEST FAIL: {r['counts']}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
