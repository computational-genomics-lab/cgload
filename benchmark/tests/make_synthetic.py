#!/usr/bin/env python3
"""Write a small RefSeq-like genome with every hard case the protocol names.

    python3 benchmark/tests/make_synthetic.py OUTDIR [--genes N] [--seed S]

Writes OUTDIR/genomic.gff and OUTDIR/assembly.fa. Used by the benchmark's own
tests and to estimate loader speeds; never part of any result.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

CHR1 = "NC_900001.1"
CHRM = "NC_900002.1"


def gene_model(rows, index, start, strand):
    gid, tid, cid = f"gene-G{index}", f"rna-NM_{index}.1", f"cds-NP_{index}.1"
    exons = [(start, start + 299), (start + 500, start + 799), (start + 1000, start + 1299)]
    end = exons[-1][1]
    dbx = f"Dbxref=GeneID:{index},HGNC:HGNC:{index}"
    rows.append((CHR1, "BestRefSeq", "gene", start, end, strand, f"ID={gid};{dbx};Name=G{index};gene_biotype=protein_coding"))
    rows.append((CHR1, "BestRefSeq", "mRNA", start, end, strand, f"ID={tid};Parent={gid};{dbx};product=protein {index}%2C isoform 1"))
    for n, (a, b) in enumerate(exons, start=1):
        rows.append((CHR1, "BestRefSeq", "exon", a, b, strand, f"ID=exon-NM_{index}.1-{n};Parent={tid};{dbx}"))
    for a, b in exons:
        rows.append((CHR1, "BestRefSeq", "CDS", a + 10 if a == exons[0][0] else a, b, strand,
                     f"ID={cid};Parent={tid};{dbx};protein_id=NP_{index}.1"))
    return end


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("--genes", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args(argv)
    rnd = random.Random(args.seed)
    args.out.mkdir(parents=True, exist_ok=True)
    rows = []
    position = 1000
    for index in range(1, args.genes + 1):
        position = gene_model(rows, index, position, rnd.choice("+-")) + 700
    chr1_length = position + 10_000
    # A pseudogene with an exon attached directly to it.
    p = position + 100
    rows.append((CHR1, "Curated Genomic", "pseudogene", p, p + 900, "+", "ID=gene-PSG1;gene_biotype=pseudogene;pseudo=true"))
    rows.append((CHR1, "Curated Genomic", "exon", p, p + 900, "+", "ID=id-PSG1;Parent=gene-PSG1;pseudo=true"))
    # A programmed frameshift: two CDS rows of one ID sharing a base.
    f = p + 2000
    rows.append((CHR1, "BestRefSeq", "gene", f, f + 1999, "+", "ID=gene-FS1;gene_biotype=protein_coding"))
    rows.append((CHR1, "BestRefSeq", "mRNA", f, f + 1999, "+", "ID=rna-FS1;Parent=gene-FS1"))
    rows.append((CHR1, "BestRefSeq", "exon", f, f + 1999, "+", "ID=exon-FS1-1;Parent=rna-FS1"))
    rows.append((CHR1, "BestRefSeq", "CDS", f + 100, f + 1000, "+", "ID=cds-FS1;Parent=rna-FS1;exception=ribosomal slippage"))
    rows.append((CHR1, "BestRefSeq", "CDS", f + 1000, f + 1800, "+", "ID=cds-FS1;Parent=rna-FS1;exception=ribosomal slippage"))
    # A non-coding gene of an INSDC ncRNA class.
    n = f + 3000
    rows.append((CHR1, "BestRefSeq", "gene", n, n + 400, "-", "ID=gene-YR1;gene_biotype=Y_RNA"))
    rows.append((CHR1, "BestRefSeq", "Y_RNA", n, n + 400, "-", "ID=rna-YR1;Parent=gene-YR1"))
    rows.append((CHR1, "BestRefSeq", "exon", n, n + 400, "-", "ID=exon-YR1-1;Parent=rna-YR1"))
    chr1_length = max(chr1_length, n + 5000)

    # The mitochondrion: circular, 20 kb.
    m = 20_000
    mito = [
        (CHRM, "RefSeq", "region", 1, m, "+", f"ID={CHRM}:1..{m};Is_circular=true;genome=mitochondrion;mol_type=genomic DNA"),
        # Crosses the origin: ends past the length, as NCBI writes it.
        (CHRM, "RefSeq", "gene", 19_800, 20_300, "+", "ID=gene-ORI1;gene_biotype=protein_coding"),
        (CHRM, "RefSeq", "CDS", 19_800, 20_300, "+", "ID=cds-ORI1;Parent=gene-ORI1"),
        # A trans-spliced gene, written as NCBI writes maize nad1: strand ?, the
        # mRNA's extent taken from one segment, so two exons lie outside it.
        (CHRM, "RefSeq", "gene", 5_000, 5_300, "?", "ID=gene-NAD1;gene=nad1;exception=trans-splicing"),
        (CHRM, "RefSeq", "mRNA", 5_000, 5_300, "?", "ID=rna-NAD1;Parent=gene-NAD1;exception=trans-splicing%2C RNA editing"),
        (CHRM, "RefSeq", "exon", 1_000, 1_200, "-", "ID=exon-NAD1-1;Parent=rna-NAD1"),
        (CHRM, "RefSeq", "exon", 5_000, 5_300, "+", "ID=exon-NAD1-2;Parent=rna-NAD1"),
        (CHRM, "RefSeq", "exon", 9_000, 9_100, "-", "ID=exon-NAD1-3;Parent=rna-NAD1"),
    ]
    with (args.out / "genomic.gff").open("w") as handle:
        handle.write("##gff-version 3\n")
        handle.write(f"##sequence-region {CHR1} 1 {chr1_length}\n")
        handle.write(f"##sequence-region {CHRM} 1 {m}\n")
        handle.write(f"{CHR1}\tRefSeq\tregion\t1\t{chr1_length}\t.\t+\t.\tID={CHR1}:1..{chr1_length};genome=chromosome;mol_type=genomic DNA\n")
        for seqid, source, ftype, start, end, strand, attrs in rows:
            phase = "0" if ftype == "CDS" else "."
            handle.write(f"{seqid}\t{source}\t{ftype}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n")
        for seqid, source, ftype, start, end, strand, attrs in mito:
            phase = "0" if ftype == "CDS" else "."
            handle.write(f"{seqid}\t{source}\t{ftype}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n")
    with (args.out / "assembly.fa").open("w") as handle:
        for name, length in ((CHR1, chr1_length), (CHRM, m)):
            handle.write(f">{name}\n")
            sequence = "".join(rnd.choice("ACGT") for _ in range(length))
            for i in range(0, length, 80):
                handle.write(sequence[i:i + 80] + "\n")
    print(f"{args.out}: {chr1_length + m:,} bp, {len(rows) + len(mito) + 1:,} rows")


if __name__ == "__main__":
    main()
