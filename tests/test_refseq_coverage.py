"""RefSeq content the benchmark's 30 genomes showed cgload refusing (D-097, D-098).

What these assert, stated before they were written:

* A `precursor_RNA`, as older NCBI GFF3 writes a microRNA precursor
  (*Xenopus tropicalis*), carries its exons and its mature miRNA, as a
  `primary_transcript` does.
* A `pseudogenic_tRNA` (rice and grape mitochondria) carries its exon, as a
  `pseudogenic_transcript` does.
* The two rows of a trans-spliced gene may overlap when they lie on opposite
  strands (*Drosophila* mod(mdg4)). The same overlap on one strand, or without
  the trans-splicing qualifier, is still refused.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import NormaliseError, normalise
from cgload.parsers.normalise import canonical_type as canonical
from cgload.parsers.tokenizer import tokenize

HEADER = "##gff-version 3\n##sequence-region chr1 1 100000\n"


def normalised(tmp_path: Path, body: str):
    path = tmp_path / "a.gff3"
    path.write_text(HEADER + body)
    tokenized = tokenize(path)
    return normalise(tokenized.features, dialects.REFSEQ, comments=tokenized.comments)


def by_id(result) -> dict:
    return {item.source_id: item for item in result.features}


def test_a_precursor_rna_carries_its_exons_and_its_mirna(tmp_path):
    result = normalised(tmp_path, (
        "chr1\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=gene-mir103-2;gene_biotype=miRNA\n"
        "chr1\tRefSeq\tprecursor_RNA\t100\t200\t.\t+\t.\tID=rna-NR_049673.1;"
        "Parent=gene-mir103-2\n"
        "chr1\tRefSeq\texon\t100\t200\t.\t+\t.\tID=exon-NR_049673.1-1;Parent=rna-NR_049673.1\n"
        "chr1\tRefSeq\tmiRNA\t120\t141\t.\t+\t.\tID=rna-mir103-2;Parent=rna-NR_049673.1\n"
        "chr1\tRefSeq\texon\t120\t141\t.\t+\t.\tID=exon-mir103-2-1;Parent=rna-mir103-2\n"
    ))
    features = by_id(result)
    assert features["rna-NR_049673.1"].feature_type == "mRNA"
    assert features["rna-NR_049673.1"].source_type == "precursor_RNA"
    assert features["rna-mir103-2"].parent_source_id == "rna-NR_049673.1"
    assert features["exon-NR_049673.1-1"].parent_source_id == "rna-NR_049673.1"


@pytest.mark.parametrize("kind", ["pseudogenic_tRNA", "pseudogenic_rRNA"])
def test_a_pseudogenic_rna_carries_its_exon(tmp_path, kind):
    result = normalised(tmp_path, (
        "chr1\tRefSeq\tpseudogene\t500\t572\t.\t+\t.\tID=gene-t02;pseudo=true\n"
        f"chr1\tRefSeq\t{kind}\t500\t572\t.\t+\t.\tID=rna-t02;Parent=gene-t02\n"
        "chr1\tRefSeq\texon\t500\t572\t.\t+\t.\tID=exon-t02-1;Parent=rna-t02\n"
    ))
    assert by_id(result)["exon-t02-1"].parent_source_id == "rna-t02"
    assert canonical(kind) == canonical("pseudogenic_transcript")


TRANS_SPLICED = (
    "chr1\tRefSeq\tgene\t21342226\t21370249\t.\t+\t.\tID=gene-m;exception=trans-splicing\n"
    "chr1\tRefSeq\tgene\t21370101\t21377401\t.\t{strand}\t.\tID=gene-m;{qualifier}\n"
)


def test_trans_spliced_rows_on_opposite_strands_may_overlap(tmp_path):
    result = normalised(tmp_path, TRANS_SPLICED.format(strand="-",
                                                       qualifier="exception=trans-splicing"))
    gene = by_id(result)["gene-m"]
    assert gene.strand == "?"
    assert [(s.start, s.end, s.strand) for s in gene.segments] == [
        (21342226, 21370249, "+"), (21370101, 21377401, "-")]


def test_trans_spliced_rows_on_one_strand_still_may_not_overlap(tmp_path):
    with pytest.raises(NormaliseError, match="overlapping"):
        normalised(tmp_path, TRANS_SPLICED.format(strand="+",
                                                  qualifier="exception=trans-splicing"))


def test_an_overlap_across_strands_without_the_qualifier_is_still_refused(tmp_path):
    with pytest.raises(NormaliseError, match="reappears with strand"):
        normalised(tmp_path, TRANS_SPLICED.format(strand="-", qualifier="note=x"))
