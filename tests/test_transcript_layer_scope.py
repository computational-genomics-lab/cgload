"""The scope of the transcript-bypass rule, pinned in both directions.

D-043 gave the rule per-gene scope. D-058 exempted `pseudogene` parents. This
file pins what remains, including the one edge that looks like a bug and is
load-bearing: `has_transcripts` tests for `mRNA` alone rather than for
TRANSCRIPT_TERMS (D-059). Widening it refuses every prokaryotic record, because
bacterial genomes carry tRNA and rRNA genes. ``test_widening_has_transcripts_to_
all_transcript_terms_would_refuse_prokaryotes`` states that cost as an
executable fact so the next reader does not have to rediscover it.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.db.vocabulary import TRANSCRIPT_TERMS
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import NormaliseError, normalise
from cgload.parsers.tokenizer import tokenize

C = "NC_050096.1"


def row(*columns: str) -> str:
    return "\t".join(columns) + "\n"


def load(tmp_path: Path, body: str):
    path = tmp_path / "sample.gff3"
    path.write_text("##gff-version 3\n" + body)
    tokenized = tokenize(path)
    detected = dialects.detect(tokenized)
    return normalise(tokenized.features, detected.profile, comments=tokenized.comments)


#: A eukaryotic gene with its transcript, so `has_transcripts` is True.
MRNA_LAYER = (
    row(C, "Gnomon", "gene", "200404", "203170", ".", "-", ".", "ID=gene-A")
    + row(C, "Gnomon", "mRNA", "200404", "203170", ".", "-", ".", "ID=rna-A;Parent=gene-A")
    + row(C, "Gnomon", "exon", "200404", "203170", ".", "-", ".", "ID=e-A;Parent=rna-A")
)

#: A gene whose transcript is absent, its exon hanging off the gene.
BARE_GENE = (
    row(C, "Gnomon", "gene", "300000", "301000", ".", "-", ".", "ID=gene-B")
    + row(C, "Gnomon", "exon", "300000", "300500", ".", "-", ".", "ID=e-B;Parent=gene-B")
)

#: A prokaryotic record. `gene`->`CDS` with no transcript is the normal shape
#: (D-036), and the tRNA and rRNA genes beside it are present in every real
#: bacterial genome -- which is exactly what makes them dangerous here.
#:
#: Column 2 is `RefSeq`, as a prokaryotic RefSeq GFF3 actually writes it. The
#: first version put `GenBank` here, which selected the *GenBank* profile for a
#: GFF3 file -- and that profile requires each computed edge to carry the rule
#: that produced it (D-032), a stamp only the GenBank tokenizer applies. The test
#: failed on a category error in its own fixture rather than on anything the rule
#: under test does.
PROKARYOTE = (
    row(C, "RefSeq", "gene", "100", "1000", ".", "+", ".", "ID=gene-b0001;locus_tag=b0001")
    + row(C, "RefSeq", "CDS", "100", "1000", ".", "+", "0", "ID=cds-b0001;Parent=gene-b0001")
    + row(C, "RefSeq", "gene", "2000", "2076", ".", "+", ".", "ID=gene-b0002;locus_tag=b0002")
    + row(C, "RefSeq", "tRNA", "2000", "2076", ".", "+", ".", "ID=rna-b0002;Parent=gene-b0002")
    + row(C, "RefSeq", "gene", "3000", "4500", ".", "+", ".", "ID=gene-b0003;locus_tag=b0003")
    + row(C, "RefSeq", "rRNA", "3000", "4500", ".", "+", ".", "ID=rna-b0003;Parent=gene-b0003")
)


def test_a_gene_missing_its_transcript_is_refused(tmp_path):
    """The rule still does its job in the case it was written for."""
    with pytest.raises(NormaliseError, match="no transcript of its own"):
        load(tmp_path, MRNA_LAYER + BARE_GENE)


def test_a_pseudogene_parent_is_exempt(tmp_path):
    """D-058. A pseudogene has no functional transcript to conceal."""
    result = load(
        tmp_path,
        MRNA_LAYER
        + row(C, "Gnomon", "pseudogene", "198423", "198910", ".", "-", ".",
              "ID=pg-1;pseudo=true")
        + row(C, "Gnomon", "exon", "198423", "198428", ".", "-", ".", "ID=pe-1;Parent=pg-1")
        + row(C, "Gnomon", "exon", "198535", "198910", ".", "-", ".", "ID=pe-2;Parent=pg-1"),
    )
    assert sum(1 for f in result.features if f.parent_source_id == "pg-1") == 2


def test_a_prokaryotic_record_loads(tmp_path):
    """D-036. No transcript layer, so the bypass rule must not fire at all."""
    result = load(tmp_path, PROKARYOTE)
    assert result.feature_count == 6


def test_widening_has_transcripts_to_all_transcript_terms_would_refuse_prokaryotes(
    tmp_path,
):
    """Why `has_transcripts` tests `mRNA` and not TRANSCRIPT_TERMS (D-059).

    The two halves of the rule look inconsistent -- `genes_with_transcripts`
    uses all four transcript terms and `has_transcripts` uses one. Widening the
    latter to match is the obvious tidy-up and it is wrong: every bacterial
    genome has tRNA and rRNA genes, so the file would read as transcript-bearing
    and every prokaryotic `gene`->`CDS` edge would be refused.

    This test asserts the *premise* of that failure rather than monkey-patching
    the module: the prokaryotic fixture contains transcript terms, contains no
    mRNA, and must load. Any change that makes `has_transcripts` true here will
    break ``test_a_prokaryotic_record_loads`` above.
    """
    result = load(tmp_path, PROKARYOTE)
    present = {f.feature_type for f in result.features}
    assert present & TRANSCRIPT_TERMS, "fixture must contain transcript terms"
    assert "mRNA" not in present, "fixture must contain no mRNA"
    assert result.feature_count == 6


@pytest.mark.xfail(
    reason="D-059, accepted gap: a eukaryotic file whose transcript layer is "
    "entirely non-coding reads as prokaryotic, so the bypass rule never runs "
    "and a genuinely missing transcript loads unchecked. Not fixed because no "
    "predicate available at this point separates it from a prokaryotic record. "
    "If this test ever passes, the predicate was improved -- update D-059.",
    strict=True,
)
def test_ncrna_only_eukaryote_with_a_missing_transcript_is_refused(tmp_path):
    ncrna_layer = (
        row(C, "Gnomon", "gene", "200404", "203170", ".", "-", ".", "ID=gene-N")
        + row(C, "Gnomon", "lnc_RNA", "200404", "203170", ".", "-", ".",
              "ID=rna-N;Parent=gene-N")
        + row(C, "Gnomon", "exon", "200404", "203170", ".", "-", ".", "ID=e-N;Parent=rna-N")
    )
    with pytest.raises(NormaliseError, match="no transcript of its own"):
        load(tmp_path, ncrna_layer + BARE_GENE)
