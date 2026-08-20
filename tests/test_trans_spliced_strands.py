"""D-068: a trans-spliced gene spans both strands.

Found by loading the maize *mitochondrion*, which the earlier scale run had
excluded. Third form of the same NCBI trans-splicing converter defect: after
too-small envelopes (D-046) and one inverted row (D-056), a gene written as two
rows sharing an ID and differing only in strand.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import NormaliseError, normalise
from cgload.parsers.tokenizer import tokenize

FIXTURE = Path(__file__).parent / "fixtures" / "broken" / "trans_spliced_strands.gff3"


def _normalise(features):
    tokenized = tokenize(FIXTURE)
    return normalise(features, dialects.REFSEQ, comments=tokenized.comments)


def _load(exclude: str | None = None):
    tokenized = tokenize(FIXTURE)
    features = [
        item
        for item in tokenized.features
        if exclude is None or exclude not in (item.id or "") + (item.first("Parent") or "")
    ]
    return normalise(features, dialects.REFSEQ, comments=tokenized.comments)


def test_a_trans_spliced_gene_merges_across_strands():
    """The segments are genuinely on opposite strands, so no single value is
    correct and the merged feature takes GFF3's `?` -- which is what NCBI itself
    writes on the mRNA row of the same gene."""
    result = _load(exclude="bad")
    gene = next(item for item in result.features if item.source_id == "gene-ZeamMp186")

    assert len(gene.segments) == 2
    assert (gene.start, gene.end) == (50490, 267232)
    assert gene.strand == "?", "a gene spanning both strands must not claim one"


def test_an_unexplained_strand_disagreement_is_still_refused():
    """The permission is gated on the qualifier, not on the dialect. Without
    `exception=trans-splicing` two rows sharing an ID and differing in strand are
    two features, and merging them would lose one."""
    with pytest.raises(NormaliseError, match="reappears with strand"):
        _load()


def test_the_qualifier_is_matched_on_its_value_not_its_key(tmp_path):
    """`exception=` covers several unrelated NCBI annotations -- ribosomal
    slippage, rearrangement required for product, annotated by transcript data.
    Only trans-splicing licenses a strand disagreement, so matching on the key
    alone would silently merge genes flagged for something else entirely."""
    path = tmp_path / "other_exception.gff3"
    path.write_text(
        "##gff-version 3\n"
        "c\tRefSeq\tgene\t100\t200\t.\t+\t.\tID=g;exception=rearrangement required for product\n"
        "c\tRefSeq\tgene\t300\t400\t.\t-\t.\tID=g;exception=rearrangement required for product\n"
    )
    tokenized = tokenize(path)
    with pytest.raises(NormaliseError, match="reappears with strand"):
        normalise(tokenized.features, dialects.REFSEQ, comments=tokenized.comments)


def test_only_the_dialects_that_need_it_carry_the_permission():
    """A permission that applies everywhere is not a permission, it is a hole."""
    assert dialects.REFSEQ.strand_exception_qualifier == "exception"
    for profile in (dialects.FUNANNOTATE, dialects.AUGUSTUS, dialects.GENBANK):
        assert profile.strand_exception_qualifier is None, profile.name
