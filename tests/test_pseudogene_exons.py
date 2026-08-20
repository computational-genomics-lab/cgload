"""RefSeq attaches exons directly to a pseudogene, and that is not a defect.

Measured on *Zea mays* B73 GCF_902167145.1: 5,195 pseudogenes and 12,757 exons
whose parent is a gene-level row. cgload refused the first of them at line 247
of chromosome 1 -- standard NCBI output, rejected by two rules that had been
written from what four fungal fixtures happened to contain.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import NormaliseError, normalise
from cgload.parsers.tokenizer import tokenize


def gff3(tmp_path: Path, body: str, name: str = "sample.gff3") -> Path:
    path = tmp_path / name
    path.write_text("##gff-version 3\n" + body)
    return path


def row(*columns: str) -> str:
    return "\t".join(columns) + "\n"


#: A gene with a proper transcript layer, so `has_transcripts` is True and the
#: rule under test is live.
WITH_TRANSCRIPT = (
    row("NC_050096.1", "Gnomon", "gene", "200404", "203170", ".", "-", ".",
        "ID=gene-A;Name=A")
    + row("NC_050096.1", "Gnomon", "mRNA", "200404", "203170", ".", "-", ".",
          "ID=rna-A;Parent=gene-A")
    + row("NC_050096.1", "Gnomon", "exon", "200404", "203170", ".", "-", ".",
          "ID=e-A;Parent=rna-A")
)

#: Verbatim shape of maize `gene-LOC103630212`, the row that stopped the load.
PSEUDOGENE_BLOCK = (
    row("NC_050096.1", "Gnomon", "pseudogene", "198423", "198910", ".", "-", ".",
        "ID=gene-LOC103630212;Name=LOC103630212;gbkey=Gene;"
        "pseudo=true;gene_biotype=pseudogene")
    + row("NC_050096.1", "Gnomon", "exon", "198423", "198428", ".", "-", ".",
          "ID=id-LOC103630212;Parent=gene-LOC103630212;gbkey=exon")
    + row("NC_050096.1", "Gnomon", "exon", "198535", "198910", ".", "-", ".",
          "ID=id-LOC103630212-2;Parent=gene-LOC103630212;gbkey=exon")
)


def _normalise(path: Path):
    tokenized = tokenize(path)
    detected = dialects.detect(tokenized)
    return normalise(tokenized.features, detected.profile, comments=tokenized.comments)


def test_pseudogene_exons_load(tmp_path):
    """A pseudogene has no functional transcript to model, so its exons attach
    to it directly. There is no missing layer being hidden."""
    result = _normalise(gff3(tmp_path, WITH_TRANSCRIPT + PSEUDOGENE_BLOCK))
    exons = [f for f in result.features if f.source_type == "exon"]
    assert len(exons) == 3
    pseudo_exons = [f for f in exons if f.parent_source_id == "gene-LOC103630212"]
    assert len(pseudo_exons) == 2


def test_a_gene_missing_its_transcript_is_still_refused(tmp_path):
    """The rule's real purpose, unchanged. A `gene` -- not a pseudogene -- whose
    transcript is absent while the rest of the file has one is an incomplete
    hierarchy dressed as a valid one."""
    body = WITH_TRANSCRIPT + (
        row("NC_050096.1", "Gnomon", "gene", "300000", "301000", ".", "-", ".",
            "ID=gene-B;Name=B")
        + row("NC_050096.1", "Gnomon", "exon", "300000", "300500", ".", "-", ".",
              "ID=e-B;Parent=gene-B")
    )
    with pytest.raises(NormaliseError, match="no transcript of its own"):
        _normalise(gff3(tmp_path, body))


def test_the_refusal_names_the_pseudogene_exemption(tmp_path):
    """Someone hitting this error should not have to read the source to learn
    that a pseudogene parent is treated differently."""
    body = WITH_TRANSCRIPT + (
        row("NC_050096.1", "Gnomon", "gene", "300000", "301000", ".", "-", ".",
            "ID=gene-B;Name=B")
        + row("NC_050096.1", "Gnomon", "exon", "300000", "300500", ".", "-", ".",
              "ID=e-B;Parent=gene-B")
    )
    with pytest.raises(NormaliseError, match="pseudogene parent is exempt"):
        _normalise(gff3(tmp_path, body))


def test_exon_may_name_a_pseudogene_parent():
    """The second rule with the same defect. `CDS` already permitted a
    pseudogene parent; `exon` did not, and the omission alone refused the file
    even after the transcript-layer rule was narrowed."""
    from cgload.db.vocabulary import ALLOWED_PARENTS

    assert "pseudogene" in ALLOWED_PARENTS["exon"]
    assert "pseudogene" in ALLOWED_PARENTS["CDS"]


def test_a_pseudogene_transcript_still_works(tmp_path):
    """Some pseudogenes do carry a transcript -- RefSeq emits `mRNA` under
    `pseudogene` for transcribed pseudogenes. Both shapes must load."""
    body = WITH_TRANSCRIPT + (
        row("NC_050096.1", "Gnomon", "pseudogene", "400000", "401000", ".", "+", ".",
            "ID=gene-P2;Name=P2")
        + row("NC_050096.1", "Gnomon", "mRNA", "400000", "401000", ".", "+", ".",
              "ID=rna-P2;Parent=gene-P2")
        + row("NC_050096.1", "Gnomon", "exon", "400000", "400500", ".", "+", ".",
              "ID=e-P2;Parent=rna-P2")
    )
    result = _normalise(gff3(tmp_path, body))
    assert any(f.source_id == "rna-P2" for f in result.features)
