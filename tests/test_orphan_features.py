"""A feature that declares no parent is top-level, whatever its type.

GFF3 nowhere requires an exon to be parented. cgload refused parentless
`exon` rows because `ALLOWED_PARENTS["exon"]` is non-empty, conflating two
different things: a feature naming a parent it may not have, and a feature
naming none at all. Only the first hides something.

Found on *Zea mays* B73 GCF_902167145.1 at line 992,090 -- two chloroplast
exons flanking a parentless intron, all carrying `number=2`, the second exon
block of a trans-spliced gene whose first exon lies elsewhere.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import NormaliseError, normalise
from cgload.parsers.tokenizer import tokenize

FIXTURE = Path(__file__).parent / "fixtures" / "refseq_maize_orphan_exons.gff3"


def _normalise(path: Path):
    tokenized = tokenize(path)
    detected = dialects.detect(tokenized)
    return normalise(tokenized.features, detected.profile, comments=tokenized.comments)


def test_parentless_exons_load_as_top_level():
    result = _normalise(FIXTURE)
    top_level = {f.source_id for f in result.features if f.parent_source_id is None}
    assert "id-NC_001666.2:129636..129867-2" in top_level
    assert "id-NC_001666.2:130408..130436-2" in top_level


def test_the_orphan_intron_loads_too():
    """`intron` is a Tier B type, so it was never subject to the rule -- but it
    belongs to the same block and must survive alongside its exons."""
    result = _normalise(FIXTURE)
    assert any(f.source_type == "intron" for f in result.features)


def test_coordinates_and_identifiers_are_preserved():
    """The identifier is NCBI's coordinate-derived fallback and must be stored
    verbatim (D-023), not normalised into something tidier."""
    result = _normalise(FIXTURE)
    orphan = next(
        f for f in result.features
        if f.source_id == "id-NC_001666.2:129636..129867-2"
    )
    assert (orphan.start, orphan.end) == (129636, 129867)
    assert orphan.attributes.get("number") == ["2"]


def test_a_disallowed_named_parent_is_still_refused(tmp_path):
    """The rule's real purpose, unchanged: naming a parent you may not have is
    an error. An exon may not hang off another exon."""
    body = (
        "NC_001666.2\tRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-A\n"
        "NC_001666.2\tRefSeq\tmRNA\t100\t900\t.\t+\t.\tID=rna-A;Parent=gene-A\n"
        "NC_001666.2\tRefSeq\texon\t100\t900\t.\t+\t.\tID=e-A;Parent=rna-A\n"
        "NC_001666.2\tRefSeq\texon\t200\t300\t.\t+\t.\tID=e-bad;Parent=e-A\n"
    )
    path = tmp_path / "bad.gff3"
    path.write_text("##gff-version 3\n" + body)
    with pytest.raises(NormaliseError, match="may only attach"):
        _normalise(path)


def test_parented_features_in_the_same_file_still_link():
    """The narrowing must not detach anything that did declare a parent."""
    result = _normalise(FIXTURE)
    exon = next(f for f in result.features if f.source_id == "exon-ZemaCp099-1")
    assert exon.parent_source_id == "rna-ZemaCp099"
