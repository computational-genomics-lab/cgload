"""Sequence for the portal: reverse complement, splicing, translation (D-113)."""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.portal.sequence import (
    GENETIC_CODES,
    AssemblySequence,
    Piece,
    coding_sequence,
    reverse_complement,
    transcript_order,
    translate,
)


def test_reverse_complement_keeps_ambiguity_codes_and_case():
    assert reverse_complement("ACGTNRYacgt") == "acgtRYNACGT"


def test_standard_code():
    assert translate("ATGGCCTAA") == "MA*"
    assert translate("ATGGCCTA") == "MA"            # a final incomplete codon is dropped
    assert translate("ATGNNNGCC") == "MXA"          # an unknown base makes an X


def test_alternative_start_codons_read_as_methionine_only_at_a_start():
    assert translate("TTGTTG") == "ML"
    assert translate("TTGTTG", from_start=False) == "LL"


def test_vertebrate_mitochondrial_code_differs_where_ncbi_says():
    assert translate("TGAAGAATA", 2, from_start=False) == "W*M"
    assert translate("TGAAGAATA", 1, from_start=False) == "*RI"


def test_every_table_has_64_codons_and_known_stops():
    for table, (amino, starts, name) in GENETIC_CODES.items():
        assert len(amino) == len(starts) == 64, table
        assert set(starts) <= {"-", "M"} and name
    assert translate("TAA", 1) == "*" and translate("TAA", 6) == "Q"


def test_an_unknown_table_is_refused():
    with pytest.raises(ValueError):
        translate("ATG", 7)


def test_transcript_order_follows_the_strand():
    pieces = [Piece(300, 400, "-"), Piece(100, 200, "-")]
    assert [piece.start for piece in transcript_order(pieces)] == [300, 100]
    plus = [Piece(300, 400, "+"), Piece(100, 200, "+")]
    assert [piece.start for piece in transcript_order(plus)] == [100, 300]


def test_trans_spliced_rows_keep_the_files_order():
    pieces = [Piece(100, 200, "+"), Piece(900, 950, "-"), Piece(500, 600, "+")]
    ordered = transcript_order(pieces, file_order=[3, 1, 2])
    assert [piece.start for piece in ordered] == [900, 500, 100]


@pytest.fixture
def assembly(tmp_path: Path) -> AssemblySequence:
    from pyfaidx import Fasta

    path = tmp_path / "a.fa"
    path.write_text(">chr\nAAACCCGGGTTTAAACCCGG\n>circ\nATGCATGCAT\n")
    Fasta(str(path)).close()  # writes the .fai, as load does
    return AssemblySequence(path, circular={"circ"})


def test_fetch_is_one_based_and_inclusive(assembly):
    assert assembly.fetch("chr", 1, 3) == "AAA"
    assert assembly.fetch("chr", 4, 9) == "CCCGGG"


def test_a_row_past_the_end_wraps_only_on_a_circular_sequence(assembly):
    assert assembly.fetch("circ", 9, 12) == "ATAT"
    assert assembly.fetch("circ", 11, 13) == "ATG"
    with pytest.raises(ValueError):
        assembly.fetch("chr", 19, 22)


def test_minus_strand_pieces_are_read_reverse_complemented(assembly):
    assert assembly.piece("chr", Piece(7, 12, "-")) == "AAACCC"


def test_phase_trims_the_first_row(assembly):
    pieces = [Piece(4, 9, "+", phase=2), Piece(13, 15, "+", phase=1)]
    assert coding_sequence(assembly, "chr", pieces) == "CGGGAAA"


def test_the_vendored_igv_is_the_one_recorded():
    """D-114: the file served is the npm release VENDORED.md names."""
    import hashlib
    import re

    import cgload.portal

    static = Path(cgload.portal.__file__).parent / "static" / "igv"
    recorded = re.search(r"SHA-256: `([0-9a-f]{64})`", (static / "VENDORED.md").read_text())
    actual = hashlib.sha256((static / "igv.min.js").read_bytes()).hexdigest()
    assert recorded and recorded.group(1) == actual
    assert "MIT License" in (static / "LICENSE").read_text()
