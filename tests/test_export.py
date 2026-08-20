"""Milestone 9: GFF3 export and the round trip.

The round-trip assertion is **on feature counts per source type, never on
identifier sets**. 175 of 350 child rows in the AUGUSTUS fixture carry no ID
(D-024), and synthesising one on export would break D-023. The stricter
assertion fails, and the obvious way to make it pass is exactly the thing these
rules exist to prevent.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest
from cgload.export import (
    EXPORT_PRAGMA,
    escape_attribute,
    export_gff3,
    is_cgload_export,
    round_trip_counts,
)
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import normalise
from cgload.parsers.tokenizer import tokenize

FIXTURES = Path(__file__).parent / "fixtures"
ALL_FIXTURES = sorted(FIXTURES.glob("*.gff3")) + sorted(FIXTURES.glob("*.gbff"))


def _load(path: Path):
    tokenized = tokenize(path)
    detected = dialects.detect(tokenized)
    result = normalise(
        tokenized.features, detected.profile, comments=tokenized.comments
    )
    return result.features, detected.profile


@pytest.mark.parametrize("path", ALL_FIXTURES, ids=lambda p: p.name)
def test_round_trip_is_exact_for_every_fixture(path):
    """Export then count independently; every source type must match.

    Verified across all five dialects: refseq, funannotate, augustus, genbank
    prokaryotic and genbank eukaryotic.
    """
    features, _profile = _load(path)
    text, summary = export_gff3(features)
    assert round_trip_counts(text) == Counter(f.source_type for f in features)
    assert summary.features == len(features)


def test_round_trip_survives_anonymous_rows():
    """AUGUSTUS emits intron, start_codon and stop_codon with only a Parent --
    175 of 350 child rows. They export without an ID and re-import as their own
    features, because both the parser and the naive counter treat an ID-less row
    that way (D-024)."""
    features, _profile = _load(FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    text, summary = export_gff3(features)
    assert summary.anonymous_features > 100
    assert round_trip_counts(text) == Counter(f.source_type for f in features)


def test_round_trip_catches_a_dropped_feature():
    """Mutation test. A check that never fires proves nothing."""
    features, _profile = _load(FIXTURES / "refseq_fgraminearum_NC_026477.gff3")
    complete = Counter(f.source_type for f in features)
    text, _summary = export_gff3(features[:-3])
    assert round_trip_counts(text) != complete


def test_column_three_is_the_source_type_not_the_canonical_term():
    """Exporting `feature_type` would rewrite AUGUSTUS's `transcript` to `mRNA`
    and make the round trip lossy on any dialect using a synonym -- the same
    reason `verify` counts `source_type`."""
    features, _profile = _load(FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    transcripts = [f for f in features if f.source_type == "transcript"]
    assert transcripts, "fixture lost its AUGUSTUS transcript rows"
    assert transcripts[0].feature_type == "mRNA", "canonical term differs, as intended"

    text, _summary = export_gff3(features)
    emitted = {
        line.split("\t")[2]
        for line in text.splitlines()
        if not line.startswith("#") and len(line.split("\t")) >= 9
    }
    assert "transcript" in emitted
    assert "mRNA" not in emitted


def test_column_two_is_per_feature_so_mixed_sources_survive():
    """An assembly loaded from several files keeps each feature's own column 2,
    with no extra machinery."""
    augustus, _ = _load(FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    refseq, _ = _load(FIXTURES / "refseq_fgraminearum_NC_026477.gff3")
    text, _summary = export_gff3(list(augustus[:5]) + list(refseq[:5]))
    programs = {
        line.split("\t")[1]
        for line in text.splitlines()
        if not line.startswith("#") and len(line.split("\t")) >= 9
    }
    assert programs == {"AUGUSTUS", "RefSeq"}


def test_nothing_is_synthesised():
    """AUGUSTUS in, no exon rows out. Stage 5 refused to invent exons because
    where AUGUSTUS predicts a UTR a synthesised exon would be wrong; export must
    not reintroduce them by the back door."""
    features, _profile = _load(FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    text, _summary = export_gff3(features)
    emitted = [
        line
        for line in text.splitlines()
        if not line.startswith("#") and len(line.split("\t")) >= 9
        and line.split("\t")[2] == "exon"
    ]
    assert not emitted


def test_absent_types_are_declared_in_the_header():
    """"The input had no exons" and "the exons were lost" must be
    distinguishable from the file alone."""
    features, _profile = _load(FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    text, summary = export_gff3(features, declared_absent=["exon"])
    assert summary.absent_types == ("exon",)
    assert "exon" in text.split("\n##")[0] or "does not emit" in text


def test_export_declares_itself():
    features, _profile = _load(FIXTURES / "refseq_mycmay_GCF_000328475.gff3")
    text, _summary = export_gff3(features, version="2.0.0")
    assert text.startswith("##gff-version 3\n")
    assert text.splitlines()[1].startswith(EXPORT_PRAGMA)
    assert is_cgload_export(text)


def test_a_source_file_is_not_mistaken_for_an_export():
    for path in ALL_FIXTURES:
        assert not is_cgload_export(path.read_text(errors="replace"))


def test_reimported_export_detects_as_cgload():
    features, _profile = _load(FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3")
    text, _summary = export_gff3(features)
    import tempfile

    path = Path(tempfile.mkdtemp()) / "reexport.gff3"
    path.write_text(text)
    assert dialects.detect(tokenize(path)).profile.name == "cgload"


def test_export_is_byte_stable():
    """Two exports of the same data must be identical, or a diff between them
    shows spurious changes and nobody trusts it."""
    features, _profile = _load(FIXTURES / "genbank_maeruginosa_NC_010296.gbff")
    first, _ = export_gff3(features)
    second, _ = export_gff3(features)
    assert first == second


@pytest.mark.parametrize(
    ("raw", "encoded"),
    [
        ("plain", "plain"),
        ("PH-1; NRRL 31084", "PH-1%3B NRRL 31084"),
        ("a=b", "a%3Db"),
        ("x,y", "x%2Cy"),
        ("100%", "100%25"),
        ("already%3B", "already%253B"),
    ],
)
def test_reserved_characters_are_percent_encoded(raw, encoded):
    """RefSeq encodes `;` as `%3B` on input. A value round-tripping through
    cgload must come back encoded the same way, or the re-import splits the
    attribute differently from the original parse."""
    assert escape_attribute(raw) == encoded
