"""Milestone 3 tests: tokenizer, dialect detection, normalisation.

What these assert, stated before they were written:

* The tokenizer reads the physical format and nothing else — nine columns,
  `#` comments anywhere, GFF3 attribute syntax with the right decode order.
* Detection identifies a registered dialect, names a recognised-but-unregistered
  one, and refuses everything else without ever falling back to a permissive
  reader.
* Normalisation groups CDS segments into one feature, keeps each anonymous row
  separate, numbers segments by coordinate, and refuses every way of quietly
  merging two features into one.
* The AUGUSTUS protein oracle catches a dropped segment.
* The three real fixtures pass end to end, with the counts they should have.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.parsers import profiles as P
from cgload.parsers.normalise import (
    NormaliseError,
    check_declared_proteins,
    counts_by_source_type,
    normalise,
)
from cgload.parsers.tokenizer import TokenizeError, parse_attributes, tokenize

FIXTURES = Path(__file__).parent / "fixtures"
REFSEQ_FILE = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
FUNANNOTATE_FILE = FIXTURES / "funannotate_aphaeospermum_pheo_arth1.gff3"
AUGUSTUS_FILE = FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3"


def gff3(tmp_path: Path, body: str, name: str = "t.gff3") -> Path:
    path = tmp_path / name
    path.write_text("##gff-version 3\n" + body)
    return path


def row(*columns: str) -> str:
    return "\t".join(columns) + "\n"


# -- tokenizer -------------------------------------------------------------


def test_nine_columns_are_required(tmp_path):
    path = gff3(tmp_path, "chr1\tRefSeq\tgene\t1\t100\t.\t+\n")
    with pytest.raises(TokenizeError, match="found 7"):
        tokenize(path)


def test_error_names_the_file_and_line(tmp_path):
    path = gff3(tmp_path, row("chr1", "RefSeq", "gene", "x", "100", ".", "+", ".", "ID=g1"))
    with pytest.raises(TokenizeError, match=r"t\.gff3:2:"):
        tokenize(path)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        (3, "0", "below 1"),
        (4, "5", "precedes start"),
        (6, "x", "not \\+, -, \\. or \\?"),
        (7, "3", "not 0, 1 or 2"),
        (5, "high", "not a number"),
    ],
)
def test_malformed_columns_are_refused(tmp_path, column, value, message):
    columns = ["chr1", "RefSeq", "CDS", "10", "100", ".", "+", ".", "ID=c1"]
    columns[column] = value
    with pytest.raises(TokenizeError, match=message):
        tokenize(gff3(tmp_path, row(*columns)))


def test_comments_are_skipped_wherever_they_appear(tmp_path):
    """AUGUSTUS writes evidence blocks between feature rows (finding 13)."""
    body = (
        "# start gene g1\n"
        + row("chr1", "AUGUSTUS", "gene", "1", "100", ".", "+", ".", "ID=g1")
        + "# protein sequence = [MKV]\n"
        + "# end gene g1\n"
        + row("chr1", "AUGUSTUS", "gene", "200", "300", ".", "+", ".", "ID=g2")
    )
    result = tokenize(gff3(tmp_path, body))
    assert result.data_line_count == 2
    assert len(result.comments) == 4


def test_percent_decoding_happens_after_splitting(tmp_path):
    """RefSeq hides a semicolon in a value as %3B (finding 9). Decoding first
    would split the attribute in two."""
    attributes = parse_attributes("strain=PH-1%3B NRRL 31084;gbkey=Src", Path("x"), 1)
    assert attributes["strain"] == ["PH-1; NRRL 31084"]
    assert set(attributes) == {"strain", "gbkey"}


def test_comma_separated_values_become_a_list(tmp_path):
    attributes = parse_attributes("Dbxref=GeneID:1,GenBank:XM_1", Path("x"), 1)
    assert attributes["Dbxref"] == ["GeneID:1", "GenBank:XM_1"]


def test_trailing_semicolon_is_not_an_empty_key(tmp_path):
    """funannotate ends every line with ';' (finding 10)."""
    attributes = parse_attributes("ID=PHE_1;Parent=PHE_0;", Path("x"), 1)
    assert set(attributes) == {"ID", "Parent"}


def test_attribute_without_equals_is_refused(tmp_path):
    with pytest.raises(TokenizeError, match="no '=' separator"):
        parse_attributes("ID=g1;broken", Path("x"), 1)


def test_attribute_lookup_is_case_insensitive(tmp_path):
    """funannotate writes DBxref, RefSeq writes Dbxref (finding 10)."""
    body = row("c", "funannotate", "gene", "1", "9", ".", "+", ".", "ID=g;DBxref=PFAM:1")
    result = tokenize(gff3(tmp_path, body))
    assert result.features[0].first("Dbxref") == "PFAM:1"


def test_a_file_with_no_feature_rows_is_refused(tmp_path):
    with pytest.raises(TokenizeError, match="no feature rows"):
        tokenize(gff3(tmp_path, "# only comments\n"))


def test_sequence_region_directive_is_captured(tmp_path):
    body = "##sequence-region chr1 1 5000\n" + row(
        "chr1", "RefSeq", "gene", "1", "9", ".", "+", ".", "ID=g"
    )
    assert tokenize(gff3(tmp_path, body)).declared_regions == {"chr1": (1, 5000)}


# -- detection -------------------------------------------------------------


def test_refseq_is_detected():
    detection = P.detect(tokenize(REFSEQ_FILE))
    assert detection.profile.name == "refseq"
    assert "RefSeq" in detection.evidence


def test_a_declared_but_unregistered_dialect_is_named_not_guessed(tmp_path, monkeypatch):
    """D-027's middle outcome. Every real dialect is registered at milestone 5,
    so this uses a synthetic profile: the mechanism must keep working, because it
    is what lets the *next* dialect be recognised before it is supported."""
    import dataclasses

    pretend = dataclasses.replace(
        P.REFSEQ, name="ensembl", sources=frozenset({"Ensembl"}), header_signature=None
    )
    monkeypatch.setattr(P, "KNOWN_UNREGISTERED", (pretend,))
    path = gff3(tmp_path, row("c", "Ensembl", "gene", "1", "9", ".", "+", ".", "ID=g"))
    with pytest.raises(P.DetectionError) as raised:
        P.detect(tokenize(path))
    assert "ensembl" in str(raised.value)
    assert "not yet enabled" in str(raised.value)


@pytest.mark.parametrize(
    ("path", "expected"),
    [(FUNANNOTATE_FILE, "funannotate"), (AUGUSTUS_FILE, "augustus")],
)
def test_the_other_two_dialects_are_now_registered(path, expected):
    """Milestone 5. They were refused at milestone 3 with their names; now they
    load, and the detection evidence is still reported."""
    detection = P.detect(tokenize(path))
    assert detection.profile.name == expected
    assert detection.evidence


def test_unknown_source_is_refused_with_its_evidence(tmp_path):
    path = gff3(tmp_path, row("c", "SomeNewTool", "gene", "1", "9", ".", "+", ".", "ID=g"))
    with pytest.raises(P.DetectionError, match="no confident match"):
        P.detect(tokenize(path))


def test_there_is_no_generic_fallback_profile():
    """The most tempting and most damaging addition (D-013). If a profile ever
    appears whose source set is empty or which matches anything, this fails."""
    for profile in P.REGISTERED + P.KNOWN_UNREGISTERED:
        assert profile.sources, f"{profile.name} claims no column-2 value"


def test_dialect_override_wins(tmp_path):
    path = gff3(tmp_path, row("c", "Mystery", "gene", "1", "9", ".", "+", ".", "ID=g"))
    detection = P.detect(tokenize(path), override="refseq")
    assert detection.confidence is P.Confidence.OVERRIDDEN


def test_unknown_override_is_refused(tmp_path):
    path = gff3(tmp_path, row("c", "RefSeq", "gene", "1", "9", ".", "+", ".", "ID=g"))
    with pytest.raises(P.DetectionError, match="unknown dialect"):
        P.detect(tokenize(path), override="ensembl")


def test_every_profile_declares_every_field():
    """Profiles fail closed (D-013): no field may be left to a default."""
    for profile in P.REGISTERED + P.KNOWN_UNREGISTERED:
        assert profile.provides, profile.name
        assert profile.parent_inference is not None, profile.name
        assert profile.stop_codon is not None, profile.name
        assert isinstance(profile.permits_anonymous_rows, bool), profile.name


def test_augustus_declares_no_exon_type():
    """Finding 5: so zero exons is a recorded expectation, not a silent drop."""
    assert "exon" not in P.AUGUSTUS.provides
    assert "exon" in P.REFSEQ.provides


# -- normalisation ---------------------------------------------------------


def normalised_from(path: Path, profile: P.Profile, **kwargs):
    tokenized = tokenize(path)
    return normalise(
        tokenized.features,
        profile,
        comments=tokenized.comments,
        unexpected=P.unexpected_types(profile, tokenized.features),
        missing=P.missing_types(profile, tokenized.features),
        **kwargs,
    )


def test_cds_segments_collapse_into_one_feature():
    result = normalised_from(REFSEQ_FILE, P.REFSEQ)
    counts = counts_by_source_type(result)
    # 56 CDS rows in the file; 10 CDS features, because segments share an ID.
    assert counts["CDS"] == 10
    assert counts["exon"] == 147


def test_segments_are_numbered_by_coordinate_not_file_order():
    """funannotate emits minus-strand segments 3'->5' (finding 7)."""
    result = normalised_from(FUNANNOTATE_FILE, P.FUNANNOTATE)
    multi = [f for f in result.features if len(f.segments) > 1]
    assert multi
    for feature in multi:
        starts = [segment.start for segment in feature.segments]
        indexes = [segment.segment_index for segment in feature.segments]
        assert starts == sorted(starts)
        assert indexes == list(range(len(indexes)))


def test_each_anonymous_row_is_its_own_feature():
    """D-024. 175 of AUGUSTUS's rows carry no ID."""
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    anonymous = [f for f in result.features if f.is_anonymous]
    assert len(anonymous) == 175
    assert {f.source_type for f in anonymous} == {"intron", "start_codon", "stop_codon"}


def test_anonymous_rows_are_refused_by_a_profile_that_forbids_them(tmp_path):
    body = row("c", "RefSeq", "gene", "1", "99", ".", "+", ".", "ID=g") + row(
        "c", "RefSeq", "exon", "1", "9", ".", "+", ".", "Parent=g"
    )
    with pytest.raises(NormaliseError, match="does not permit anonymous rows"):
        normalised_from(gff3(tmp_path, body), P.REFSEQ)


def test_codon_rows_stay_inert_and_do_not_inflate_cds():
    """D-022. 49 CDS features and 49+49 codon features, not 147 CDS."""
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    counts = counts_by_source_type(result)
    assert counts["CDS"] == 49
    assert counts["start_codon"] == 49
    assert counts["stop_codon"] == 49
    assert {f.feature_type for f in result.features if f.source_type == "stop_codon"} == {
        "stop_codon"
    }


def test_augustus_transcript_maps_to_mrna_while_source_type_is_kept():
    """D-017's two columns, on real input."""
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    transcripts = [f for f in result.features if f.source_type == "transcript"]
    assert transcripts
    assert all(f.feature_type == "mRNA" for f in transcripts)


def test_ids_are_verbatim_with_no_suffix_stripping():
    """D-023. Stripping -T1 would collide with the parent gene."""
    result = normalised_from(FUNANNOTATE_FILE, P.FUNANNOTATE)
    mrnas = [f for f in result.features if f.feature_type == "mRNA"]
    assert mrnas
    for feature in mrnas:
        assert str(feature.source_id).endswith("-T1")
        assert feature.source_id != feature.parent_source_id


# -- the ways two features could silently become one ----------------------


def test_repeated_id_on_a_non_segment_type_is_an_error(tmp_path):
    """Exons never share an ID (finding 6)."""
    body = (
        row("c", "RefSeq", "gene", "1", "999", ".", "+", ".", "ID=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m;Parent=g")
        + row("c", "RefSeq", "exon", "1", "99", ".", "+", ".", "ID=e;Parent=m")
        + row("c", "RefSeq", "exon", "200", "299", ".", "+", ".", "ID=e;Parent=m")
    )
    with pytest.raises(NormaliseError, match="repeats on a exon row"):
        normalised_from(gff3(tmp_path, body), P.REFSEQ)


@pytest.mark.parametrize(
    ("second", "message"),
    [
        (("c2", "RefSeq", "CDS", "200", "299", ".", "+", "0", "ID=c;Parent=m"), "seqid"),
        (("c", "RefSeq", "CDS", "200", "299", ".", "-", "0", "ID=c;Parent=m"), "strand"),
        (("c", "RefSeq", "CDS", "200", "299", ".", "+", "0", "ID=c;Parent=m2"), "parent"),
        (("c", "RefSeq", "CDS", "50", "150", ".", "+", "0", "ID=c;Parent=m"), "overlapping"),
    ],
)
def test_a_repeated_id_that_is_not_a_discontinuous_feature_is_an_error(
    tmp_path, second, message
):
    """The D-014 exemption, one violation at a time. Each of these would
    otherwise merge two features into one, silently."""
    body = (
        row("c", "RefSeq", "gene", "1", "999", ".", "+", ".", "ID=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m;Parent=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m2;Parent=g")
        + row("c2", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m3;Parent=g")
        + row("c", "RefSeq", "CDS", "1", "99", ".", "+", "0", "ID=c;Parent=m")
        + row(*second)
    )
    with pytest.raises(NormaliseError, match=message):
        normalised_from(gff3(tmp_path, body), P.REFSEQ)


def test_phase_and_score_may_differ_between_segments(tmp_path):
    """Excluded from the sameness check by design: phase varies across CDS
    segments, and AUGUSTUS gives each segment its own score."""
    body = (
        row("c", "RefSeq", "gene", "1", "999", ".", "+", ".", "ID=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m;Parent=g")
        + row("c", "RefSeq", "CDS", "1", "99", "0.9", "+", "0", "ID=c;Parent=m")
        + row("c", "RefSeq", "CDS", "200", "299", "0.1", "+", "2", "ID=c;Parent=m")
    )
    result = normalised_from(gff3(tmp_path, body), P.REFSEQ)
    cds = next(f for f in result.features if f.source_type == "CDS")
    assert [s.phase for s in cds.segments] == [0, 2]


def test_a_dangling_parent_is_an_error(tmp_path):
    body = row("c", "RefSeq", "exon", "1", "99", ".", "+", ".", "ID=e;Parent=absent")
    with pytest.raises(NormaliseError, match="not defined anywhere"):
        normalised_from(gff3(tmp_path, body), P.REFSEQ)


def test_multiple_parents_are_refused(tmp_path):
    body = (
        row("c", "RefSeq", "gene", "1", "999", ".", "+", ".", "ID=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m1;Parent=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m2;Parent=g")
        + row("c", "RefSeq", "exon", "1", "99", ".", "+", ".", "ID=e;Parent=m1,m2")
    )
    with pytest.raises(NormaliseError, match="2 parents"):
        normalised_from(gff3(tmp_path, body), P.REFSEQ)


def test_a_tier_a_child_cannot_hang_off_an_inert_parent(tmp_path):
    """Linkage rules assume Tier A parents (D-017)."""
    body = (
        row("c", "RefSeq", "sequence_feature", "1", "999", ".", "+", ".", "ID=x")
        + row("c", "RefSeq", "exon", "1", "99", ".", "+", ".", "ID=e;Parent=x")
    )
    with pytest.raises(NormaliseError, match="may only attach to"):
        normalised_from(gff3(tmp_path, body), P.REFSEQ)


def test_refseq_pseudogene_parents_are_accepted():
    """Finding 8: mRNA and tRNA off a pseudogene is real and correct."""
    result = normalised_from(REFSEQ_FILE, P.REFSEQ)
    pseudo = {
        f.source_id for f in result.features if f.feature_type == "pseudogene"
    }
    children = [f for f in result.features if f.parent_source_id in pseudo]
    assert children
    assert {f.feature_type for f in children} <= {"mRNA", "tRNA"}


def test_id_prefix_is_applied_to_children_too(tmp_path):
    """D-014. A prefix that renamed parents but not Parent references would
    dangle every child."""
    body = (
        row("c", "RefSeq", "gene", "1", "999", ".", "+", ".", "ID=g")
        + row("c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m;Parent=g")
    )
    result = normalised_from(gff3(tmp_path, body), P.REFSEQ, id_prefix="FUN_")
    assert {f.source_id for f in result.features} == {"FUN_g", "FUN_m"}
    assert next(f for f in result.features if f.source_id == "FUN_m").parent_source_id == "FUN_g"


def test_unimplemented_parent_inference_refuses_rather_than_falls_back(tmp_path):
    """A profile asking for synthesised inference must error, not quietly use
    explicit Parent. Structural is implemented as of milestone 6, so this now
    pins the remaining unimplemented strategy (GTF, deferred to v2.1)."""
    import dataclasses

    profile = dataclasses.replace(
        P.REFSEQ, parent_inference=P.ParentInference.SYNTHESISED, name="pretend"
    )
    body = row("c", "RefSeq", "gene", "1", "99", ".", "+", ".", "ID=g")
    with pytest.raises(NormaliseError, match="not implemented"):
        normalised_from(gff3(tmp_path, body), profile)


def test_a_computed_edge_without_a_stamped_method_is_refused(tmp_path):
    """D-032. A profile that does not state its own parents must have the rule
    recorded per feature; a silent None would become `explicit_parent` and make
    a false provenance claim about the containment algorithm."""
    import dataclasses

    profile = dataclasses.replace(P.REFSEQ, name="pretend", states_own_parents=False)
    body = row("c", "RefSeq", "gene", "1", "999", ".", "+", ".", "ID=g") + row(
        "c", "RefSeq", "mRNA", "1", "999", ".", "+", ".", "ID=m;Parent=g"
    )
    with pytest.raises(NormaliseError, match="must carry the rule"):
        normalised_from(gff3(tmp_path, body), profile)


# -- the AUGUSTUS protein oracle (D-025) ----------------------------------


def test_declared_proteins_agree_with_the_cds_spans():
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    assert len(result.proteins) == 49
    assert check_declared_proteins(result) == []


def test_the_oracle_catches_a_dropped_segment():
    """The mutation D-008 asks for, on the coordinate arithmetic. Delete one CDS
    segment and the declared protein length no longer follows."""
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    victim = next(f for f in result.features if f.feature_type == "CDS" and len(f.segments) > 1)
    victim.segments.pop()
    problems = check_declared_proteins(result)
    assert problems, "dropping a CDS segment went undetected"


def test_the_oracle_catches_an_off_by_one_span():
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    victim = next(f for f in result.features if f.feature_type == "CDS")
    victim.segments[0].end += 3
    assert check_declared_proteins(result)


def test_the_oracle_is_silent_for_dialects_without_declared_proteins():
    result = normalised_from(REFSEQ_FILE, P.REFSEQ)
    assert result.proteins == {}
    assert check_declared_proteins(result) == []


# -- end to end on the real fixtures --------------------------------------


def test_refseq_fixture_normalises_to_the_expected_shape():
    result = normalised_from(REFSEQ_FILE, P.REFSEQ)
    assert counts_by_source_type(result) == {
        "region": 1,
        "gene": 95,
        "pseudogene": 3,
        "mRNA": 11,
        "tRNA": 61,
        "rRNA": 26,
        "exon": 147,
        "CDS": 10,
    }
    assert result.unexpected_types == set()


def test_funannotate_fixture_normalises_to_the_expected_shape():
    result = normalised_from(FUNANNOTATE_FILE, P.FUNANNOTATE)
    assert counts_by_source_type(result) == {
        "gene": 46,
        "mRNA": 21,
        "tRNA": 25,
        "exon": 165,
        "CDS": 21,
    }


def test_augustus_fixture_normalises_to_the_expected_shape():
    result = normalised_from(AUGUSTUS_FILE, P.AUGUSTUS)
    assert counts_by_source_type(result) == {
        "gene": 49,
        "transcript": 49,
        "CDS": 49,
        "intron": 77,
        "start_codon": 49,
        "stop_codon": 49,
    }
    assert result.missing_types == set()


def test_refseq_region_row_carries_the_decoded_strain():
    result = normalised_from(REFSEQ_FILE, P.REFSEQ)
    region = next(f for f in result.features if f.source_type == "region")
    assert region.attributes["strain"] == ["PH-1; NRRL 31084"]
    assert region.attributes["Dbxref"] == ["taxon:229533"]


def test_a_gzipped_gff3_is_read(tmp_path):
    """NCBI serves annotation gzipped, and the uploaded F. graminearum files
    arrived named `.gff` while being gzip. Deciding by magic bytes rather than by
    name (D-035) has to apply to both tokenizers, not only GenBank's."""
    import gzip

    path = tmp_path / "compressed.gff3"
    with gzip.open(path, "wt") as handle:
        handle.write("##gff-version 3\n")
        handle.write("c\tRefSeq\tgene\t1\t99\t.\t+\t.\tID=g1\n")
    result = tokenize(path)
    assert result.data_line_count == 1
    assert result.features[0].id == "g1"
