"""Milestone 6: GenBank flat files.

Every test here corresponds to a fact measured on a real file, named in the
docstring. Where a number appears, it was counted, not assumed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from cgload.parsers import profiles as dialects
from cgload.parsers.genbank import (
    GenBankError,
    parse_location,
    tokenize_genbank,
)
from cgload.parsers.tokenizer import is_genbank, tokenize

FIXTURES = Path(__file__).parent / "fixtures"
EUKARYOTE = FIXTURES / "genbank_dscam1_CG12164.gbff"
PROKARYOTE = FIXTURES / "genbank_maeruginosa_NC_010296.gbff"


# ---------------------------------------------------------------------------
# Location grammar
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("expression", "spans", "strand"),
    [
        ("1..96", [(1, 96)], "+"),
        ("complement(27220..27435)", [(27220, 27435)], "-"),
        ("join(1..96,200..300)", [(1, 96), (200, 300)], "+"),
        ("complement(join(1..96,200..300))", [(1, 96), (200, 300)], "-"),
        ("<1..422", [(1, 422)], "+"),
        ("complement(<27220..>27435)", [(27220, 27435)], "-"),
        ("order(1..10,20..30)", [(1, 10), (20, 30)], "+"),
        ("451", [(451, 451)], "+"),
    ],
)
def test_location_grammar(expression, spans, strand):
    assert parse_location(expression, Path("x"), 1) == (spans, strand)


def test_mixed_strand_join_is_refused():
    """One strand is stored per feature, so a mixed-strand join cannot be
    represented without discarding information. Refuse rather than pick one."""
    with pytest.raises(GenBankError, match="mixed-strand"):
        parse_location("join(1..96,complement(200..300))", Path("x"), 1)


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def test_genbank_detected_by_content_not_extension(tmp_path):
    """NCBI serves .gbff, .gb, .gbk and gzipped variants; users rename them."""
    renamed = tmp_path / "no_useful_extension.txt"
    renamed.write_bytes(PROKARYOTE.read_bytes())
    assert is_genbank(renamed)
    assert dialects.detect(tokenize(renamed)).profile.name == "genbank"


def test_gzipped_genbank_is_read(tmp_path):
    """The M. aeruginosa download arrives gzipped with a .gbff name."""
    import gzip

    packed = tmp_path / "packed.gbff"
    packed.write_bytes(gzip.compress(PROKARYOTE.read_bytes()))
    assert is_genbank(packed)
    assert tokenize_genbank(packed).data_line_count > 0


def test_gff3_dialects_are_unaffected():
    """The GenBank branch must not change how the three GFF3 dialects detect."""
    for name in ("refseq", "funannotate", "augustus"):
        matches = [p for p in FIXTURES.glob("*.gff3") if name.split("_")[0] in p.name]
        for path in matches:
            assert dialects.detect(tokenize(path)).profile.name == name


# ---------------------------------------------------------------------------
# Prokaryotic linkage: locus_tag
# ---------------------------------------------------------------------------


def test_prokaryote_selects_locus_tag_strategy():
    """M. aeruginosa NIES-843 has 5,808 genes and zero mRNA records.

    The discriminator is mRNA specifically. An earlier version tested for
    'any transcript type', which selects the eukaryotic strategy for every
    prokaryote -- tRNA and rRNA are present -- leaving every CDS unlinked
    because containment needs a transcript layer that is not there.
    """
    tokenized = tokenize_genbank(PROKARYOTE)
    strategies = [c for _, c in tokenized.comments if "linkage strategy" in c]
    assert strategies == ["# cgload linkage strategy: locus_tag"]
    assert not any(f.type == "mRNA" for f in tokenized.features)


def test_prokaryote_links_every_child_to_its_gene():
    """Every non-gene record carries a locus_tag matching exactly one gene."""
    tokenized = tokenize_genbank(PROKARYOTE)
    genes = {
        f.attributes["ID"][0]
        for f in tokenized.features
        if f.type in ("gene", "pseudogene")
    }
    children = [f for f in tokenized.features if f.type not in ("gene", "pseudogene")]
    unlinked = [f for f in children if "Parent" not in f.attributes]
    assert not unlinked, f"{len(unlinked)} child rows have no parent"
    for child in children:
        assert child.attributes["Parent"][0] in genes


def test_pseudogene_cds_does_not_collide_with_its_own_gene():
    """526 of 5,808 M. aeruginosa loci are pseudogenes whose CDS has no
    protein_id. Identified by locus_tag alone, such a CDS would take its own
    gene's identifier and the two rows would merge into one feature (D-030)."""
    tokenized = tokenize_genbank(PROKARYOTE)
    pseudo_cds = [
        f
        for f in tokenized.features
        if f.type == "CDS" and "pseudo" in f.attributes
    ]
    assert pseudo_cds, "fixture carries no pseudogene CDS"
    for cds in pseudo_cds:
        identifier = cds.attributes["ID"][0]
        parent = cds.attributes["Parent"][0]
        assert identifier != parent
        assert identifier.endswith(":CDS")


def test_pseudogene_rows_are_flagged_for_the_completeness_check():
    """A pseudogene CDS has no /translation, so any check computing protein
    length from CDS spans must skip it rather than report a mismatch."""
    tokenized = tokenize_genbank(PROKARYOTE)
    flagged = [f for f in tokenized.features if "pseudo" in f.attributes]
    assert flagged
    assert all("translation" not in f.attributes for f in flagged)


def test_multi_segment_cds_becomes_one_feature_with_several_segments():
    """join() spans a programmed frameshift; 37 CDSs in the full genome."""
    tokenized = tokenize_genbank(PROKARYOTE)
    by_id: dict[str, int] = {}
    for feature in tokenized.features:
        if feature.type == "CDS":
            by_id[feature.attributes["ID"][0]] = by_id.get(feature.attributes["ID"][0], 0) + 1
    assert any(count > 1 for count in by_id.values()), "fixture carries no join() CDS"


# ---------------------------------------------------------------------------
# Eukaryotic linkage: containment, then the declared note, then refusal
# ---------------------------------------------------------------------------


def test_eukaryote_selects_containment_strategy():
    tokenized = tokenize_genbank(EUKARYOTE)
    strategies = [c for _, c in tokenized.comments if "linkage strategy" in c]
    assert strategies == ["# cgload linkage strategy: containment"]


def test_every_cds_resolves_to_exactly_one_transcript():
    tokenized = tokenize_genbank(EUKARYOTE)
    coding = [f for f in tokenized.features if f.type == "CDS"]
    assert coding
    assert all("Parent" in f.attributes for f in coding)


def test_cds_pairing_matches_the_files_own_ground_truth():
    """The CDS /note names its transcript: 'from transcript CG12164-RA'.

    Measured on the full Dscam1 region (79 transcripts, 75 isoforms on one
    gene): containment alone resolves 75 of 79 uniquely; the note resolves the
    remaining 4.
    """
    text = EUKARYOTE.read_text()
    table = text.split("FEATURES")[1].split("\nORIGIN")[0]

    truth: dict[str, str] = {}
    for record in re.split(r"\n(?=     \S)", table):
        if not re.match(r"\s+CDS\s", record):
            continue
        protein = re.search(r'/protein_id="([^"]+)"', record)
        note = re.search(r'/note="([^"]*)"', record, re.S)
        if not (protein and note):
            continue
        named = re.search(r"from transcript (\S+?);", re.sub(r"\s+", " ", note.group(1)))
        if named:
            truth[protein.group(1)] = named.group(1)
    assert truth, "fixture lost its /note ground truth"

    tokenized = tokenize_genbank(EUKARYOTE)
    transcript_note = {
        f.attributes["ID"][0]: " ".join(f.attributes.get("note", [])).split(";")[0].strip()
        for f in tokenized.features
        if f.type == "mRNA"
    }

    def suffix(text_: str) -> str | None:
        found = re.search(r"-([A-Z]+)\b", text_ or "")
        return found.group(1) if found else None

    checked = 0
    for feature in tokenized.features:
        if feature.type != "CDS":
            continue
        identifier = feature.attributes["ID"][0]
        if identifier not in truth:
            continue
        parent = feature.attributes["Parent"][0]
        assert suffix(transcript_note[parent]) == suffix(truth[identifier]), (
            f"{identifier} paired with {parent}, but the file says "
            f"{truth[identifier]}"
        )
        checked += 1
    assert checked, "no CDS could be checked against ground truth"


def test_positional_pairing_is_not_used():
    """The specification originally said to link the nth CDS to the nth
    transcript, 'which GenBank guarantees'. It does not: transcripts carry UTRs
    and therefore sort differently from their own coding sequences. Measured on
    Dscam1, positional pairing is wrong for 2 of 79 -- silently.

    CG12164 is the counterexample and is in this fixture: its mRNAs appear in
    the order RB, RA while its CDSs appear RA, RB.
    """
    tokenized = tokenize_genbank(EUKARYOTE)
    order: list[tuple[str, str]] = []
    seen: set[str] = set()
    for feature in tokenized.features:
        if feature.type not in ("mRNA", "CDS"):
            continue
        identifier = feature.attributes["ID"][0]
        if identifier in seen:
            continue
        seen.add(identifier)
        gene = feature.attributes.get("gene", [""])[0]
        if gene == "CG12164":
            order.append((feature.type, identifier))

    transcripts = [i for t, i in order if t == "mRNA"]
    coding = [i for t, i in order if t == "CDS"]
    assert len(transcripts) == len(coding) == 2

    tokenized_by_id = {
        f.attributes["ID"][0]: f for f in tokenized.features if "ID" in f.attributes
    }
    positional = dict(zip(coding, transcripts, strict=True))
    actual = {c: tokenized_by_id[c].attributes["Parent"][0] for c in coding}
    assert actual != positional, (
        "linkage agrees with file order for CG12164, which is the case known to "
        "be wrong -- positional pairing has been reintroduced"
    )


def test_ambiguity_without_a_tiebreak_is_refused(tmp_path):
    """Mutation test. Strip the /note qualifiers and the two CG12164 isoforms
    become indistinguishable: both transcripts contain both CDSs. The loader
    must name the gene and both candidates rather than choose."""
    stripped = tmp_path / "no_notes.gbff"
    stripped.write_text(
        re.sub(r'\n\s{21}/note="[^"]*"', "", EUKARYOTE.read_text(), flags=re.S)
    )
    with pytest.raises(GenBankError) as excinfo:
        tokenize_genbank(stripped)
    message = str(excinfo.value)
    assert "CG12164" in message
    assert "will not guess" in message


def test_identifiers_come_from_the_file():
    """No identifier is invented. Every ID is a protein_id, transcript_id, or a
    locus_tag (with the record type appended where the tag alone would collide).
    """
    for path in (EUKARYOTE, PROKARYOTE):
        text = path.read_text()
        for feature in tokenize_genbank(path).features:
            identifier = feature.attributes.get("ID", [None])[0]
            if identifier is None:
                continue
            base = identifier.split(":")[0]
            assert base in text, f"{identifier!r} does not appear in {path.name}"


def test_partial_features_are_flagged_not_silently_trimmed():
    """'<27220..27435' means the feature runs off the assembled sequence. The
    coordinate is real; it is the completeness claim that changes."""
    tokenized = tokenize_genbank(PROKARYOTE)
    partial = [f for f in tokenized.features if "partial" in f.attributes]
    assert partial, "fixture carries no partial feature"
    assert all(f.start > 0 and f.end >= f.start for f in partial)


# =========================================================================
# Milestone 6b: linkage provenance, slippage, and the end-to-end load.
# =========================================================================

import sqlalchemy as _sa  # noqa: E402
from click.testing import CliRunner  # noqa: E402
from cgload.cli.main import cli  # noqa: E402
from cgload.config import load_config  # noqa: E402
from cgload.db.engine import build_engine  # noqa: E402
from cgload.db.init import create_schema  # noqa: E402
from cgload.db.schema import feature as _feature  # noqa: E402
from cgload.db.vocabulary import COMPUTED_LINKAGE, LINKAGE_METHODS  # noqa: E402
from cgload.loader import LoadError, load_organism  # noqa: E402
from cgload.parsers import profiles as _P  # noqa: E402
from cgload.parsers.normalise import normalise as _normalise  # noqa: E402


def _normalised(path):
    tokenized = tokenize_genbank(path)
    detection = _P.detect(tokenized)
    return _normalise(
        tokenized.features,
        detection.profile,
        comments=tokenized.comments,
        unexpected=_P.unexpected_types(detection.profile, tokenized.features),
        missing=_P.missing_types(detection.profile, tokenized.features),
    )


# -- linkage provenance (D-032) -------------------------------------------


def test_a_genbank_feature_is_never_recorded_as_explicit_parent():
    """The regression guard for the provenance claim.

    GenBank carries no parent pointers, so `explicit_parent` would be a false
    claim about provenance -- and it is the claim the containment algorithm is
    published on. This is the kind of bug that shows up in the paper rather than
    in a stack trace, so it gets its own test.
    """
    for path in (PROKARYOTE, EUKARYOTE):
        methods = {
            item.linkage_method for item in _normalised(path).features
            if item.parent_source_id is not None
        }
        assert "explicit_parent" not in methods, path.name
        assert methods <= set(LINKAGE_METHODS)


def test_prokaryotic_edges_are_recorded_as_locus_tag():
    """Bacteria: the qualifier is shared and unique, so nothing is computed."""
    linked = [
        item for item in _normalised(PROKARYOTE).features
        if item.parent_source_id is not None
    ]
    assert linked
    assert {item.linkage_method for item in linked} == {"locus_tag"}


def test_eukaryotic_cds_edges_are_structural_or_note_never_locus_tag():
    """The fly: CDS-to-transcript is the ambiguous edge and must record which of
    the two rules resolved it. Transcript-to-gene is a shared qualifier and is
    correctly `locus_tag`, so the assertion is scoped to CDS."""
    features = _normalised(EUKARYOTE).features
    coding = [
        item for item in features
        if item.feature_type == "CDS" and item.parent_source_id is not None
    ]
    assert coding
    methods = {item.linkage_method for item in coding}
    assert methods <= {"structural", "note"}
    assert "structural" in methods, "containment resolved nothing; fixture changed"
    assert "note" in methods, "the ambiguous pair is gone; fixture no longer exercises it"


def test_computed_and_stated_linkage_are_distinguishable():
    """`structural` must be classed as computed, or a reader cannot tell an
    edge cgload inferred from one the file stated."""
    assert "structural" in COMPUTED_LINKAGE
    assert "explicit_parent" not in COMPUTED_LINKAGE
    assert "locus_tag" not in COMPUTED_LINKAGE


def test_the_linkage_attribute_is_not_also_stored_as_an_attribute():
    """One source of truth. A second copy in `feature_attribute` could disagree
    with the column."""
    from cgload.parsers.profiles import LINKAGE_ATTRIBUTE

    for path in (PROKARYOTE, EUKARYOTE):
        for item in _normalised(path).features:
            assert LINKAGE_ATTRIBUTE not in item.attributes


# -- ribosomal slippage (D-033) -------------------------------------------


def test_a_ribosomal_slippage_cds_loads_with_overlapping_segments():
    """join(57926..58406,58406..59148) reads base 58406 twice. That is what a
    programmed frameshift is, and the GFF3 non-overlap rule would refuse a valid
    record."""
    overlapping = []
    for item in _normalised(PROKARYOTE).features:
        spans = [(segment.start, segment.end) for segment in item.segments]
        if any(a[1] >= b[0] for a, b in zip(spans, spans[1:], strict=False)):
            overlapping.append(item)
    assert overlapping, "the fixture no longer contains a slippage record"
    for item in overlapping:
        assert item.attributes.get("ribosomal_slippage") is not None


def test_an_unexplained_overlap_is_still_an_error(tmp_path):
    """The permission is gated on the qualifier, not on the format. Without
    /ribosomal_slippage an overlap in a GenBank file is a duplicate."""
    import dataclasses

    from cgload.parsers.normalise import NormaliseError
    from cgload.parsers.tokenizer import RawFeature

    profile = dataclasses.replace(_P.GENBANK, overlap_qualifier="ribosomal_slippage")
    rows = [
        RawFeature(
            line_number=n,
            seqid="c",
            source="GenBank",
            type="CDS",
            start=start,
            end=end,
            score=None,
            strand="+",
            phase=0,
            attributes={"ID": ["shared"]},
        )
        for n, (start, end) in enumerate([(100, 300), (250, 500)], start=1)
    ]
    with pytest.raises(NormaliseError, match="overlapping"):
        _normalise(rows, profile)


# -- the load path (the gap milestone 6 left) ------------------------------


def _fasta_for(path, tmp_path):
    """A FASTA long enough for the record's coordinates, named for its
    accession."""
    import random
    import re

    text = path.read_text()
    version = re.search(r"^VERSION\s+(\S+)", text, re.M)
    accession = re.search(r"^ACCESSION\s+(\S+)", text, re.M)
    seqid = (version or accession).group(1)
    coordinates = [
        int(value)
        for line in text.splitlines()
        if line.startswith("     ") and ".." in line
        for value in re.findall(r"\b(\d{3,})\b", line)
    ]
    length = max(coordinates) + 1000
    random.seed(5)
    target = tmp_path / f"{seqid}.fa"
    with target.open("w") as handle:
        handle.write(f">{seqid}\n")
        sequence = "".join(random.choice("ACGT") for _ in range(length))
        for index in range(0, length, 60):
            handle.write(sequence[index : index + 60] + "\n")
    return target


def _config_for(path, tmp_path, name):
    fasta = _fasta_for(path, tmp_path)
    annotation = tmp_path / path.name
    if not annotation.exists():
        annotation.write_bytes(path.read_bytes())
    config = tmp_path / f"{name}.toml"
    config.write_text(
        f'[organism]\ngenus = "Testus"\nspecies = "{name}"\nstrain = "S"\n'
        f'ncbi_taxon_id = 1\n\n[assembly]\nname = "{name}"\n'
        f'fasta = "{fasta.name}"\n\n[[files.annotation]]\n'
        f'path = "{annotation.name}"\n'
    )
    return config


@pytest.fixture
def engine(tmp_path):
    built = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(built)
    return built


@pytest.mark.parametrize(
    ("fixture", "name", "expected"),
    [
        (PROKARYOTE, "prokaryote",
         {"gene": 28, "CDS": 21, "tRNA": 2, "rRNA": 2, "ncRNA": 2, "tmRNA": 1}),
        (EUKARYOTE, "eukaryote", {"gene": 2, "mRNA": 5, "CDS": 5}),
    ],
)
def test_a_genbank_file_loads_end_to_end(engine, tmp_path, fixture, name, expected):
    """Never run before 6b: the format is dispatched by content, parsed, linked
    and stored, in one transaction with the assembly."""
    config = load_config(_config_for(fixture, tmp_path, name))
    report = load_organism(engine, config, tmp_path / "data")
    assert report.files[0].dialect == "genbank"
    assert report.files[0].counts_by_type == expected


def test_the_stored_linkage_method_survives_to_the_database(engine, tmp_path):
    """The provenance must be queryable, not just computed: `explain` reads this
    column, and so will the paper's accuracy table."""
    config = load_config(_config_for(EUKARYOTE, tmp_path, "eukaryote"))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        rows = dict(
            connection.execute(
                _sa.select(_feature.c.linkage_method, _sa.func.count())
                .where(_feature.c.source_type == "CDS")
                .group_by(_feature.c.linkage_method)
            ).all()
        )
    assert set(rows) <= {"structural", "note"}
    assert sum(rows.values()) > 0


def test_the_slippage_cds_reaches_the_database_as_two_segments(engine, tmp_path):
    config = load_config(_config_for(PROKARYOTE, tmp_path, "prokaryote"))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        rows = connection.execute(
            _sa.select(_feature.c.source_id, _feature.c.start, _feature.c.end)
            .where(_feature.c.source_id == "WP_148204717.1")
            .order_by(_feature.c.segment_index)
        ).all()
    assert [(row.start, row.end) for row in rows] == [(57926, 58406), (58406, 59148)]


def test_a_file_named_gbff_that_is_not_genbank_is_refused(engine, tmp_path):
    """The name lies both ways (D-035). A GFF3 renamed .gbff must not be read as
    GFF3 by accident, because the dispatch is by content."""
    liar = tmp_path / "liar.gbff"
    liar.write_text("##gff-version 3\nc\tRefSeq\tgene\t1\t99\t.\t+\t.\tID=g\n")
    config = _config_for(EUKARYOTE, tmp_path, "eukaryote")
    config.write_text(config.read_text().replace(EUKARYOTE.name, liar.name))
    with pytest.raises(LoadError, match="LOCUS"):
        load_organism(engine, load_config(config), tmp_path / "data")


def test_the_load_command_reports_the_genbank_dialect(engine, tmp_path):
    config = _config_for(PROKARYOTE, tmp_path, "prokaryote")
    result = CliRunner().invoke(
        cli,
        ["load", "--config", str(config), "--data-dir", str(tmp_path / "data"),
         "--db-url", str(engine.url)],
    )
    assert result.exit_code == 0, result.output
    assert "genbank" in result.output


def test_no_spurious_missing_type_warning_for_genbank(engine, tmp_path):
    """One profile serves both kingdoms, so `provides` is a union: a bacterial
    record legitimately has no mRNA. Warning on that trains users to ignore
    warnings, which costs more than it gains."""
    config = load_config(_config_for(PROKARYOTE, tmp_path, "prokaryote"))
    report = load_organism(engine, config, tmp_path / "data")
    assert not [w for w in report.warnings if "expects" in w]


# -- the shared non-conformance guard -------------------------------------


def test_pseudogenes_and_partials_share_one_exemption_predicate():
    """Three unrelated things produce a CDS that will not divide by three.
    Milestone 8's completeness check must consult the same predicate as the
    AUGUSTUS oracle, or the two will disagree about which features are exempt."""
    from cgload.conformance import exemption_for, is_exempt

    assert is_exempt({"pseudo": ["true"]})
    assert is_exempt({"pseudo": []})  # GenBank /pseudo carries no value
    assert is_exempt({"partial": ["true"]})
    assert is_exempt({"artificial_location": ["heterogeneous population"]})
    assert is_exempt({"ribosomal_slippage": []})
    assert not is_exempt({"product": ["hypothetical protein"]})
    assert not is_exempt({"partial": ["false"]}), "an explicit false is not an exemption"
    assert exemption_for({"pseudo": ["true"]}).flag == "pseudo"


def test_the_prokaryotic_fixture_exercises_the_exemption():
    """4 of its 28 loci are pseudogenes, so this is not a hypothetical path."""
    from cgload.conformance import is_exempt

    exempt = [
        item for item in _normalised(PROKARYOTE).features if is_exempt(item.attributes)
    ]
    assert len(exempt) >= 4
