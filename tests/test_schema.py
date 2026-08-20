"""Milestone 1 tests.

What these assert, stated before they were written (working agreement §7):

* The schema can be created and the vocabulary is seeded exactly once.
* The two type columns are genuinely two columns, and nothing in the schema
  lets `verify` read the normalised one by accident.
* A discontinuous feature -- one source_id, several segments -- stores as
  several rows and still counts as one feature.
* Two features that differ only by a normalised ID cannot silently merge.
* The vocabulary hash changes when linkage semantics change, and a database
  seeded under one hash refuses a later one.
* init does not overwrite an existing database, and never accepts a password
  in an argument.
"""

from __future__ import annotations

import pytest
from cgload.db import vocabulary as vocab
from cgload.db.engine import CredentialInArgumentError, build_engine, validated_url
from cgload.db.init import (
    SchemaAlreadyPresentError,
    VocabularyMismatchError,
    assert_vocabulary_current,
    create_schema,
)
from cgload.db.schema import (
    assembly,
    feature,
    metadata,
    organism,
    sequence_region,
    source_file,
    strain,
    vocabulary,
    vocabulary_meta,
)
from sqlalchemy import func, insert, select, text, update
from sqlalchemy.exc import IntegrityError


@pytest.fixture
def engine(tmp_path):
    eng = build_engine(f"sqlite:///{tmp_path / 'test.db'}")
    create_schema(eng)
    return eng


@pytest.fixture
def region(engine):
    """A minimal organism/strain/assembly/region/file chain."""
    with engine.begin() as conn:
        org = conn.execute(
            insert(organism).values(genus="Phytophthora", species="megasperma")
        ).inserted_primary_key[0]
        st = conn.execute(
            insert(strain).values(organism_id=org, name="reference")
        ).inserted_primary_key[0]
        asm = conn.execute(
            insert(assembly).values(
                strain_id=st,
                name="CJ26",
                vocabulary_version=vocab.VOCABULARY_VERSION,
                vocabulary_hash=vocab.content_hash(),
                cgload_version="2.0.0.dev0",
            )
        ).inserted_primary_key[0]
        sf = conn.execute(
            insert(source_file).values(
                assembly_id=asm, path="a.gff3", sha256="0" * 64, file_format="gff3"
            )
        ).inserted_primary_key[0]
        reg = conn.execute(
            insert(sequence_region).values(assembly_id=asm, source_id="scaffold_12", length=500000)
        ).inserted_primary_key[0]
    return {"assembly_id": asm, "sequence_region_id": reg, "source_file_id": sf}


def _seg(region, **kw):
    base = {
        "assembly_id": region["assembly_id"],
        "sequence_region_id": region["sequence_region_id"],
        "source_file_id": region["source_file_id"],
        "strand": "+",
    }
    return {**base, **kw}


# -- schema shape ----------------------------------------------------------


def test_init_creates_every_table(engine):
    from sqlalchemy import inspect

    assert set(inspect(engine).get_table_names()) == set(metadata.tables)


def test_vocabulary_meta_is_a_singleton(engine):
    with engine.begin() as conn, pytest.raises(IntegrityError):
        conn.execute(insert(vocabulary_meta).values(vocabulary_meta_id=2, version=1,
                                                    content_hash="x" * 64))


def test_source_type_and_feature_type_are_independent_columns(engine, region):
    """The whole point of D-017. An AUGUSTUS transcript stores both terms."""
    with engine.begin() as conn:
        conn.execute(insert(feature), _seg(region, source_id="g1.t1",
                                           source_type="transcript", feature_type="mRNA",
                                           start=100, end=900))
        row = conn.execute(
            select(feature.c.source_type, feature.c.feature_type)
        ).one()
    assert (row.source_type, row.feature_type) == ("transcript", "mRNA")


def test_verify_counts_source_type_not_feature_type(engine, region):
    """Two source terms collapsing to one canonical term must remain two on
    the counting side, or the primary completeness check desynchronises."""
    with engine.begin() as conn:
        conn.execute(insert(feature), [
            _seg(region, source_id="a", source_type="transcript", feature_type="mRNA",
                 start=1, end=10),
            _seg(region, source_id="b", source_type="mRNA", feature_type="mRNA",
                 start=20, end=30),
        ])
        by_source = dict(conn.execute(
            select(feature.c.source_type, func.count(func.distinct(feature.c.source_id)))
            .group_by(feature.c.source_type)
        ).all())
    assert by_source == {"transcript": 1, "mRNA": 1}


# -- discontinuous features ------------------------------------------------


def test_multi_segment_feature_is_many_rows_but_one_feature(engine, region):
    """RefSeq repeats one ID across every CDS segment. This must load."""
    spans = [(100, 200), (300, 400), (500, 600)]
    with engine.begin() as conn:
        conn.execute(insert(feature), [
            _seg(region, source_id="cds-XP_1.1", source_type="CDS", feature_type="CDS",
                 start=s, end=e, phase=0, segment_index=i)
            for i, (s, e) in enumerate(spans)
        ])
        rows = conn.execute(select(func.count()).select_from(feature)).scalar_one()
        features = conn.execute(
            select(func.count(func.distinct(feature.c.source_id)))
        ).scalar_one()
    assert (rows, features) == (3, 1)


def test_identical_segment_is_rejected(engine, region):
    """A malformed file repeating one segment must not double-count it."""
    with engine.begin() as conn:
        conn.execute(insert(feature), _seg(region, source_id="cds-1", source_type="CDS",
                                           feature_type="CDS", start=100, end=200))
        with pytest.raises(IntegrityError):
            conn.execute(insert(feature), _seg(region, source_id="cds-1", source_type="CDS",
                                               feature_type="CDS", start=100, end=200))


def test_same_id_on_two_regions_is_allowed(engine, region):
    """The uniqueness constraint is per region, not per assembly: IDs are only
    required to be unique within a sequence region."""
    with engine.begin() as conn:
        other = conn.execute(
            insert(sequence_region).values(
                assembly_id=region["assembly_id"], source_id="scaffold_13", length=1000
            )
        ).inserted_primary_key[0]
        conn.execute(insert(feature), _seg(region, source_id="x", source_type="gene",
                                           feature_type="gene", start=1, end=9))
        conn.execute(insert(feature), {
            **_seg(region, source_id="x", source_type="gene", feature_type="gene",
                   start=1, end=9),
            "sequence_region_id": other,
        })
        assert conn.execute(select(func.count()).select_from(feature)).scalar_one() == 2


# -- constraints that catch malformed input --------------------------------


@pytest.mark.parametrize(
    "bad",
    [
        {"start": 0, "end": 10},          # coordinates are 1-based
        {"start": 500, "end": 100},       # end before start
        {"strand": "x"},                  # column 7 vocabulary
        {"phase": 3},                     # column 8 vocabulary
        {"segment_index": -1},
        {"linkage_method": "fuzzy_product_name"},  # name-similarity matching, barred
    ],
)
def test_malformed_values_are_rejected(engine, region, bad):
    values = _seg(region, source_id="f", source_type="CDS", feature_type="CDS",
                  start=100, end=200)
    values.update(bad)
    with engine.begin() as conn, pytest.raises(IntegrityError):
        conn.execute(insert(feature), values)


def test_foreign_keys_are_enforced_on_sqlite(engine, region):
    """SQLite disables FK enforcement by default; build_engine turns it on."""
    with engine.begin() as conn, pytest.raises(IntegrityError):
        conn.execute(insert(feature), _seg(region, source_id="orphan", source_type="CDS",
                                           feature_type="CDS", start=1, end=2,
                                           parent_id=999999))


# -- vocabulary versioning -------------------------------------------------


def test_hash_covers_hierarchy_rules_not_just_terms(monkeypatch):
    """Editing ALLOWED_PARENTS changes linkage semantics, so it must change
    the hash even though no term was added."""
    before = vocab.content_hash()
    monkeypatch.setitem(vocab.ALLOWED_PARENTS, "CDS", frozenset({"mRNA", "gene"}))
    assert vocab.content_hash() != before


def test_hash_covers_linkage_methods(monkeypatch):
    before = vocab.content_hash()
    monkeypatch.setattr(vocab, "LINKAGE_METHODS", (*vocab.LINKAGE_METHODS, "guessed"))
    assert vocab.content_hash() != before


def test_stored_hash_guards_later_commands(engine, monkeypatch):
    assert_vocabulary_current(engine)  # seeded database agrees with itself
    monkeypatch.setitem(vocab.SYNONYMS, "novel_transcript_term", "mRNA")
    with pytest.raises(VocabularyMismatchError, match="remap"):
        assert_vocabulary_current(engine)


def test_every_synonym_resolves_into_tier_a():
    assert set(vocab.SYNONYMS.values()) <= set(vocab.HIERARCHY_TERMS)


def test_no_synonym_is_also_a_canonical_term():
    """A term cannot be both, or its canonical_term is ambiguous at seed time."""
    assert not set(vocab.SYNONYMS) & set(vocab.HIERARCHY_TERMS)


def test_allowed_parents_covers_exactly_tier_a():
    assert set(vocab.ALLOWED_PARENTS) == set(vocab.HIERARCHY_TERMS)


def test_allowed_parents_only_names_tier_a_terms():
    """A Tier A child may not hang off a term the loader has no rules for."""
    for parents in vocab.ALLOWED_PARENTS.values():
        assert parents <= set(vocab.HIERARCHY_TERMS)


def test_seeded_rows_match_the_module(engine):
    with engine.connect() as conn:
        stored = {r.term: r.canonical_term for r in conn.execute(select(vocabulary)).all()}
    assert stored == {r["term"]: r["canonical_term"] for r in vocab.seed_rows()}


# -- init behaviour --------------------------------------------------------


def test_init_refuses_an_existing_database(engine):
    with pytest.raises(SchemaAlreadyPresentError):
        create_schema(engine)


def test_init_force_recreates(engine):
    with engine.begin() as conn:
        conn.execute(insert(organism).values(genus="G", species="s"))
    create_schema(engine, force=True)
    with engine.connect() as conn:
        assert conn.execute(select(func.count()).select_from(organism)).scalar_one() == 0


def test_password_in_url_is_refused():
    with pytest.raises(CredentialInArgumentError):
        validated_url("postgresql://user:secret@host/db")


def test_password_from_environment_is_allowed():
    url = validated_url("postgresql://user:secret@host/db", from_environment=True)
    assert "secret" not in url.render_as_string(hide_password=True)


def test_sqlite_needs_no_credentials(tmp_path):
    build_engine(f"sqlite:///{tmp_path / 'x.db'}")


def test_no_rtree_virtual_table(engine):
    """D-002: one code path across SQLite, MySQL and PostgreSQL. RTREE has no
    equivalent on the servers, so its presence would void that claim."""
    with engine.connect() as conn:
        sql = conn.execute(text("SELECT sql FROM sqlite_master WHERE sql IS NOT NULL")).scalars()
        assert not any("VIRTUAL TABLE" in s.upper() for s in sql)


def test_schema_is_idempotent_under_repeated_create_all(engine):
    metadata.create_all(engine)  # must not raise


def test_update_cannot_smuggle_a_bad_strand(engine, region):
    with engine.begin() as conn:
        conn.execute(insert(feature), _seg(region, source_id="f", source_type="gene",
                                           feature_type="gene", start=1, end=2))
        with pytest.raises(IntegrityError):
            conn.execute(update(feature).values(strand="x"))


# -- regression: codon rows must not inflate CDS counts ---------------------


def test_codon_terms_are_not_mapped_to_cds():
    """AUGUSTUS and funannotate include start and stop codons inside their
    terminal CDS rows. Mapping the codon rows to canonical CDS would count the
    same bases twice, so `query --type CDS` would return more CDS features
    than the gene has."""
    for term in ("start_codon", "stop_codon"):
        assert term not in vocab.SYNONYMS
        assert term not in vocab.HIERARCHY_TERMS
        assert term in vocab.NOT_CANONICALISED


def test_augustus_gene_counts_are_not_inflated_by_codon_rows(engine, region):
    """The real AUGUSTUS shape: one transcript, four CDS segments sharing an
    ID, plus a start_codon and a stop_codon row."""
    rows = [
        _seg(region, source_id="g1", source_type="gene", feature_type="gene",
             start=644, end=1840, strand="-"),
        _seg(region, source_id="g1.t1", source_type="transcript", feature_type="mRNA",
             start=644, end=1840, strand="-"),
    ]
    spans = [(644, 739), (801, 1348), (1395, 1601), (1651, 1840)]
    rows += [
        _seg(region, source_id="g1.t1.cds", source_type="CDS", feature_type="CDS",
             start=s, end=e, strand="-", segment_index=i)
        for i, (s, e) in enumerate(spans)
    ]
    rows += [
        _seg(region, source_id="g1.t1.stop_codon.1", source_type="stop_codon",
             feature_type="stop_codon", start=644, end=646, strand="-"),
        _seg(region, source_id="g1.t1.start_codon.1", source_type="start_codon",
             feature_type="start_codon", start=1838, end=1840, strand="-"),
    ]
    with engine.begin() as conn:
        conn.execute(insert(feature), rows)
        canonical = dict(conn.execute(
            select(feature.c.feature_type, func.count(func.distinct(feature.c.source_id)))
            .group_by(feature.c.feature_type)
        ).all())
        source = dict(conn.execute(
            select(feature.c.source_type, func.count(func.distinct(feature.c.source_id)))
            .group_by(feature.c.source_type)
        ).all())

    assert canonical["CDS"] == 1, "codon rows must not be counted as CDS features"
    assert canonical["mRNA"] == 1
    # And the counting side still reports the file's own words, one each.
    assert source == {"gene": 1, "transcript": 1, "CDS": 1,
                      "start_codon": 1, "stop_codon": 1}


def test_no_canonical_term_is_also_marked_uncanonicalised():
    assert not set(vocab.NOT_CANONICALISED) & set(vocab.HIERARCHY_TERMS)
    assert not set(vocab.NOT_CANONICALISED) & set(vocab.SYNONYMS)


# -- anonymous rows (AUGUSTUS intron / start_codon / stop_codon) ------------


def test_a_row_with_no_id_is_stored_and_counted_as_its_own_feature(engine, region):
    """AUGUSTUS gives intron, start_codon and stop_codon rows a Parent and no
    ID. Each is its own feature; collapsing them by a shared NULL would report
    one where the file has many, silently, in the headline check."""
    from cgload.db.counts import feature_count_expression

    with engine.begin() as conn:
        conn.execute(insert(feature), [
            _seg(region, source_id=None, source_type="stop_codon",
                 feature_type="stop_codon", start=s, end=s + 2, strand="-")
            for s in (644, 2820, 7664)
        ])
        counted = conn.execute(
            select(feature_count_expression()).where(
                feature.c.assembly_id == region["assembly_id"]
            )
        ).scalar_one()
    assert counted == 3


def test_named_and_anonymous_rows_count_together(engine, region):
    """One four-segment CDS (one feature) plus two anonymous introns (two)."""
    from cgload.db.counts import feature_count_expression

    rows = [
        _seg(region, source_id="g1.t1.cds", source_type="CDS", feature_type="CDS",
             start=s, end=e, segment_index=i)
        for i, (s, e) in enumerate([(644, 739), (801, 1348), (1395, 1601), (1651, 1840)])
    ]
    # segment_index given explicitly: a batch insert whose dicts do not all
    # carry the same keys is rejected by SQLAlchemy rather than falling back to
    # the server default, so the loader must emit uniform rows. Worth pinning.
    rows += [
        _seg(region, source_id=None, source_type="intron", feature_type="intron",
             start=740, end=800, segment_index=0),
        _seg(region, source_id=None, source_type="intron", feature_type="intron",
             start=1349, end=1394, segment_index=0),
    ]
    with engine.begin() as conn:
        conn.execute(insert(feature), rows)
        counted = conn.execute(
            select(feature_count_expression()).where(
                feature.c.assembly_id == region["assembly_id"]
            )
        ).scalar_one()
    assert counted == 3
