"""Milestone 7 tests: eggNOG and InterProScan.

What these assert, stated before they were written:

* Both formats are recognised by content and parsed into uniform hits.
* The multi-valued cells are split correctly, including COG's separator-free
  letter run, which a comma split gets wrong.
* Every query protein must match a feature by default, and a file that
  half-matches is refused with a count rather than loaded quietly.
* A query matching two features is refused outright: attaching a protein's
  annotation to the wrong isoform is invisible and unrecoverable.
* A repeated domain at two positions survives as two rows; an exact duplicate
  collapses to one.
* Loading functional annotation does not change any feature count, because it
  attaches to the graph rather than becoming part of it.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from cgload.config import load_config
from cgload.db.counts import feature_count_expression
from cgload.db.engine import build_engine
from cgload.db.init import create_schema
from cgload.db.schema import annotation_run, feature, protein_annotation
from cgload.functional import FunctionalLoadError, load_functional
from cgload.loader import load_organism
from cgload.parsers.functional import (
    FunctionalParseError,
    detect_functional_format,
    parse_functional,
    parse_interproscan,
)
from sqlalchemy import func, select

FIXTURES = Path(__file__).parent / "fixtures"
REFSEQ_GFF = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
EGGNOG = FIXTURES / "eggnog_fgraminearum.emapper.annotations"
INTERPRO = FIXTURES / "interproscan_fgraminearum.tsv"


# -- parsing ---------------------------------------------------------------


def test_both_formats_are_recognised_by_content():
    assert detect_functional_format(EGGNOG) == "eggnog"
    assert detect_functional_format(INTERPRO) == "interproscan"


def test_an_unrecognised_file_is_refused(tmp_path):
    path = tmp_path / "mystery.tsv"
    path.write_text("some\tshort\tline\n")
    with pytest.raises(FunctionalParseError, match="not recognised"):
        detect_functional_format(path)


def test_eggnog_version_is_captured():
    assert parse_functional(EGGNOG).tool_version == "emapper-2.1.12"


def test_eggnog_columns_are_read_from_the_header_not_by_position(tmp_path):
    """Column order has changed between eggNOG versions, so a positional reader
    silently mis-assigns every field after the change."""
    path = tmp_path / "reordered.emapper.annotations"
    path.write_text(
        "#query\tPFAMs\tDescription\tevalue\n"
        "XP_1.1\tPF00067\tCytochrome P450\t1e-40\n"
    )
    hits = parse_functional(path).hits
    pfam = [hit for hit in hits if hit.analysis == "Pfam"]
    assert [hit.accession for hit in pfam] == ["PF00067"]


def test_eggnog_data_before_the_header_is_refused(tmp_path):
    path = tmp_path / "headerless.emapper.annotations"
    path.write_text("XP_1.1\tPF00067\n#query\tPFAMs\n")
    with pytest.raises(FunctionalParseError, match="before the '#query' header"):
        parse_functional(path, tool="eggnog")


def test_cog_categories_are_letters_not_a_comma_list():
    """`EGP` is three COG categories with no separator. Splitting on commas
    yields one accession named 'EGP', and a query for COG 'E' then misses it."""
    hits = parse_functional(EGGNOG).hits
    cog = {hit.accession for hit in hits if hit.analysis == "COG"}
    assert cog, "the fixture no longer has COG categories"
    assert all(len(accession) == 1 for accession in cog)
    assert {"E", "G", "P"} <= cog


def test_placeholder_values_do_not_become_accessions():
    """eggNOG writes '-' for absent. A row for a term named '-' would appear on
    thousands of proteins and mean nothing."""
    for hit in parse_functional(EGGNOG).hits:
        assert hit.accession not in ("-", "", "NA", ".")


def test_multi_valued_cells_become_several_hits():
    hits = parse_functional(EGGNOG).hits
    go = [hit for hit in hits if hit.analysis == "GO"]
    assert len(go) > len({hit.query_id for hit in go})


def test_interproscan_column_count_is_validated(tmp_path):
    path = tmp_path / "short.tsv"
    path.write_text("XP_1.1\tmd5\t100\tPfam\tPF00067\n")
    with pytest.raises(FunctionalParseError, match="at least 11"):
        parse_interproscan(path)


def test_interproscan_reversed_span_is_refused(tmp_path):
    path = tmp_path / "reversed.tsv"
    path.write_text(
        "XP_1.1\tmd5\t100\tPfam\tPF00067\tCytochrome\t90\t20\t1.0\tT\t12-08-2026\n"
    )
    with pytest.raises(FunctionalParseError, match="precedes start"):
        parse_interproscan(path)


def test_an_evalue_in_the_score_column_is_recorded_as_an_evalue():
    """InterProScan puts the e-value in the score column for HMM analyses.
    Reporting 1.2e-40 as a 'score' inverts the ordering a user would sort by."""
    hits = [hit for hit in parse_functional(INTERPRO).hits if hit.evalue is not None]
    assert hits
    assert all(hit.evalue < 0.01 for hit in hits)
    assert all(hit.score is None for hit in hits)


def test_the_interpro_entry_is_a_separate_assignment():
    """Two signatures from different analyses commonly map to one InterPro entry,
    so 'which proteins have IPR017853' must give one row per protein."""
    hits = parse_functional(INTERPRO).hits
    interpro = [hit for hit in hits if hit.analysis == "InterPro"]
    assert interpro
    assert all(hit.accession.startswith("IPR") for hit in interpro)


def test_queries_with_no_hits_still_count_as_queries(tmp_path):
    """The match rate must be computed against every protein the tool was asked
    about, not only those that produced an annotation."""
    path = tmp_path / "sparse.emapper.annotations"
    path.write_text(
        "#query\tPFAMs\tGOs\n" "XP_1.1\tPF00067\t-\n" "XP_2.1\t-\t-\n"
    )
    parsed = parse_functional(path)
    assert parsed.queries == {"XP_1.1", "XP_2.1"}
    assert len(parsed.hits) == 1


# -- loading ---------------------------------------------------------------


@pytest.fixture
def loaded(tmp_path):
    """A real RefSeq assembly, loaded, ready to annotate."""
    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(engine)

    coordinates = [
        int(column)
        for line in REFSEQ_GFF.read_text().splitlines()
        if not line.startswith("#") and line.count("\t") >= 8
        for column in line.split("\t")[3:5]
    ]
    length = max(coordinates or [10000]) + 1000
    # Written in blocks rather than character by character: the fixture needs a
    # FASTA long enough for the record's coordinates (~5 Mb) and generating that
    # one base at a time dominated the test run.
    fasta = tmp_path / "NC_026477.1.fa"
    block = ("ACGT" * 15) + "\n"
    with fasta.open("w") as handle:
        handle.write(">NC_026477.1\n")
        handle.write(block * (length // 60 + 1))

    annotation = tmp_path / REFSEQ_GFF.name
    annotation.write_text(REFSEQ_GFF.read_text())
    config = tmp_path / "organism.toml"
    config.write_text(
        f'[organism]\ngenus = "Fusarium"\nspecies = "graminearum"\nstrain = "PH-1"\n'
        f'ncbi_taxon_id = 229533\n\n[assembly]\nname = "ASM24013v3"\n'
        f'fasta = "{fasta.name}"\n\n[[files.annotation]]\npath = "{annotation.name}"\n'
    )
    report = load_organism(engine, load_config(config), tmp_path / "data")
    return engine, report.assembly_id, tmp_path


@pytest.mark.parametrize("fixture", [EGGNOG, INTERPRO])
def test_a_functional_file_attaches_to_real_features(loaded, fixture):
    engine, assembly_id, _ = loaded
    report = load_functional(engine, assembly_id, fixture)

    assert report.query_count == 20
    assert report.matched_count == 20
    assert report.unmatched_count == 0
    assert report.match_rate == 1.0
    assert report.hit_count > 0

    with engine.connect() as connection:
        stored = connection.execute(
            select(func.count()).select_from(protein_annotation)
        ).scalar_one()
        run = connection.execute(select(annotation_run)).one()
    assert stored == report.hit_count
    assert run.matched_count == 20


def test_the_join_goes_through_protein_id(loaded):
    """The queries are RefSeq protein accessions, which appear as the
    `protein_id` attribute rather than as any feature's ID."""
    engine, assembly_id, _ = loaded
    load_functional(engine, assembly_id, EGGNOG)
    with engine.connect() as connection:
        pairs = connection.execute(
            select(protein_annotation.c.query_id, feature.c.source_id)
            .select_from(
                protein_annotation.join(
                    feature, feature.c.feature_id == protein_annotation.c.feature_id
                )
            )
            .distinct()
        ).all()
    assert pairs
    assert all(pair.query_id != pair.source_id for pair in pairs)


def test_loading_functional_annotation_changes_no_feature_count(loaded):
    """D-013's scope limit made checkable: functional annotation attaches to the
    graph rather than becoming part of it, so `verify` is untouched by it."""
    engine, assembly_id, _ = loaded
    with engine.connect() as connection:
        before = dict(
            connection.execute(
                select(feature.c.source_type, feature_count_expression()).group_by(
                    feature.c.source_type
                )
            ).all()
        )
    load_functional(engine, assembly_id, EGGNOG)
    load_functional(engine, assembly_id, INTERPRO)
    with engine.connect() as connection:
        after = dict(
            connection.execute(
                select(feature.c.source_type, feature_count_expression()).group_by(
                    feature.c.source_type
                )
            ).all()
        )
    assert before == after


def test_a_repeated_domain_at_two_positions_survives_as_two_rows(loaded):
    """InterProScan reports a repeated zinc finger twice, at different positions.
    Both are real; collapsing them would lose the architecture."""
    engine, assembly_id, _ = loaded
    load_functional(engine, assembly_id, INTERPRO)
    with engine.connect() as connection:
        rows = connection.execute(
            select(protein_annotation.c.protein_start)
            .where(protein_annotation.c.accession == "PF00096")
            .order_by(protein_annotation.c.protein_start)
        ).scalars().all()
    assert len(rows) >= 2
    assert len(set(rows)) == len(rows)


def test_an_exact_duplicate_collapses(loaded, tmp_path):
    """eggNOG lists a term once per source column, so the same GO term can arrive
    twice for one protein. That is one fact, not two."""
    engine, assembly_id, _ = loaded
    queries = re.findall(r"protein_id=([^;\n]+)", REFSEQ_GFF.read_text())
    path = tmp_path / "dupes.emapper.annotations"
    path.write_text(
        "#query\tGOs\n"
        f"{queries[0]}\tGO:0005515,GO:0005515,GO:0005515\n"
    )
    report = load_functional(engine, assembly_id, path)
    assert report.hit_count == 1


# -- refusals --------------------------------------------------------------


def test_a_file_that_matches_nothing_is_refused_with_a_count(loaded, tmp_path):
    """The classic silent partial success. A functional file produced from a
    different annotation must not load quietly."""
    engine, assembly_id, _ = loaded
    path = tmp_path / "wrong.emapper.annotations"
    path.write_text("#query\tPFAMs\nNOT_A_REAL_PROTEIN.1\tPF00067\n")
    with pytest.raises(FunctionalLoadError, match="match no feature"):
        load_functional(engine, assembly_id, path)


def test_a_partial_match_is_refused_unless_asked_for(loaded, tmp_path):
    engine, assembly_id, _ = loaded
    queries = re.findall(r"protein_id=([^;\n]+)", REFSEQ_GFF.read_text())
    path = tmp_path / "half.emapper.annotations"
    path.write_text(
        "#query\tPFAMs\n"
        f"{queries[0]}\tPF00067\n"
        "NOT_REAL_1.1\tPF00067\n"
    )
    with pytest.raises(FunctionalLoadError, match="50.0%"):
        load_functional(engine, assembly_id, path)

    report = load_functional(engine, assembly_id, path, minimum_match_rate=0.5)
    assert report.matched_count == 1
    assert report.unmatched_count == 1


def test_the_unmatched_count_is_stored_not_only_reported(loaded, tmp_path):
    """An unmatched query leaves no row behind, so without this column the fact
    that half a file went nowhere is unrecoverable after the load."""
    engine, assembly_id, _ = loaded
    queries = re.findall(r"protein_id=([^;\n]+)", REFSEQ_GFF.read_text())
    path = tmp_path / "half.emapper.annotations"
    path.write_text(f"#query\tPFAMs\n{queries[0]}\tPF00067\nNOT_REAL_1.1\tPF00067\n")
    load_functional(engine, assembly_id, path, minimum_match_rate=0.5)
    with engine.connect() as connection:
        run = connection.execute(select(annotation_run)).one()
    assert (run.query_count, run.matched_count, run.unmatched_count) == (2, 1, 1)


def test_a_query_matching_two_features_is_refused(loaded, tmp_path):
    """Attaching a protein's annotation to the wrong isoform is invisible and
    unrecoverable, so an ambiguous identifier stops the load."""
    engine, assembly_id, _ = loaded
    from cgload.db.schema import feature_attribute
    from sqlalchemy import insert

    with engine.begin() as connection:
        targets = connection.execute(
            select(feature.c.feature_id)
            .where(feature.c.feature_type == "CDS", feature.c.segment_index == 0)
            .limit(2)
        ).scalars().all()
        connection.execute(
            insert(feature_attribute),
            [
                {"feature_id": target, "key": "protein_id", "value": "AMBIGUOUS.1"}
                for target in targets
            ],
        )
    path = tmp_path / "ambiguous.emapper.annotations"
    path.write_text("#query\tPFAMs\nAMBIGUOUS.1\tPF00067\n")
    with pytest.raises(FunctionalLoadError, match="matches 2 features"):
        load_functional(engine, assembly_id, path)


def test_attaching_to_a_missing_assembly_is_refused(loaded):
    engine, _assembly_id, _ = loaded
    with pytest.raises(FunctionalLoadError, match="no assembly with id"):
        load_functional(engine, 999, EGGNOG)


def test_the_same_file_cannot_be_attached_twice(loaded):
    engine, assembly_id, _ = loaded
    load_functional(engine, assembly_id, EGGNOG)
    with pytest.raises(FunctionalLoadError, match="already attached"):
        load_functional(engine, assembly_id, EGGNOG)


# -- through the config and the load command ------------------------------


def test_functional_files_load_in_the_same_transaction_as_the_genome(tmp_path):
    """A failure attaching annotation must not leave a genome loaded with its
    annotation half-attached."""
    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(engine)

    fasta = tmp_path / "NC_026477.1.fa"
    length = 6_100_000
    with fasta.open("w") as handle:
        handle.write(">NC_026477.1\n")
        handle.write((("ACGT" * 15) + "\n") * (length // 60))

    (tmp_path / REFSEQ_GFF.name).write_text(REFSEQ_GFF.read_text())
    broken = tmp_path / "wrong.emapper.annotations"
    broken.write_text("#query\tPFAMs\nNOT_A_REAL_PROTEIN.1\tPF00067\n")
    config = tmp_path / "organism.toml"
    config.write_text(
        f'[organism]\ngenus = "Fusarium"\nspecies = "graminearum"\nstrain = "PH-1"\n\n'
        f'[assembly]\nname = "ASM24013v3"\nfasta = "{fasta.name}"\n\n'
        f'[[files.annotation]]\npath = "{REFSEQ_GFF.name}"\n\n'
        f'[[files.functional]]\npath = "{broken.name}"\n'
    )

    with pytest.raises(FunctionalLoadError, match="match no feature"):
        load_organism(engine, load_config(config), tmp_path / "data")

    with engine.connect() as connection:
        for table in (feature, annotation_run, protein_annotation):
            assert (
                connection.execute(select(func.count()).select_from(table)).scalar_one() == 0
            )
    assert [p for p in (tmp_path / "data").rglob("*") if p.is_file()] == []


def test_config_reads_functional_entries(tmp_path):
    config = tmp_path / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "G"\nspecies = "s"\n\n[assembly]\nname = "A"\n'
        'fasta = "x.fa"\n\n[[files.functional]]\npath = "a.annotations"\n\n'
        '[[files.functional]]\npath = "b.tsv"\ntool = "interproscan"\n'
        "minimum_match_rate = 0.8\n"
    )
    parsed = load_config(config)
    assert len(parsed.functional) == 2
    assert parsed.functional[0].minimum_match_rate == 1.0
    assert parsed.functional[1].tool == "interproscan"
    assert parsed.functional[1].minimum_match_rate == 0.8
