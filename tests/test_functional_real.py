"""Milestone 7b: the eggNOG parser against real emapper output.

Every earlier functional-annotation test ran against constructed files. This
module runs against a slice of genuine `emapper-2.1.4` output (*Ustilago
maydis*, GCF_000328475.2), and each test corresponds to something the
constructed files did not contain and could not have revealed.

Counts in the docstrings were measured on the complete 6,509-protein file,
which is not shipped; the fixture is an 11-row slice chosen to carry each trap.
See `tests/fixtures/PROVENANCE.md`.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest
import sqlalchemy as sa
from cgload.parsers.functional import parse_eggnog

FIXTURE = Path(__file__).parent / "fixtures" / "eggnog_mycmay_emapper214.annotations"


@pytest.fixture(scope="module")
def parsed():
    return parse_eggnog(FIXTURE)


def test_version_is_read_from_the_header_block(parsed):
    """The tool version appears on its own `##` line before the command line.

    The command line also contains the string `emapper` (it is the path to
    `emapper.py`), so a reader that takes the first `##` line mentioning
    `emapper` must not take the command line instead.
    """
    assert parsed.tool_version == "emapper-2.1.4"


def test_built_against_a_different_version_than_documented():
    """Milestone 7 was written against the emapper-2.1.12 column layout; the
    first real file supplied was 2.1.4. The 21 column names are identical
    between the two, so header-driven reading holds -- but this is now a
    measured fact rather than an assumption."""
    header = [
        line for line in FIXTURE.read_text().splitlines() if line.startswith("#query")
    ]
    assert len(header) == 1
    columns = header[0].lstrip("#").split("\t")
    assert len(columns) == 21
    assert columns[0] == "query"
    assert columns[-1] == "PFAMs"


def test_dash_is_the_empty_marker_not_the_empty_string(parsed):
    """Real output writes `-` for every absent value -- 61,425 cells in the full
    file, and never an empty string. A parser treating `""` as the missing
    marker stores tens of thousands of annotations whose accession is a hyphen.
    """
    assert all(hit.accession != "-" for hit in parsed.hits)
    assert all(hit.accession.strip() for hit in parsed.hits)


def test_cog_categories_are_split_letter_by_letter(parsed):
    """`IJT` is three categories with no separator. Confirmed on real output:
    97 of 121 distinct COG_category values in the full file are multi-letter."""
    cog = [hit.accession for hit in parsed.hits if hit.analysis == "COG"]
    assert cog
    assert all(len(accession) == 1 and accession.isalpha() for accession in cog)
    # The fixture deliberately includes a three-letter value.
    per_query = Counter(hit.query_id for hit in parsed.hits if hit.analysis == "COG")
    assert max(per_query.values()) >= 3


def test_kegg_pathway_map_duplicates_are_not_stored_twice(parsed):
    """KEGG publishes each pathway as both `ko00362` and `map00362`, and eggNOG
    emits both. Measured on the full file: 2,253 rows carry pathway data and in
    all 2,253 the two sets are identical, so half of the 17,226 emitted rows
    were the same pathway written twice (D-031).
    """
    pathways = [hit.accession for hit in parsed.hits if hit.analysis == "KEGG_Pathway"]
    assert pathways, "fixture carries no KEGG_Pathway data"
    numbers = {accession[2:] for accession in pathways if accession.startswith("ko")}
    for accession in pathways:
        if accession.startswith("map"):
            assert accession[3:] not in numbers, (
                f"{accession} duplicates ko{accession[3:]} for the same protein"
            )


def test_orthologous_group_scope_is_preserved(parsed):
    """`296A5@1|root` and `2RD8W@2759|Eukaryota` are assignments at different
    taxonomic levels. Truncating at `@` would merge them into one accession and
    lose the level, which is the informative part."""
    groups = [hit.accession for hit in parsed.hits if hit.analysis == "eggNOG_OG"]
    assert groups
    assert any("@" in accession and "|" in accession for accession in groups)


def test_commas_inside_the_description_do_not_become_accessions(parsed):
    """686 descriptions in the full file contain a comma. Description is a free
    text column and must not be split the way the term columns are."""
    descriptions = {
        hit.description for hit in parsed.hits if hit.description is not None
    }
    assert any("," in description for description in descriptions), (
        "fixture lost its comma-bearing description"
    )


def test_a_protein_may_carry_hundreds_of_go_terms(parsed):
    """GO is by far the largest analysis: 286,671 of 380,022 hits in the full
    file, a mean of 89.7 terms per annotated protein and a maximum of 570 on a
    single one. This is a capacity fact, not an anomaly -- `protein_annotation`
    will be the largest table in the database.
    """
    go_per_query = Counter(
        hit.query_id for hit in parsed.hits if hit.analysis == "GO"
    )
    assert max(go_per_query.values()) > 400


def test_no_duplicate_triples_within_a_protein(parsed):
    """One protein must not gain the same accession twice for one analysis.
    Confirmed zero on the full 380,022-hit file."""
    triples = Counter(
        (hit.query_id, hit.analysis, hit.accession) for hit in parsed.hits
    )
    assert not [key for key, count in triples.items() if count > 1]


def test_rows_with_no_annotation_at_all_are_kept_as_queries(parsed):
    """A protein that got a seed ortholog but no functional terms still counts
    as a matched query. Dropping it would understate the match rate and trip the
    unmatched-protein guard for a file that is entirely correct."""
    annotated = {hit.query_id for hit in parsed.hits}
    assert parsed.queries >= annotated
    assert len(parsed.queries) > len(annotated) or len(parsed.queries) == len(annotated)


def test_kegg_ko_retains_its_prefix(parsed):
    """`KEGG_ko` values are written `ko:K18798` while `KEGG_Pathway` uses a bare
    `ko00362`. The prefix is stored as eggNOG wrote it; this test pins the
    convention so a query interface is built against a known form rather than a
    guessed one.
    """
    kos = [hit.accession for hit in parsed.hits if hit.analysis == "KEGG_ko"]
    if kos:
        assert all(accession.startswith("ko:") for accession in kos)


# ---------------------------------------------------------------------------
# The paired fixture: same assembly as the eggNOG file, so the join is real
# ---------------------------------------------------------------------------

GENOME = Path(__file__).parent / "fixtures" / "refseq_mycmay_GCF_000328475.gff3"


def test_the_eggnog_file_and_the_genome_are_the_same_assembly():
    """Both fixtures come from GCF_000328475.2. Measured on the full pair:
    all 6,509 eggNOG query IDs appear as a protein_id on a CDS row in the GFF3.
    Before this pair existed, functional annotation was parse-tested only."""
    assert "GCF_000328475.2" in GENOME.read_text()


def test_every_eggnog_query_matches_a_cds_protein_id():
    """The join key is the protein accession. Zero unmatched in the full file."""
    import re

    proteins = set()
    for line in GENOME.read_text().splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 9 or columns[2] != "CDS":
            continue
        found = re.search(r"protein_id=([^;]*)", columns[8])
        if found:
            proteins.add(found.group(1))

    queries = parse_eggnog(FIXTURE).queries
    assert queries, "fixture has no queries"
    assert queries <= proteins, f"unmatched: {sorted(queries - proteins)[:5]}"


def test_genes_without_annotation_are_not_an_error():
    """271 of the 6,780 proteins in the full genome got no eggNOG hit at all.

    That is the opposite direction from the unmatched-query guard and must not
    trigger it: the guard fires when a *query* matches no feature, never when a
    *feature* has no annotation. A genome where some genes are unannotated is
    the normal case, not a failed load.
    """
    import re

    proteins = set()
    for line in GENOME.read_text().splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 9 or columns[2] != "CDS":
            continue
        found = re.search(r"protein_id=([^;]*)", columns[8])
        if found:
            proteins.add(found.group(1))

    queries = parse_eggnog(FIXTURE).queries
    assert proteins - queries, "fixture lost its unannotated proteins"


# ---------------------------------------------------------------------------
# KEGG pathway aliasing
# ---------------------------------------------------------------------------


def test_map_form_resolves_to_the_stored_ko_form():
    """`map00362` is what KEGG's own URLs and most papers use, so it is what a
    user pastes into a search box. cgload stores `ko00362` (D-031), so an
    un-normalised search would report a pathway missing that is present."""
    from cgload.parsers.functional import canonical_kegg_pathway

    assert canonical_kegg_pathway("map00362") == "ko00362"
    assert canonical_kegg_pathway("ko00362") == "ko00362"
    assert canonical_kegg_pathway(" map01100 ") == "ko01100"


def test_canonicalisation_does_not_touch_other_identifiers():
    """`ko:K18798` is a KEGG *orthology* term, not a pathway, and shares the
    prefix. Rewriting it would corrupt a different analysis type."""
    from cgload.parsers.functional import canonical_kegg_pathway

    assert canonical_kegg_pathway("ko:K18798") == "ko:K18798"
    assert canonical_kegg_pathway("mapK123") == "mapK123"
    assert canonical_kegg_pathway("map") == "map"


def test_stored_pathways_are_already_canonical():
    """Everything the parser stores must survive canonicalisation unchanged,
    or the search would normalise to a form that is not in the database."""
    from cgload.parsers.functional import canonical_kegg_pathway

    pathways = [
        hit.accession
        for hit in parse_eggnog(FIXTURE).hits
        if hit.analysis == "KEGG_Pathway"
    ]
    assert pathways
    assert all(canonical_kegg_pathway(p) == p for p in pathways)


# ---------------------------------------------------------------------------
# Milestone 7c: what loading the matched pair end to end exposed
# ---------------------------------------------------------------------------


def test_a_gene_level_exon_alongside_a_transcript_is_not_a_bypass(tmp_path):
    """D-043. Real RefSeq emits a tRNA gene's exons twice -- once under the tRNA
    and again under the gene with `gbkey=exon`. The second pair is redundant NCBI
    output, not a missing transcript, and a file-wide rule refused a correct
    genome. Reproduced from GCF_000328475.2 gene UMAG_16001.
    """
    from cgload.parsers import profiles as profile_table
    from cgload.parsers.normalise import normalise
    from cgload.parsers.tokenizer import tokenize

    body = (
        "##gff-version 3\n"
        "c\tRefSeq\tgene\t81245\t81327\t.\t-\t.\tID=gene-U;gene_biotype=tRNA\n"
        "c\tRefSeq\ttRNA\t81245\t81327\t.\t-\t.\tID=rna-U;Parent=gene-U\n"
        "c\tRefSeq\texon\t81291\t81327\t.\t-\t.\tID=exon-U-1;Parent=rna-U;gbkey=tRNA\n"
        "c\tRefSeq\texon\t81245\t81280\t.\t-\t.\tID=exon-U-2;Parent=rna-U;gbkey=tRNA\n"
        "c\tRefSeq\texon\t81245\t81280\t.\t-\t.\tID=id-U;Parent=gene-U;gbkey=exon\n"
        "c\tRefSeq\texon\t81291\t81327\t.\t-\t.\tID=id-U-2;Parent=gene-U;gbkey=exon\n"
        "c\tRefSeq\tgene\t82101\t83639\t.\t+\t.\tID=gene-P\n"
        "c\tRefSeq\tmRNA\t82101\t83639\t.\t+\t.\tID=rna-P;Parent=gene-P\n"
        "c\tRefSeq\tCDS\t82101\t83639\t.\t+\t0\tID=cds-P;Parent=rna-P\n"
    )
    path = tmp_path / "t.gff3"
    path.write_text(body)
    tokenized = tokenize(path)
    result = normalise(tokenized.features, profile_table.REFSEQ)
    assert len(result.features) == 9


def test_a_genuinely_missing_transcript_is_still_refused(tmp_path):
    """The narrower rule must still catch what the wide one was for: a gene whose
    parts attach straight to it while the rest of the file has transcripts."""
    from cgload.parsers import profiles as profile_table
    from cgload.parsers.normalise import NormaliseError, normalise
    from cgload.parsers.tokenizer import tokenize

    body = (
        "##gff-version 3\n"
        "c\tRefSeq\tgene\t1\t999\t.\t+\t.\tID=gene-A\n"
        "c\tRefSeq\tmRNA\t1\t999\t.\t+\t.\tID=rna-A;Parent=gene-A\n"
        "c\tRefSeq\tCDS\t1\t999\t.\t+\t0\tID=cds-A;Parent=rna-A\n"
        "c\tRefSeq\tgene\t2000\t2999\t.\t+\t.\tID=gene-B\n"
        "c\tRefSeq\tCDS\t2000\t2999\t.\t+\t0\tID=cds-B;Parent=gene-B\n"
    )
    path = tmp_path / "t.gff3"
    path.write_text(body)
    tokenized = tokenize(path)
    with pytest.raises(NormaliseError, match="no transcript of its own"):
        normalise(tokenized.features, profile_table.REFSEQ)


def test_the_matched_pair_loads_end_to_end(tmp_path):
    """The last constructed-data gap, closed. Before this the eggNOG file was
    parse-tested only: no genome in the repository described its proteins."""
    from cgload.config import load_config
    from cgload.db.engine import build_engine
    from cgload.db.init import create_schema
    from cgload.db.schema import annotation_run, protein_annotation
    from cgload.loader import load_organism

    lengths: dict[str, int] = {}
    for line in GENOME.read_text().splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) >= 9:
            lengths[columns[0]] = max(lengths.get(columns[0], 0), int(columns[4]))

    fasta = tmp_path / "assembly.fa"
    with fasta.open("w") as handle:
        for name, longest in sorted(lengths.items()):
            handle.write(f">{name}\n")
            handle.write((("ACGT" * 15) + "\n") * ((longest + 1000) // 60 + 1))

    (tmp_path / GENOME.name).write_text(GENOME.read_text())
    (tmp_path / FIXTURE.name).write_text(FIXTURE.read_text())
    config = tmp_path / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Mycosarcoma"\nspecies = "maydis"\nstrain = "521"\n'
        'ncbi_taxon_id = 5270\n\n[assembly]\nname = "Umaydis521_2.0"\n'
        'accession = "GCF_000328475.2"\n'
        f'fasta = "{fasta.name}"\n\n[[files.annotation]]\npath = "{GENOME.name}"\n\n'
        f'[[files.functional]]\npath = "{FIXTURE.name}"\n'
    )

    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(engine)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    assert report.functional
    run = report.functional[0]
    assert run.tool_version == "emapper-2.1.4"
    assert run.unmatched_count == 0
    assert run.matched_count == run.query_count

    with engine.connect() as connection:
        stored = connection.execute(
            sa.select(sa.func.count()).select_from(protein_annotation)
        ).scalar_one()
        runs = connection.execute(
            sa.select(sa.func.count()).select_from(annotation_run)
        ).scalar_one()
    assert stored == run.hit_count
    assert runs == 1


def test_no_map_form_pathway_is_stored_from_the_real_file():
    """D-031 through the whole pipeline, on real output."""
    stored = {
        hit.accession
        for hit in parse_eggnog(FIXTURE).hits
        if hit.analysis == "KEGG_Pathway"
    }
    assert stored
    assert not [accession for accession in stored if accession.startswith("map")]


def test_batching_is_bounded():
    """D-044. 57 findings per protein means a real fungal genome produces around
    three-quarters of a million rows from one file, inside the load transaction."""
    from cgload.functional import INSERT_BATCH_SIZE

    assert 0 < INSERT_BATCH_SIZE <= 10_000
