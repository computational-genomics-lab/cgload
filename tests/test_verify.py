"""Milestone 8: completeness checking.

Every number in these tests was measured on real NCBI files. Where a test names
a full-record count, the file is not shipped -- see `tests/fixtures/PROVENANCE.md`
for accessions and the command to reproduce.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest
from cgload.verify import (
    Outcome,
    VerifyError,
    compare_across_formats,
    compare_with_reference,
    count_source_naively,
    exit_code,
    parse_ncbi_feature_counts,
    render,
)

FIXTURES = Path(__file__).parent / "fixtures"
NCBI_COUNTS = FIXTURES / "ncbi_feature_count_GCF_000328475.txt"
MYCMAY_GFF = FIXTURES / "refseq_mycmay_GCF_000328475.gff3"
GENBANK_PRO = FIXTURES / "genbank_maeruginosa_NC_010296.gbff"
GENBANK_EUK = FIXTURES / "genbank_dscam1_CG12164.gbff"


# ---------------------------------------------------------------------------
# The naive pass
# ---------------------------------------------------------------------------


def test_naive_count_groups_multi_segment_features():
    """A CDS spanning several exons is one feature written on several rows.

    Measured on the full *U. maydis* GFF3: 9,739 CDS rows, 6,780 CDS features.
    NCBI publishes 6,780. Counting rows would overstate by 43%.
    """
    counts = count_source_naively(MYCMAY_GFF)
    rows = sum(
        1
        for line in MYCMAY_GFF.read_text().splitlines()
        if not line.startswith("#") and len(line.split("\t")) >= 9
        and line.split("\t")[2] == "CDS"
    )
    assert counts["CDS"] < rows, "segments were not grouped"


def test_naive_count_reads_genbank_records_as_whole_features():
    """A GenBank record already is the feature: `join(1..96,200..300)` is one
    CDS on one record, so no grouping is needed and none is applied."""
    counts = count_source_naively(GENBANK_PRO)
    assert counts["gene"] > 0
    assert counts["CDS"] > 0
    assert "source" not in counts, "the source record is sequence metadata"


def test_naive_pass_imports_no_parser_code():
    """D-008: the check must not share code with the thing it checks.

    A verifier that reused the parser could only catch storage bugs. The bugs
    that matter most -- a `return` indented one level too deep, a silent
    `else: pass` -- happen during parsing and are invisible to any check that
    parses the file the same way.

    Scanned as code, not as text. A plain substring search also matched the words
    `tokenize` and `normalise` inside explanatory comments, so writing a comment
    about why the naive pass exists broke the guard protecting it. Comments
    cannot import anything; the AST is what decides.
    """
    import ast

    source = (Path(__file__).parents[1] / "src" / "cgload" / "verify.py").read_text()
    tree = ast.parse(source)

    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
            imported.update(f"{node.module}.{alias.name}" for alias in node.names)

    for module in sorted(imported):
        assert not module.startswith("cgload.parsers"), (
            f"verify.py imports {module!r}; the naive pass must share no code "
            f"with the parser (D-008)"
        )

    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    } | {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    for forbidden in ("tokenize", "normalise", "detect", "canonical_type"):
        assert forbidden not in called, f"verify.py calls {forbidden!r}"


# ---------------------------------------------------------------------------
# NCBI published counts: three traps in a seven-line file
# ---------------------------------------------------------------------------


def test_counts_are_summed_across_class_within_a_feature():
    """`gene` appears three times -- protein_coding 6,762, rRNA 34, tRNA 111.
    Only the sum, 6,907, equals the GFF3's gene count. Reading one row gives
    6,762 and reports 145 missing genes on a correct load."""
    counts = parse_ncbi_feature_counts(NCBI_COUNTS, "GCF_000328475.2")
    assert counts["gene"] == 6907


def test_na_in_unique_ids_falls_back_to_placements():
    """The `tRNA` row of this real file has `na` in Unique Ids and 111 in
    Placements. A reader taking only the first column drops the type entirely
    and then reports 111 missing tRNAs."""
    counts = parse_ncbi_feature_counts(NCBI_COUNTS, "GCF_000328475.2")
    assert counts["tRNA"] == 111


def test_assembly_is_matched_on_the_full_assembly_column():
    """`Assembly-unit accession` is GCF_000328485.2 while the assembly is
    GCF_000328475.2 -- one digit apart. Matching the wrong column selects
    nothing, and a verifier that silently found nothing would report every
    feature as missing."""
    with pytest.raises(VerifyError, match="Full Assembly"):
        parse_ncbi_feature_counts(NCBI_COUNTS, "GCF_000328485.2")


# ---------------------------------------------------------------------------
# Tier 3 fails only when it can mean something
# ---------------------------------------------------------------------------


def test_reference_check_is_skipped_for_a_partial_assembly():
    """Published totals describe the whole assembly. Every fixture in this
    repository is one or two sequences of many, so on fixtures this tier must
    skip -- and say so, rather than reporting a discrepancy it cannot stand
    behind (D-033)."""
    result = compare_with_reference(
        Counter({"gene": 15}),
        Counter({"gene": 6907}),
        sequences_loaded=2,
        sequences_in_assembly=23,
    )
    assert result.outcome is Outcome.SKIPPED
    assert "2 of 23" in result.detail


def test_reference_check_passes_on_the_whole_real_assembly():
    """Measured: the naive pass over the complete *U. maydis* GFF3 equals NCBI's
    published totals on all five shared types, including CDS at feature level."""
    source = Counter(
        {"gene": 6907, "mRNA": 6780, "CDS": 6780, "tRNA": 111, "rRNA": 34}
    )
    reference = parse_ncbi_feature_counts(NCBI_COUNTS, "GCF_000328475.2")
    result = compare_with_reference(
        source, reference, sequences_loaded=23, sequences_in_assembly=23
    )
    assert result.outcome is Outcome.PASS


def test_a_skipped_tier_is_not_a_failure():
    skipped = compare_with_reference(
        Counter(), Counter(), sequences_loaded=1, sequences_in_assembly=9
    )
    assert exit_code([skipped]) == 0


# ---------------------------------------------------------------------------
# Tier 2: the cross-parser oracle
# ---------------------------------------------------------------------------


def test_pseudogene_naming_difference_is_not_a_discrepancy():
    """GFF3 gives a pseudogene its own column-3 type; GenBank writes it as a
    `gene` record with `/pseudo`. Measured on GCF_000240135.3: GFF3 has 13,714
    gene + 11 pseudogene, GenBank has 13,725 gene. 13,714 + 11 = 13,725, so the
    two agree exactly once the mapping is applied.
    """
    from cgload.verify import GENBANK_TYPE_EQUIVALENCE

    assert GENBANK_TYPE_EQUIVALENCE["gene"] == ("gene", "pseudogene")


def test_cross_parser_catches_a_deleted_feature(tmp_path):
    """Mutation test. Delete one gene from the GFF3 and the GenBank count no
    longer agrees. A GFF3 parser that loses a feature cannot make the GenBank
    count lose it too, which is why this check cannot be faked.
    """
    original = GENBANK_EUK

    # The agreeing file has to agree on *every* non-structural type, not only on
    # `gene`. The tiers compare the union now (D-072), so a GFF3 carrying genes
    # alone no longer matches a GenBank file that also has transcripts and coding
    # sequences: it has lost all of both, which is precisely the total-loss case
    # the union was introduced to catch. Under the old intersection semantics this
    # file passed, which is why the test had to change with the rule.
    counts = count_source_naively(original)
    assert counts["gene"], "fixture no longer has genes"

    rows = ["##gff-version 3\n"]
    position = 1
    for feature_type, total in sorted(counts.items()):
        for index in range(total):
            rows.append(
                f"chr1\ttest\t{feature_type}\t{position}\t{position + 8}\t.\t+\t.\t"
                f"ID={feature_type}{index}\n"
            )
            position += 10

    agreeing = tmp_path / "agree.gff3"
    agreeing.write_text("".join(rows))
    agreed = compare_across_formats(agreeing, original)
    assert agreed.outcome is Outcome.PASS, agreed.mismatches

    # Drop exactly one gene row.
    victim = next(row for row in rows if "\tgene\t" in row)
    disagreeing = tmp_path / "disagree.gff3"
    disagreeing.write_text("".join(row for row in rows if row is not victim))
    result = compare_across_formats(disagreeing, original)
    assert result.outcome is Outcome.FAIL
    assert "gene" in result.mismatches


def test_types_absent_by_construction_are_not_reported_as_loss():
    """GFF3 states exons explicitly; GenBank implies them from the location
    expression. Zero GenBank exons is a fact about the format, not a loss --
    the same distinction the `provides` profile field exists to record."""
    result = compare_across_formats(MYCMAY_GFF, GENBANK_PRO)
    assert "exon" not in result.expected
    assert "exon" not in result.actual


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def test_failure_exits_two_and_names_the_delta():
    failing = compare_with_reference(
        Counter({"gene": 6900}),
        Counter({"gene": 6907}),
        sequences_loaded=23,
        sequences_in_assembly=23,
    )
    assert failing.outcome is Outcome.FAIL
    assert exit_code([failing]) == 2
    assert failing.mismatches["gene"] == (6907, 6900)
    assert "-7" in render([failing])


def test_every_tier_appears_in_the_report_including_skips():
    """A silently omitted tier is indistinguishable from a passing one."""
    results = [
        compare_with_reference(
            Counter({"gene": 1}), Counter({"gene": 1}),
            sequences_loaded=1, sequences_in_assembly=1,
        ),
        compare_with_reference(
            Counter(), Counter(), sequences_loaded=1, sequences_in_assembly=9
        ),
    ]
    text = render(results)
    assert "PASS" in text
    assert "----" in text


def test_one_feature_under_its_two_names_is_not_a_discrepancy(tmp_path):
    """D-102. NCBI's GFF3 writes a Sequence Ontology term where its GenBank file
    writes an INSDC key, and for `ncRNA` and `repeat_region` the term depends on
    a qualifier. Measured on the benchmark's five fungi with organelle-free
    differences: each family agrees once folded. A real loss still fails."""
    from cgload.verify import compare_across_formats

    gff = tmp_path / "a.gff3"
    rows = ["##gff-version 3\n"]
    for index, kind in enumerate(("gene", "pseudogene", "snoRNA", "snRNA", "ncRNA",
                                  "transcript", "origin_of_replication",
                                  "mobile_genetic_element", "long_terminal_repeat",
                                  "repeat_region", "regulatory_region")):
        rows.append(f"chr1\tRefSeq\t{kind}\t{10 * index + 1}\t{10 * index + 5}\t.\t+\t.\t"
                    f"ID=f{index}\n")
    gff.write_text("".join(rows))
    keys = ("gene", "gene", "ncRNA", "ncRNA", "ncRNA", "misc_RNA", "rep_origin",
            "mobile_element", "repeat_region", "repeat_region", "regulatory",
            "assembly_gap")
    genbank = tmp_path / "a.gbff"
    genbank.write_text(
        "LOCUS       chr1                     200 bp    DNA     linear   CON 01-JAN-2020\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..200\n"
        + "".join(f"     {key:<16}{10 * index + 1}..{10 * index + 5}\n"
                  for index, key in enumerate(keys))
        + "ORIGIN\n//\n")
    assert compare_across_formats(gff, genbank).outcome is Outcome.PASS

    gff.write_text("".join(row for row in rows if "\tsnRNA\t" not in row))
    lost = compare_across_formats(gff, genbank)
    assert lost.outcome is Outcome.FAIL
    assert lost.mismatches == {"ncRNA": (2, 3)}



#: A nuclear record with one ordinary gene, one pseudogene stated as a CDS with
#: no mRNA, and a gene sharing the pseudogene's symbol; and a mitochondrion whose
#: GenBank record states CDS alone. Modelled on yeast (GCF_000146045.2): YDR134C
#: and NC_001224.1.
SPLIT_GENBANK = (
    "LOCUS       NC_000001                900 bp    DNA     linear   CON 01-JAN-2020\n"
    "VERSION     NC_000001.1\n"
    "FEATURES             Location/Qualifiers\n"
    "     source          1..900\n"
    "     gene            10..90\n"
    '                     /locus_tag="YA1"\n'
    '                     /db_xref="GeneID:101"\n'
    "     mRNA            10..90\n"
    '                     /locus_tag="YA1"\n'
    '                     /db_xref="GeneID:101"\n'
    "     CDS             20..80\n"
    '                     /locus_tag="YA1"\n'
    '                     /db_xref="GeneID:101"\n'
    "     gene            110..190\n"
    '                     /gene="ABC1"\n'
    '                     /locus_tag="YB1"\n'
    "                     /pseudo\n"
    '                     /db_xref="GeneID:102"\n'
    "     CDS             110..190\n"
    '                     /gene="ABC1"\n'
    '                     /locus_tag="YB1"\n'
    '                     /note="Blocked reading frame"\n'
    "                     /pseudo\n"
    '                     /db_xref="GeneID:102"\n'
    "     gene            210..290\n"
    '                     /gene="ABC1"\n'
    '                     /locus_tag="YC1"\n'
    "     mRNA            210..290\n"
    '                     /gene="ABC1"\n'
    '                     /locus_tag="YC1"\n'
    "     CDS             220..280\n"
    '                     /gene="ABC1"\n'
    '                     /locus_tag="YC1"\n'
    "ORIGIN\n//\n"
    "LOCUS       NC_000002                400 bp    DNA     circular CON 01-JAN-2020\n"
    "VERSION     NC_000002.1\n"
    "FEATURES             Location/Qualifiers\n"
    "     source          1..400\n"
    '                     /organelle="mitochondrion"\n'
    "     gene            10..90\n"
    '                     /locus_tag="Q1"\n'
    "     CDS             10..90\n"
    '                     /locus_tag="Q1"\n'
    "     gene            110..190\n"
    '                     /locus_tag="Q2"\n'
    "     CDS             110..190\n"
    '                     /locus_tag="Q2"\n'
    "     rRNA            210..290\n"
    '                     /locus_tag="Q3"\n'
    "ORIGIN\n//\n"
)


def split_gff(*, drop: str = "", extra: str = "") -> str:
    """The same genes as NCBI's GFF3 writes them: an mRNA for every one."""
    rows = [
        "##gff-version 3\n",
        "NC_000001.1\tRefSeq\tgene\t10\t90\t.\t+\t.\tID=gene-YA1;locus_tag=YA1\n",
        "NC_000001.1\tRefSeq\tmRNA\t10\t90\t.\t+\t.\tID=rna-NM_1.1;Parent=gene-YA1;"
        "Dbxref=GeneID:101;locus_tag=YA1;transcript_id=NM_1.1\n",
        "NC_000001.1\tRefSeq\tCDS\t20\t80\t.\t+\t0\tID=cds-1;Parent=rna-NM_1.1\n",
        "NC_000001.1\tRefSeq\tpseudogene\t110\t190\t.\t+\t.\tID=gene-YB1;locus_tag=YB1\n",
        "NC_000001.1\tRefSeq\tmRNA\t110\t190\t.\t+\t.\tID=rna-YB1;Parent=gene-YB1;"
        "Dbxref=GeneID:102;gene=ABC1;locus_tag=YB1;pseudo=true\n",
        "NC_000001.1\tRefSeq\tCDS\t110\t190\t.\t+\t0\tID=cds-YB1;Parent=rna-YB1\n",
        "NC_000001.1\tRefSeq\tgene\t210\t290\t.\t+\t.\tID=gene-YC1;locus_tag=YC1\n",
        "NC_000001.1\tRefSeq\tmRNA\t210\t290\t.\t+\t.\tID=rna-NM_3.1;Parent=gene-YC1;"
        "gene=ABC1;locus_tag=YC1;transcript_id=NM_3.1\n",
        "NC_000001.1\tRefSeq\tCDS\t220\t280\t.\t+\t0\tID=cds-3;Parent=rna-NM_3.1\n",
        "NC_000002.1\tRefSeq\tgene\t10\t90\t.\t+\t.\tID=gene-Q1;locus_tag=Q1\n",
        "NC_000002.1\tRefSeq\tmRNA\t10\t90\t.\t+\t.\tID=rna-Q1;Parent=gene-Q1;locus_tag=Q1\n",
        "NC_000002.1\tRefSeq\tCDS\t10\t90\t.\t+\t0\tID=cds-Q1;Parent=rna-Q1\n",
        "NC_000002.1\tRefSeq\tgene\t110\t190\t.\t+\t.\tID=gene-Q2;locus_tag=Q2\n",
        "NC_000002.1\tRefSeq\tmRNA\t110\t190\t.\t+\t.\tID=rna-Q2;Parent=gene-Q2;locus_tag=Q2\n",
        "NC_000002.1\tRefSeq\tCDS\t110\t190\t.\t+\t0\tID=cds-Q2;Parent=rna-Q2\n",
        "NC_000002.1\tRefSeq\trRNA\t210\t290\t.\t+\t.\tID=rna-Q3;locus_tag=Q3\n",
    ]
    return "".join(row for row in rows if not (drop and f"ID={drop};" in row)) + extra


def test_mrna_is_not_compared_where_genbank_states_cds_without_mrna(tmp_path):
    """D-103, D-105. NCBI's GFF3 gives an mRNA to every gene whose CDS its
    GenBank file states with no mRNA: yeast's 19 mitochondrial genes and its six
    pseudogenes with a blocked reading frame. Those are set aside, and the
    report says how many."""
    genbank = tmp_path / "a.gbff"
    genbank.write_text(SPLIT_GENBANK)
    gff = tmp_path / "a.gff3"
    gff.write_text(split_gff())
    result = compare_across_formats(gff, genbank)
    assert result.outcome is Outcome.PASS, result.mismatches
    assert "3 GFF3 mRNA" in result.detail
    assert result.expected["mRNA"] == result.actual["mRNA"] == 2


def test_mrna_is_still_compared_for_every_other_gene(tmp_path):
    """Losing an ordinary gene's mRNA still fails, including the gene that
    shares a symbol with the pseudogene, and so does an mRNA the GFF3 gives a
    gene whose GenBank record has no CDS to account for it."""
    genbank = tmp_path / "a.gbff"
    genbank.write_text(SPLIT_GENBANK)
    gff = tmp_path / "a.gff3"
    for lost in ("rna-NM_1.1", "rna-NM_3.1"):
        gff.write_text(split_gff(drop=lost))
        assert compare_across_formats(gff, genbank).mismatches == {"mRNA": (1, 2)}
    gff.write_text(split_gff(extra="NC_000002.1\tRefSeq\tmRNA\t210\t290\t.\t+\t.\t"
                                   "ID=rna-Q3b;locus_tag=Q3\n"))
    assert compare_across_formats(gff, genbank).mismatches == {"mRNA": (3, 2)}


def test_a_family_member_is_matched_without_regard_to_case(tmp_path):
    """D-106. INSDC qualifier values are case-insensitive, and NCBI's GFF3
    copies /rpt_type into column 3 as written: *S. pombe*'s LTRs are
    `LONG_TERMINAL_REPEAT` rows, its GenBank file's `repeat_region` records."""
    gff = tmp_path / "a.gff3"
    gff.write_text("##gff-version 3\n"
                   "chr1\tRefSeq\tLONG_TERMINAL_REPEAT\t1\t300\t.\t+\t.\tID=id-1\n")
    genbank = tmp_path / "a.gbff"
    genbank.write_text(
        "LOCUS       chr1                    1000 bp    DNA     linear   CON 01-JAN-2020\n"
        "FEATURES             Location/Qualifiers\n"
        "     source          1..1000\n"
        "     repeat_region   1..300\n"
        '                     /rpt_type=LONG_TERMINAL_REPEAT\n'
        "ORIGIN\n//\n")
    assert compare_across_formats(gff, genbank).outcome is Outcome.PASS
