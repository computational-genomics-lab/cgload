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
