"""Regressions for four gaps found by reading the verification code.

Each corresponds to a case where a wrong load reported success.
"""

from __future__ import annotations

from collections import Counter

from cgload.verify import (
    Outcome,
    check_parent_edges,
    compare_with_reference,
    exit_code,
)

# --- D-072: total loss of a type must not be suppressed by an intersection ---


def test_total_loss_of_a_type_fails():
    """An intersection made the most alarming failure the one that could not
    fail: a type absent from one side was dropped from the comparison. Losing
    321 of 322 tRNAs failed; losing all 322 passed."""
    result = compare_with_reference(
        Counter({"gene": 100, "mRNA": 100}),
        Counter({"gene": 100, "mRNA": 100, "tRNA": 322}),
        sequences_loaded=1,
        sequences_in_assembly=1,
    )
    assert result.outcome is Outcome.FAIL
    assert result.mismatches["tRNA"] == (322, 0)


def test_partial_loss_still_fails():
    result = compare_with_reference(
        Counter({"tRNA": 321}),
        Counter({"tRNA": 322}),
        sequences_loaded=1,
        sequences_in_assembly=1,
    )
    assert result.outcome is Outcome.FAIL


def test_a_correct_load_still_passes():
    result = compare_with_reference(
        Counter({"gene": 100, "tRNA": 322}),
        Counter({"gene": 100, "tRNA": 322}),
        sequences_loaded=1,
        sequences_in_assembly=1,
    )
    assert result.outcome is Outcome.PASS


# --- D-074: parent edges are computed and were never verified ---


def test_zero_parent_edges_fails():
    """A GenBank load producing no hierarchy at all passed every tier: each
    feature was present and counted, and nothing looked at the edges."""
    result = check_parent_edges({None: 5793}, dialect="genbank")
    assert result.outcome is Outcome.FAIL
    assert "hierarchy is absent" in result.detail


def test_some_unlinked_children_fails():
    result = check_parent_edges({"containment": 100, None: 7}, dialect="genbank")
    assert result.outcome is Outcome.FAIL
    assert "7" in result.detail


def test_a_fully_linked_load_passes():
    result = check_parent_edges(
        {"containment": 13312, "locus_tag": 322}, dialect="genbank"
    )
    assert result.outcome is Outcome.PASS
    assert "13,634" in result.detail


def test_genbank_may_not_claim_explicit_parents():
    """GenBank carries no parent pointers, so an edge recorded as stated by the
    file means the tokenizer's synthesised attribute was read back as the
    file's own -- the provenance bug, arriving from the other direction."""
    result = check_parent_edges(
        {"explicit_parent": 86, "containment": 62}, dialect="genbank"
    )
    assert result.outcome is Outcome.FAIL
    assert "no parent pointers" in result.detail


def test_gff3_explicit_parents_are_expected():
    result = check_parent_edges({"explicit_parent": 354}, dialect="refseq")
    assert result.outcome is Outcome.PASS


def test_no_parented_types_skips_rather_than_passing():
    result = check_parent_edges({}, dialect="genbank")
    assert result.outcome is Outcome.SKIPPED
    assert exit_code([result]) == 0
