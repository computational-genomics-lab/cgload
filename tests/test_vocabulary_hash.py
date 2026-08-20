"""The content hash must move when linkage meaning moves (D-017, D-061)."""

from __future__ import annotations

from cgload.db import vocabulary
from cgload.db.vocabulary import content_hash


def test_hash_covers_allowed_parents(monkeypatch):
    before = content_hash()
    widened = dict(vocabulary.ALLOWED_PARENTS)
    widened["exon"] = widened["exon"] | {"repeat_region"}
    monkeypatch.setattr(vocabulary, "ALLOWED_PARENTS", widened)
    assert content_hash() != before


def test_hash_covers_transcript_terms(monkeypatch):
    """D-061. Adding a transcript term changes which files load."""
    before = content_hash()
    monkeypatch.setattr(
        vocabulary, "TRANSCRIPT_TERMS", vocabulary.TRANSCRIPT_TERMS | {"miRNA"}
    )
    assert content_hash() != before


def test_hash_covers_transcript_layer_required(monkeypatch):
    """D-061. Dropping `exon` here silently permits every concealed transcript."""
    before = content_hash()
    monkeypatch.setattr(
        vocabulary, "TRANSCRIPT_LAYER_REQUIRED", frozenset({"CDS"})
    )
    assert content_hash() != before


def test_hash_covers_the_synonym_table(monkeypatch):
    """The sentinel must be a term not already mapped.

    The first version used `miRNA`, which D-062 had since added to SYNONYMS, so
    the 'extended' table was identical to the real one and the hash correctly did
    not move. The test failed while the code was right -- a test asserting a
    property of its own fixture rather than of the system.
    """
    sentinel = "cgload_test_sentinel_term"
    assert sentinel not in vocabulary.SYNONYMS

    before = content_hash()
    extended = dict(vocabulary.SYNONYMS)
    extended[sentinel] = "ncRNA"
    monkeypatch.setattr(vocabulary, "SYNONYMS", extended)
    assert content_hash() != before


def test_hash_is_stable_across_calls():
    assert content_hash() == content_hash()
