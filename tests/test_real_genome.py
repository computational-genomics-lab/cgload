# tests/test_real_genome.py
"""Nightly tier: the cross-parser check on a full published assembly.

Not a duplicate of the bundled-pair test. That pair is 43 genes sliced from one
region; this is 13,725 genes across 19 sequences, and it is the only place the
two parsers are compared at the scale the paper reports.
"""
import os
from pathlib import Path

import pytest

from cgload.verify import Outcome, compare_across_formats

GENOMES = Path(os.environ.get("CGLOAD_REAL_GENOMES", "/tmp/real-genomes"))


@pytest.mark.slow
def test_fusarium_agrees_across_formats():
    gff3 = GENOMES / "fusarium_genomic.gff"
    gbff = GENOMES / "fusarium_genomic.gbff"
    # Fail rather than skip. A nightly tier that quietly skips when its inputs
    # are missing reports success for a check that never ran.
    assert gff3.exists() and gbff.exists(), (
        f"real genome files not found under {GENOMES}; the workflow must fetch "
        f"them before running this tier"
    )
    result = compare_across_formats(gff3, gbff)
    assert result.outcome is Outcome.PASS, result.mismatches
