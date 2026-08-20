"""Every valid fixture must normalise without a violation.

A check that fires on correct data is worse than no check: it teaches whoever
hits it that the checks are noise. This sweep is the standing guard against
that, and it is the test the last five instances of the permitted-set bug class
would each have caught before a real genome did.

`tests/fixtures/broken/` is excluded by design -- those files are deliberately
invalid and are asserted on individually elsewhere.

NOTE ON SCOPE: the delivered archive contained no `tests/fixtures/` directory,
so this sweep has never actually run. It is written to *fail loudly rather than
skip silently* once merged into a tree that has fixtures: if the directory is
absent the test skips with a message saying so, and if it is present but the
sweep finds nothing to parse the test fails. Do not read a green run of this
file in isolation as evidence that the sweep passed.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import NormaliseError, normalise
from cgload.parsers.tokenizer import tokenize

FIXTURES = Path(__file__).parent / "fixtures"
BROKEN = FIXTURES / "broken"

#: Extensions the tokenizer is expected to handle. Widen as profiles are added.
SUFFIXES = {".gff", ".gff3", ".gtf", ".gbk", ".gb", ".genbank", ".embl"}


def valid_fixtures() -> list[Path]:
    if not FIXTURES.is_dir():
        return []
    return sorted(
        path
        for path in FIXTURES.rglob("*")
        if path.is_file()
        and path.suffix.lower() in SUFFIXES
        and BROKEN not in path.parents
    )


def test_the_fixture_directory_exists():
    """Guard against this whole file skipping into meaninglessness."""
    if not FIXTURES.is_dir():
        pytest.skip(
            f"{FIXTURES} does not exist. The sweep cannot run. This is a skip, "
            f"not a pass -- see the module docstring."
        )
    assert valid_fixtures(), (
        f"{FIXTURES} exists but contains no parseable fixture with a suffix in "
        f"{sorted(SUFFIXES)}. Either the suffix list is stale or the fixtures moved; "
        f"either way the sweep is silently checking nothing."
    )


@pytest.mark.parametrize("fixture", valid_fixtures(), ids=lambda p: p.name)
def test_valid_fixture_normalises_without_violation(fixture: Path):
    """Zero false positives. A refusal here means a rule is refusing correct data.

    If a fixture legitimately fails, it is not a valid fixture -- move it to
    `tests/fixtures/broken/` and assert on its specific error there. Do not
    relax the rule to make this pass.
    """
    tokenized = tokenize(fixture)
    detected = dialects.detect(tokenized)
    try:
        result = normalise(
            tokenized.features, detected.profile, comments=tokenized.comments
        )
    except NormaliseError as exc:
        pytest.fail(
            f"{fixture.relative_to(FIXTURES)} is a valid fixture but was refused:\n"
            f"    {exc}\n"
            f"Either a rule is refusing correct data, or this fixture belongs in "
            f"tests/fixtures/broken/."
        )
    assert result.feature_count > 0, f"{fixture.name} normalised to zero features"
