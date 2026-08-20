"""D-011 claims the roadmap and the CLI cannot drift. This is the test that
makes that true: it reads docs/roadmap.md rather than a hardcoded list.
"""

import re
from pathlib import Path

from cgload.cli.main import PENDING, cli

ROADMAP = Path(__file__).resolve().parents[1] / "docs" / "roadmap.md"


def _roadmap_milestones() -> set[int]:
    rows = re.findall(r"^\|\s*(\d+)\s*\|", ROADMAP.read_text(), flags=re.MULTILINE)
    return {int(n) for n in rows}


def test_every_pending_command_targets_a_real_milestone():
    milestones = _roadmap_milestones()
    assert milestones, "no milestone table parsed from docs/roadmap.md"
    for name, (milestone, _summary) in PENDING.items():
        assert milestone in milestones, f"{name} points at milestone {milestone}, not in roadmap"


def test_pending_table_and_registered_commands_agree():
    """Once a command is implemented it must be deleted from PENDING in the
    same commit, so anything registered and not pending is real code."""
    assert set(PENDING) <= set(cli.commands)
