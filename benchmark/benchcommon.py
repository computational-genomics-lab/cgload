"""Small helpers the benchmark's runners share. Standard library only.

Nothing here judges a result; it runs commands, measures them, and appends rows
to the result tables, so the runners stay short and every table is written the
same way: tab-separated, one header, one row per genome and tool, appended as
each one finishes, so an interrupted run keeps what it finished.
"""

from __future__ import annotations

import contextlib
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MEASURE = HERE / "measure.py"
#: The genome list; BENCH_GENOMES points elsewhere for the benchmark's own tests.
GENOMES = Path(os.environ.get("BENCH_GENOMES", HERE / "genomes.tsv"))


def genomes() -> list[dict]:
    lines = [line for line in GENOMES.read_text().splitlines()
             if line.strip() and not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def genome(accession: str) -> dict:
    for row in genomes():
        if row["accession"] == accession:
            return row
    raise SystemExit(f"{accession} is not in {GENOMES}")


def measured(name: str, command: list[str], steps: Path, *, timeout: float | None = None,
             env: dict | None = None) -> dict:
    """Run ``command`` under measure.py; its stdout, stderr and measurement go to
    steps/NAME.out, .err and .json. Returns the measurement."""
    steps.mkdir(parents=True, exist_ok=True)
    record = steps / f"{name}.json"
    wrapper = [sys.executable, str(MEASURE), "--json", str(record),
               "--stdout", str(steps / f"{name}.out"), "--stderr", str(steps / f"{name}.err")]
    if timeout:
        wrapper += ["--timeout", str(timeout)]
    subprocess.run([*wrapper, "--", *command], check=False,
                   env={**os.environ, **(env or {})})
    if not record.exists():
        return {"command": command, "exit_code": -1, "wall_s": None, "peak_rss_mb": None,
                "timed_out": False}
    return json.loads(record.read_text())


def text_of(steps: Path, name: str) -> str:
    parts = []
    for suffix in (".out", ".err"):
        path = steps / f"{name}{suffix}"
        if path.exists():
            parts.append(path.read_text(errors="replace"))
    return "\n".join(parts)


def first_lines(text: str, count: int = 3) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return " | ".join(lines[:count])[:600]


def append_row(table: Path, row: dict) -> None:
    """Append one row; write the header first if the table is new. Locked, as
    several runners may append to one table at once."""
    import gffmodel

    gffmodel.locked_append(table, "\t".join(row),
                           ["\t".join(_cell(value) for value in row.values())])


def _cell(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value).replace("\t", " ").replace("\n", " ")


def done(table: Path, **match) -> bool:
    """Whether the table already has a row with these values (for resuming)."""
    if not table.exists():
        return False
    with table.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if all(row.get(key) == str(value) for key, value in match.items()):
                return True
    return False


def tree_bytes(path: Path) -> int:
    total = 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            with contextlib.suppress(OSError):
                total += (Path(root) / name).stat().st_size
    return total


def safe_name(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


#: A verify tier line: two spaces, the name, the outcome, the detail.
TIER = re.compile(r"^  (\S.*?)\s{2,}(PASS|FAIL|----)\s{2}(.*)$")


def tiers(text: str) -> dict[str, str]:
    found = {}
    for line in text.splitlines():
        match = TIER.match(line)
        if match:
            name = re.sub(r"\s*\[.*\]$", "", match.group(1)).strip()
            outcome = {"----": "skipped"}.get(match.group(2), match.group(2).lower())
            # Several files can each have a tier; one failure fails the tier.
            if found.get(name) != "fail":
                found[name] = outcome
    return found


def git_commit(repo: Path) -> str:
    done_ = subprocess.run(["git", "-C", str(repo), "describe", "--always", "--dirty",
                            "--tags"], capture_output=True, text=True, check=False)
    return done_.stdout.strip() or "unknown"
