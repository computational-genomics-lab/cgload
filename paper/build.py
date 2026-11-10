#!/usr/bin/env python3
"""Build the manuscript for review: paper/build/manuscript.{md,docx}.

    python3 paper/build.py

Each line of manuscript.md of the form <!-- include tables/NAME.md --> is
replaced by that table, as paper/figures.py wrote it, and pandoc makes a Word
file with the figures in it. Run paper/figures.py first.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INCLUDE = re.compile(r"^<!-- include (tables/[\w.-]+\.md) -->$")
# Notes to the authors, which pandoc would otherwise drop with the comment.
NOTE = re.compile(r"<!--\s*(NOT FOR SUBMISSION.*?|AUTHORS:.*?)-->", re.S)


def assemble(source: Path) -> str:
    lines = []
    for line in source.read_text().splitlines():
        match = INCLUDE.match(line.strip())
        if match:
            table = HERE / match.group(1)
            if not table.exists():
                raise SystemExit(f"{table} is missing: run paper/figures.py first")
            lines.append(table.read_text().rstrip("\n"))
        else:
            lines.append(line)
    text = "\n".join(lines) + "\n"

    def note(match: re.Match) -> str:
        said = " ".join(match.group(1).removeprefix("AUTHORS:").split())
        return f"\n**[TO THE AUTHORS: {said}]**\n"

    return NOTE.sub(note, text)


def main() -> int:
    out = HERE / "build"
    out.mkdir(exist_ok=True)
    text = assemble(HERE / "manuscript.md")
    built = out / "manuscript.md"
    built.write_text(text)
    pandoc = shutil.which("pandoc")
    if pandoc is None:
        print(f"wrote {built}; pandoc is not installed, so no Word file", file=sys.stderr)
        return 0
    # Fixed arguments: the program found above and paths under paper/.
    subprocess.run([pandoc, str(built), "--from", "markdown", "--to", "docx",  # noqa: S603
                    "--resource-path", str(HERE), "--output", str(out / "manuscript.docx")],
                   check=True)
    print(f"wrote {built} and {out / 'manuscript.docx'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
