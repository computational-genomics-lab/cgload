#!/usr/bin/env python3
"""Build the manuscript for review: paper/build/manuscript.{md,docx}.

    python3 paper/build.py

manuscript.md is the source. This script

* replaces each line `<!-- include tables/NAME.md -->` by that table, as
  paper/figures.py wrote it (run that first);
* numbers the citations. The text cites by key, `[@gff3]` or
  `[@tabix; @htslib]`, and the reference list at the end defines each key on a
  line `[@key]: reference`. References are numbered in the order they are
  first cited, as the journal asks, so a number can never point at the wrong
  reference. A key that is cited and not defined, or defined and never cited,
  stops the build;
* shows each note to the authors (`<!-- AUTHORS: ... -->`) in bold, where
  pandoc would drop the comment;
* makes a Word file with pandoc, with typographic substitution off, so that
  `--force` stays two hyphens.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
INCLUDE = re.compile(r"^<!-- include (tables/[\w.-]+\.md) -->$")
NOTE = re.compile(r"<!--\s*AUTHORS:(.*?)-->", re.S)
COMMENT = re.compile(r"<!--.*?-->", re.S)
DEFINITION = re.compile(r"^\[@([\w-]+)\]:\s+(.+)$")
CITATION = re.compile(r"\[(@[\w-]+(?:\s*;\s*@[\w-]+)*)\]")


def include_tables(text: str) -> str:
    lines = []
    for line in text.splitlines():
        match = INCLUDE.match(line.strip())
        if match:
            table = HERE / match.group(1)
            if not table.exists():
                raise SystemExit(f"{table} is missing: run paper/figures.py first")
            lines.append(table.read_text().rstrip("\n"))
        else:
            lines.append(line)
    return "\n".join(lines) + "\n"


def number_citations(text: str) -> str:
    """Replace keys by numbers in order of first citation, and write the list."""
    references: dict[str, str] = {}
    body = []
    for line in text.splitlines():
        match = DEFINITION.match(line)
        if match:
            if match.group(1) in references:
                raise SystemExit(f"reference [@{match.group(1)}] is defined twice")
            references[match.group(1)] = match.group(2).strip()
        else:
            body.append(line)
    order: list[str] = []

    def numbered(match: re.Match) -> str:
        keys = [key.strip().lstrip("@") for key in match.group(1).split(";")]
        for key in keys:
            if key not in references:
                raise SystemExit(f"[@{key}] is cited but not defined")
            if key not in order:
                order.append(key)
        return "[" + ", ".join(str(order.index(key) + 1) for key in keys) + "]"

    text = CITATION.sub(numbered, "\n".join(body))
    unused = sorted(set(references) - set(order))
    if unused:
        raise SystemExit("defined but never cited: " + ", ".join(f"[@{k}]" for k in unused))
    listed = "\n".join(f"{i}. {references[key]}" for i, key in enumerate(order, 1))
    marker = "{{references}}"
    if text.count(marker) != 1:
        raise SystemExit(f"manuscript.md needs exactly one line {marker}")
    return text.replace(marker, listed) + "\n"


def assemble(source: Path) -> str:
    text = include_tables(source.read_text())

    def note(match: re.Match) -> str:
        return f"\n**[TO THE AUTHORS: {' '.join(match.group(1).split())}]**\n"

    text = NOTE.sub(note, text)
    text = COMMENT.sub("", text)
    return number_citations(text)


def main() -> int:
    out = HERE / "build"
    out.mkdir(exist_ok=True)
    built = out / "manuscript.md"
    built.write_text(assemble(HERE / "manuscript.md"))
    pandoc = shutil.which("pandoc")
    if pandoc is None:
        print(f"wrote {built}; pandoc is not installed, so no Word file", file=sys.stderr)
        return 0
    # Fixed arguments: the program found above and paths under paper/.
    subprocess.run([pandoc, str(built), "--from", "markdown-smart", "--to", "docx",  # noqa: S603
                    "--resource-path", str(HERE), "--output", str(out / "manuscript.docx")],
                   check=True)
    # The additional files, under the names the manuscript gives them.
    extras = {"tables/tableS1_genomes.csv": "additional_file_1_table_S1.csv",
              "tables/tableS2_faults.csv": "additional_file_2_table_S2.csv",
              "tables/tableS3_timing.csv": "additional_file_3_table_S3.csv",
              "tables/tableS4_sequences.csv": "additional_file_4_table_S4.csv",
              "llm_use.md": "additional_file_5_llm_use.md"}
    for source, name in extras.items():
        shutil.copyfile(HERE / source, out / name)
    for figure in sorted((HERE / "figures").glob("fig*.pdf")):
        shutil.copyfile(figure, out / figure.name)
    print(f"wrote {built}, {out / 'manuscript.docx'}, {len(extras)} additional files "
          "and the figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
