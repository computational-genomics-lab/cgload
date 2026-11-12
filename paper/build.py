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
  `--force` stays two hyphens, and gives it what BMC Bioinformatics asks of a
  submission: A4 pages, double line spacing in the text, line numbers and page
  numbers. Figures are not embedded: the journal takes them as separate files,
  and their legends are a section of the manuscript.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import zipfile
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


W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
FOOTER = (f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:ftr xmlns:w="{W}">'
          '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:fldChar w:fldCharType="begin"/></w:r>'
          '<w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
          '<w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>')
SECTION = ('<w:sectPr><w:footerReference w:type="default" r:id="rIdPageFooter"/>'
           '<w:pgSz w:w="11906" w:h="16838"/>'
           '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
           'w:header="720" w:footer="720" w:gutter="0"/>'
           '<w:lnNumType w:countBy="1" w:restart="continuous"/></w:sectPr>')


def journal_layout(docx: Path) -> None:
    """Double spacing in body text, line numbers and page numbers, on A4."""
    with zipfile.ZipFile(docx) as archive:
        parts = {name: archive.read(name) for name in archive.namelist()}

    def edit(name: str, old: str, new: str) -> None:
        text = parts[name].decode("utf-8")
        if text.count(old) != 1:
            raise SystemExit(f"{name}: expected one {old!r} in pandoc's output")
        parts[name] = text.replace(old, new).encode("utf-8")

    document = re.sub(r"<w:sectPr\b.*?</w:sectPr>|<w:sectPr\b[^>]*/>", "",
                      parts["word/document.xml"].decode("utf-8"), flags=re.S)
    parts["word/document.xml"] = document.encode("utf-8")
    edit("word/document.xml", "</w:body>", SECTION + "</w:body>")
    # Body Text carries the prose; tables, lists and references are Compact.
    edit("word/styles.xml", '<w:spacing w:after="180" w:before="180" />',
         '<w:spacing w:after="180" w:before="180" w:line="480" w:lineRule="auto" />')
    edit("word/styles.xml", '<w:spacing w:after="36" w:before="36" />',
         '<w:spacing w:after="36" w:before="36" w:line="240" w:lineRule="auto" />')
    parts["word/footer1.xml"] = FOOTER.encode("utf-8")
    edit("word/_rels/document.xml.rels", "</Relationships>",
         '<Relationship Id="rIdPageFooter" Type="http://schemas.openxmlformats.org/officeDocument/'
         '2006/relationships/footer" Target="footer1.xml"/></Relationships>')
    edit("[Content_Types].xml", "</Types>",
         '<Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-'
         'officedocument.wordprocessingml.footer+xml"/></Types>')
    with zipfile.ZipFile(docx, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, data in parts.items():
            archive.writestr(name, data)


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
    journal_layout(out / "manuscript.docx")
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
