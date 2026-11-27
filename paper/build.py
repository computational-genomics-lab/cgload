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
* puts the measured numbers in place of each `{{portal:KEY}}` and
  `{{blast:KEY}}`, from tables/portal_numbers.json and
  tables/blast_numbers.json, which paper/figures.py writes from the results.
  Until those results exist a file is empty, and each number is shown as
  **[portal:KEY]**; once it has numbers, a KEY it lacks stops the build;
* makes a Word file with pandoc, with typographic substitution off, so that
  `--force` stays two hyphens, and gives it what BMC Bioinformatics asks of a
  submission: A4 pages, double line spacing in the text, line numbers and page
  numbers. Figures are not embedded: the journal takes them as separate files,
  and their legends are a section of the manuscript.
"""

from __future__ import annotations

import html
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
from xml.dom import minidom

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


MEASURED = ("portal", "blast", "blast1")


def measured_numbers(text: str) -> str:
    for kind in MEASURED:
        pattern = re.compile(r"\{\{" + kind + r":([a-z0-9_]+)\}\}")
        path = HERE / "tables" / f"{kind}_numbers.json"
        numbers = json.loads(path.read_text()) if path.exists() else {}
        wanted = set(pattern.findall(text))
        if numbers and wanted - numbers.keys():
            raise SystemExit(f"no {kind} number for "
                             f"{', '.join(sorted(wanted - numbers.keys()))}")
        if not numbers and wanted:
            print(f"the {kind} results are pending: {len(wanted)} numbers shown by name",
                  file=sys.stderr)
        text = pattern.sub(lambda m, numbers=numbers, kind=kind:
                           numbers.get(m.group(1), f"**[{kind}:{m.group(1)}]**"), text)
    return text


def assemble(source: Path) -> str:
    text = measured_numbers(include_tables(source.read_text()))

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


TEXT_WIDTH = 11906 - 2 * 1440   # twips: A4 less the margins of SECTION
TABLE_POINTS = 10               # type size of text in tables; code is set a point smaller
CELL_ROOM = 2 * 108 + 60        # twips: the Table style's cell margins, and some slack
SHORT = 1080                    # twips: a cell whose text fits in this stays on one line
ZWSP = chr(0x200B)              # a break opportunity without a space (Table 1's headers)
BREAKS = re.compile(f"([\\s{ZWSP}]+)")
TABLE_TEXT_STYLE = ('<w:style w:type="paragraph" w:customStyle="1" w:styleId="TableText">'
                    '<w:name w:val="Table Text"/><w:basedOn w:val="Compact"/><w:qFormat/>'
                    f'<w:rPr><w:sz w:val="{2 * TABLE_POINTS}"/>'
                    f'<w:szCs w:val="{2 * TABLE_POINTS}"/></w:rPr></w:style>')


class Measure:
    """Width in twips of text as Word sets it in a table: Cambria, bold in the header
    row, and Consolas for code. Metric-compatible fonts (Caladea) or wider ones
    (DejaVu) stand in where those are missing, and without Pillow and matplotlib a
    width per character does."""

    SIZE = 200   # pixels at which the fonts are loaded

    def __init__(self) -> None:
        self.cache: dict[tuple[str, str], float] = {}
        try:
            from matplotlib import font_manager
            from PIL import ImageFont
        except ImportError:
            self.fonts = None
            return

        def find(families: list[str], weight: str):
            for family in families:
                properties = font_manager.FontProperties(family=family, weight=weight)
                try:
                    path = font_manager.findfont(properties, fallback_to_default=False)
                except ValueError:
                    continue
                return ImageFont.truetype(path, self.SIZE)
            return None

        serif = ["Cambria", "Caladea", "DejaVu Serif"]
        mono = ["Consolas", "Inconsolata", "DejaVu Sans Mono"]
        self.fonts = {"text": find(serif, "normal"), "bold": find(serif, "bold"),
                      "code": find(mono, "normal")}

    def __call__(self, text: str, kind: str) -> float:
        key = (text, kind)
        if key not in self.cache:
            points = TABLE_POINTS - 1 if kind == "code" else TABLE_POINTS
            font = self.fonts.get(kind) if self.fonts else None
            if font is None:
                ems = len(text) * (0.6 if kind != "text" else 0.55)
            else:
                ems = font.getlength(text) / self.SIZE
            self.cache[key] = ems * points * 20
        return self.cache[key]


def cell_widths(cell: str, header: bool, measure: Measure) -> tuple[float, float]:
    """The widest unbreakable piece of a cell's text, and the whole text on one line."""
    longest = line = token = 0.0
    for run in re.findall(r"<w:r>(.*?)</w:r>", cell, re.S):
        kind = ("code" if "VerbatimChar" in run
                else "bold" if header or re.search(r"<w:b ?/>", run) else "text")
        text = html.unescape("".join(re.findall(r"<w:t(?: [^>]*)?>(.*?)</w:t>", run, re.S)))
        for piece in BREAKS.split(text):
            if not piece:
                continue
            if BREAKS.fullmatch(piece):
                longest, token = max(longest, token), 0.0
                line += measure(" ", kind) if piece.strip(ZWSP) else 0.0
            else:
                width = measure(piece, kind)
                token += width
                line += width
    return max(longest, token), line


def fit_table(table: str, measure: Measure) -> str:
    """Column widths for one table, as Word's automatic layout would choose them:
    each column at least as wide as its widest word or number, and the rest of
    the text width shared out by how much wider each column's text would like to
    be. Text in the table is set at TABLE_POINTS."""
    words: list[float] = []    # the widest word or number of each column
    short: list[float] = []    # the same, or a short cell's whole text if wider
    most: list[float] = []     # the widest cell, on one line
    for number, row in enumerate(re.findall(r"<w:tr\b.*?</w:tr>", table, re.S)):
        for column, cell in enumerate(re.findall(r"<w:tc>.*?</w:tc>", row, re.S)):
            low, high = cell_widths(cell, number == 0, measure)
            if column == len(words):
                words.append(0.0)
                short.append(0.0)
                most.append(0.0)
            words[column] = max(words[column], low + CELL_ROOM)
            short[column] = max(short[column], max(low, min(high, SHORT)) + CELL_ROOM)
            most[column] = max(most[column], high + CELL_ROOM)
    least = short if sum(short) < TEXT_WIDTH else words
    if sum(most) <= TEXT_WIDTH:
        widths = [w * TEXT_WIDTH / sum(most) for w in most]
    elif sum(least) >= TEXT_WIDTH:
        print("note: a table is too wide for its words; Word will break some", file=sys.stderr)
        widths = [w * TEXT_WIDTH / sum(least) for w in least]
    else:
        spare = TEXT_WIDTH - sum(least)
        wants = [high - low for low, high in zip(least, most, strict=True)]
        widths = [low + spare * want / sum(wants)
                  for low, want in zip(least, wants, strict=True)]
    grid = "".join(f'<w:gridCol w:w="{round(w)}" />' for w in widths)
    table = re.sub(r"<w:tblGrid>.*?</w:tblGrid>", f"<w:tblGrid>{grid}</w:tblGrid>", table,
                   count=1, flags=re.S)
    table = table.replace('<w:pStyle w:val="Compact" />', '<w:pStyle w:val="TableText" />')
    code = f'<w:rStyle w:val="VerbatimChar" /><w:sz w:val="{2 * TABLE_POINTS - 2}" />'
    return table.replace('<w:rStyle w:val="VerbatimChar" />', code)


def fit_tables(document: str) -> str:
    measure = Measure()
    return re.sub(r"<w:tbl>.*?</w:tbl>", lambda m: fit_table(m.group(0), measure), document,
                  flags=re.S)


#: The order of child elements that the Office Open XML schema requires, for the
#: elements whose children pandoc writes in another order. Word opens such files
#: anyway; a strict reader need not.
SCHEMA_ORDER = {
    "w:settings": [
        "w:writeProtection", "w:view", "w:zoom", "w:removePersonalInformation",
        "w:removeDateAndTime", "w:doNotDisplayPageBoundaries", "w:displayBackgroundShape",
        "w:printPostScriptOverText", "w:printFractionalCharacterWidth", "w:printFormsData",
        "w:embedTrueTypeFonts", "w:embedSystemFonts", "w:saveSubsetFonts", "w:saveFormsData",
        "w:mirrorMargins", "w:alignBordersAndEdges", "w:bordersDoNotSurroundHeader",
        "w:bordersDoNotSurroundFooter", "w:gutterAtTop", "w:hideSpellingErrors",
        "w:hideGrammaticalErrors", "w:activeWritingStyle", "w:proofState", "w:formsDesign",
        "w:attachedTemplate", "w:linkStyles", "w:stylePaneFormatFilter",
        "w:stylePaneSortMethod", "w:documentType", "w:mailMerge", "w:revisionView",
        "w:trackRevisions", "w:doNotTrackMoves", "w:doNotTrackFormatting",
        "w:documentProtection", "w:autoFormatOverride", "w:styleLockTheme",
        "w:styleLockQFSet", "w:defaultTabStop", "w:autoHyphenation",
        "w:consecutiveHyphenLimit", "w:hyphenationZone", "w:doNotHyphenateCaps",
        "w:showEnvelope", "w:summaryLength", "w:clickAndTypeStyle", "w:defaultTableStyle",
        "w:evenAndOddHeaders", "w:bookFoldRevPrinting", "w:bookFoldPrinting",
        "w:bookFoldPrintingSheets", "w:drawingGridHorizontalSpacing",
        "w:drawingGridVerticalSpacing", "w:displayHorizontalDrawingGridEvery",
        "w:displayVerticalDrawingGridEvery", "w:doNotUseMarginsForDrawingGridOrigin",
        "w:drawingGridHorizontalOrigin", "w:drawingGridVerticalOrigin",
        "w:doNotShadeFormData", "w:noPunctuationKerning", "w:characterSpacingControl",
        "w:printTwoOnOne", "w:strictFirstAndLastChars", "w:noLineBreaksAfter",
        "w:noLineBreaksBefore", "w:savePreviewPicture", "w:doNotValidateAgainstSchema",
        "w:saveInvalidXml", "w:ignoreMixedContent", "w:alwaysShowPlaceholderText",
        "w:doNotDemarcateInvalidXml", "w:saveXmlDataOnly", "w:useXSLTWhenSaving",
        "w:saveThroughXslt", "w:showXMLTags", "w:alwaysMergeEmptyNamespace",
        "w:updateFields", "w:hdrShapeDefaults", "w:footnotePr", "w:endnotePr", "w:compat",
        "w:docVars", "w:rsids", "m:mathPr", "w:attachedSchema", "w:themeFontLang",
        "w:clrSchemeMapping", "w:doNotIncludeSubdocsInStats", "w:doNotAutoCompressPictures",
        "w:forceUpgrade", "w:captions", "w:readModeInkLockDown", "w:smartTagType",
        "sl:schemaLibrary", "w:shapeDefaults", "w:doNotEmbedSmartTags", "w:decimalSymbol",
        "w:listSeparator"],
    "w:style": [
        "w:name", "w:aliases", "w:basedOn", "w:next", "w:link", "w:autoRedefine",
        "w:hidden", "w:uiPriority", "w:semiHidden", "w:unhideWhenUsed", "w:qFormat",
        "w:locked", "w:personal", "w:personalCompose", "w:personalReply", "w:rsid",
        "w:pPr", "w:rPr", "w:tblPr", "w:trPr", "w:tcPr", "w:tblStylePr"],
    "w:pPr": [
        "w:pStyle", "w:keepNext", "w:keepLines", "w:pageBreakBefore", "w:framePr",
        "w:widowControl", "w:numPr", "w:suppressLineNumbers", "w:pBdr", "w:shd", "w:tabs",
        "w:suppressAutoHyphens", "w:kinsoku", "w:wordWrap", "w:overflowPunct",
        "w:topLinePunct", "w:autoSpaceDE", "w:autoSpaceDN", "w:bidi", "w:adjustRightInd",
        "w:snapToGrid", "w:spacing", "w:ind", "w:contextualSpacing", "w:mirrorIndents",
        "w:suppressOverlap", "w:jc", "w:textDirection", "w:textAlignment",
        "w:textboxTightWrap", "w:outlineLvl", "w:divId", "w:cnfStyle", "w:rPr", "w:sectPr",
        "w:pPrChange"],
    "w:rPr": [
        "w:ins", "w:del", "w:moveFrom", "w:moveTo", "w:rStyle", "w:rFonts", "w:b", "w:bCs",
        "w:i", "w:iCs", "w:caps", "w:smallCaps", "w:strike", "w:dstrike", "w:outline",
        "w:shadow", "w:emboss", "w:imprint", "w:noProof", "w:snapToGrid", "w:vanish",
        "w:webHidden", "w:color", "w:spacing", "w:w", "w:kern", "w:position", "w:sz",
        "w:szCs", "w:highlight", "w:u", "w:effect", "w:bdr", "w:shd", "w:fitText",
        "w:vertAlign", "w:rtl", "w:cs", "w:em", "w:lang", "w:eastAsianLayout",
        "w:specVanish", "w:oMath", "w:rPrChange"],
    "w:tblPr": [
        "w:tblStyle", "w:tblpPr", "w:tblOverlap", "w:bidiVisual", "w:tblStyleRowBandSize",
        "w:tblStyleColBandSize", "w:tblW", "w:jc", "w:tblCellSpacing", "w:tblInd",
        "w:tblBorders", "w:shd", "w:tblLayout", "w:tblCellMar", "w:tblLook", "w:tblCaption",
        "w:tblDescription", "w:tblPrChange"],
    "w:tcPr": [
        "w:cnfStyle", "w:tcW", "w:gridSpan", "w:hMerge", "w:vMerge", "w:tcBorders", "w:shd",
        "w:noWrap", "w:tcMar", "w:textDirection", "w:tcFitText", "w:vAlign", "w:hideMark",
        "w:headers", "w:tcPrChange"],
}


def schema_order(xml: str) -> str:
    """Put the children of SCHEMA_ORDER's elements in the order the schema requires,
    drop stray text between them, and give list identifiers their eight hex digits."""
    # The XML is pandoc's output from our own manuscript, not outside input.
    dom = minidom.parseString(xml.encode("utf-8"))  # noqa: S318
    for tag, order in SCHEMA_ORDER.items():
        rank = {name: i for i, name in enumerate(order)}
        for element in dom.getElementsByTagName(tag):
            children = [node for node in element.childNodes
                        if node.nodeType == node.ELEMENT_NODE]
            for node in list(element.childNodes):
                element.removeChild(node)
            for node in sorted(children, key=lambda node: rank.get(node.tagName, len(order))):
                element.appendChild(node)
    for nsid in dom.getElementsByTagName("w:nsid"):
        nsid.setAttribute("w:val", nsid.getAttribute("w:val").zfill(8))
    return dom.toxml(encoding="UTF-8").decode("utf-8")


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
    parts["word/document.xml"] = fit_tables(document).encode("utf-8")
    edit("word/document.xml", "</w:body>", SECTION + "</w:body>")
    # Body Text carries the prose; lists and references are Compact, and tables
    # Table Text, a smaller Compact.
    edit("word/styles.xml", '<w:spacing w:after="180" w:before="180" />',
         '<w:spacing w:after="180" w:before="180" w:line="480" w:lineRule="auto" />')
    edit("word/styles.xml", '<w:spacing w:after="36" w:before="36" />',
         '<w:spacing w:after="36" w:before="36" w:line="240" w:lineRule="auto" />')
    edit("word/styles.xml", "</w:styles>", TABLE_TEXT_STYLE + "</w:styles>")
    parts["word/footer1.xml"] = FOOTER.encode("utf-8")
    edit("word/_rels/document.xml.rels", "</Relationships>",
         '<Relationship Id="rIdPageFooter" Type="http://schemas.openxmlformats.org/officeDocument/'
         '2006/relationships/footer" Target="footer1.xml"/></Relationships>')
    edit("[Content_Types].xml", "</Types>",
         '<Override PartName="/word/footer1.xml" ContentType="application/vnd.openxmlformats-'
         'officedocument.wordprocessingml.footer+xml"/></Types>')
    for name in ("word/document.xml", "word/styles.xml", "word/settings.xml",
                 "word/numbering.xml"):
        if name in parts:
            parts[name] = schema_order(parts[name].decode("utf-8")).encode("utf-8")
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
              "llm_use.md": "additional_file_5_llm_use.md",
              "tables/tableS5_portal.csv": "additional_file_6_table_S5.csv",
              "tables/tableS6_blast.csv": "additional_file_8_table_S6.csv"}
    for source, name in extras.items():
        if source.endswith(".md"):
            # Notes to the authors are HTML comments; the submitted file leaves them
            # out, and the build says which are still open.
            text = (HERE / source).read_text()
            for note in re.findall(r"<!--\s*(AUTHORS.*?)-->", text, re.S):
                print(f"note to the authors in {source}: {' '.join(note.split())}",
                      file=sys.stderr)
            (out / name).write_text(re.sub(r"\n{3,}", "\n\n", COMMENT.sub("", text)))
        else:
            shutil.copyfile(HERE / source, out / name)
    for figure in sorted((HERE / "figures").glob("fig*.pdf")):
        shutil.copyfile(figure, out / figure.name)
    print(f"wrote {built}, {out / 'manuscript.docx'}, {len(extras)} additional files "
          "and the figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
