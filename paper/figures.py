#!/usr/bin/env python3
"""Make the paper's figures and tables from the benchmark results.

    python3 paper/figures.py [RESULTS] [--out paper]

RESULTS is a results folder as benchmark/run.sh leaves it; by default the one
the paper reports, paper/data/benchmark-run-2.2. The portal's measurements
(PROTOCOL.md section 13) are read from --portal, by default
paper/data/portal-v0.3.0; while that folder does not exist, their tables say
"pending". Figures go to OUT/figures/ (PDF for the journal, PNG to look at),
tables to OUT/tables/ as Markdown, and the portal's numbers for the text to
OUT/tables/portal_numbers.json, which paper/build.py puts in place of each
{{portal:KEY}}. Every number is read from the results; none is typed in here.
Figure 3, the portal's pages, is made by paper/screenshots.py. Needs matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sqlite3
import statistics
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "src"))
from cgload.db.schema import metadata as SCHEMA  # noqa: E402

# Colours: the light-mode reference palette of the dataviz method (slots 1 to 3,
# checked together on white: worst colour-blind distance 9.2, normal 24.0).
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
NOTHING, FAILED = "#f0efec", "#d6d5cf"
RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
        "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

FULL_WIDTH = 170 / 25.4  # inches: a full page width of BMC Bioinformatics, 170 mm
SMALL = 7  # points: the journal's smallest text

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": AXIS, "ytick.color": AXIS,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5, "grid.linestyle": "-",
    "axes.axisbelow": True, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
})

TOOL_NAMES = {
    "cgload-sqlite": "cgload (SQLite)", "cgload-postgres": "cgload (PostgreSQL)",
    "gffutils-default": "gffutils, default", "gffutils": "gffutils, create_unique",
    "txdbmaker": "txdbmaker", "chado": "Chado (python-chado)", "tabix": "tabix",
}
CLADES = [("vertebrate", "Vertebrates"), ("insect", "Insects"), ("plant", "Plants"),
          ("fungus", "Fungi")]


# --- reading ---------------------------------------------------------------------------

def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(text: str | None) -> int:
    return int(text) if text not in (None, "", "NA") else 0


def species(organism: str) -> str:
    """Genus, species and a subspecies (Canis lupus familiaris); the strain is in
    the genome table."""
    words = organism.split()
    keep = 3 if len(words) > 2 and words[2].isalpha() and words[2].islower() else 2
    return " ".join(words[:keep])


def short_species(organism: str) -> str:
    """H. sapiens, for labels where the full name does not fit."""
    words = organism.split()
    return f"{words[0][0]}. {words[1]}"


class Results:
    def __init__(self, folder: Path):
        self.folder = folder
        self.genomes = read_tsv(folder / "tables" / "scale.tsv")
        self.order = sorted(self.genomes, key=lambda g: (
            [c for c, _ in CLADES].index(g["clade"]), -number(g["rows"])))
        self.loaders = read_tsv(folder / "compare_loaders.tsv")
        self.tools = {row["tool"]: row for row in read_tsv(folder / "tables" / "tools.tsv")}
        self.cases = read_tsv(folder / "tables" / "cases.tsv")
        self.faults = read_tsv(folder / "tables" / "faults.tsv")
        self.fault_trials = read_tsv(folder / "faults.tsv")
        self.kills = read_tsv(folder / "tables" / "kills.tsv")
        self.agreement = read_tsv(folder / "tables" / "agreement.tsv")
        self.timing = read_tsv(folder / "tables" / "timing.tsv")
        self.scale = {(row["accession"], row["backend"]): row
                      for row in read_tsv(folder / "cgload_scale.tsv")}
        self.headline = json.loads((folder / "tables" / "headline.json").read_text())
        self.name = {g["accession"]: g["organism"] for g in self.genomes}

    def loader(self, accession: str, tool: str) -> dict:
        """One loader's comparison for one genome, in compare_loaders.tsv's columns.

        cgload's own comparison is its round trip, kept in cgload_scale.tsv; there
        "changed" counts only the changes it did not declare."""
        if tool.startswith("cgload-"):
            row = self.scale[(accession, tool.removeprefix("cgload-"))]
            return {"outcome": row["outcome"], "truth": row["rt_truth"],
                    "kept": row["rt_kept"], "changed": row["rt_changed"],
                    "changed_declared": row["rt_changed_declared"],
                    "changed_silent": row["rt_changed_silent"], "lost": row["rt_lost"],
                    "lost_silent": row["rt_lost_silent"]}
        found = [r for r in self.loaders if r["accession"] == accession and r["tool"] == tool]
        if len(found) != 1:
            raise SystemExit(f"{accession} {tool}: {len(found)} rows in compare_loaders.tsv")
        return found[0]


# --- figure helpers --------------------------------------------------------------------

def save(figure, out: Path, name: str) -> None:
    out.mkdir(parents=True, exist_ok=True)
    # No dates in the files, so that drawing again changes nothing.
    figure.savefig(out / f"{name}.pdf", metadata={"CreationDate": None})
    figure.savefig(out / f"{name}.png", dpi=300, metadata={"Software": None})
    plt.close(figure)


def hex_rgb(colour: str) -> tuple[float, float, float]:
    colour = colour.lstrip("#")
    return tuple(int(colour[i:i + 2], 16) / 255 for i in (0, 2, 4))


def luminance(colour: str) -> float:
    def channel(c):
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (channel(c) for c in hex_rgb(colour))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def ink_on(fill: str) -> str:
    """White or near-black, whichever reads better on the fill."""
    light = luminance(fill)
    return INK if (light + 0.05) / 0.05 > 1.05 / (light + 0.05) else "#ffffff"


def ramp_colour(fraction: float) -> str:
    """A step of the one-hue ramp, light for little, dark for much."""
    fraction = min(max(fraction, 0.0), 1.0)
    position = fraction * (len(RAMP) - 1)
    low = int(position)
    high = min(low + 1, len(RAMP) - 1)
    weight = position - low
    a, b = hex_rgb(RAMP[low]), hex_rgb(RAMP[high])
    mixed = [round((x + (y - x) * weight) * 255) for x, y in zip(a, b, strict=True)]
    return "#" + "".join(f"{c:02x}" for c in mixed)


def log_seconds(axis, which: str = "x") -> None:
    def label(value, _position):
        if value < 1:
            milliseconds = value * 1000
            return f"{milliseconds:g} ms"
        return f"{value:g} s"
    formatter = matplotlib.ticker.FuncFormatter(label)
    target = axis.xaxis if which == "x" else axis.yaxis
    target.set_major_formatter(formatter)


def thousands(value, _position) -> str:
    if value >= 1e6:
        return f"{value / 1e6:g} M"
    if value >= 1e3:
        return f"{value / 1e3:g} k"
    return f"{value:g}"


# --- Figure 1: how cgload loads and verifies -----------------------------------------

CHECK = RAMP[0]   # boxes that check, in the one-hue ramp's lightest step
STEP = NOTHING    # boxes that load


def box(axis, x, y, w, h, title, body, fill):
    from matplotlib.patches import FancyBboxPatch
    axis.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.05",
                                  facecolor=fill, edgecolor="none"))
    axis.text(x + 0.08, y + h - 0.08, title, ha="left", va="top", fontsize=7,
              fontweight="bold", color=INK)
    axis.text(x + 0.08, y + h - 0.27, body, ha="left", va="top", fontsize=SMALL,
              color=INK2, linespacing=1.35)


def arrow(axis, start, end, colour=MUTED):
    axis.annotate("", xy=end, xytext=start,
                  arrowprops={"arrowstyle": "-|>", "color": colour, "linewidth": 0.8,
                              "shrinkA": 0, "shrinkB": 0, "mutation_scale": 7})


def figure_design(out: Path) -> None:
    """Figure 1: every step of a load and of `verify`, with what stops each and
    what each writes. The step numbers 1 to 9 are those of the paper's list."""
    width, height = FULL_WIDTH, 7.75
    figure = plt.figure(figsize=(width, height))
    axis = figure.add_axes((0, 0, 1, 1))
    axis.set_xlim(0, width)
    axis.set_ylim(0, height)
    axis.axis("off")
    line = 0.106                       # inches per line: 7-point text, 109% leading
    pad_top, pad_bottom = 0.046, 0.036
    stop_fill = "#fdeee8"           # a light tint of ORANGE, for what stops a load
    left, right = 0.04, width - 0.17   # the right margin carries one connector
    corridor = width - 0.085

    renderer = figure.canvas.get_renderer()
    pad = 0.07                         # inches between a box's edge and its text

    def wrap(text: str, room: float, weight: str = "normal") -> list[str]:
        """Break ``text`` into lines no wider than ``room`` inches, as drawn."""
        font = FontProperties(family="sans-serif", size=SMALL, weight=weight)
        lines: list[str] = []
        for word in text.split():
            trial = f"{lines[-1]} {word}" if lines else word
            wide, _, _ = renderer.get_text_width_height_descent(trial, font, ismath=False)
            if lines and wide / figure.dpi <= room:
                lines[-1] = trial
            else:
                lines.append(word)
        return lines

    def rect(x, top, w, h, fill, edge=None):
        axis.add_patch(Rectangle((x, top - h), w, h, facecolor=fill,
                                 edgecolor=edge or "none", linewidth=0.6))

    def write(x, top, rows):
        """rows: (text, colour, weight); one line each, from the top down."""
        y = top - pad_top
        for text, colour, weight in rows:
            if text:
                axis.text(x, y, text, ha="left", va="top", fontsize=SMALL, color=colour,
                          fontweight=weight)
            y -= line
        return y

    def heading(y, label, title):
        axis.text(left, y, label, ha="left", va="top", fontsize=9, fontweight="bold",
                  color=INK)
        axis.text(left + 0.2, y - 0.01, title, ha="left", va="top", fontsize=8,
                  fontweight="bold", color=INK)

    def flow(x, top, bottom):
        axis.annotate("", xy=(x, bottom), xytext=(x, top),
                      arrowprops={"arrowstyle": "-|>", "color": MUTED, "linewidth": 0.8,
                                  "shrinkA": 0, "shrinkB": 0, "mutation_scale": 7})

    # --- a: the load ---------------------------------------------------------------
    y = height - 0.05
    heading(y, "a", "cgload load: from the input files to a committed genome")
    y -= 0.25

    inputs = [
        ("Configuration (TOML)", "organism, strain, assembly and the files to load"),
        ("Assembly", "FASTA, one record per sequence"),
        ("Functional annotation", "eggNOG-mapper or InterProScan output (optional)"),
        ("Annotation", "GFF3 (RefSeq, funannotate, AUGUSTUS) or GenBank"),
    ]
    gap = 0.08
    w_in = (right - left - 3 * gap) / 4
    h_in = pad_top + pad_bottom + 3 * line
    for i, (title, body) in enumerate(inputs):
        x = left + i * (w_in + gap)
        rect(x, y, w_in, h_in, "white", AXIS)
        write(x + pad, y, [(title, INK, "bold")]
              + [(text, INK2, "normal") for text in wrap(body, w_in - 2 * pad)])
    annotation_box = (left + 3 * (w_in + gap), y, w_in, h_in)
    under = y - h_in - 0.06
    axis.plot([left + w_in / 2, left + w_in / 2, right - w_in / 2, right - w_in / 2],
              [y - h_in, under, under, y - h_in], color=MUTED, linewidth=0.8)
    flow(width / 2, under, under - 0.12)
    y -= h_in + 0.21

    # Phase 1: five checks that read the files and write nothing.
    axis.text(left, y, "Read and check every file. Nothing is written yet.", ha="left",
              va="top", fontsize=7.5, fontweight="bold", color=INK)
    y -= 0.18
    reads = [
        ("Tokenise", "Split each line into nine columns and attributes.",
         "a line without nine columns; a bad coordinate, strand or phase",
         "declared lengths"),
        ("Detect the dialect", "Match column 2 and the header with the five dialect "
         "profiles.", "no profile matches, or several do", "the dialect and its evidence"),
        ("Normalise", "Join rows that share an ID; order segments; map types; find "
         "parents.", "two parents; a missing or wrong-typed parent; one ID on two features",
         "how parents were found"),
        ("Check parent spans", "A child must lie inside its parent; no row may run "
         "backwards.", "a violation, unless --repair-envelopes fixes the parent",
         "each repair, old and new"),
        ("Check across files", "Compare the files of a genome; check AUGUSTUS proteins.",
         "an ID in two files; a file listed twice; a protein that does not fit its CDS",
         "types no profile declares"),
    ]
    w_r = (right - left - 4 * gap) / 5
    room = w_r - 2 * pad
    blocks = []
    for title, body, stops, records in reads:
        rows = [(text, INK, "bold") for text in wrap(title, room, "bold")]
        rows += [(text, INK, "normal") for text in wrap(body, room)]
        rows += [("Stops the load if", ORANGE, "bold")]
        rows += [(text, INK, "normal") for text in wrap(stops, room)]
        rows += [("Records", INK2, "bold")]
        rows += [(text, INK2, "normal") for text in wrap(records, room)]
        blocks.append(rows)
    h_r = pad_top + pad_bottom + max(len(rows) for rows in blocks) * line
    for i, rows in enumerate(blocks):
        x = left + i * (w_r + gap)
        rect(x, y, w_r, h_r, CHECK)
        write(x + pad, y, rows)
        if i:
            arrow(axis, (x - gap, y - pad_top - line / 2), (x, y - pad_top - line / 2))
    y -= h_r + 0.06
    axis.text(left, y, "A stop names the file, the line and, where one exists, the option "
              "that would help. All files are read and checked before anything is "
              "written.", ha="left", va="top", fontsize=SMALL, color=INK2)
    flow(width / 2, y - 0.12, y - 0.23)
    y -= 0.25

    # Phase 2: the nine steps of the transaction.
    axis.text(left, y, "One transaction: rows into the database, files into a staging "
              "folder", ha="left", va="top", fontsize=7.5, fontweight="bold", color=INK)
    y -= 0.18
    w_stop, w_write = 1.78, 1.45
    x_step = left + w_stop + 0.1
    x_write = right - w_write
    w_step = x_write - 0.1 - x_step
    rail = x_step + 0.11                       # the line that joins the step numbers
    x_box = x_step + 0.24
    w_box = x_step + w_step - x_box
    for x, label, colour in ((left, "Stops the load if", ORANGE), (x_box, "Step", INK),
                             (x_write, "Writes", INK)):
        axis.text(x, y, label, ha="left", va="top", fontsize=SMALL, fontweight="bold",
                  color=colour)
    y -= 0.15
    steps = [
        ("Stage the assembly", "Copy the FASTA to a staging folder; index it with pyfaidx.",
         "a duplicate or empty sequence", "assembly.fa, assembly.fa.fai"),
        ("Write the genome rows", "Find or add the organism and strain; add the assembly "
         "with its checksum, and a row per file and per sequence.",
         "the assembly is already loaded, unless --force",
         "organism, strain, assembly, source_file, sequence_region"),
        ("Write the features", "Bulk insert: one row per segment, each with its bin.",
         "a sequence the FASTA lacks; a row past the end of a linear sequence",
         "feature"),
        ("Link the parents", "Second pass: point each child at its parent's first row.",
         "a stated parent with no stored row", "feature.parent_id"),
        ("Write the attributes", "Column 9 except ID and Parent, one row per value.",
         "", "feature_attribute; proteins.faa (AUGUSTUS)"),
        ("Attach functional annotation", "Match each protein name to one feature: by ID, "
         "or by protein_id, transcript_id, locus_tag and others.",
         "a name fitting two features; too few matches (default: all)",
         "annotation_run, protein_annotation"),
        ("End-of-load checks", "Compare the stored counts with the independent count: "
         "features by type, rows, and Parent= links.", "any count that differs", None),
        ("Write out", "PostgreSQL's planner statistics; the annotation as GFF3, "
         "BGZF-compressed and tabix-indexed by cgload; the manifest.", "",
         "annotation.gff3.gz, .tbi, manifest.json"),
        ("Commit", "Rename the staging folder into place, then commit.",
         "the commit fails (the folder is then moved back)",
         "the assembly's folder in the data folder"),
    ]
    first_rail = None
    for number, (title, body, stops, writes) in enumerate(steps, 1):
        body_rows = wrap(body, w_box - 2 * pad)
        stop_rows = wrap(stops, w_stop - 2 * pad) if stops else []
        write_rows = wrap(writes, w_write - 2 * pad) if writes else []
        n = max(1 + len(body_rows), len(stop_rows), len(write_rows), 2)
        h = pad_top + pad_bottom + n * line
        check = number == 7
        rect(x_box, y, w_box, h, CHECK if check else STEP)
        write(x_box + pad, y, [(title, INK, "bold")]
              + [(text, INK, "normal") for text in body_rows])
        mid = y - pad_top - line / 2
        if first_rail is None:
            first_rail = mid
        axis.add_patch(matplotlib.patches.Circle((rail, mid), 0.085, facecolor="white",
                                                 edgecolor=INK2, linewidth=0.6, zorder=3))
        axis.text(rail, mid, str(number), ha="center", va="center", fontsize=SMALL,
                  fontweight="bold", color=INK, zorder=4)
        last_rail = mid
        if stop_rows:
            rect(left, y, w_stop, h, stop_fill)
            write(left + pad, y, [(text, INK, "normal") for text in stop_rows])
            arrow(axis, (rail - 0.085, mid), (left + w_stop, mid), colour=ORANGE)
        if write_rows:
            rect(x_write, y, w_write, h, "white", AXIS)
            write(x_write + pad, y, [(text, INK2, "normal") for text in write_rows])
            arrow(axis, (x_box + w_box, mid), (x_write, mid))
        if check:
            # The independent count: separate code, reading the annotation file
            # itself, delivers the numbers that step 7 compares.
            rect(x_write, y, w_write, h, RAMP[2])
            write(x_write + pad, y, [("Independent count", INK, "bold")]
                  + [(text, INK, "normal")
                     for text in wrap("separate code reads each annotation file again",
                                      w_write - 2 * pad)])
            arrow(axis, (x_write, mid), (x_box + w_box, mid))
            ax, atop, aw, ah = annotation_box
            start = atop - ah / 2
            axis.plot([ax + aw, corridor, corridor], [start, start, y - 0.05],
                      color=MUTED, linewidth=0.8)
            arrow(axis, (corridor, y - 0.05), (x_write + w_write, y - 0.05))
        y -= h + 0.035
    axis.plot([rail, rail], [first_rail, last_rail], color=INK2, linewidth=0.6, zorder=2)
    rect(left, y, right - left, pad_top + pad_bottom + 2 * line, stop_fill)
    write(left + pad, y, [
        ("Any stop in steps 1 to 9: the transaction is rolled back and the staging folder "
         "deleted, so nothing changes.", INK, "normal"),
        ("A killed process: the database rolls back too, but the staging folder stays on "
         "disk, outside any assembly, until deleted by hand.", INK, "normal")])
    y -= pad_top + pad_bottom + 2 * line + 0.12

    # --- b: verify -----------------------------------------------------------------
    heading(y, "b", "cgload verify, at any later time: five tiers")
    y -= 0.22
    tiers = [
        ("Source count", "The independent count of the file equals the database's."),
        ("Parent links", "Every feature the file gives a parent has a stored link."),
        ("Stored assembly", "FASTA and index present, SHA-256 as at the load."),
        ("Cross-parser (optional)", "Counts per type in the GFF3 equal those of a "
         "GenBank file."),
        ("Reference (optional)", "Counts per type equal NCBI's published totals."),
    ]
    rows = [[(title, INK, "bold")] + [(text, INK, "normal") for text in wrap(body, room)]
            for title, body in tiers]
    # (the same column width and wrap as the five checks of panel a)
    h_t = pad_top + pad_bottom + max(len(r) for r in rows) * line
    for i, block in enumerate(rows):
        x = left + i * (w_r + gap)
        rect(x, y, w_r, h_t, CHECK)
        write(x + pad, y, block)
    y -= h_t + 0.06
    axis.text(left, y, "The tiers compare counts. They cannot tell which feature is "
              "which: a renamed feature, a wrong parent or a changed value keeps every "
              "count.", ha="left", va="top", fontsize=SMALL, color=INK2)
    y -= line
    axis.text(left, y, "Each tier prints PASS, FAIL, or why it could not run; verify exits "
              "with status 2 if a tier fails, and 0 otherwise, also when a tier is skipped.",
              ha="left", va="top", fontsize=SMALL, color=INK2)
    print(f"Figure 1: {y - line - 0.05:.2f} in to spare")
    save(figure, out, "fig1")


# --- Figure 2: the schema, drawn from cgload.db.schema ---------------------------------

#: Columns shown for each table, in the box's order. The figure fails if one of
#: them is not in the schema; columns left out are named in the legend's pointer
#: to src/cgload/db/schema.py.
SCHEMA_COLUMNS = {
    "organism": ["organism_id", "ncbi_taxon_id", "genus", "species"],
    "strain": ["strain_id", "organism_id", "name"],
    "assembly": ["assembly_id", "strain_id", "name", "accession", "fasta_path",
                 "fasta_sha256", "vocabulary_hash", "cgload_version"],
    "source_file": ["source_file_id", "assembly_id", "path", "sha256", "detected_dialect"],
    "sequence_region": ["sequence_region_id", "assembly_id", "source_id", "length",
                        "circular"],
    "feature": ["feature_id", "assembly_id", "sequence_region_id", "source_file_id",
                "source_id", "source_type", "feature_type", "start", "end", "strand",
                "phase", "segment_index", "parent_id", "linkage_method", "source_line",
                "bin"],
    "feature_attribute": ["feature_attribute_id", "feature_id", "key", "value"],
    "annotation_run": ["annotation_run_id", "assembly_id", "tool", "sha256", "query_count",
                       "matched_count"],
    "protein_annotation": ["protein_annotation_id", "annotation_run_id", "feature_id",
                           "analysis", "accession", "protein_start", "protein_end"],
    "vocabulary": ["term", "canonical_term", "tier"],
    "vocabulary_meta": ["version", "content_hash", "database_id"],
}
#: What one row of each table stands for, printed under its name.
SCHEMA_GRAIN = {
    "organism": "one row per species", "strain": "one row per strain",
    "assembly": "one row per loaded assembly", "source_file": "one row per input file",
    "sequence_region": "one row per FASTA sequence", "feature": "one row per segment",
    "feature_attribute": "one row per attribute value",
    "annotation_run": "one row per functional file",
    "protein_annotation": "one row per assignment",
    "vocabulary": "one row per type term", "vocabulary_meta": "a single row",
}
SCHEMA_GROUPS = {   # fill of each box: genome and provenance, features, function, seed
    "organism": NOTHING, "strain": NOTHING, "assembly": NOTHING, "source_file": NOTHING,
    "sequence_region": NOTHING, "feature": RAMP[2], "feature_attribute": RAMP[2],
    "annotation_run": RAMP[0], "protein_annotation": RAMP[0],
    "vocabulary": "white", "vocabulary_meta": "white",
}
ROW = 0.118   # inches per column line at 7 pt


def figure_schema(out: Path) -> None:
    if set(SCHEMA_COLUMNS) != set(SCHEMA.tables):
        raise SystemExit(f"the schema has {sorted(SCHEMA.tables)}; the figure draws "
                         f"{sorted(SCHEMA_COLUMNS)}")
    width, height = FULL_WIDTH, 4.6
    figure = plt.figure(figsize=(width, height))
    axis = figure.add_axes((0, 0, 1, 1))
    axis.set_xlim(0, width)
    axis.set_ylim(0, height)
    axis.axis("off")
    w = 1.5
    col = [0.04, 1.70, 3.36, 5.12]
    last = width - col[3] - 0.02   # the right column holds the longest names
    place = {   # table: (column, top)
        "organism": (0, 4.52), "strain": (0, 3.58), "vocabulary": (0, 1.62),
        "vocabulary_meta": (0, 0.80),
        "assembly": (1, 4.52), "source_file": (1, 3.10), "sequence_region": (1, 2.04),
        "feature": (2, 4.52),
        "feature_attribute": (3, 4.52), "protein_annotation": (3, 3.58),
        "annotation_run": (3, 2.28),
    }
    boxes = {}
    for name, columns in SCHEMA_COLUMNS.items():
        table = SCHEMA.tables[name]
        missing = [c for c in columns if c not in table.c]
        if missing:
            raise SystemExit(f"{name} has no column {missing}")
        c, top = place[name]
        x = col[c]
        bw = last if c == 3 else w
        h = (len(columns) + 2) * ROW + 0.1
        y = top - h
        boxes[name] = (x, y, bw, h)
        fill = SCHEMA_GROUPS[name]
        axis.add_patch(Rectangle((x, y), bw, h, facecolor=fill,
                                 edgecolor=AXIS if fill == "white" else "none",
                                 linewidth=0.6))
        axis.text(x + 0.06, top - 0.05, name, ha="left", va="top", fontsize=7,
                  fontweight="bold", color=INK)
        axis.text(x + 0.06, top - 0.05 - ROW, SCHEMA_GRAIN[name], ha="left", va="top",
                  fontsize=SMALL, style="italic", color=INK2)
        for i, column_name in enumerate(columns):
            column = table.c[column_name]
            mark = "PK" if column.primary_key else ("FK" if column.foreign_keys else "")
            yy = top - 0.05 - (i + 2) * ROW
            axis.text(x + 0.1, yy, column_name, ha="left", va="top", fontsize=SMALL,
                      color=INK if mark else INK2, family="monospace")
            if mark:
                axis.text(x + bw - 0.03, yy, mark, ha="right", va="top", fontsize=SMALL,
                          color=MUTED)

    def row_y(name, column_name):
        x, y, w_, h = boxes[name]
        top = y + h
        index = SCHEMA_COLUMNS[name].index(column_name)
        return top - 0.05 - (index + 2) * ROW - ROW / 2

    def link(child, column_name, parent, route=None):
        """child.column -> parent's key, from the side facing the parent."""
        cx, _, cw, _ = boxes[child]
        px, _, pw, _ = boxes[parent]
        key = next(c.name for c in SCHEMA.tables[parent].primary_key.columns)
        y0, y1 = row_y(child, column_name), row_y(parent, key)
        if px > cx:
            x0, x1 = cx + cw, px
        elif px < cx:
            x0, x1 = cx, px + pw
        else:   # same column: go round the left side
            x0 = x1 = cx
        if route is None:
            route = (x0 + x1) / 2 if x0 != x1 else x0 - 0.04
        axis.plot([x0, route, route, x1], [y0, y0, y1, y1], color=MUTED, linewidth=0.6,
                  solid_joinstyle="miter")
        axis.plot([x1], [y1], marker="o", markersize=2.2, color=MUTED)

    link("strain", "organism_id", "organism")
    link("assembly", "strain_id", "strain")
    link("source_file", "assembly_id", "assembly", route=1.665)
    link("sequence_region", "assembly_id", "assembly", route=1.665)
    link("feature", "assembly_id", "assembly", route=3.245)
    link("feature", "source_file_id", "source_file", route=3.285)
    link("feature", "sequence_region_id", "sequence_region", route=3.325)
    link("feature_attribute", "feature_id", "feature", route=5.0)
    link("protein_annotation", "feature_id", "feature", route=5.0)
    link("protein_annotation", "annotation_run_id", "annotation_run", route=5.09)
    link("annotation_run", "assembly_id", "assembly", route=3.245)
    # parent_id points at a row of the same table
    fx, _, fw, _ = boxes["feature"]
    y0, y1 = row_y("feature", "parent_id"), row_y("feature", "feature_id")
    right = fx + fw + 0.06
    axis.plot([fx + fw, right, right, fx + fw], [y0, y0, y1, y1], color=MUTED,
              linewidth=0.6)
    axis.plot([fx + fw], [y1], marker="o", markersize=2.2, color=MUTED)

    legend = [(NOTHING, "Genome and provenance"), (RAMP[2], "Annotation features"),
              (RAMP[0], "Functional annotation"), ("white", "Written by cgload init")]
    lx, ly = col[1], 0.62
    for i, (fill, label) in enumerate(legend):
        x = lx + (i % 2) * 2.0
        y = ly - 0.17 - (i // 2) * 0.19
        axis.add_patch(Rectangle((x, y), 0.16, 0.11, facecolor=fill,
                                 edgecolor=AXIS if fill == "white" else "none",
                                 linewidth=0.6))
        axis.text(x + 0.22, y + 0.055, label, ha="left", va="center", fontsize=SMALL,
                  color=INK2)
    axis.text(lx, ly + 0.08, "PK: primary key. FK: foreign key; a line joins it to the "
              "key it names, marked with a dot.", ha="left", va="center", fontsize=SMALL,
              color=INK2)
    save(figure, out, "fig2")


# --- Table 1: how the example file is stored, by running cgload on it -----------------

def table_example() -> str:
    """Load examples/minimal into a new SQLite database and show its first gene."""
    example = HERE.parent / "examples" / "minimal"
    with tempfile.TemporaryDirectory() as scratch:
        database = Path(scratch) / "example.db"
        environment = {"CGLOAD_DB_URL": f"sqlite:///{database}", "PATH": "/usr/bin:/bin",
                       "PYTHONPATH": str(HERE.parent / "src")}
        command = [sys.executable, "-m", "cgload"]
        for args in (["init"], ["load", "--config", str(example / "organism.toml"),
                                "--data-dir", str(Path(scratch) / "data")]):
            subprocess.run(command + args, check=True, cwd=example, env=environment,  # noqa: S603
                           capture_output=True)
        with sqlite3.connect(database) as connection:
            rows = connection.execute(
                "SELECT feature_id, source_id, source_type, start, \"end\", phase, "
                "segment_index, parent_id FROM feature "
                "WHERE source_id LIKE 'DEMO_000001%' ORDER BY feature_id").fetchall()
    body = []
    for row in rows:
        cells = ["" if value is None else str(value) for value in row]
        cells[1] = f"`{cells[1]}`"
        body.append(cells)
    header = ["feature_id", "source_id", "source_type", "start", "end", "phase",
              "segment_index", "parent_id"]
    # Break the long names after the underscore, so Word does not break them
    # inside a word.
    header = [name.replace("_", "_\u200b") for name in header]
    return md_table(header, body, "rllrrrrr", widths=[9, 22, 11, 9, 9, 7, 11, 9])


# --- Figure 4: the cost of a load (Q1) ------------------------------------------------

def figure_scale(results: Results, out: Path) -> None:
    panels = [
        ("load_wall_s", 1, "Load time (s)", "s"),
        ("verify_wall_s", 1, "Verify time (s)", "s"),
        ("load_peak_mb", 1 / 1024, "Peak memory of the load (GiB)", "GiB"),
        ("db_bytes", 1 / 1024 ** 3, "Database size (GiB)", "GiB"),
    ]
    backends = [("sqlite", "SQLite", BLUE), ("postgres", "PostgreSQL", ORANGE)]
    figure, axes = plt.subplots(2, 2, figsize=(FULL_WIDTH, 4.5))
    largest = max(results.genomes, key=lambda g: number(g["rows"]))
    for axis, (column, factor, title, unit), letter in zip(axes.flat, panels, "abcd",
                                                            strict=True):
        for backend, label, colour in backends:
            xs, ys = [], []
            for genome in results.genomes:
                row = results.scale[(genome["accession"], backend)]
                xs.append(number(genome["rows"]))
                ys.append(float(row[column]) * factor)
            axis.plot(xs, ys, linestyle="none", marker="o", markersize=4.6,
                      markerfacecolor=colour, markeredgecolor="white",
                      markeredgewidth=0.7, label=label, zorder=3)
        # The largest genome's two values, beside its points: the higher one
        # nudged up and the lower one down, so that close values do not collide.
        values = {backend: float(results.scale[(largest["accession"], backend)][column])
                  * factor for backend, _, _ in backends}
        higher = max(values, key=values.get)
        for backend, value in values.items():
            # Seconds from 10 up are whole, as the text gives them; GiB keep one
            # decimal (29.5 GiB, not 29), as in the text and Table S1.
            shown = (f"{value:,.0f} {unit}" if value >= 10 and unit == "s"
                     else f"{value:,.1f} {unit}")
            up = backend == higher
            axis.annotate(shown, (number(largest["rows"]), value),
                          xytext=(6, 2.5 if up else -2.5), textcoords="offset points",
                          ha="left", va="bottom" if up else "top", fontsize=SMALL, color=INK2)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(thousands))
        axis.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _p: f"{v:g}"))
        axis.set_xlim(1e4, 2e7)
        axis.grid(True, which="major")
        axis.grid(False, which="minor")
        axis.set_title(title, loc="left", color=INK, fontweight="bold")
        axis.text(-0.16, 1.06, letter, transform=axis.transAxes, fontsize=9,
                  fontweight="bold", color=INK, va="bottom")
        axis.set_xlabel("Rows in the GFF3 file")
    handles = [Line2D([], [], linestyle="none", marker="o", markersize=5,
                      markerfacecolor=colour, markeredgecolor="white", label=label)
               for _, label, colour in backends]
    figure.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
                  bbox_to_anchor=(0.5, 1.0), handletextpad=0.3, columnspacing=1.5)
    figure.tight_layout(rect=(0, 0, 1, 0.95), h_pad=1.6, w_pad=2.0)
    save(figure, out, "fig5")


# --- Figure 3: what each loader did with each genome (Q3) -----------------------------

LOADER_COLUMNS = [
    ("cgload-sqlite", "cgload", "SQLite and PostgreSQL"),
    ("gffutils-default", "gffutils", "default"),
    ("gffutils", "gffutils", "create_unique"),
    ("txdbmaker", "txdbmaker", "transcripts"),
    ("chado", "Chado", "python-chado"),
]


def silent_share(row: dict) -> float | None:
    """Percent of the file's features (transcripts) changed or lost without a word."""
    if row["outcome"] not in ("loaded", "loaded_with_repair"):
        return None
    truth = number(row["truth"])
    return 100 * (number(row["changed_silent"]) + number(row["lost_silent"])) / truth


def figure_loaders(results: Results, out: Path) -> None:
    for accession in (g["accession"] for g in results.genomes):
        sqlite = results.loader(accession, "cgload-sqlite")
        postgres = results.loader(accession, "cgload-postgres")
        same = ("outcome", "kept", "changed", "changed_declared", "lost")
        if any(sqlite[key] != postgres[key] for key in same):
            raise SystemExit(f"{accession}: the two cgload back ends differ; one column "
                             "cannot stand for both")
    shares = [silent_share(r) for r in results.loaders if r["tool"] != "cgload-postgres"]
    top = max(s for s in shares if s is not None)
    scale_max = 5 * (int(top // 5) + 1)

    rows_per_clade = [(clade, title, [g for g in results.order if g["clade"] == clade])
                      for clade, title in CLADES]
    row_h, gap_h, head_h = 0.165, 0.09, 0.42
    n_rows = len(results.order)
    body_h = n_rows * row_h + (len(CLADES) - 1) * gap_h
    legend_h = 0.62
    height = head_h + body_h + legend_h
    figure = plt.figure(figsize=(FULL_WIDTH, height))
    left, names_w = 0.28, 1.62
    col_w = (FULL_WIDTH - left - names_w - 0.08) / len(LOADER_COLUMNS)
    axis = figure.add_axes((0, 0, 1, 1))
    axis.set_xlim(0, FULL_WIDTH)
    axis.set_ylim(height, 0)
    axis.axis("off")
    x0 = left + names_w
    pad = 0.012  # inches of white between cells (about 2 px at screen size)

    for i, (_tool, name, setting) in enumerate(LOADER_COLUMNS):
        cx = x0 + (i + 0.5) * col_w
        axis.text(cx, head_h - 0.20, name, ha="center", va="bottom", fontsize=7.5,
                  fontweight="bold", color=INK)
        axis.text(cx, head_h - 0.06, setting, ha="center", va="bottom", fontsize=SMALL,
                  color=INK2)

    y = head_h
    for _clade, title, genomes in rows_per_clade:
        top_y = y
        for genome in genomes:
            accession = genome["accession"]
            repaired = genome["sqlite_outcome"] == "loaded_with_repair"
            label = species(genome["organism"]) + (" †" if repaired else "")
            axis.text(x0 - 0.06, y + row_h / 2, label, ha="right", va="center",
                      fontsize=SMALL, style="italic", color=INK)
            for i, (tool, _name, _setting) in enumerate(LOADER_COLUMNS):
                row = results.loader(accession, tool)
                share = silent_share(row)
                if share is None:
                    fill = FAILED
                    text = "timed out" if row["outcome"] == "timed_out" else "not loaded"
                    colour = INK2
                elif share == 0:
                    fill, text, colour = NOTHING, "0", INK2
                else:
                    fill = ramp_colour(share / scale_max)
                    text = "<0.1" if share < 0.1 else f"{share:.1f}"
                    colour = ink_on(fill)
                axis.add_patch(Rectangle((x0 + i * col_w + pad, y + pad),
                                         col_w - 2 * pad, row_h - 2 * pad,
                                         facecolor=fill, edgecolor="none"))
                axis.text(x0 + (i + 0.5) * col_w, y + row_h / 2, text, ha="center",
                          va="center", fontsize=SMALL, color=colour)
            y += row_h
        axis.text(left * 0.55, (top_y + y) / 2, title, rotation=90, ha="center",
                  va="center", fontsize=7, color=INK, fontweight="bold")
        axis.plot([left * 0.95, left * 0.95], [top_y + 0.02, y - 0.02], color=AXIS,
                  linewidth=0.6)
        y += gap_h

    # Key: the ramp, then "nothing" and "failed" chips.
    key_y = head_h + body_h + 0.20
    bar_x, bar_w, bar_h = x0, 2.2, 0.10
    steps = 60
    for k in range(steps):
        axis.add_patch(Rectangle((bar_x + k * bar_w / steps, key_y), bar_w / steps + 0.002,
                                 bar_h, facecolor=ramp_colour((k + 0.5) / steps),
                                 edgecolor="none"))
    for tick in range(0, scale_max + 1, 5):
        tx = bar_x + bar_w * tick / scale_max
        axis.plot([tx, tx], [key_y + bar_h, key_y + bar_h + 0.04], color=AXIS, linewidth=0.6)
        axis.text(tx, key_y + bar_h + 0.06, f"{tick}%", ha="center", va="top",
                  fontsize=SMALL, color=INK2)
    axis.text(bar_x, key_y - 0.05, "Features changed or lost without a message "
              "(% of the file's features)", ha="left", va="bottom", fontsize=SMALL, color=INK)
    chip_x = bar_x + bar_w + 0.35
    for fill, text in ((NOTHING, "none"), (FAILED, "not loaded")):
        axis.add_patch(Rectangle((chip_x, key_y), 0.18, bar_h, facecolor=fill,
                                 edgecolor="none"))
        axis.text(chip_x + 0.24, key_y + bar_h / 2, text, ha="left", va="center",
                  fontsize=SMALL, color=INK2)
        chip_x += 0.24 + 0.08 * len(text) + 0.25
    save(figure, out, "fig4")


# --- Figure 5: query waits (Q6) -------------------------------------------------------

TASKS = [("gene", "T1  A gene and its parts"),
         ("region", "T2  A 1 Mb window"),
         ("count", "T3  Counts by type"),
         ("sequence", "T4  Largest sequence as GFF3")]
QUERY_TOOLS = [  # (tool, marker, face, edge)
    ("cgload-sqlite", "o", BLUE, BLUE),
    ("cgload-postgres", "o", "white", BLUE),
    ("gffutils", "s", ORANGE, ORANGE),
    ("tabix", "D", AQUA, AQUA),
]
RESOLUTION = 0.001  # the clock's resolution, in seconds


def sequence_times(results: Results) -> dict:
    """T4 per sequence: rows written, and each tool's median over the repeated asks."""
    asked = defaultdict(lambda: defaultdict(list))
    rows = {}
    for r in read_tsv(results.folder / "query_timing.tsv"):
        if r["task"] != "sequence" or r["exit_code"] != "0":
            continue
        key = (r["accession"], r["target"])
        rows.setdefault(key, set()).add(number(r["answer"]))
        if r["first"] == "False":
            asked[key][r["tool"]].append(float(r["wall_s"]))
    found = {}
    for key, tools in asked.items():
        if len(rows[key]) != 1:
            raise SystemExit(f"{key}: the tools wrote different numbers of rows, {rows[key]}")
        found[key] = {"rows": next(iter(rows[key])),
                      "median_s": {tool: statistics.median(v) for tool, v in tools.items()}}
    return found


def largest_sequences(results: Results) -> dict:
    """Of each genome's T4 sequences, the one with most rows."""
    largest = {}
    for (accession, target), value in sequence_times(results).items():
        if accession not in largest or value["rows"] > largest[accession][1]["rows"]:
            largest[accession] = ((accession, target), value)
    return dict(largest.values())


def figure_queries(results: Results, out: Path) -> None:
    accessions = sorted({r["accession"] for r in results.timing},
                        key=lambda a: number(next(g["rows"] for g in results.genomes
                                                  if g["accession"] == a)))
    median = {(r["accession"], r["view"], r["tool"], r["task"]): float(r["repeat_median_s"])
              for r in results.timing}
    # T4's three sequences differ a hundredfold in size (the one with most
    # features and two drawn from those with at least 1,000), so their median
    # says little. The figure shows the largest; the table has all three.
    for (accession, _target), tools in largest_sequences(results).items():
        for tool, seconds in tools["median_s"].items():
            median[(accession, "command", tool, "sequence")] = seconds
    views = [("command", "One command per answer"), ("library", "One process, asked in turn")]
    figure, axes = plt.subplots(2, 4, figsize=(FULL_WIDTH, 3.9), sharey=True)
    for row_index, (view, _view_title) in enumerate(views):
        for col_index, (task, task_title) in enumerate(TASKS):
            axis = axes[row_index][col_index]
            asked = [t for t in QUERY_TOOLS
                     if any((a, view, t[0], task) in median for a in accessions)]
            if not asked:
                axis.axis("off")
                continue
            offsets = {t[0]: (k - (len(QUERY_TOOLS) - 1) / 2) * 0.17
                       for k, t in enumerate(QUERY_TOOLS)}
            for tool, marker, face, edge in asked:
                for y, accession in enumerate(accessions):
                    value = median.get((accession, view, tool, task))
                    if value is None:
                        continue
                    axis.plot(max(value, RESOLUTION), y + offsets[tool], linestyle="none",
                              marker=marker, markersize=4.2 if marker != "D" else 3.8,
                              markerfacecolor=face, markeredgecolor=edge,
                              markeredgewidth=0.9 if face == "white" else 0.0, zorder=3)
            axis.set_xscale("log")
            axis.set_xlim(RESOLUTION / 1.6, 60)
            # A gridline at every decade, a label at every other one.
            axis.set_xticks([1e-3, 1e-2, 1e-1, 1, 10])
            axis.set_xticklabels(["1 ms", "", "0.1 s", "", "10 s"])
            axis.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            axis.grid(True, axis="x")
            axis.grid(False, axis="y")
            axis.set_yticks(range(len(accessions)))
            axis.set_yticklabels([short_species(results.name[a]) for a in accessions],
                                 style="italic")
            axis.set_ylim(len(accessions) - 0.5, -0.5)
            axis.tick_params(axis="y", length=0)
            if row_index == 0:
                axis.set_title(task_title, loc="left", fontsize=SMALL, color=INK,
                               fontweight="bold")
    handles = [Line2D([], [], linestyle="none", marker=marker, markersize=5,
                      markerfacecolor=face, markeredgecolor=edge,
                      markeredgewidth=0.9 if face == "white" else 0.0,
                      label=TOOL_NAMES[tool])
               for tool, marker, face, edge in QUERY_TOOLS]
    legend_axis = axes[1][3]
    legend_axis.axis("off")
    legend_axis.legend(handles=handles, loc="center left", frameon=False,
                       title="Median of repeated asks", title_fontsize=SMALL,
                       handletextpad=0.4, labelspacing=0.7)
    figure.tight_layout(w_pad=0.8, h_pad=1.4)
    figure.subplots_adjust(left=0.155)
    for row_index, (_view, view_title) in enumerate(views):
        box = axes[row_index][0].get_position()
        figure.text(0.012, (box.y0 + box.y1) / 2, view_title, rotation=90, ha="left",
                    va="center", fontsize=7, fontweight="bold", color=INK)
    save(figure, out, "fig6")


# --- tables ----------------------------------------------------------------------------

def md_table(header: list[str], rows: list[list[str]], align: str | None = None,
             widths: list[int] | None = None) -> str:
    """A pipe table. ``widths`` sets the relative column widths in the Word file:
    pandoc sizes the columns of a wide table by the dashes under each header."""
    align = align or "l" + "r" * (len(header) - 1)
    widths = widths or [1] * len(header)

    def mark(a: str, w: int) -> str:
        dashes = "-" * (2 * w)
        return {"l": ":" + dashes, "r": dashes + ":", "c": ":" + dashes + ":"}[a]

    lines = ["| " + " | ".join(header) + " |",
             "| " + " | ".join(mark(a, w) for a, w in zip(align, widths, strict=True)) + " |"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def write_csv(table: str, path: Path) -> None:
    """A Markdown table as CSV: emphasis marks off, thousands separators off."""
    lines = [line for line in table.splitlines() if line.startswith("|")]
    # With a byte-order mark Excel reads the file as UTF-8, so a dash survives,
    # and does not take a file that starts with "ID" for a SYLK file.
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        for index, line in enumerate(lines):
            if index == 1:              # the alignment row
                continue
            cells = [cell.strip().replace("*", "").replace("`", "")
                     for cell in line.strip().strip("|").split("|")]
            writer.writerow([cell.replace(",", "") if re.fullmatch(r"[\d,]+(\.\d+)?", cell)
                             else cell for cell in cells])


def n(value: int) -> str:
    return f"{value:,}"


def table_loaders(results: Results) -> str:
    """One column per loader: what it did with the genomes it loaded. A load that
    did not finish stored nothing, so its genomes are counted as not loaded and
    not compared; the split into announced and silent means nothing for a load
    that reports its own failure. cgload has one column, since its two database
    systems gave the same result (figure_loaders checks that)."""
    keys = ("truth", "kept", "changed", "changed_declared", "changed_silent", "lost",
            "lost_silent")
    accessions = [g["accession"] for g in results.genomes]
    columns = [("cgload-sqlite", "cgload"), ("gffutils-default", "gffutils, default"),
               ("gffutils", "gffutils, create_unique"), ("txdbmaker", "txdbmaker"),
               ("chado", "Chado (python-chado)")]
    measures = ["Genomes loaded", "Genomes not loaded", "Compared", "Kept",
                "Changed, declared", "Changed, announced", "Changed, silent",
                "Lost, announced", "Lost, silent"]
    cells: dict[str, list[str]] = {}
    for tool, label in columns:
        each = [results.loader(accession, tool) for accession in accessions]
        done = [r for r in each if r["outcome"] in ("loaded", "loaded_with_repair")]
        total = {k: sum(number(r[k]) for r in done) for k in keys}
        summary = results.tools[tool]  # the benchmark's own totals, as a cross-check
        if tool != "chado" and total["truth"] != number(summary["truth"]):
            raise SystemExit(f"{tool}: {total['truth']} here, {summary['truth']} in tools.tsv")
        counts = [total["truth"], total["kept"], total["changed_declared"],
                  total["changed"] - total["changed_silent"], total["changed_silent"],
                  total["lost"] - total["lost_silent"], total["lost_silent"]]
        shown = [n(c) for c in counts] if total["truth"] else ["–"] * len(counts)
        if tool == "txdbmaker":
            shown[0] += " transcripts"
        cells[label] = [str(len(done)), str(len(each) - len(done)), *shown]
    header = ["", *(label for _tool, label in columns)]
    rows = [[measure, *(cells[label][i] for _tool, label in columns)]
            for i, measure in enumerate(measures)]
    return md_table(header, rows, "l" + "r" * len(columns))


CASE_NAMES = {
    "multi_row": "Written on several rows", "multi_valued": "Attribute with several values",
    "rows_disagree": "Rows carry different attributes", "overlapping_rows": "Rows overlap",
    "trans_spliced": "Trans-spliced", "organelle": "On an organelle",
    "exon_under_pseudogene": "Exon under a pseudogene", "escapes_parent": "Outside its parent",
    "crosses_origin": "Crosses the origin", "strand_unknown": "Strand unknown (?)",
    "past_end": "Past the end of its sequence", "inverted": "Start after end",
}


def chado_cases_finished(results: Results) -> tuple[int, dict[str, Counter]]:
    """Chado's hard cases in the loads that finished only: a failed load stored
    nothing, so counting its genomes would put every crash under every case."""
    totals: dict[str, Counter] = defaultdict(Counter)
    finished = [r for r in results.loaders if r["tool"] == "chado" and r["outcome"] == "loaded"]
    for row in finished:
        comparison = json.loads(
            (results.folder / "compare" / f"{row['accession']}.chado.json").read_text())
        for case, counts in comparison["by_case"].items():
            totals[case].update(counts)
    return len(finished), totals


def table_cases(results: Results) -> str:
    """Per hard case, what each loader altered. txdbmaker counts transcripts in
    the genomes it loaded, and Chado features in its loads that finished; both
    carry their own totals."""
    by = {(r["tool"], r["case"]): r for r in results.cases}
    finished, chado = chado_cases_finished(results)
    # Most features first; ties by name, so that drawing again changes nothing.
    cases = sorted({r["case"] for r in results.cases},
                   key=lambda c: (-number(by[("cgload-sqlite", c)]["truth"])
                                  if ("cgload-sqlite", c) in by else 0, c))

    def counts(tool, case) -> dict | None:
        if tool == "chado":
            return {k: str(v) for k, v in chado[case].items()} if case in chado else None
        return by.get((tool, case))

    def cell(tool, case):
        r = counts(tool, case)
        if r is None:
            return "–"
        altered = number(r.get("changed")) + number(r.get("lost"))
        silent = number(r.get("changed_silent")) + number(r.get("lost_silent"))
        declared = number(r.get("changed_declared"))
        own_total = tool in ("txdbmaker", "chado")
        of = f" of {n(number(r['truth']))}" if own_total else ""
        if declared and not altered:
            return f"{n(declared)} repaired, declared{of}"
        if not altered:
            return f"0{of}" if own_total else "none"
        if silent == altered:
            return f"{n(altered)}{of}, all silent"
        if silent:
            return f"{n(altered)}{of}, {n(silent)} silent"
        return f"{n(altered)}{of}, all announced"

    tools = ["cgload-sqlite", "gffutils", "chado", "txdbmaker"]
    header = ["Hard case", "Features", "cgload", "gffutils, create_unique",
              f"Chado, {finished} finished loads", "txdbmaker (transcripts)"]
    rows = [[CASE_NAMES.get(c, c), n(number(by[("cgload-sqlite", c)]["truth"]))]
            + [cell(t, c) for t in tools] for c in cases]
    return md_table(header, rows, "lrllll")


FAULT_GROUPS = [
    ("Features deleted, duplicated or retyped in the database",
     {"F01", "F02", "F03", "F09", "F11", "F19"}),
    ("Parent links removed from the database", {"F04", "F05", "F18"}),
    ("Stored assembly damaged or moved", {"F13", "F14", "F15", "F17"}),
    ("Input file edited after the load", {"F16"}),
    ("Parser drops, invents or retypes rows", {"P01", "P02", "P03", "P04", "P06"}),
    # Strand, start, an attribute, an ID or a parent changed: every count stays.
    ("Values, IDs or parents changed, counts kept (database)",
     {"F06", "F07", "F08", "F10", "F12"}),
    ("Values or IDs changed, counts kept (parser)", {"P05", "P07", "P08"}),
]


def caught_by(trials: list[dict]) -> str:
    """The checks that caught a fault, read from each trial's tiers and message."""
    found = []

    def add(label):
        if label not in found:
            found.append(label)

    for trial in trials:
        message = trial["message"]
        if trial["tier_source_count"] == "fail":
            add("source count")
        if trial["tier_source_count"] == "skipped" and "has changed since" in message:
            add("input's SHA-256 (count skipped, exit 0)")
        if trial["tier_parent_edges"] == "fail":
            add("parent links")
        if trial["tier_stored_assembly"] == "fail":
            add("stored assembly")
        if trial["load_refused"] == "True":
            for text, label in (("which is not defined", "load refused: parent not defined"),
                                ("did not store the parent links",
                                 "load rolled back: links short"),
                                ("did not store what the file contains",
                                 "load rolled back: counts short"),
                                ("two stored rows share", "load rolled back: row stored twice")):
                if text in message:
                    add(label)
    return ", ".join(found)


def table_faults(results: Results) -> str:
    by = {r["fault_id"]: r for r in results.faults}
    rows = []
    for title, ids in FAULT_GROUPS:
        members = [by[i] for i in sorted(ids)]
        trials = sum(number(r["genomes"]) for r in members)
        detected = sum(number(r["detected"]) for r in members)
        as_expected = sum(number(r["as_expected"]) for r in members)
        expected = {r["expected"] for r in members}
        caught = caught_by([r for r in results.fault_trials
                            if r["fault_id"] in ids and r["observed"] == "detected"])
        failed = failing_exit(results, ids)
        rows.append([title, ", ".join(sorted(ids)), " / ".join(sorted(expected)),
                     f"{detected} of {trials}", f"{failed} of {trials}",
                     f"{as_expected} of {trials}", caught or "–"])
    header = ["Fault", "Faults", "Expected", "Reported", "Exit 2", "As expected",
              "Caught by"]
    return md_table(header, rows, "lllrrrl")


def failing_exit(results: Results, ids) -> int:
    """Trials of these faults whose command exited with a failing status. A fault
    can be reported and still exit 0: F16, where verify notices the edited input,
    skips the source count and passes."""
    return sum(1 for r in results.fault_trials
               if r["fault_id"] in ids and r["exit_code"] not in ("0", "NA", ""))


def table_faults_full(results: Results) -> str:
    rows = [[r["fault_id"], r["fault"].replace("_", " "), r["layer"], r["expected"],
             f"{number(r['detected'])} of {number(r['genomes'])}",
             f"{failing_exit(results, {r['fault_id']})} of {number(r['genomes'])}",
             f"{number(r['as_expected'])} of {number(r['genomes'])}",
             caught_by([x for x in results.fault_trials if x["fault_id"] == r["fault_id"]
                        and x["observed"] == "detected"]) or "–"] for r in results.faults]
    header = ["Fault ID", "Fault", "Layer", "Expected", "Reported", "Exit 2", "As expected",
              "Caught by"]
    return md_table(header, rows, "llllrrrl")


def table_kills(results: Results) -> str:
    """What each interrupted load left, from the trials themselves."""
    trials = read_tsv(results.folder / "kills.tsv")
    labels = {"cgload": "cgload, new load", "cgload-force": "cgload, `--force` reload",
              "gffutils": "gffutils, create_unique", "chado": "Chado (python-chado)"}
    groups = defaultdict(list)
    for trial in trials:
        groups[(trial["tool"], trial["accession"])].append(trial)
    rows = []
    for (tool, accession), group in groups.items():
        states = defaultdict(int)
        for trial in group:
            states[trial["state"]] += 1
        notes = []
        if tool == "gffutils":
            opens = sum(trial["opens"] == "True" for trial in group)
            notes.append(f"gffutils opens {opens} of the {len(group)} files left")
        if tool == "cgload":
            again = sum(trial["reload_ok"] == "True" for trial in group)
            notes.append(f"the repeated load succeeded {again} times")
        if tool == "cgload-force":
            passed = sum(trial["verify_after"] == "pass" for trial in group)
            notes.append(f"the earlier load passed verify {passed} times")
        if tool.startswith("cgload"):
            left = [number(trial["leftover_bytes"]) for trial in group]
            if any(left):
                notes.append(f"{sum(1 for b in left if b)} left staging files, up to "
                             f"{max(left) / 1024 ** 2:,.0f} MiB")
        rows.append([labels.get(tool, tool), species(results.name[accession]),
                     str(len(group)),
                     ", ".join(f"{v} {k}" for k, v in sorted(states.items())),
                     "; ".join(notes) or "–"])
    header = ["Loader", "Genome", "Kills", "State afterwards", "Notes"]
    return md_table(header, rows, "llrll")


def table_agreement(results: Results) -> str:
    """cgload's structural checks and the census, row by row of the census."""
    rows_of = read_tsv(results.folder / "agreement.tsv")
    census = [r for r in rows_of if r["side"] == "survey"]
    only = [r for r in rows_of if r["side"] == "cgload"]
    by_agreement = Counter(r["agreement"] for r in census)
    by_repair = Counter(r["cgload_repair"] for r in census)
    parents = {(r["accession"], r["detail"]) for r in census
               if r["cgload_repair"] == "accepted once its parent is repaired"}
    expected = {"both accept as valid", "both object", "survey: documented; cgload objects"}
    if set(by_agreement) - expected:
        raise SystemExit(f"agreement.tsv has classes this table does not show: "
                         f"{set(by_agreement) - expected}")
    totals = {k: sum(number(r[k]) for r in results.agreement)
              for k in ("survey_rows", "survey_error_cgload_accepts", "cgload_only")}
    if totals["survey_rows"] != len(census) or totals["cgload_only"] != len(only):
        raise SystemExit("agreement.tsv and tables/agreement.tsv disagree")
    rows = [
        ["Rows listed by the census", n(len(census))],
        ["*Without `--repair-envelopes`*", ""],
        ["Both accept: the only finding is the strand `?`",
         n(by_agreement["both accept as valid"])],
        ["Both object: the census classes the row as an error or as unexplained",
         n(by_agreement["both object"])],
        ["The census classes the row as documented; cgload refuses the file",
         n(by_agreement["survey: documented; cgload objects"])],
        ["The census classes the row as an error; cgload accepts it",
         n(totals["survey_error_cgload_accepts"])],
        ["*With `--repair-envelopes`*", ""],
        ["Accepted as written", n(by_repair["accepted"])],
        [f"Accepted once the row's parent is repaired ({len(parents)} parents)",
         n(by_repair["accepted once its parent is repaired"])],
        ["The row's own feature is repaired", n(by_repair["repaired"])],
        ["Stored as written, with a warning", n(by_repair["stored as written, with a warning"])],
        ["Stored as written: the disagreement exists only across an origin",
         n(by_repair["stored as written (D-091)"])],
        ["*Features only cgload objects to, each the parent of a listed row*", n(len(only))],
    ]
    if sum(number(r[1]) for r in rows[7:12]) != len(census):
        raise SystemExit("the rows under --repair-envelopes do not add up to the census's rows")
    return md_table(["", "Rows"], rows, "lr")


def name_md(organism: str) -> str:
    """Genus and species in italics, the strain upright."""
    kept = species(organism)
    rest = organism[len(kept):].strip()
    return f"*{kept}*" + (f" {rest}" if rest else "")


def input_files(results: Results) -> dict[str, tuple[str, str]]:
    """Each input GFF3 file: its annotation release, from the census, and its MD5.

    An accession does not pin a RefSeq annotation: NCBI replaces some files in
    place. The release named in the file's own header and the MD5 the benchmark
    checked before each run do."""
    survey = HERE.parent / "survey" / "results" / "survey.tsv"
    census = {row["accession"]: row for row in read_tsv(survey) if row["arm"] == "refseq"}
    files = {}
    for row in read_tsv(results.folder / "inputs.tsv"):
        accession = row["accession"]
        if accession not in census:
            raise SystemExit(f"{accession} is not in the census")
        if census[accession]["file_md5"] != row["gff_md5"]:
            raise SystemExit(f"{accession}: the census read another file than the benchmark")
        files[accession] = (census[accession]["annotation_release"], row["gff_md5"])
    return files


def table_genomes(results: Results) -> str:
    files = input_files(results)
    rows = []
    for g in results.order:
        a = results.scale[(g["accession"], "sqlite")]
        b = results.scale[(g["accession"], "postgres")]
        outcome = "repaired" if g["sqlite_outcome"] == "loaded_with_repair" else "loaded"
        release, md5 = files[g["accession"]]
        rows.append([name_md(g["organism"]), g["accession"], release, md5, g["clade"],
                     n(number(g["rows"])),
                     n(number(g["features"])), outcome,
                     f"{float(a['load_wall_s']):,.1f}", f"{float(b['load_wall_s']):,.1f}",
                     f"{float(a['load_peak_mb']) / 1024:.1f}",
                     f"{float(b['load_peak_mb']) / 1024:.1f}",
                     f"{float(a['verify_wall_s']):.1f}", f"{float(b['verify_wall_s']):.1f}",
                     f"{number(a['db_bytes']) / 1024 ** 3:.2f}",
                     f"{number(b['db_bytes']) / 1024 ** 3:.2f}",
                     f"{number(a['data_bytes']) / 1024 ** 3:.2f}"])
    header = ["Organism", "Accession", "Annotation release", "GFF3 MD5", "Clade", "Rows",
              "Features", "cgload",
              "Load s, SQLite", "Load s, PostgreSQL", "Peak GiB, SQLite",
              "Peak GiB, PostgreSQL", "Verify s, SQLite", "Verify s, PostgreSQL",
              "Database GiB, SQLite", "Database GiB, PostgreSQL", "Data folder GiB"]
    return md_table(header, rows, "lllllrrl" + "r" * 9)


def table_sequences(results: Results) -> str:
    """T4 per sequence: the median over three unlike sequences hides the largest."""
    size = {g["accession"]: number(g["rows"]) for g in results.genomes}
    tools = [tool for tool, *_ in QUERY_TOOLS]
    rows = []
    for (accession, target), value in sorted(
            sequence_times(results).items(),
            key=lambda item: (size[item[0][0]], -item[1]["rows"])):
        rows.append([f"*{species(results.name[accession])}*", target, n(value["rows"]),
                     *[f"{value['median_s'][tool]:.2f}" if tool in value["median_s"] else "–"
                       for tool in tools]])
    header = ["Genome", "Sequence", "Rows written", *[TOOL_NAMES[tool] + ", s" for tool in tools]]
    return md_table(header, rows, "llr" + "r" * len(tools))


# --- The portal (Q9; PROTOCOL.md section 13) --------------------------------------------

PENDING = "pending"
REQUEST_KINDS = [  # (kind in portal_requests.tsv, what the paper calls it)
    ("search", "Search by the gene's name"),
    ("gene_page", "Gene page"),
    ("region_page", "Region page, gene and 5 kb either side"),
    ("protein", "Protein, as FASTA"),
    ("range_annotation", "64 KiB of the annotation file"),
    ("range_assembly", "64 KiB of the assembly file"),
]
TRANSLATION_CLASSES = ["no_internal_stop", "stop_declared_exception", "stop_partial",
                       "stop_pseudo", "stop_unexplained", "failed"]


class Portal:
    """What benchmark/run_portal.py left in its results folder."""

    def __init__(self, folder: Path):
        self.folder = folder
        self.loads = read_tsv(folder / "portal_loads.tsv")
        self.index = read_tsv(folder / "portal_index.tsv")
        self.requests = read_tsv(folder / "portal_requests.tsv")
        self.found = read_tsv(folder / "portal_found.tsv")
        self.concurrency = read_tsv(folder / "portal_concurrency.tsv")
        self.translation = read_tsv(folder / "portal_translation.tsv")
        self.security = read_tsv(folder / "portal_security.tsv")
        self.readonly = read_tsv(folder / "portal_readonly.tsv")
        self.environment = json.loads((folder / "portal_environment.json").read_text())

    def times(self, backend: str, kind: str, accession: str | None = None) -> list[float]:
        return [float(r["ms"]) for r in self.requests
                if r["backend"] == backend and r["kind"] == kind
                and (accession is None or r["accession"] == accession)]


def nearest_rank(values: list[float], share: float) -> float:
    """The percentile run_portal.py reports: the value at that rank, no interpolation."""
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(round(share * (len(ordered) - 1))))]


def ms(value: float) -> str:
    return f"{value:,.0f}" if value >= 100 else f"{value:.1f}"


def table_portal(portal: Portal | None) -> str:
    """Table 8: how long each kind of request took, on SQLite with all 30 genomes in
    one database, and on PostgreSQL with human."""
    header = ["Request", "Requests", "Median, ms", "95th percentile, ms", "Slowest, ms",
              "Median, ms (PostgreSQL, human)", "95th percentile, ms (PostgreSQL, human)"]
    rows = []
    for kind, label in REQUEST_KINDS:
        if portal is None:
            rows.append([label] + [PENDING] * 6)
            continue
        times = portal.times("sqlite", kind)
        pg = portal.times("postgresql", kind)
        rows.append([label, n(len(times)), ms(statistics.median(times)),
                     ms(nearest_rank(times, 0.95)), ms(max(times)),
                     ms(statistics.median(pg)) if pg else "–",
                     ms(nearest_rank(pg, 0.95)) if pg else "–"])
    return md_table(header, rows, "lrrrrrr", [5, 2, 2, 2, 2, 3, 3])


def table_portal_genomes(results: Results, portal: Portal | None) -> str:
    """Table S5: the portal per genome, in the order of Table S1."""
    header = ["Organism", "Accession", "Loaded", "Index terms", "Index s", "Index check",
              *[label.split(",")[0] + ", median ms" for _, label in REQUEST_KINDS[:4]],
              "Genes searched", "Found", "Found on page 1", "CDS translated",
              "No stop inside", "Stop, exception declared", "Stop, partial CDS",
              "Stop, pseudogene", "Stop, unexplained"]
    rows = []
    for g in results.order:
        accession = g["accession"]
        if portal is None:
            rows.append([name_md(g["organism"]), accession] + [PENDING] * (len(header) - 2))
            continue
        load = next((r for r in portal.loads if r["accession"] == accession
                     and r["backend"] == "sqlite"), None)
        index = next((r for r in portal.index if r["accession"] == accession
                      and r["backend"] == "sqlite"), None)
        found = [r for r in portal.found if r["accession"] == accession]
        tr = next((r for r in portal.translation if r["accession"] == accession), None)
        medians = []
        for kind, _ in REQUEST_KINDS[:4]:
            times = portal.times("sqlite", kind, accession)
            medians.append(ms(statistics.median(times)) if times else "–")
        rows.append([
            name_md(g["organism"]), accession,
            {"loaded": "yes", "loaded_with_repair": "repaired"}.get(
                load["outcome"], load["outcome"]) if load else "–",
            n(number(index["terms"])) if index else "–",
            f"{float(index['seconds']):,.1f}" if index else "–",
            index["check"] if index else "–",
            *medians, n(len(found)), n(sum(r["found"] == "True" for r in found)),
            n(sum(r["page"] == "1" for r in found)),
            *([n(number(tr["checked"]))] + [n(number(tr[key])) for key in
                                             TRANSLATION_CLASSES[:5]] if tr else ["–"] * 6)])
    return md_table(header, rows, "ll" + "l" + "rr" + "l" + "r" * (len(header) - 6))


def portal_versions(portal: Portal | None) -> list[list[str]]:
    """Table 2's rows for the portal's run, as its environment recorded them."""
    from_package = HERE.parent / "src" / "cgload" / "portal" / "static" / "igv" / "VENDORED.md"
    igv = re.search(r"igv\.js (\d+\.\d+\.\d+)", from_package.read_text()).group(1)
    if portal is None:
        return [["cgload", "0.3.0", "the portal (Q9; run pending)"],
                ["Flask, Werkzeug, Jinja2", PENDING, "the portal"],
                ["waitress", PENDING, "the portal's web server"],
                ["igv.js", igv, "the portal's genome browser, served with cgload"]]
    env = portal.environment
    packages = env["packages"]
    when = ""
    if env.get("started"):
        import datetime

        day = datetime.date.fromisoformat(env["started"])
        when = f", run on {day.day} {day:%B %Y}"
    return [["cgload", env["cgload"], f"the portal (Q9), at commit {env['commit']}{when}"],
            ["Flask, Werkzeug, Jinja2",
             f"{packages['flask']}, {packages['werkzeug']}, {packages['jinja2']}", "the portal"],
            ["waitress", packages["waitress"], "the portal's web server"],
            ["igv.js", igv, "the portal's genome browser, served with cgload"]]


def portal_numbers(portal: Portal) -> dict[str, str]:
    """The portal's numbers the text quotes, as the text writes them."""
    out: dict[str, str] = {}
    loads = [r for r in portal.loads if r["backend"] == "sqlite"]
    out["genomes"] = n(len(loads))
    out["loaded"] = n(sum(r["outcome"] in ("loaded", "loaded_with_repair") for r in loads))
    out["repaired"] = n(sum(r["outcome"] == "loaded_with_repair" for r in loads))
    out["load_hours"] = f"{sum(float(r['wall_s']) for r in loads) / 3600:.1f}"
    index = [r for r in portal.index if r["backend"] == "sqlite"]
    out["index_pass"] = n(sum(r["check"] == "PASS" for r in index))
    out["index_terms_million"] = f"{sum(number(r['terms']) for r in index) / 1e6:.1f}"
    out["index_minutes"] = f"{sum(float(r['seconds']) for r in index) / 60:.0f}"
    human = next((r for r in index if r["accession"] == "GCF_000001405.40"), None)
    if human:
        out["index_human_minutes"] = f"{float(human['seconds']) / 60:.1f}"
    db = portal.readonly[0]
    out["database_gib"] = f"{number(db['bytes']) / 1024 ** 3:.1f}"
    out["database_unchanged"] = "unchanged" if db["unchanged"] == "True" else "CHANGED"
    sqlite = [r for r in portal.requests if r["backend"] == "sqlite"]
    out["requests"] = n(len(sqlite))
    out["requests_not_ok"] = n(sum(not 200 <= int(r["status"]) < 400 for r in sqlite))
    for kind, _ in REQUEST_KINDS:
        times = portal.times("sqlite", kind)
        out[f"{kind}_median"] = ms(statistics.median(times))
        out[f"{kind}_p95"] = ms(nearest_rank(times, 0.95))
        out[f"{kind}_max"] = ms(max(times))
        pg = portal.times("postgresql", kind)
        if pg:
            out[f"{kind}_pg_median"] = ms(statistics.median(pg))
            out[f"{kind}_pg_p95"] = ms(nearest_rank(pg, 0.95))
    every = [float(r["ms"]) for r in sqlite]
    out["all_p95"] = ms(nearest_rank(every, 0.95))
    out["all_max_s"] = f"{max(every) / 1000:.1f}"
    out["under_1s_percent"] = f"{100 * sum(t < 1000 for t in every) / len(every):.1f}"
    found = portal.found
    out["genes"] = n(len(found))
    out["found"] = n(sum(r["found"] == "True" for r in found))
    out["found_percent"] = f"{100 * sum(r['found'] == 'True' for r in found) / len(found):.1f}"
    out["found_first_page"] = n(sum(r["page"] == "1" for r in found))
    missed = [r for r in found if r["found"] != "True"]
    out["missed"] = n(len(missed))
    out["missed_few"] = n(sum(1 for r in missed if r["matches"] and int(r["matches"]) <= 1000))
    for row in portal.concurrency:
        c = row["clients"]
        out[f"c{c}_per_second"] = f"{float(row['per_second']):,.0f}"
        out[f"c{c}_p50"] = ms(float(row["p50_ms"]))
        out[f"c{c}_p95"] = ms(float(row["p95_ms"]))
        out[f"c{c}_p99"] = ms(float(row["p99_ms"]))
        out[f"c{c}_requests"] = n(number(row["requests"]))
    out["load_test_errors"] = n(sum(number(r["errors"]) for r in portal.concurrency))
    total = {key: sum(number(r[key]) for r in portal.translation)
             for key in ["checked", *TRANSLATION_CLASSES]}
    for key, value in total.items():
        out[f"cds_{key}"] = n(value)
    out["cds_clean_percent"] = f"{100 * total['no_internal_stop'] / total['checked']:.2f}"
    out["cds_explained"] = n(total["stop_declared_exception"] + total["stop_partial"]
                             + total["stop_pseudo"])
    probe = portal.security[0]
    statuses = json.loads(probe["statuses"])
    out["probe_requests"] = n(number(probe["requests"]))
    out["probe_server_errors"] = n(number(probe["server_errors"]))
    out["probe_missing_headers"] = n(number(probe["missing_security_headers"]))
    out["probe_failures"] = "none" if not probe["failures"] else probe["failures"]
    for status in ("200", "400", "404", "405"):
        out[f"probe_{status}"] = n(statuses.get(status, 0))
    env = portal.environment
    out["cgload"] = env["cgload"]
    out["cpus"] = n(env["cpus"])
    return out


VERSIONS = [  # (what the paper calls it, file, kind, package name, what it was used for)
    ("Python", "env.txt", "conda", "python", "cgload, gffutils and the benchmark scripts"),
    ("SQLite", "env.txt", "conda", "libsqlite", "cgload and gffutils databases"),
    ("PostgreSQL", "env.txt", "conda", "postgresql", "cgload and Chado databases"),
    ("SQLAlchemy", "env-pip.txt", "pip", "SQLAlchemy", "cgload"),
    ("pyfaidx", "env-pip.txt", "pip", "pyfaidx", "cgload"),
    ("click", "env-pip.txt", "pip", "click", "cgload"),
    ("psycopg", "env-pip.txt", "pip", "psycopg", "cgload on PostgreSQL"),
    ("HTSlib (tabix, bgzip)", "env.txt", "conda", "htslib", "query timing"),
    ("gffutils", "env-pip.txt", "pip", "gffutils", "loader compared"),
    ("R", "env-r.txt", "conda", "r-base", "txdbmaker"),
    ("txdbmaker", "env-r.txt", "conda", "bioconductor-txdbmaker", "loader compared"),
    ("GenomicFeatures", "env-r.txt", "conda", "bioconductor-genomicfeatures", "txdbmaker"),
    ("rtracklayer", "env-r.txt", "conda", "bioconductor-rtracklayer", "txdbmaker"),
    ("python-chado", "env-chado-pip.txt", "pip", "chado", "loader compared"),
    ("Python, for python-chado", "env-chado.txt", "conda", "python", "python-chado"),
    ("SQLAlchemy, for python-chado", "env-chado-pip.txt", "pip", "SQLAlchemy", "python-chado"),
    ("Biopython", "env-chado-pip.txt", "pip", "biopython", "python-chado"),
    ("bcbio-gff", "env-chado-pip.txt", "pip", "bcbio-gff", "python-chado"),
]


def version_of(text: str, kind: str, package: str) -> str:
    if kind == "pip":
        found = re.search(rf"^{re.escape(package)}==(\S+)$", text, re.M | re.I)
    else:
        found = re.search(rf"/{re.escape(package)}-([0-9][^-/]*)-[^-/]+\.(?:conda|tar\.bz2)\s*$",
                          text, re.M)
    if not found:
        raise SystemExit(f"no version of {package} in the recorded environment")
    return found.group(1)


#: The release whose code the benchmark ran (commits 8c999e2 and 9b97178).
BENCHMARKED_CGLOAD = "0.2.1"


def table_versions(results: Results, portal: Portal | None) -> str:
    """Versions as the server's environments recorded them at the time of the run."""
    folder = results.folder / "environment"
    # The version the benchmark ran, not the one installed now: 0.2.2 changed
    # how export writes columns 2 and 6 after the runs (D-109, D-110).
    rows = [["cgload", BENCHMARKED_CGLOAD,
             "the loader presented; benchmarked at commits 8c999e2 and 9b97178, the same "
             "code apart from its version string and comments"],
            ["cgload", "0.2.2", "the check of GFF3 columns 2, 6 and 8, run on 20 November "
             "2026 (tag v0.2.2, commit dbb51ce)"]]
    for label, file, kind, package, used in VERSIONS:
        rows.append([label, version_of((folder / file).read_text(), kind, package), used])
    chado = {r["tool_version"] for r in results.loaders if r["tool"] == "chado"}
    if chado != {"python-chado 2.3.9, Chado 1.31"}:
        raise SystemExit(f"Chado's recorded version is {chado}")
    # The schema and the ontology are files that benchmark/run.sh fetches at fixed
    # commits; chado-inputs.log records their MD5. The lnc_RNA result depends on the
    # ontology's version, so both are named.
    at = [row[0] for row in rows].index("python-chado") + 1
    rows[at:at] = [["Chado schema", "1.31 (GMOD/Chado commit ef04588)", "python-chado's database"],
                   ["Sequence Ontology", "SO-Ontologies commit 4340c14",
                    "python-chado's database (so.obo)"]]
    rows[2:2] = portal_versions(portal)
    return md_table(["Software", "Version", "Used for"], rows, "lll")


def table_timing(results: Results) -> str:
    rows = []
    names = dict(TASKS)
    size = {g["accession"]: number(g["rows"]) for g in results.genomes}
    for r in sorted(results.timing,
                    key=lambda r: (r["view"], r["task"], size[r["accession"]], r["tool"])):
        rows.append([r["view"], names.get(r["task"], r["task"]).split("  ")[0],
                     species(results.name[r["accession"]]), TOOL_NAMES.get(r["tool"], r["tool"]),
                     r["targets"], r["first_run_median_s"], r["repeat_median_s"]])
    header = ["View", "Task", "Genome", "Tool", "Targets", "First ask, median s",
              "Repeated asks, median s"]
    return md_table(header, rows, "lllllrr")


def write_tables(results: Results, portal: Portal | None, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    numbers = portal_numbers(portal) if portal else {}
    (out / "portal_numbers.json").write_text(json.dumps(numbers, indent=1, sort_keys=True) + "\n")
    tables = {
        "table7_portal.md": table_portal(portal),
        "tableS5_portal.md": table_portal_genomes(results, portal),
        "table0_example.md": table_example(),
        "table1_loaders.md": table_loaders(results),
        "table2_cases.md": table_cases(results),
        "table3_faults.md": table_faults(results),
        "table4_kills.md": table_kills(results),
        "table5_agreement.md": table_agreement(results),
        "tableS1_genomes.md": table_genomes(results),
        "tableS2_faults.md": table_faults_full(results),
        "tableS3_timing.md": table_timing(results),
        "tableS4_sequences.md": table_sequences(results),
        "table6_versions.md": table_versions(results, portal),
    }
    for name, text in tables.items():
        (out / name).write_text(text)
        if name.startswith("tableS"):   # the additional files, for a spreadsheet
            write_csv(text, out / name.replace(".md", ".csv"))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", nargs="?", type=Path,
                        default=HERE / "data" / "benchmark-run-2.2")
    parser.add_argument("--portal", type=Path, default=HERE / "data" / "portal-v0.3.0")
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args(argv)
    results = Results(args.results)
    portal = Portal(args.portal) if (args.portal / "portal_summary.txt").exists() else None
    if portal is None:
        print(f"no portal results in {args.portal}: its tables say pending", file=sys.stderr)
    figure_design(args.out / "figures")
    figure_schema(args.out / "figures")
    figure_scale(results, args.out / "figures")
    figure_loaders(results, args.out / "figures")
    figure_queries(results, args.out / "figures")
    write_tables(results, portal, args.out / "tables")
    print(f"figures in {args.out / 'figures'}, tables in {args.out / 'tables'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
