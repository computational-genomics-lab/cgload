#!/usr/bin/env python3
"""Make the paper's figures and tables from the benchmark results.

    python3 paper/figures.py [RESULTS] [--out paper]

RESULTS is a results folder as benchmark/run.sh leaves it; by default the one
the paper reports, paper/data/benchmark-run-2.2. Figures go to OUT/figures/
(PDF for the journal, PNG to look at), tables to OUT/tables/ as Markdown.
Every number is read from RESULTS; none is typed in here. Needs matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent

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
    width, height = FULL_WIDTH, 3.55
    figure = plt.figure(figsize=(width, height))
    axis = figure.add_axes((0, 0, 1, 1))
    axis.set_xlim(0, width)
    axis.set_ylim(0, height)
    axis.axis("off")
    margin, gap = 0.06, 0.155
    w = (width - 2 * margin - 4 * gap) / 5
    xs = [margin + i * (w + gap) for i in range(5)]
    top, h = 2.47, 1.03
    box(axis, xs[0], top, w, h, "Input files",
        "GFF3 and FASTA\nGenBank (optional)\neggNOG-mapper,\nInterProScan (optional)", STEP)
    box(axis, xs[1], top, w, h, "Read, normalise",
        "Detect the dialect:\nRefSeq, funannotate,\nAUGUSTUS, GenBank.\nRefuse an unknown file.",
        STEP)
    box(axis, xs[2], top, w, h, "Structural checks",
        "Child inside its parent,\nstart not after end,\nrow inside its sequence.\n"
        "Refuse, or repair and\nrecord each repair.", CHECK)
    box(axis, xs[3], top, w, h, "One transaction",
        "Rows into SQLite or\nPostgreSQL. FASTA\nbgzipped and indexed in\na staging folder, "
        "moved\nin before the commit.", STEP)
    box(axis, xs[4], top, w, h, "End-of-load check",
        "Stored counts against\nthe independent count.\nShort: roll back, and\nthe database is as\n"
        "it was before.", CHECK)
    for a, b in zip(xs, xs[1:], strict=False):
        arrow(axis, (a + w, top + h / 2), (b, top + h / 2))

    mid, mh = 1.42, 0.72
    box(axis, xs[1], mid, 2 * w + gap, mh, "Independent count of the file",
        "Separate code that imports nothing from cgload\ncounts features by type, rows "
        "of multi-row\nfeatures and Parent= statements.", CHECK)
    box(axis, xs[4], mid, w, mh, "Database", "and data folder,\nchecksums recorded", STEP)
    arrow(axis, (xs[0] + w / 2, top), (xs[1], mid + mh / 2))
    arrow(axis, (xs[2] + w, mid + mh / 2 + 0.12), (xs[4] + 0.25, top))
    arrow(axis, (xs[4] + w / 2, top), (xs[4] + w / 2, mid + mh))

    low, lh = 0.08, 0.74
    rule = low + lh + 0.08
    axis.plot([xs[0], xs[4] + w], [rule, rule], color=AXIS, linewidth=0.6)
    axis.text(xs[1] + 0.08, rule + 0.05, "cgload verify, at any later time", ha="left",
              va="bottom", fontsize=7, fontweight="bold", color=INK)
    tiers = [("Source count", "Independent count\nagainst the database"),
             ("Parent links", "Every Parent= is a\nstored link"),
             ("Stored assembly", "FASTA and index\nagainst their checksums"),
             ("Cross-parser", "GFF3 against GenBank\n(optional)"),
             ("Reference", "File against NCBI's\npublished totals\n(optional)")]
    for x, (title, body) in zip(xs, tiers, strict=True):
        box(axis, x, low, w, lh, title, body, CHECK)
    arrow(axis, (xs[1] + 0.3, mid), (xs[0] + w / 2, low + lh))
    arrow(axis, (xs[4] + w / 2, mid), (xs[4] + w / 2, rule))
    save(figure, out, "fig1")


# --- Figure 2: the cost of a load (Q1) ------------------------------------------------

def figure_scale(results: Results, out: Path) -> None:
    panels = [
        ("load_wall_s", 1, "Load time (s)", "s"),
        ("verify_wall_s", 1, "Verify time (s)", "s"),
        ("load_peak_mb", 1 / 1024, "Peak memory of the load (GB)", "GB"),
        ("db_bytes", 1 / 1024 ** 3, "Database size (GB)", "GB"),
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
            shown = f"{value:,.0f} {unit}" if value >= 10 else f"{value:.1f} {unit}"
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
    save(figure, out, "fig2")


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
    save(figure, out, "fig3")


# --- Figure 4: query waits (Q6) -------------------------------------------------------

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
    save(figure, out, "fig4")


# --- tables ----------------------------------------------------------------------------

def md_table(header: list[str], rows: list[list[str]], align: str | None = None) -> str:
    align = align or "l" + "r" * (len(header) - 1)
    marks = {"l": ":--", "r": "--:", "c": ":-:"}
    lines = ["| " + " | ".join(header) + " |",
             "| " + " | ".join(marks[a] for a in align) + " |"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines) + "\n"


def write_csv(table: str, path: Path) -> None:
    """A Markdown table as CSV: emphasis marks off, thousands separators off."""
    lines = [line for line in table.splitlines() if line.startswith("|")]
    with path.open("w", newline="") as handle:
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
    cases = sorted({r["case"] for r in results.cases},
                   key=lambda c: -number(by[("cgload-sqlite", c)]["truth"])
                   if ("cgload-sqlite", c) in by else 0)

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
    ("Values changed, counts kept (database)", {"F06", "F07", "F08", "F10", "F12"}),
    ("Values changed, counts kept (parser)", {"P05", "P07", "P08"}),
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
        rows.append([title, ", ".join(sorted(ids)), " / ".join(sorted(expected)),
                     f"{detected} of {trials}", f"{as_expected} of {trials}", caught or "–"])
    header = ["Fault", "Faults", "Expected", "Detected", "As expected", "Caught by"]
    return md_table(header, rows, "lllrrl")


def table_faults_full(results: Results) -> str:
    rows = [[r["fault_id"], r["fault"].replace("_", " "), r["layer"], r["expected"],
             f"{number(r['detected'])} of {number(r['genomes'])}",
             f"{number(r['as_expected'])} of {number(r['genomes'])}",
             caught_by([x for x in results.fault_trials if x["fault_id"] == r["fault_id"]
                        and x["observed"] == "detected"]) or "–"] for r in results.faults]
    header = ["ID", "Fault", "Layer", "Expected", "Detected", "As expected", "Caught by"]
    return md_table(header, rows, "llllrrl")


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
                             f"{max(left) / 1024 ** 2:,.0f} MB")
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


def table_genomes(results: Results) -> str:
    rows = []
    for g in results.order:
        a = results.scale[(g["accession"], "sqlite")]
        b = results.scale[(g["accession"], "postgres")]
        outcome = "repaired" if g["sqlite_outcome"] == "loaded_with_repair" else "loaded"
        rows.append([name_md(g["organism"]), g["accession"], g["clade"], n(number(g["rows"])),
                     n(number(g["features"])), outcome,
                     f"{float(a['load_wall_s']):,.1f}", f"{float(b['load_wall_s']):,.1f}",
                     f"{float(a['load_peak_mb']) / 1024:.1f}",
                     f"{float(b['load_peak_mb']) / 1024:.1f}",
                     f"{float(a['verify_wall_s']):.1f}", f"{float(b['verify_wall_s']):.1f}",
                     f"{number(a['db_bytes']) / 1024 ** 3:.2f}",
                     f"{number(b['db_bytes']) / 1024 ** 3:.2f}",
                     f"{number(a['data_bytes']) / 1024 ** 3:.2f}"])
    header = ["Organism", "Accession", "Clade", "Rows", "Features", "cgload",
              "Load s, SQLite", "Load s, PostgreSQL", "Peak GB, SQLite",
              "Peak GB, PostgreSQL", "Verify s, SQLite", "Verify s, PostgreSQL",
              "Database GB, SQLite", "Database GB, PostgreSQL", "Data folder GB"]
    return md_table(header, rows, "lllrrl" + "r" * 9)


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


def table_versions(results: Results) -> str:
    """Versions as the server's environments recorded them at the time of the run."""
    folder = results.folder / "environment"
    released = re.search(r'__version__ = "([^"]+)"',
                         (HERE.parent / "src" / "cgload" / "__init__.py").read_text()).group(1)
    rows = [["cgload", released, "the loader presented; benchmarked at commits 8c999e2 and "
             "9b97178, the same code apart from its version string and comments"]]
    for label, file, kind, package, used in VERSIONS:
        rows.append([label, version_of((folder / file).read_text(), kind, package), used])
    chado = {r["tool_version"] for r in results.loaders if r["tool"] == "chado"}
    if chado != {"python-chado 2.3.9, Chado 1.31"}:
        raise SystemExit(f"Chado's recorded version is {chado}")
    rows.insert(-3, ["Chado schema", "1.31", "python-chado's database"])
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


def write_tables(results: Results, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    tables = {
        "table1_loaders.md": table_loaders(results),
        "table2_cases.md": table_cases(results),
        "table3_faults.md": table_faults(results),
        "table4_kills.md": table_kills(results),
        "table5_agreement.md": table_agreement(results),
        "tableS1_genomes.md": table_genomes(results),
        "tableS2_faults.md": table_faults_full(results),
        "tableS3_timing.md": table_timing(results),
        "tableS4_sequences.md": table_sequences(results),
        "table6_versions.md": table_versions(results),
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
    parser.add_argument("--out", type=Path, default=HERE)
    args = parser.parse_args(argv)
    results = Results(args.results)
    figure_design(args.out / "figures")
    figure_scale(results, args.out / "figures")
    figure_loaders(results, args.out / "figures")
    figure_queries(results, args.out / "figures")
    write_tables(results, args.out / "tables")
    print(f"figures in {args.out / 'figures'}, tables in {args.out / 'tables'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
