"""The gene model drawing on a feature's page: an SVG made on the server.

A drawing made here needs no script in the browser and is the same picture
whatever browser reads it. Every label is escaped, since a feature's name is
whatever its file said.
"""

from __future__ import annotations

from html import escape

from cgload.db.export_source import StoredFeature
from cgload.portal.data import EXON_TYPES, FamilyView, display_name

WIDTH = 960
LEFT = 170          # room for a row's label
RIGHT = 16
ROW = 22
TOP = 40
#: Rows drawn at most; the table below the drawing lists every member.
MAX_ROWS = 150


def _ticks(start: int, end: int, count: int = 6) -> list[int]:
    span = max(end - start, 1)
    raw = span / count
    magnitude = 10 ** (len(str(int(raw))) - 1)
    step = next(step * magnitude for step in (1, 2, 5, 10) if step * magnitude >= raw)
    first = ((start + step - 1) // step) * step
    return list(range(first, end + 1, step))


def _label(item: StoredFeature) -> str:
    text = display_name(item)
    return text if len(text) <= 24 else text[:23] + "…"


def gene_model_svg(family: FamilyView, link) -> str:
    """The family drawn to scale: one row per transcript (or per child), exons
    as boxes, coding parts taller, introns as lines; the selected feature's
    row is marked. ``link(item)`` gives the URL of a member's page."""
    top = family.top
    start, end = top.start, top.end
    for item in family.members:
        start, end = min(start, item.start), max(end, item.end)
    span = max(end - start + 1, 1)
    inner = WIDTH - LEFT - RIGHT

    def x(position: int) -> float:
        return LEFT + (position - start) / span * inner

    def box(segment_start: int, segment_end: int, height: int, y: float, css: str) -> str:
        left = x(segment_start)
        width = max(x(segment_end + 1) - left, 1.0)
        return (f'<rect class="{css}" x="{left:.1f}" y="{y - height / 2:.1f}" '
                f'width="{width:.1f}" height="{height}"/>')

    children = family.children
    rows: list[tuple[StoredFeature, list[StoredFeature]]] = []
    direct_parts = [child for child in children.get(top.row_id, [])
                    if child.feature_type in EXON_TYPES or child.feature_type == "CDS"]
    transcripts = [child for child in children.get(top.row_id, [])
                   if child not in direct_parts]
    for transcript in transcripts:
        rows.append((transcript, _parts(transcript, children)))
    if direct_parts or not transcripts:
        rows.insert(0, (top, direct_parts))
    shown = rows[:MAX_ROWS]
    height = TOP + ROW * (len(shown) + 1) + (16 if len(rows) > MAX_ROWS else 0)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" class="model" viewBox="0 0 {WIDTH} '
             f'{height}" role="img" aria-label="Gene model of {escape(_label(top))}">']
    # Scale
    parts.append(f'<line class="axis" x1="{LEFT}" y1="22" x2="{WIDTH - RIGHT}" y2="22"/>')
    for tick in _ticks(start, end):
        tx = x(tick)
        parts.append(f'<line class="axis" x1="{tx:.1f}" y1="18" x2="{tx:.1f}" y2="26"/>')
        parts.append(f'<text class="tick" x="{tx:.1f}" y="14" text-anchor="middle">'
                     f'{tick:,}</text>')
    parts.append(f'<text class="tick" x="4" y="14">{escape(top.seqid)}</text>')
    arrow = {"+": "→", "-": "←"}.get(top.strand, "")
    selected = family.selected.row_id

    for index, (owner, pieces) in enumerate(shown):
        y = TOP + ROW * index + ROW / 2
        marked = owner.row_id == selected or any(piece.row_id == selected for piece in pieces)
        css = "row selected" if marked else "row"
        parts.append(f'<a href="{escape(link(owner))}"><g class="{css}">')
        parts.append(f'<rect class="hit" x="0" y="{y - ROW / 2:.1f}" width="{WIDTH}" '
                     f'height="{ROW}"/>')
        parts.append(f'<text class="label" x="4" y="{y + 4:.1f}">'
                     f'{escape(_label(owner))} {arrow}</text>')
        parts.append(f'<line class="intron" x1="{x(owner.start):.1f}" y1="{y:.1f}" '
                     f'x2="{x(owner.end + 1):.1f}" y2="{y:.1f}"/>')
        exons = [piece for piece in pieces if piece.feature_type in EXON_TYPES]
        coding = [piece for piece in pieces if piece.feature_type == "CDS"]
        if not exons:
            # A transcript whose exons are its own rows, as GenBank and GTF
            # write them, or a feature with no parts at all.
            kind = "cds" if owner.feature_type == "CDS" else "exon"
            for segment in owner.segments:
                parts.append(box(segment.start, segment.end, 14 if kind == "cds" else 10, y,
                                 kind))
        for exon in exons:
            css_part = "utr" if exon.feature_type != "exon" else "exon"
            for segment in exon.segments:
                parts.append(box(segment.start, segment.end, 10, y, css_part))
        for cds in coding:
            for segment in cds.segments:
                parts.append(box(segment.start, segment.end, 14, y, "cds"))
        parts.append("</g></a>")
    if len(rows) > MAX_ROWS:
        parts.append(f'<text class="tick" x="{LEFT}" y="{height - 6}">'
                     f'{len(rows) - MAX_ROWS:,} more rows are listed below</text>')
    parts.append("</svg>")
    return "".join(parts)


def _parts(owner: StoredFeature, children: dict[int, list[StoredFeature]]
           ) -> list[StoredFeature]:
    """The exons, UTRs and CDS below a transcript, at any depth."""
    found, stack = [], list(children.get(owner.row_id, []))
    while stack:
        item = stack.pop()
        if item.feature_type in EXON_TYPES or item.feature_type == "CDS":
            found.append(item)
        stack.extend(children.get(item.row_id, []))
    return found
