"""Layer 1: the tab-delimited tokenizer (D-013).

This module knows the *physical* format — nine tab-separated columns, `#`
comments, GFF3 attribute syntax — and nothing about pipeline conventions. It
does not know what RefSeq or AUGUSTUS are, does not map feature types, and does
not resolve parents. Everything dialect-specific lives in a profile.

Deliberately kept dumb, because a tokenizer that started making decisions would
give the dialect profiles somewhere to hide.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import unquote

#: A GFF3 line has exactly nine columns.
COLUMNS = 9

#: Directives worth keeping. `##sequence-region` gives a declared length that
#: the RefSeq profile cross-checks against the FASTA (finding 11).
SEQUENCE_REGION = re.compile(r"^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)\s*$")


class TokenizeError(ValueError):
    """A line cannot be read as GFF3. Always names the file and line number."""

    def __init__(self, path: Path, line_number: int, message: str) -> None:
        super().__init__(f"{path.name}:{line_number}: {message}")
        self.path = path
        self.line_number = line_number


@dataclass(frozen=True)
class RawFeature:
    """One data line, split and decoded. No interpretation applied.

    ``attributes`` maps a key to a list of values, because GFF3 comma-separates
    multiple values within one key (finding 10) and a single-value assumption
    silently discards all but the first.
    """

    line_number: int
    seqid: str
    source: str
    type: str
    start: int
    end: int
    score: float | None
    strand: str
    phase: int | None
    attributes: dict[str, list[str]]
    #: The file wrote start > end. Only ever True when tokenize() was called
    #: with tolerate_inverted=True, i.e. under --repair-envelopes. The
    #: coordinates are stored exactly as written; nothing is swapped here,
    #: because a swap would invent an extent the file never claimed. The
    #: envelope repair recomputes it from the feature's children instead
    #: (D-056).
    inverted: bool = False

    def first(self, key: str) -> str | None:
        """The first value for ``key``, matched case-insensitively.

        funannotate writes ``DBxref`` where RefSeq writes ``Dbxref``
        (finding 10), so a case-sensitive lookup silently misses one of them.
        """
        lowered = key.lower()
        for name, values in self.attributes.items():
            if name.lower() == lowered and values:
                return values[0]
        return None

    def all(self, key: str) -> list[str]:
        lowered = key.lower()
        for name, values in self.attributes.items():
            if name.lower() == lowered:
                return list(values)
        return []

    @property
    def id(self) -> str | None:
        """``ID``, or None. None is meaningful: AUGUSTUS emits intron and codon
        rows with only a Parent, and each such row is its own feature (D-024)."""
        return self.first("ID")

    @property
    def parents(self) -> list[str]:
        """Every parent named. GFF3 permits `Parent=a,b`, and the comma split
        already happened in ``parse_attributes``, so this must read the whole
        list -- taking only the first value would silently discard a parent and
        attach the feature to one arbitrary branch of the graph."""
        return self.all("Parent")


@dataclass
class TokenizedFile:
    features: list[RawFeature] = field(default_factory=list)
    #: `##sequence-region` declarations: source_id -> (start, end).
    declared_regions: dict[str, tuple[int, int]] = field(default_factory=dict)
    #: Every `#` line, in order, with its line number. The AUGUSTUS profile
    #: reads protein blocks out of these (D-025); other profiles ignore them.
    comments: list[tuple[int, str]] = field(default_factory=list)
    #: Distinct column-2 values. Detection reads this (D-013).
    sources: set[str] = field(default_factory=set)
    data_line_count: int = 0


def _open(path: Path):
    """Open plain or gzipped, decided by magic bytes rather than by name.

    NCBI serves annotation gzipped, and the uploaded *F. graminearum* files
    arrived named `.gff` and `.gbff` while both were gzip. The GenBank tokenizer
    already read the magic bytes; this one did not, so a straight NCBI download
    failed with a decoding error naming no cause (D-035 applied consistently).
    """
    with open(path, "rb") as probe:
        if probe.read(2) == b"\x1f\x8b":
            import gzip

            return gzip.open(path, "rt", encoding="utf-8", errors="strict")
    return path.open("r", encoding="utf-8", errors="strict")


def parse_attributes(raw: str, path: Path, line_number: int) -> dict[str, list[str]]:
    """Split column 9, then percent-decode.

    Order matters and is the whole point (finding 9): RefSeq writes
    ``strain=PH-1%3B NRRL 31084``, where ``%3B`` is a semicolon. Decoding before
    splitting would break the field apart; not decoding at all stores a
    corrupted value.
    """
    attributes: dict[str, list[str]] = {}
    if raw.strip() in ("", "."):
        return attributes

    for chunk in raw.split(";"):
        if not chunk.strip():
            # funannotate terminates every line with a trailing ';' (finding 10).
            continue
        if "=" not in chunk:
            raise TokenizeError(
                path, line_number, f"attribute {chunk.strip()!r} has no '=' separator"
            )
        key, _, value = chunk.partition("=")
        key = unquote(key.strip())
        values = [unquote(part) for part in value.split(",")]
        attributes.setdefault(key, []).extend(values)
    return attributes


def _integer(text: str, path: Path, line_number: int, what: str) -> int:
    try:
        return int(text)
    except ValueError as exc:
        raise TokenizeError(path, line_number, f"{what} {text!r} is not an integer") from exc


def is_genbank(path: Path) -> bool:
    """Does this file begin with a GenBank LOCUS line?

    Checked by content, not by extension: NCBI serves `.gbff`, `.gb`, `.gbk`
    and gzipped variants of each, and users rename them.
    """
    import gzip as _gzip

    with open(path, "rb") as probe:
        magic = probe.read(2)
    opener = _gzip.open if magic == b"\x1f\x8b" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                return stripped.startswith("LOCUS")
    return False


def tokenize(path: Path, *, tolerate_inverted: bool = False) -> TokenizedFile:
    """Read an annotation file into raw features.

    Dispatches on content: a GenBank flat file goes to the GenBank tokenizer,
    which returns the same RawFeature shape so everything downstream is shared.

    The comment rule is "first character is ``#``", applied wherever the line
    appears — AUGUSTUS writes whole evidence blocks between feature rows
    (finding 13). D-008's naive counting pass must use the identical rule.
    """
    if is_genbank(path):
        from cgload.parsers.genbank import tokenize_genbank

        return tokenize_genbank(path)

    path = Path(path)
    result = TokenizedFile()

    with _open(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.startswith("#"):
                result.comments.append((line_number, line))
                match = SEQUENCE_REGION.match(line)
                if match:
                    result.declared_regions[match.group(1)] = (
                        int(match.group(2)),
                        int(match.group(3)),
                    )
                continue

            columns = line.split("\t")
            if len(columns) != COLUMNS:
                raise TokenizeError(
                    path,
                    line_number,
                    f"expected {COLUMNS} tab-separated columns, found {len(columns)}",
                )

            start = _integer(columns[3], path, line_number, "start")
            end = _integer(columns[4], path, line_number, "end")
            if start < 1:
                raise TokenizeError(path, line_number, f"start {start} is below 1")
            inverted = end < start
            if inverted and not tolerate_inverted:
                raise TokenizeError(
                    path, line_number, f"end {end} precedes start {start}"
                )

            strand = columns[6]
            if strand not in ("+", "-", ".", "?"):
                # '?' is valid GFF3: stranded, but the strand is unknown. Maize's
                # trans-spliced mitochondrial mRNAs use it, because their exons lie
                # on opposite strands and no single value is correct.
                raise TokenizeError(
                    path, line_number, f"strand {strand!r} is not +, -, . or ?"
                )

            phase: int | None = None
            if columns[7] != ".":
                phase = _integer(columns[7], path, line_number, "phase")
                if phase not in (0, 1, 2):
                    raise TokenizeError(path, line_number, f"phase {phase} is not 0, 1 or 2")

            score: float | None = None
            if columns[5] != ".":
                try:
                    score = float(columns[5])
                except ValueError as exc:
                    raise TokenizeError(
                        path, line_number, f"score {columns[5]!r} is not a number"
                    ) from exc

            result.features.append(
                RawFeature(
                    line_number=line_number,
                    seqid=unquote(columns[0]),
                    source=columns[1],
                    type=columns[2],
                    start=start,
                    end=end,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes=parse_attributes(columns[8], path, line_number),
                    inverted=inverted,
                )
            )
            result.sources.add(columns[1])
            result.data_line_count += 1

    if not result.features:
        raise TokenizeError(path, 0, "contains no feature rows")

    return result


def iter_data_lines(path: Path) -> Iterator[tuple[int, list[str]]]:
    """Minimal independent pass over the data lines.

    Used by D-008's naive counting pass. Shares the comment rule with
    ``tokenize`` and nothing else: no attribute parsing, no validation, no
    profile. Kept in this module so the two comment rules cannot drift, which is
    the single point of contact D-008 permits.
    """
    with _open(Path(path)) as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n")
            if not stripped.strip() or stripped.startswith("#"):
                continue
            yield line_number, stripped.split("\t")
