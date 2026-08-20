"""Milestone 9: GFF3 export.

cgload exports **its own dialect**, not the one the data came in as. Emitting
AUGUSTUS would mean synthesising the exon rows stage 5 deliberately refused to
invent; emitting RefSeq would mean fabricating `gbkey` and `Dbxref` values that
were never in the input. A file claiming to be RefSeq that is not is the same
class of dishonesty as a load reporting a success it cannot back up (D-035).

Four rules, and each one exists to keep the round trip exact:

1. **Column 2 is the stored `source_program`, verbatim.** It is per-feature
   already, so an assembly loaded from several files exports correctly with no
   extra machinery -- an AUGUSTUS gene and a RefSeq gene keep their own column 2
   in the same output file.
2. **Column 3 is `source_type`, not `feature_type`.** The same reason `verify`
   counts `source_type`: exporting the canonical term would silently rewrite
   `transcript` to `mRNA` and make the round trip lossy on any dialect that uses
   a synonym.
3. **Nothing is synthesised.** AUGUSTUS in, no exon rows out. A header comment
   records which types are absent and why, so a downstream operator learns it
   from the file rather than from a support thread.
4. **The header declares `##cgload-export`**, the version and the vocabulary
   hash, so a re-imported file is recognised as cgload output rather than
   misdetected as the dialect it was derived from.

**Anonymous features are the one place the round trip cannot be exact by
identifier.** A row with no ID -- AUGUSTUS `intron`, `start_codon` and
`stop_codon`, 175 of 350 child rows in one fixture -- has nothing to write in
column 9. Synthesising an identifier would break D-023's rule that identifiers
are verbatim from the file, so such rows are emitted without one. They re-import
correctly because both the parser and the naive counter treat an ID-less row as
its own feature (D-024).

The consequence is that **the round-trip assertion is on feature counts per
source type, never on identifier sets.** A stricter assertion fails, and the
obvious way to make it pass is to synthesise IDs, which is the thing this rule
exists to prevent.
"""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

#: Emitted as the first line after `##gff-version 3`. Detection keys on this,
#: and it must outrank column 2: a file exported from an AUGUSTUS load still
#: says `AUGUSTUS` in column 2 and would otherwise be misdetected as AUGUSTUS.
EXPORT_PRAGMA = "##cgload-export"

#: The attribute carrying how each parent edge was determined. Duplicated as a
#: literal from ``cgload.parsers.profiles`` rather than imported, to keep the
#: exporter free of parser imports -- a test asserts the two agree.
LINKAGE_ATTRIBUTE = "cgload_linkage_method"

#: Characters that must be escaped in a GFF3 attribute value. The GFF3 spec
#: requires percent-encoding of these; a raw `;` or `=` would split the field
#: and a raw tab would split the row. RefSeq already encodes `;` as `%3B` in
#: input, so a value round-tripping through cgload must come back encoded the
#: same way or the re-import splits differently from the original parse.
_MUST_ESCAPE = {
    "%": "%25",  # first, or it double-encodes the others
    ";": "%3B",
    "=": "%3D",
    "&": "%26",
    ",": "%2C",
    "\t": "%09",
    "\n": "%0A",
    "\r": "%0D",
}

_RESERVED_ORDER = ("ID", "Parent", "Name", "Alias", "Target", "Gap", "Note")


def escape_attribute(value: str) -> str:
    """Percent-encode the characters GFF3 reserves in an attribute value."""
    for character, replacement in _MUST_ESCAPE.items():
        value = value.replace(character, replacement)
    return value


@dataclass(frozen=True)
class ExportSummary:
    """What was written, so the caller can assert on it without re-reading."""

    rows: int
    features: int
    per_source_type: Counter
    absent_types: tuple[str, ...]
    anonymous_features: int


def _format_score(score: float | None) -> str:
    if score is None:
        return "."
    # AUGUSTUS writes real scores; everything else writes '.'. Trailing zeros
    # are stripped so a re-export of a re-import is byte-identical.
    return f"{score:g}"


def _attribute_column(feature, *, include_linkage: bool) -> str:
    """Build column 9, reserved keys first, then the rest alphabetically.

    Deterministic ordering matters: exporting twice must give the same bytes, or
    a diff between two exports shows spurious changes and nobody trusts it.
    """
    pairs: list[tuple[str, str]] = []
    if feature.source_id is not None:
        pairs.append(("ID", feature.source_id))
    if feature.parent_source_id is not None:
        pairs.append(("Parent", feature.parent_source_id))

    extra = {
        key: values
        for key, values in feature.attributes.items()
        if key not in {"ID", "Parent"}
    }
    for key in _RESERVED_ORDER:
        if key in extra:
            pairs.append((key, ",".join(extra.pop(key))))
    for key in sorted(extra):
        pairs.append((key, ",".join(extra[key])))

    if include_linkage and feature.linkage_method:
        # Not from the source file, so it is namespaced. A consumer can tell
        # cgload's own annotation from the file's at a glance, and a re-import
        # keeps it as an ordinary attribute rather than mistaking it for one of
        # the file's own.
        # Must be the name the normaliser reads (D-053). Writing a different
        # spelling silently discards the provenance on re-import: a containment
        # edge comes back as `explicit_parent`, which is the false claim D-032
        # was written to prevent, arriving by the export path instead.
        pairs.append((LINKAGE_ATTRIBUTE, feature.linkage_method))

    return ";".join(f"{key}={escape_attribute(value)}" for key, value in pairs)


def export_gff3(
    features: Sequence,
    *,
    version: str = "2.0.0",
    vocabulary_hash: str = "",
    declared_absent: Iterable[str] = (),
    sequence_regions: dict[str, int] | None = None,
    include_linkage: bool = True,
) -> tuple[str, ExportSummary]:
    """Serialise features to cgload-dialect GFF3.

    ``declared_absent`` names types the source dialect never provides -- passed
    through from the profile's ``provides`` declaration -- so the file records
    the difference between "the input had no exons" and "the exons were lost".
    """
    lines: list[str] = ["##gff-version 3", f"{EXPORT_PRAGMA} {version}"]
    if vocabulary_hash:
        lines.append(f"##cgload-vocabulary {vocabulary_hash}")

    for name, length in sorted((sequence_regions or {}).items()):
        lines.append(f"##sequence-region {name} 1 {length}")

    absent = tuple(sorted(set(declared_absent)))
    if absent:
        lines.append(
            "# The source dialect does not emit these types, so their absence "
            "is a property of the input and not a loss during loading:"
        )
        lines.append("#   " + ", ".join(absent))

    programs = sorted({f.source_program for f in features})
    if programs:
        lines.append("# source programs present: " + ", ".join(programs))

    rows = 0
    anonymous = 0
    per_type: Counter = Counter()

    # Stable order: by sequence, then start, then source_id so two exports of
    # the same data are byte-identical.
    ordered = sorted(
        features,
        key=lambda f: (f.seqid, f.start, f.end, f.source_type, f.source_id or ""),
    )
    for feature in ordered:
        per_type[feature.source_type] += 1
        if feature.source_id is None:
            anonymous += 1
        attributes = _attribute_column(feature, include_linkage=include_linkage)
        for segment in sorted(feature.segments, key=lambda s: s.start):
            rows += 1
            lines.append(
                "\t".join(
                    (
                        feature.seqid,
                        feature.source_program,
                        feature.source_type,
                        str(segment.start),
                        str(segment.end),
                        _format_score(segment.score),
                        feature.strand,
                        "." if segment.phase is None else str(segment.phase),
                        attributes,
                    )
                )
            )

    text = "\n".join(lines) + "\n"
    summary = ExportSummary(
        rows=rows,
        features=len(ordered),
        per_source_type=per_type,
        absent_types=absent,
        anonymous_features=anonymous,
    )
    return text, summary


def is_cgload_export(text: str) -> bool:
    """Does this text declare itself as cgload output?

    Detection must consult this **before** column 2. An export of an AUGUSTUS
    load still carries `AUGUSTUS` in column 2, and the augustus profile would
    otherwise claim it -- then apply AUGUSTUS's rules to a file that no longer
    follows them.
    """
    for line in text.splitlines():
        if line.startswith(EXPORT_PRAGMA):
            return True
        if not line.startswith("#"):
            return False
    return False


def round_trip_counts(text: str) -> Counter:
    """Count features per source type from exported text, naively.

    Deliberately independent of the parser, for the same reason `verify` is
    (D-008): a round-trip check that used the parser on both sides could only
    prove the parser is self-consistent, not that the export is complete.
    """
    seen: set[tuple[str, str]] = set()
    anonymous: Counter = Counter()
    for line in text.splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) < 9:
            continue
        found = re.search(r"(?:^|;)ID=([^;]*)", columns[8])
        if found:
            seen.add((columns[2], found.group(1)))
        else:
            anonymous[columns[2]] += 1
    counts = Counter(source_type for source_type, _identifier in seen)
    counts.update(anonymous)
    return counts
