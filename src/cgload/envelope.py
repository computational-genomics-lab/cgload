"""Milestone 10: parent envelopes that contradict their own children.

A parent feature's coordinates must span every child it claims. When they do
not, one of the two is wrong -- and cgload can tell which, because the child
coordinates are the ones the annotation pipeline computed directly while the
parent envelope is derived.

**Why this check exists.** Every other structural rule in cgload is per-row:
start before end, strand in the vocabulary, coordinates inside the sequence.
A row can satisfy all of them and still be wrong, and nothing downstream can
tell. Measured on *Zea mays* B73 (GCF_902167145.1, 1.2 M features): 36
structural problems, of which **1 fails a per-row check and 34 do not**. The 34
load cleanly with envelopes wrong by up to 281 kb. Only the parent-child
relationship exposes them.

**The provider bug.** All 34 are organellar -- mitochondrial (`ZeamMp*`) and
chloroplast (`ZemaCp*`) -- across seven trans-spliced genes. NCBI's GFF3
converter emits the parent envelope from **one segment** of the trans-spliced
set rather than from their union, and which segment varies: for
`gene-ZeamMp186` it is the last, for `gene-ZeamMp016` the first, for
`gene-ZeamMp071` a middle one. The same assembly's GenBank file has the
coordinates right, which is how we know the children are correct.

Note that the affected parents carry strand `-` (six), `?` (two) and `+` (one).
An early hypothesis that `?` marked the affected rows was wrong and would have
scoped a fix to a third of the cases.

**Repair runs bottom-up, and this is not optional.** For `ZeamMp186` both the
mRNA and the gene envelope are wrong. Repairing the gene from its immediate
child -- the still-broken mRNA -- yields `50,490..267,232`, which is plausible
and wrong. Repairing leaves-upward yields `50,490..548,772`, which matches the
GenBank join exactly. A single top-down pass produces a number that looks fixed
and is not, which is worse than refusing.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace


class EnvelopeError(ValueError):
    """A parent envelope contradicts its children and was not repaired."""


@dataclass(frozen=True)
class Violation:
    """One child lying outside the parent that claims it."""

    parent_id: str
    parent_type: str
    parent_strand: str
    parent_start: int
    parent_end: int
    child_id: str | None
    child_type: str
    child_start: int
    child_end: int

    #: Set when the feature contradicts itself rather than its parent, i.e. the
    #: source row wrote start > end (D-056). Rendered differently, because the
    #: parent/child wording describes nothing in that case.
    @property
    def is_self_contradiction(self) -> bool:
        return self.child_type == "(self)"

    def describe(self) -> str:
        if self.is_self_contradiction:
            return (
                f"{self.parent_type} {self.parent_id} is written "
                f"{self.parent_start:,}..{self.parent_end:,}, which runs backwards, and "
                f"has no children to recompute an extent from"
            )
        child = self.child_id or "<anonymous>"
        return (
            f"{self.child_type} {child} {self.child_start:,}..{self.child_end:,} "
            f"lies outside {self.parent_type} {self.parent_id} "
            f"{self.parent_start:,}..{self.parent_end:,} (strand {self.parent_strand!r})"
        )


def _is_inverted(feature) -> bool:
    """Did the source file write start > end on any of this feature's rows?

    Only ever True when the file was tokenized with ``tolerate_inverted``, i.e.
    under ``--repair-envelopes``. Such a row is the same provider bug as a
    too-small envelope, arriving in a form that a per-row check happens to
    catch: *Zea mays* `ZeamMp017` is written `691,776..267,232`, and both
    numbers are segment coordinates from the trans-spliced set rather than its
    extent. Its children span 122,146..548,772, which is the correct answer and
    is recoverable exactly as for `ZeamMp186` (D-056).
    """
    return any(getattr(segment, "inverted", False) for segment in feature.segments)


@dataclass(frozen=True)
class Repair:
    """One envelope recomputed from the features it contains."""

    feature_id: str
    feature_type: str
    seqid: str
    old_start: int
    old_end: int
    new_start: int
    new_end: int

    @property
    def widened_by(self) -> int:
        return (self.old_start - self.new_start) + (self.new_end - self.old_end)

    def describe(self) -> str:
        return (
            f"{self.feature_type} {self.feature_id} on {self.seqid}: "
            f"{self.old_start:,}..{self.old_end:,} -> "
            f"{self.new_start:,}..{self.new_end:,} "
            f"(+{self.widened_by:,} bp, recomputed from children)"
        )


def _children_by_parent(features: Sequence) -> dict[str, list]:
    index: dict[str, list] = defaultdict(list)
    for feature in features:
        if feature.parent_source_id:
            index[feature.parent_source_id].append(feature)
    return index


def _depth_order(features: Sequence) -> list:
    """Features ordered leaves-first.

    Depth is distance from a root. Sorting descending processes the deepest
    features first, so by the time a parent is examined every descendant has
    already been repaired. A cycle -- which the parser should already have
    refused -- degrades to input order rather than looping.
    """
    by_id = {f.source_id: f for f in features if f.source_id}
    depth_cache: dict[str, int] = {}

    def depth(feature, seen: frozenset = frozenset()) -> int:
        key = feature.source_id
        if key in depth_cache:
            return depth_cache[key]
        parent_id = feature.parent_source_id
        if not parent_id or parent_id not in by_id or parent_id in seen:
            value = 0
        else:
            value = 1 + depth(by_id[parent_id], seen | {key})
        if key:
            depth_cache[key] = value
        return value

    return sorted(features, key=depth, reverse=True)


def find_violations(features: Sequence) -> list[Violation]:
    """Every child lying outside the parent that claims it."""
    by_id = {f.source_id: f for f in features if f.source_id}
    violations: list[Violation] = []

    for feature in features:
        # An inverted row contradicts itself before any parent is consulted.
        if _is_inverted(feature):
            violations.append(
                Violation(
                    parent_id=feature.source_id or "<anonymous>",
                    parent_type=feature.source_type,
                    parent_strand=feature.strand,
                    parent_start=feature.segments[0].start,
                    parent_end=feature.segments[0].end,
                    child_id=None,
                    child_type="(self)",
                    child_start=feature.segments[0].start,
                    child_end=feature.segments[0].end,
                )
            )

    for feature in features:
        parent_id = feature.parent_source_id
        if not parent_id or parent_id not in by_id:
            continue
        parent = by_id[parent_id]
        if feature.start < parent.start or feature.end > parent.end:
            violations.append(
                Violation(
                    parent_id=parent_id,
                    parent_type=parent.source_type,
                    parent_strand=parent.strand,
                    parent_start=parent.start,
                    parent_end=parent.end,
                    child_id=feature.source_id,
                    child_type=feature.source_type,
                    child_start=feature.start,
                    child_end=feature.end,
                )
            )
    return violations


def repair_envelopes(features: Sequence) -> tuple[list, list[Repair]]:
    """Recompute contradicted envelopes from their children, leaves upward.

    Returns the repaired features and a record of every change. Features whose
    envelope already contains their children are returned untouched and are not
    reported.

    **A multi-segment parent is never repaired.** Widening it would mean picking
    which segment to extend, and the file gives no basis for that choice --
    exactly the kind of guess this project refuses. Such a case remains a
    violation and, under the default policy, refuses the load.
    """
    working = list(features)
    by_id = {f.source_id: f for f in working if f.source_id}
    children = _children_by_parent(working)
    repairs: list[Repair] = []

    for feature in _depth_order(working):
        key = feature.source_id
        if not key:
            continue
        kids = children.get(key)
        if not kids:
            continue

        required_start = min(child.start for child in kids)
        required_end = max(child.end for child in kids)
        inverted = _is_inverted(feature)
        if (
            not inverted
            and feature.start <= required_start
            and feature.end >= required_end
        ):
            continue

        if len(feature.segments) != 1:
            # Ambiguous: no basis for choosing which segment to widen.
            continue

        segment = feature.segments[0]
        old_start, old_end = segment.start, segment.end
        if inverted:
            # The written pair is not an extent at all -- both numbers are
            # segment coordinates -- so widening from them is meaningless. Take
            # the children alone.
            new_start, new_end = required_start, required_end
        else:
            new_start = min(feature.start, required_start)
            new_end = max(feature.end, required_end)

        repaired = replace(
            feature,
            segments=[
                replace(segment, start=new_start, end=new_end, inverted=False)
            ],
        )
        # Rebind so a shallower feature examined later sees the repaired child.
        by_id[key] = repaired
        working[working.index(feature)] = repaired
        for _parent_id, group in children.items():
            for index, child in enumerate(group):
                if child is feature:
                    group[index] = repaired

        repairs.append(
            Repair(
                feature_id=key,
                feature_type=feature.source_type,
                seqid=feature.seqid,
                old_start=old_start,
                old_end=old_end,
                new_start=new_start,
                new_end=new_end,
            )
        )

    return working, repairs


def apply_policy(
    features: Sequence, *, repair: bool
) -> tuple[list, list[Repair], list[Violation]]:
    """Detect, optionally repair, then re-check.

    The re-check is the point. A repair that does not clear the violation means
    something other than a truncated envelope is wrong, and proceeding on a
    feature that still contradicts its children would be exactly the silent
    partial success this project exists to prevent.
    """
    violations = find_violations(features)
    if not violations:
        return list(features), [], []

    if not repair:
        return list(features), [], violations

    repaired, repairs = repair_envelopes(features)
    remaining = find_violations(repaired)
    return repaired, repairs, remaining


def format_report(repairs: Iterable[Repair], remaining: Iterable[Violation]) -> str:
    lines: list[str] = []
    repairs = list(repairs)
    remaining = list(remaining)
    if repairs:
        lines.append(f"{len(repairs)} envelope(s) recomputed from their children:")
        lines.extend(f"    {repair.describe()}" for repair in repairs)
    if remaining:
        lines.append(f"{len(remaining)} envelope violation(s) could not be repaired:")
        lines.extend(f"    {violation.describe()}" for violation in remaining)
    return "\n".join(lines)
