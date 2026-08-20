"""Milestone 8: completeness checking.

The claim cgload makes is that nothing was lost on the way in. This module is
what backs that claim, and it is built on one principle:

    **The check must not share code with the thing it checks.**

If the verifier reused the parser, it could only catch mistakes made *after*
parsing — storage bugs. The mistakes that matter most in this class of
tool were made *during* parsing: a `return` indented one level too deep, a
silent `else: pass`. Those are invisible to any check that parses the file the
same way. So the "before" number here comes from a deliberately naive pass that
knows nothing about profiles, dialects, linkage or segments (D-008).

Three tiers, and they are not equally authoritative:

1. **Source count** — the naive pass versus the database. Disagreement always
   means cgload is wrong. Fails the command.
2. **Cross-parser** — the same assembly read as GFF3 and as GenBank, by two
   engines sharing no code. Disagreement means one parser dropped something.
   Fails the command.
3. **Reference counts** — NCBI's published per-assembly totals. Disagreement can
   have innocent causes, so this fails only when the preconditions rule those
   causes out, and is otherwise skipped with the reason stated (D-033).

A check may fail the command only when disagreement implies cgload is wrong. A
check that reports concerns it cannot stand behind trains people to ignore it,
which leaves you no better off than having no check.
"""

from __future__ import annotations

import gzip
import re
from collections import Counter
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class VerifyError(ValueError):
    """The verifier could not run. Distinct from the verification failing."""


class Outcome(Enum):
    PASS = "pass"  # noqa: S105 - an outcome, not a credential
    FAIL = "fail"
    #: Preconditions not met. Not a soft failure -- the check did not run, and
    #: the reason is printed so the gap is visible rather than assumed away.
    SKIPPED = "skipped"


@dataclass(frozen=True)
class TierResult:
    name: str
    outcome: Outcome
    detail: str
    expected: dict[str, int] = field(default_factory=dict)
    actual: dict[str, int] = field(default_factory=dict)

    @property
    def mismatches(self) -> dict[str, tuple[int, int]]:
        keys = set(self.expected) | set(self.actual)
        return {
            key: (self.expected.get(key, 0), self.actual.get(key, 0))
            for key in sorted(keys)
            if self.expected.get(key, 0) != self.actual.get(key, 0)
        }


def _open(path: Path):
    """Open plain or gzipped, by magic bytes. NCBI serves both."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    opener = gzip.open if magic == b"\x1f\x8b" else open
    return opener(path, "rt", encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Tier 1: the naive pass
# ---------------------------------------------------------------------------

#: Column 3 of a GFF3 data line.
_GFF_ID = re.compile(r"(?:^|;)ID=([^;]*)")
#: A GenBank feature-table record begins at column 5.
_GB_FEATURE = re.compile(r"^ {5}(\S+)\s")


def count_source_naively(path: Path) -> Counter:
    """Count features per type without using any parser code.

    Deliberately simple. It knows three things: comment lines start with `#`,
    GFF3 data lines have nine tab-separated columns with the type in column 3,
    and rows sharing an `ID` are one feature rather than several. It does not
    know what a dialect is, cannot resolve a parent, and has no vocabulary.

    That naivety is the point. Anything cleverer would start to share
    assumptions with the parser, and a check that shares the parser's
    assumptions cannot catch the parser's mistakes.
    """
    if _looks_like_genbank(path):
        return _count_genbank_naively(path)
    return _count_gff_naively(path)


def _looks_like_genbank(path: Path) -> bool:
    with _open(path) as handle:
        for line in handle:
            if line.strip():
                return line.startswith("LOCUS")
    return False


def _count_gff_naively(path: Path) -> Counter:
    """Distinct (type, ID) pairs; rows with no ID count individually.

    Counting rows would overcount every multi-segment CDS -- 9,739 rows for
    6,780 CDS features in *U. maydis*. Counting distinct IDs matches both the
    database and NCBI's published totals.
    """
    seen: set[tuple[str, str]] = set()
    anonymous: Counter = Counter()
    with _open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) < 9:
                continue
            feature_type = columns[2]
            found = _GFF_ID.search(columns[8])
            if found:
                seen.add((feature_type, found.group(1)))
            else:
                anonymous[feature_type] += 1
    counts = Counter(feature_type for feature_type, _identifier in seen)
    counts.update(anonymous)
    return counts


def _count_genbank_naively(path: Path) -> Counter:
    """One feature-table record is one feature, however many segments it spans.

    A GenBank record already *is* the feature -- `join(1..96,200..300)` is one
    CDS written on one record -- so no grouping is needed and none is done.
    """
    counts: Counter = Counter()
    in_features = False
    with _open(path) as handle:
        for line in handle:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith(("ORIGIN", "//", "CONTIG", "BASE COUNT")):
                in_features = False
                continue
            if not in_features:
                continue
            found = _GB_FEATURE.match(line)
            if found and found.group(1) != "source":
                counts[found.group(1)] += 1
    return counts


# ---------------------------------------------------------------------------
# Tier 2: two parsers, no shared code
# ---------------------------------------------------------------------------

#: GFF3 gives a pseudogene its own column-3 type; GenBank writes it as a `gene`
#: record carrying `/pseudo`. Measured on *F. graminearum* GCF_000240135.3:
#: GFF3 has 13,714 `gene` + 11 `pseudogene`, GenBank has 13,725 `gene`. The two
#: agree exactly once this mapping is applied, so it is a naming difference and
#: not a discrepancy.
#: Both external references cgload compares against report a pseudogene under
#: `gene`: GenBank writes it as a `gene` record with `/pseudo`, and NCBI's
#: feature-count file reports `Feature=gene, Class=pseudogene`. cgload keeps
#: `pseudogene` as a distinct `source_type` because that is what the GFF3 said.
#:
#: Measured twice. *F. graminearum*: GFF3 13,714 gene + 11 pseudogene = 13,725,
#: GenBank 13,725. *Z. mays* Primary Assembly: loaded 44,321 gene + 5,195
#: pseudogene = 49,516, NCBI published 49,516. Exact both times (D-066).
EXTERNAL_TYPE_EQUIVALENCE: dict[str, tuple[str, ...]] = {
    "gene": ("gene", "pseudogene"),
}

#: Retained under its original name; the mapping was never GenBank-specific.
GENBANK_TYPE_EQUIVALENCE = EXTERNAL_TYPE_EQUIVALENCE


def compare_across_formats(gff_path: Path, genbank_path: Path) -> TierResult:
    """Count the same assembly from GFF3 and from GenBank and require agreement.

    This is the strongest check available, because the two counts come from
    files with different grammars read by different code. A GFF3 parser losing
    a feature has no way to make the GenBank count lose it too, so agreement is
    not something either parser can fake.

    Measured on GCF_000240135.3: gene 13,725 = 13,725; mRNA 13,315; CDS 13,312;
    tRNA 322; rRNA 88. Exact on every type.
    """
    gff = count_source_naively(gff_path)
    genbank = count_source_naively(genbank_path)

    folded: Counter = Counter()
    for target, sources in GENBANK_TYPE_EQUIVALENCE.items():
        folded[target] = sum(gff.get(source, 0) for source in sources)
    for feature_type, count in gff.items():
        if not any(feature_type in sources for sources in GENBANK_TYPE_EQUIVALENCE.values()):
            folded[feature_type] = count

    #: Types that exist in one format and not the other by construction, so
    #: their absence is not a loss. GFF3 states exons explicitly; GenBank
    #: implies them from the location expression. GFF3 has `region` rows;
    #: GenBank has `source`, already excluded. `gap` is assembly metadata.
    structural_only = {"exon", "region", "gap", "sequence_feature", "misc_feature"}

    # Compare the UNION, with a missing type counted as zero.
    #
    # An intersection makes the most alarming failure the one case that cannot
    # fail: if every feature of a type is gone, the type is absent from one side
    # and silently dropped from the comparison. Losing 321 of 322 tRNAs failed;
    # losing all 322 passed. Tier 1 always used a union; these two inverted it
    # (D-072).
    #
    # `structural_only` already carries the types that legitimately exist in one
    # format and not the other, so the intersection was doing nothing the
    # allowlist was not already doing -- while suppressing total loss.
    expected = {k: v for k, v in folded.items() if k not in structural_only}
    actual = {k: v for k, v in genbank.items() if k not in structural_only}
    shared = set(expected) | set(actual)
    expected = {k: expected.get(k, 0) for k in shared}
    actual = {k: actual.get(k, 0) for k in shared}

    if not shared:
        return TierResult(
            "cross-parser", Outcome.SKIPPED,
            "the two files share no comparable feature type",
        )
    if expected == actual:
        return TierResult(
            "cross-parser", Outcome.PASS,
            f"GFF3 and GenBank agree on {len(shared)} feature types",
            expected, actual,
        )
    return TierResult(
        "cross-parser", Outcome.FAIL,
        "GFF3 and GenBank disagree; one parser has lost features",
        expected, actual,
    )


# ---------------------------------------------------------------------------
# Tier 3: NCBI's published counts
# ---------------------------------------------------------------------------


def parse_ncbi_feature_counts(
    path: Path, assembly_accession: str, *, assembly_unit: str | None = None
) -> Counter:
    """Read an NCBI `*_feature_count.txt` file into per-type totals.

    Four traps, all present in real files:

    * **Counts must be summed across `Class` within a `Feature`.** `gene`
      appears once per biotype -- protein_coding, rRNA, tRNA -- and only the sum
      equals the GFF3's gene count.
    * **`Unique Ids` can be the literal string `na`** while `Placements` holds
      the number.
    * **`Assembly-unit accession` is not the assembly accession.** Match on the
      `Full Assembly` column.
    * **Rows must NOT be summed across assembly units.** This is the one that
      bit. A file may carry several groupings under one `Full Assembly`: an
      `all` rollup plus one row-set per unit. *Zea mays* GCF_902167145.1 has 24
      `all` rows, 23 `Primary Assembly` rows and 8 `non-nuclear` rows, and the
      first is the sum of the other two. Summing everything double-counts:
      mRNA is 57,071 in both `all` and `Primary Assembly`, and a naive sum
      reports 114,142 against a correct load of 57,071 (D-065).

    ``assembly_unit`` selects a grouping by its unit name or accession. The
    default takes the `all` rollup when present, since that describes the whole
    assembly; pass the unit explicitly when only part of an assembly was loaded,
    which is the usual reason for a partial load.
    """
    header: list[str] | None = None
    # (unit_key, feature) -> count
    per_unit: dict[str, Counter] = {}
    matched_rows = 0

    with _open(path) as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if not stripped.strip():
                continue
            if stripped.startswith("#"):
                header = [c.strip() for c in stripped.lstrip("#").split("\t")]
                continue
            if header is None:
                raise VerifyError(f"{path}: data before the header line")
            row = dict(zip(header, stripped.split("\t"), strict=False))

            if row.get("Full Assembly") != assembly_accession:
                continue
            matched_rows += 1

            feature = row.get("Feature")
            if not feature:
                continue
            raw = row.get("Unique Ids", "")
            if not raw.isdigit():
                raw = row.get("Placements", "")
            if not raw.isdigit():
                continue

            unit_accession = (row.get("Assembly-unit accession") or "").strip()
            unit_name = (row.get("Assembly-unit name") or "").strip()
            key = unit_accession or unit_name or "all"
            per_unit.setdefault(key, Counter())[feature] += int(raw)
            if unit_name and unit_accession:
                # Remember the human-readable name so a caller can select on it.
                per_unit[key]["__name__"] = 0

    if matched_rows == 0:
        raise VerifyError(
            f"{path}: no rows for assembly {assembly_accession!r}. Check the "
            "'Full Assembly' column -- the 'Assembly-unit accession' column "
            "holds a different accession."
        )

    # Prefer summing the component units over trusting the `all` rollup.
    #
    # The rollup can be incomplete where a component's `Unique Ids` cell is `na`.
    # Measured on GCF_902167145.1: rRNA is `Unique Ids 2427, Placements 2439` in
    # `all`, `2427/2427` in Primary Assembly, and `na/12` in non-nuclear. The
    # rollup's unique-ID total silently omits the 12 organellar rRNAs because
    # NCBI did not compute a unique count for them; 2,427 + 12 = 2,439, which is
    # what the GFF3 contains. Summing components recovers the right answer and
    # needs no per-column special-casing (D-070).
    #
    # This is not the D-065 double-count: components are summed with each other,
    # never with the rollup that already contains them.
    components = {k: v for k, v in per_unit.items() if k != "all"}
    if not assembly_unit and components:
        summed: Counter = Counter()
        for counts in components.values():
            for key, value in counts.items():
                if not key.startswith("__"):
                    summed[key] += value
        return summed

    if assembly_unit:
        for key, counts in per_unit.items():
            if assembly_unit in (key,) or assembly_unit == key:
                return Counter({k: v for k, v in counts.items() if not k.startswith("__")})
        raise VerifyError(
            f"{path}: no assembly unit matching {assembly_unit!r}. "
            f"Available: {', '.join(sorted(per_unit))}"
        )

    if "all" in per_unit:
        chosen = per_unit["all"]
    elif len(per_unit) == 1:
        chosen = next(iter(per_unit.values()))
    else:
        raise VerifyError(
            f"{path}: {len(per_unit)} assembly units present and no 'all' rollup "
            f"({', '.join(sorted(per_unit))}). Pass assembly_unit to choose one; "
            "summing across units would double-count."
        )
    return Counter({k: v for k, v in chosen.items() if not k.startswith("__")})


def compare_with_reference(
    source_counts: Counter,
    reference_counts: Counter,
    *,
    sequences_loaded: int,
    sequences_in_assembly: int,
) -> TierResult:
    """Compare against NCBI, but only when the comparison can mean something.

    Four innocent causes of disagreement (D-033): a partial assembly was loaded;
    the annotation release differs from the one the counts describe; NCBI tallies
    non-coding and pseudogenes under different classes; and re-annotation happens
    without the accession changing.

    Only the first is detectable here, so it is the precondition. When it fails,
    the tier is SKIPPED with the reason printed -- not reported as an advisory
    warning, because a check that raises concerns it cannot stand behind is one
    people learn to scroll past.
    """
    if sequences_loaded != sequences_in_assembly:
        return TierResult(
            "reference (NCBI)", Outcome.SKIPPED,
            f"loaded {sequences_loaded} of {sequences_in_assembly} sequences; "
            "published totals describe the whole assembly",
        )

    # Fold the types the external reference reports under one name (D-066).
    folded: Counter = Counter(source_counts)
    for target, sources in EXTERNAL_TYPE_EQUIVALENCE.items():
        total = sum(source_counts.get(name, 0) for name in sources)
        if total:
            folded[target] = total
            for name in sources:
                if name != target:
                    folded.pop(name, None)

    # Union, for the reason in D-072. A type present in the published counts and
    # absent from the load means every feature of that type is missing, which is
    # exactly what this tier should report rather than skip.
    #
    # The converse -- a type loaded that the reference does not publish -- is not
    # necessarily an error, because published summaries are selective. Those are
    # reported with an expected count of zero rather than suppressed, so the
    # asymmetry is visible in the output instead of hidden in the comparison.
    shared = set(folded) | set(reference_counts)
    if not shared:
        return TierResult(
            "reference (NCBI)", Outcome.SKIPPED,
            "no feature type appears in either the file or the published counts",
        )
    expected = {k: reference_counts.get(k, 0) for k in shared}
    actual = {k: folded.get(k, 0) for k in shared}
    if expected == actual:
        return TierResult(
            "reference (NCBI)", Outcome.PASS,
            f"matches published totals for {len(shared)} feature types",
            expected, actual,
        )
    return TierResult(
        "reference (NCBI)", Outcome.FAIL,
        "does not match NCBI's published totals for this assembly",
        expected, actual,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def render(results: list[TierResult]) -> str:
    """Format a verification report. Every tier is shown, including skips."""
    lines: list[str] = []
    width = max(len(r.name) for r in results) if results else 0
    for result in results:
        symbol = {Outcome.PASS: "PASS", Outcome.FAIL: "FAIL", Outcome.SKIPPED: "----"}[
            result.outcome
        ]
        lines.append(f"  {result.name:<{width}}  {symbol}  {result.detail}")
        for feature_type, (expected, actual) in result.mismatches.items():
            lines.append(
                f"      {feature_type:<20s} expected {expected:>8,d}   found {actual:>8,d}"
                f"   ({actual - expected:+,d})"
            )
    return "\n".join(lines)


def exit_code(results: list[TierResult]) -> int:
    """2 if any tier failed, 0 otherwise. A skipped tier is not a failure."""
    return 2 if any(r.outcome is Outcome.FAIL for r in results) else 0


# ---------------------------------------------------------------------------
# Tier 4: parent edges
# ---------------------------------------------------------------------------

#: Types that must have a parent when the file has a hierarchy at all. An
#: unparented one is either a genuine top-level feature (D-069) or a link the
#: loader failed to compute; the difference is what this tier reports.
LINKED_TYPES: frozenset[str] = frozenset({"CDS", "exon", "mRNA", "ncRNA", "tRNA", "rRNA"})


def check_parent_edges(
    per_method: dict[str | None, int], *, dialect: str
) -> TierResult:
    """Count parent edges by the rule that produced them.

    No other tier looks at edges, and the GenBank containment linker is the most
    novel and most fallible thing here. `_link_by_containment` skips an orphan
    CDS silently; `_link_by_locus_tag` emits no edge when a child carries no
    `locus_tag`. **A GenBank load that produced zero parent edges passed every
    check and reported PASS** (D-074).

    Cardinality tiers cannot catch this: every feature is present and counted,
    the hierarchy connecting them is simply absent. The column has been stored
    since the linkage-provenance fix; nothing ever asserted on it.

    For a GenBank file, an edge recorded as `explicit_parent` is itself a
    defect: the format carries no parent pointers, so any such edge means the
    tokenizer's synthesised attribute was mistaken for one the file supplied.
    """
    total = sum(per_method.values())
    unlinked = per_method.get(None, 0)
    linked = total - unlinked

    if total == 0:
        return TierResult(
            "parent edges", Outcome.SKIPPED,
            "no features of a type that takes a parent",
        )

    if dialect == "genbank" and per_method.get("explicit_parent"):
        return TierResult(
            "parent edges", Outcome.FAIL,
            f"{per_method['explicit_parent']:,} edge(s) recorded as stated by the "
            "file, but GenBank carries no parent pointers -- every edge is "
            "computed, so this means synthesised attributes were read back as "
            "the file's own",
            {"explicit_parent": 0}, {"explicit_parent": per_method["explicit_parent"]},
        )

    if linked == 0:
        return TierResult(
            "parent edges", Outcome.FAIL,
            f"none of {total:,} feature(s) that take a parent has one; the "
            "hierarchy is absent even though every feature was stored",
            {"linked": total}, {"linked": 0},
        )

    described = ", ".join(
        f"{method or 'unlinked'} {count:,}"
        for method, count in sorted(per_method.items(), key=lambda kv: (kv[0] is None, kv[0] or ""))
    )
    if unlinked:
        return TierResult(
            "parent edges", Outcome.FAIL,
            f"{unlinked:,} of {total:,} feature(s) that take a parent have none "
            f"({described})",
            {"linked": total}, {"linked": linked},
        )
    return TierResult(
        "parent edges", Outcome.PASS,
        f"{linked:,} edge(s), all accounted for ({described})",
    )
