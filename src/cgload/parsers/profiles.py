"""Layer 2: dialect profiles, expressed as data (D-013).

A profile declares seven things and contains no parsing logic. Adding a dialect
is a profile plus a fixture, never a new parser — which is what keeps D-008's
verification and its mutation test applicable to one engine rather than five.

**Profiles fail closed.** Every field is explicitly declared; there are no
defaults. A dialect that does not state its stop-codon convention is a
programming error caught at import, not a silent assumption at load time. Two
semantic gaps in this model surfaced during review (coordinate reconciliation,
provided-types) and two more from reading real files, so the enumeration is
demonstrably not finished; failing closed converts the next omission from silent
wrongness into a loud failure.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum

from cgload.parsers.tokenizer import RawFeature, TokenizedFile

#: Attribute the GenBank tokenizer uses to hand the normaliser the rule it used
#: to compute each parent edge. Stripped before attributes are stored: the value
#: lives in ``feature.linkage_method``, and a second copy could disagree.
LINKAGE_ATTRIBUTE = "cgload_linkage_method"


class StopCodonConvention(Enum):
    """Where the stop codon lives relative to the CDS."""

    #: Inside the terminal CDS. GFF3 dialects. Codon rows, where present, are
    #: inert features and must not be mapped to CDS (D-022).
    INSIDE_CDS = "inside_cds"
    #: Outside the CDS, as a separate row to be absorbed into the terminal CDS
    #: span. Ensembl-style GTF. Not yet exercised; arrives at milestone 5.
    OUTSIDE_CDS = "outside_cds"


class ParentInference(Enum):
    """How a child finds its parent."""

    #: An explicit `Parent` attribute. All three GFF3 dialects.
    EXPLICIT = "explicit"
    #: gene_id/transcript_id, with the gene and transcript rows synthesised
    #: because GTF often omits them. Milestone 5.
    SYNTHESISED = "synthesised"
    #: Exon-interval containment. GenBank, milestone 6 (D-007).
    STRUCTURAL = "structural"


class LinkageStrategy(Enum):
    """How a GenBank child finds its parent (D-029).

    GenBank carries no ``Parent``, so the edge is computed. Which computation is
    correct depends on the annotation, and getting it wrong is silent, so it is
    declared rather than inferred at load time.
    """

    #: ``locus_tag`` is an explicit unique key shared by the gene row and its
    #: single child. Prokaryotic annotation. No ambiguity is possible.
    LOCUS_TAG = "locus_tag"
    #: CDS-to-transcript by exon-interval containment, with the declared
    #: ``from transcript`` note as a tiebreak, refusing when both fail.
    #: Eukaryotic annotation.
    CONTAINMENT = "containment"
    #: Choose by whether the file contains transcript records at all.
    AUTO = "auto"
    #: Not a GenBank profile; the field is unused.
    NOT_APPLICABLE = "n/a"


class Confidence(Enum):
    HIGH = "high"
    LOW = "low"
    OVERRIDDEN = "overridden"


@dataclass(frozen=True)
class CommentExtraction:
    """Dialect-specific comment handling (D-025). AUGUSTUS only, so far."""

    name: str
    #: Reads the comment list, returns identifier -> extracted payload.
    extract: Callable[[list[tuple[int, str]]], dict[str, str]]


@dataclass(frozen=True)
class Profile:
    """A dialect. Data, not code."""

    name: str
    #: Column 2 values that identify this dialect.
    sources: frozenset[str]
    #: Optional header corroboration, matched against comment lines.
    header_signature: re.Pattern[str] | None
    parent_inference: ParentInference
    #: Only meaningful when ``parent_inference`` is STRUCTURAL (D-029).
    linkage_strategy: LinkageStrategy
    stop_codon: StopCodonConvention
    #: Feature types this dialect is expected to emit. Absence of a declared
    #: type is a warning; absence of an undeclared type is expected (finding 5).
    provides: frozenset[str]
    #: Types whose absence earns a warning. A subset of ``provides``: some
    #: declared types are optional in real files, and warning about those trains
    #: users to ignore warnings, which costs more than the warning gains.
    warns_if_missing: frozenset[str]
    #: Types where a repeated ID is *always* an error, never a discontinuous
    #: feature. Stated as a deny-list rather than an allow-list (D-049): real
    #: RefSeq gives `sequence_feature` a join() location and shares one ID across
    #: its segments, exactly as CDS does, so an allow-list of {"CDS"} refused a
    #: valid record. What is actually observed is narrower and more useful --
    #: exons always get their own IDs in every dialect -- so that is what is
    #: declared, and everything else is judged by the D-014 conditions alone.
    unique_id_types: frozenset[str]
    #: Whether rows may omit `ID` entirely (D-024).
    permits_anonymous_rows: bool
    comment_extraction: CommentExtraction | None
    #: A qualifier whose presence permits the segments of one feature to
    #: overlap. GenBank's /ribosomal_slippage is the only known case: a
    #: programmed frameshift deliberately reads one base twice, so
    #: ``join(57926..58406,58406..59148)`` is a valid record, not a duplicate.
    #: None means overlap is always an error, which is correct for GFF3
    #: (D-033).
    overlap_qualifier: str | None = None
    #: A qualifier whose presence permits the segments of one feature to sit on
    #: *different strands*. NCBI's `exception=trans-splicing` is the only known
    #: case: a trans-spliced gene is assembled from segments on opposite strands,
    #: so no single strand is correct for the feature as a whole and GFF3's answer
    #: is `?` (D-068). None means differing strands are always an error, which is
    #: right for every other dialect.
    strand_exception_qualifier: str | None = None
    #: True when the format states parents itself. False when cgload computes
    #: them, in which case the tokenizer must stamp the method per feature and
    #: the normaliser must not overwrite it (D-032).
    states_own_parents: bool = True

    def matches(self, tokenized: TokenizedFile) -> tuple[bool, Confidence, str]:
        """Does this profile describe the file? Always returns its evidence."""
        overlap = tokenized.sources & self.sources
        if not overlap:
            return False, Confidence.LOW, f"column 2 is {sorted(tokenized.sources)}"

        evidence = f"column 2 is {sorted(overlap)}"
        if self.header_signature is not None:
            for _, comment in tokenized.comments:
                if self.header_signature.search(comment):
                    return True, Confidence.HIGH, f"{evidence}; header matched"
            return True, Confidence.LOW, f"{evidence}; no header corroboration"
        return True, Confidence.HIGH, evidence


# --------------------------------------------------------------------------
# AUGUSTUS comment extraction (D-025)
# --------------------------------------------------------------------------

PROTEIN_OPEN = re.compile(r"^#\s*protein sequence = \[(.*)$")
GENE_MARKER = re.compile(r"^#\s*start gene (\S+)\s*$")


def extract_augustus_proteins(comments: list[tuple[int, str]]) -> dict[str, str]:
    """Recover the declared translations.

    The only place AUGUSTUS emits protein sequence: a block opening with
    ``# protein sequence = [`` and continuing over ``#`` lines to a closing
    ``]``. Keyed by the gene named in the preceding ``# start gene`` marker.
    """
    proteins: dict[str, str] = {}
    gene: str | None = None
    buffer: list[str] | None = None

    for _, line in comments:
        marker = GENE_MARKER.match(line)
        if marker:
            gene = marker.group(1)
            continue

        opening = PROTEIN_OPEN.match(line)
        if opening:
            buffer = [opening.group(1).strip()]
        elif buffer is not None:
            buffer.append(line.lstrip("#").strip())
        else:
            continue

        if buffer[-1].endswith("]"):
            buffer[-1] = buffer[-1][:-1]
            if gene is not None:
                proteins[gene] = "".join(buffer)
            buffer = None

    return proteins


AUGUSTUS_PROTEINS = CommentExtraction(
    name="augustus_protein_blocks", extract=extract_augustus_proteins
)


# --------------------------------------------------------------------------
# The registry. Milestone 3 registers RefSeq only; the other two are declared
# here so detection can *recognise* them and say so, while loading them stays
# an error until their fixtures pass (D-013's fail-closed rule).
# --------------------------------------------------------------------------

REFSEQ = Profile(
    name="refseq",
    sources=frozenset({"RefSeq", "Gnomon", "BestRefSeq", "tRNAscan-SE", "cmsearch"}),
    header_signature=re.compile(r"#!processor NCBI annotwriter|#!genome-build-accession"),
    parent_inference=ParentInference.EXPLICIT,
    linkage_strategy=LinkageStrategy.NOT_APPLICABLE,
    stop_codon=StopCodonConvention.INSIDE_CDS,
    provides=frozenset(
        {"region", "gene", "pseudogene", "mRNA", "tRNA", "rRNA", "ncRNA", "exon", "CDS"}
    ),
    warns_if_missing=frozenset({"gene", "mRNA", "exon", "CDS"}),
    unique_id_types=frozenset({"exon"}),
    permits_anonymous_rows=False,
    comment_extraction=None,
    overlap_qualifier=None,
    # RefSeq writes trans-spliced organellar genes as several rows sharing one ID
    # with differing strands, flagged `exception=trans-splicing` (D-068).
    strand_exception_qualifier="exception",
    states_own_parents=True,
)

FUNANNOTATE = Profile(
    name="funannotate",
    sources=frozenset({"funannotate"}),
    header_signature=None,
    parent_inference=ParentInference.EXPLICIT,
    linkage_strategy=LinkageStrategy.NOT_APPLICABLE,
    stop_codon=StopCodonConvention.INSIDE_CDS,
    provides=frozenset({"gene", "mRNA", "tRNA", "exon", "CDS"}),
    warns_if_missing=frozenset({"gene", "mRNA", "exon", "CDS"}),
    unique_id_types=frozenset({"exon"}),
    permits_anonymous_rows=False,
    comment_extraction=None,
    overlap_qualifier=None,
    states_own_parents=True,
)

AUGUSTUS = Profile(
    name="augustus",
    sources=frozenset({"AUGUSTUS"}),
    header_signature=re.compile(r"generated with AUGUSTUS"),
    parent_inference=ParentInference.EXPLICIT,
    linkage_strategy=LinkageStrategy.NOT_APPLICABLE,
    stop_codon=StopCodonConvention.INSIDE_CDS,
    # No `exon`: the structure is implicit in the CDS rows (finding 5). Zero
    # exons from an AUGUSTUS load is therefore expected, not a silent drop.
    provides=frozenset(
        {"gene", "transcript", "CDS", "intron", "start_codon", "stop_codon"}
    ),
    warns_if_missing=frozenset({"gene", "transcript", "CDS"}),
    unique_id_types=frozenset({"exon"}),
    permits_anonymous_rows=True,
    comment_extraction=AUGUSTUS_PROTEINS,
    overlap_qualifier=None,
    states_own_parents=True,
)

GENBANK = Profile(
    name="genbank",
    # Synthesised by the GenBank tokenizer: a flat file has no column 2.
    sources=frozenset({"GenBank"}),
    header_signature=re.compile(r"cgload linkage strategy"),
    parent_inference=ParentInference.STRUCTURAL,
    linkage_strategy=LinkageStrategy.AUTO,
    stop_codon=StopCodonConvention.INSIDE_CDS,
    provides=frozenset(
        {
            "gene", "pseudogene", "mRNA", "tRNA", "rRNA", "ncRNA", "tmRNA",
            "CDS", "regulatory", "repeat_region", "mobile_element", "misc_RNA",
        }
    ),
    # Empty, deliberately. One profile serves both kingdoms, so `provides` is a
    # union across them: a bacterial record legitimately has no mRNA and a
    # eukaryotic slice legitimately has no tRNA. Warning on either would fire on
    # every correct file.
    warns_if_missing=frozenset(),
    # A GenBank record IS the feature: join(1..96,200..300) is one CDS in one
    # record, emitted as several segments sharing the identifier the file gave it,
    # and transcripts are multi-exon the same way. Nothing here has a repeated ID
    # that is not a discontinuous feature, so the deny-list is empty.
    unique_id_types=frozenset(),
    # source, regulatory and repeat_region records carry no usable identifier.
    permits_anonymous_rows=True,
    comment_extraction=None,
    # A programmed frameshift reads one base twice, so the segments of one CDS
    # legitimately overlap. Gated on the qualifier: without it, an overlap in a
    # GenBank file is still an error (D-033).
    overlap_qualifier="ribosomal_slippage",
    # GenBank writes the same genes as one record with a join() across strands,
    # handled by the location parser rather than by segment merging.
    strand_exception_qualifier=None,
    # The format carries no parent pointers. Every edge is computed, so the
    # tokenizer stamps the method it used and the normaliser honours it (D-032).
    states_own_parents=False,
)

#: Header pragma a cgload export writes to identify itself. Duplicated from
#: ``cgload.export`` as a literal rather than imported, because ``profiles`` must
#: not depend on the exporter: the profile table is a set of declarations, and a
#: cycle through the writer would make loading depend on writing.
EXPORT_PRAGMA = "##cgload-export"

CGLOAD = Profile(
    name="cgload",
    # Column 2 of an export is the original producing program, preserved per
    # feature (D-052), so this profile claims no source value at all. It is
    # selected by the header pragma alone -- which is why `detect` checks the
    # pragma before consulting column 2.
    sources=frozenset({"cgload"}),
    header_signature=re.compile(re.escape(EXPORT_PRAGMA)),
    parent_inference=ParentInference.EXPLICIT,
    linkage_strategy=LinkageStrategy.NOT_APPLICABLE,
    stop_codon=StopCodonConvention.INSIDE_CDS,
    # An export may carry any type the source dialect held, and the vocabulary is
    # open (D-017), so nothing is declared and nothing is expected.
    provides=frozenset(
        {
            "region", "gene", "pseudogene", "mRNA", "transcript", "tRNA", "rRNA",
            "ncRNA", "tmRNA", "exon", "CDS", "intron", "start_codon", "stop_codon",
            "five_prime_UTR", "three_prime_UTR", "repeat_region", "sequence_feature",
        }
    ),
    warns_if_missing=frozenset(),
    # Whatever the source shared, the export writes back the same way. An export
    # of AUGUSTUS shares CDS ids; an export of RefSeq shares them for
    # sequence_feature too (D-049), so the deny-list stays as narrow as the
    # dialects it may have come from.
    unique_id_types=frozenset({"exon"}),
    # Preserved from the source: an AUGUSTUS export still has 175 ID-less rows,
    # because inventing identifiers on export would breach D-023.
    permits_anonymous_rows=True,
    comment_extraction=None,
    overlap_qualifier="ribosomal_slippage",
    # False, even though every exported feature does carry a `Parent`. The export
    # writes the *recorded* linkage method per feature, and that stamp is the
    # truth: an export of a GenBank load has Parent attributes whose edges were
    # computed by containment, and treating the format as self-stating would
    # relabel all of them `explicit_parent` on re-import (D-053). That is D-032's
    # false-provenance failure arriving by the export path.
    states_own_parents=False,
)

#: Profiles cgload will *load*. Each was added with a real fixture and its own
#: end-to-end test, per milestone 5: a dialect is a profile plus a fixture, never
#: a parser.
REGISTERED: tuple[Profile, ...] = (REFSEQ, FUNANNOTATE, AUGUSTUS, GENBANK, CGLOAD)

#: Profiles cgload can *name* but not yet load. Detection reports these as a
#: recognised-but-unregistered dialect, which is a clearer error than "no
#: confident match" and still refuses to guess. Empty is the correct value here:
#: the mechanism stays, because it is what lets the next dialect be recognised
#: before it is supported, and D-027's three-outcome rule is still tested using a
#: synthetic profile rather than by keeping a real one disabled.
KNOWN_UNREGISTERED: tuple[Profile, ...] = ()


class DetectionError(RuntimeError):
    """No confident match, an ambiguous match, or a recognised-but-unregistered
    dialect. Never falls back to a permissive generic profile — that fallback is
    the most tempting and most damaging thing that could be added here (D-013).
    """


@dataclass(frozen=True)
class Detection:
    profile: Profile
    confidence: Confidence
    evidence: str


def detect(tokenized: TokenizedFile, *, override: str | None = None) -> Detection:
    """Identify the dialect, or refuse. Always reports the evidence."""
    if override is not None:
        for profile in REGISTERED:
            if profile.name == override:
                return Detection(profile, Confidence.OVERRIDDEN, f"--dialect {override}")
        names = ", ".join(p.name for p in REGISTERED)
        raise DetectionError(f"unknown dialect {override!r}; registered profiles: {names}")

    # A cgload export declares itself in a header pragma, and that declaration
    # outranks column 2 (D-052). Column 2 preserves the *original* producing
    # program -- an export of an AUGUSTUS load still says AUGUSTUS there, and
    # correctly so -- but the file no longer follows AUGUSTUS's conventions, so
    # letting the augustus profile claim it would apply the wrong rules to a file
    # that says plainly what it is.
    for _, comment in tokenized.comments:
        if comment.startswith(EXPORT_PRAGMA):
            return Detection(CGLOAD, Confidence.HIGH, f"header declares {EXPORT_PRAGMA}")

    matches = []
    for profile in REGISTERED:
        matched, confidence, evidence = profile.matches(tokenized)
        if matched:
            matches.append(Detection(profile, confidence, evidence))

    if len(matches) == 1:
        return matches[0]

    if len(matches) > 1:
        names = ", ".join(m.profile.name for m in matches)
        raise DetectionError(
            f"ambiguous: {names} all match (column 2 is "
            f"{sorted(tokenized.sources)}). Pass --dialect to choose."
        )

    for profile in KNOWN_UNREGISTERED:
        matched, _, evidence = profile.matches(tokenized)
        if matched:
            raise DetectionError(
                f"recognised as {profile.name} ({evidence}), which is declared but not "
                f"yet enabled for loading. cgload will not fall back to a generic "
                f"reader, because a misread annotation produces a plausible and wrong "
                f"gene graph with no error."
            )

    registered = ", ".join(p.name for p in REGISTERED)
    raise DetectionError(
        f"no confident match: column 2 is {sorted(tokenized.sources)}, and no "
        f"registered profile claims it (registered: {registered}). Pass --dialect "
        f"to force one."
    )


def unexpected_types(profile: Profile, features: list[RawFeature]) -> set[str]:
    """Types present in the file that the profile does not declare.

    Not an error: the vocabulary is open (D-017), so an undeclared type loads
    verbatim. Surfaced in the load summary so it is visible rather than silent.
    """
    return {feature.type for feature in features} - profile.provides


def missing_types(profile: Profile, features: list[RawFeature]) -> set[str]:
    """Types the profile declares but the file does not contain.

    This is what makes "the file had no exons" distinguishable from "we dropped
    the exons" (finding 5). AUGUSTUS declares no `exon`, so its absence is
    expected; a RefSeq file with no `exon` rows earns a warning.
    """
    present = {feature.type for feature in features}
    return {declared for declared in profile.warns_if_missing if declared not in present}
