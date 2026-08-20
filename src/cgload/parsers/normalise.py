"""Raw rows to a validated feature graph.

This is the one place that turns tokenizer output into something the loader can
insert. Everything it enforces is a decision written down elsewhere; the
docstrings name which, so a reader can check the code against the decision
rather than against my intent.

Nothing here touches the database. The loader (milestone 4) consumes the result.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from cgload.conformance import is_exempt
from cgload.db.vocabulary import (
    ALLOWED_PARENTS,
    HIERARCHY_TERMS,
    LINKAGE_METHODS,
    SYNONYMS,
    TRANSCRIPT_LAYER_REQUIRED,
    TRANSCRIPT_TERMS,
)
from cgload.parsers.profiles import LINKAGE_ATTRIBUTE, ParentInference, Profile
from cgload.parsers.tokenizer import RawFeature


class NormaliseError(ValueError):
    """A file that cannot be loaded without guessing. Always names the evidence."""


@dataclass
class Segment:
    """One row. A feature is one or more segments sharing an ID."""

    line_number: int
    start: int
    end: int
    phase: int | None
    score: float | None
    segment_index: int = 0
    #: The source row wrote start > end (D-056). Carried through so envelope
    #: repair can recompute the extent from children; never silently swapped.
    inverted: bool = False


@dataclass
class NormalisedFeature:
    """One feature, ready to insert.

    ``source_id`` is verbatim from the file (D-023) and may be None for an
    anonymous row (D-024). ``source_type`` is column 3 verbatim and is what
    ``verify`` counts; ``feature_type`` is the canonical term used for linkage
    and queries (D-017).
    """

    source_id: str | None
    source_type: str
    feature_type: str
    seqid: str
    strand: str
    source_program: str
    segments: list[Segment]
    parent_source_id: str | None = None
    linkage_method: str | None = None
    attributes: dict[str, list[str]] = field(default_factory=dict)

    @property
    def start(self) -> int:
        return min(segment.start for segment in self.segments)

    @property
    def end(self) -> int:
        return max(segment.end for segment in self.segments)

    @property
    def is_anonymous(self) -> bool:
        return self.source_id is None


@dataclass
class NormalisedFile:
    features: list[NormalisedFeature]
    #: Types present that the profile does not declare (open vocabulary, D-017).
    unexpected_types: set[str]
    #: Declared types absent from the file (finding 5).
    missing_types: set[str]
    #: Proteins recovered from comment blocks, keyed by gene ID (D-025).
    proteins: dict[str, str] = field(default_factory=dict)

    @property
    def feature_count(self) -> int:
        """The number the primary completeness check compares (D-008/D-024)."""
        return len(self.features)


def canonical_type(source_type: str) -> str:
    """Synonym mapping. Equal to the input for any unmapped term (D-017)."""
    return SYNONYMS.get(source_type, source_type)


def normalise(
    features: list[RawFeature],
    profile: Profile,
    *,
    id_prefix: str | None = None,
    comments: list[tuple[int, str]] | None = None,
    unexpected: set[str] | None = None,
    missing: set[str] | None = None,
) -> NormalisedFile:
    """Group segments, resolve parents, and validate the hierarchy."""
    grouped = _group_segments(features, profile, id_prefix=id_prefix)
    _order_segments(grouped)
    _resolve_parents(grouped, profile, id_prefix=id_prefix)
    _check_hierarchy(grouped)

    proteins: dict[str, str] = {}
    if profile.comment_extraction is not None and comments:
        extracted = profile.comment_extraction.extract(comments)
        # Comment blocks name the gene by its file ID, so a per-file id_prefix
        # must be applied here too. Without it the oracle looks up `g1` while the
        # features are stored as `AUG_g1`, every gene appears to have no CDS, and
        # the load fails on a correct file -- caught by running a real prefixed
        # load, not by any unit test.
        proteins = {
            str(_prefixed(gene, id_prefix)): sequence for gene, sequence in extracted.items()
        }

    return NormalisedFile(
        features=grouped,
        unexpected_types=unexpected or set(),
        missing_types=missing or set(),
        proteins=proteins,
    )


def _prefixed(identifier: str | None, id_prefix: str | None) -> str | None:
    """Apply the per-file prefix (D-014). The only transformation applied to an
    ID anywhere in cgload, and it is user-supplied rather than dialect-derived,
    which is why D-008's naive pass is permitted to apply it too."""
    if identifier is None:
        return None
    return f"{id_prefix}{identifier}" if id_prefix else identifier


def _group_segments(
    features: list[RawFeature], profile: Profile, *, id_prefix: str | None
) -> list[NormalisedFeature]:
    """One feature per distinct ID; one feature per anonymous row.

    Enforces the D-014 exemption: repeated IDs are permitted only when they form
    a discontinuous feature — same seqid, type, strand and parent, with
    non-overlapping spans — and only for types the profile says share IDs
    (finding 6). Everything else is an error naming both line numbers.
    """
    by_id: dict[str, NormalisedFeature] = {}
    ordered: list[NormalisedFeature] = []
    first_line: dict[str, int] = {}

    for raw in features:
        identifier = _prefixed(raw.id, id_prefix)
        segment = Segment(
            line_number=raw.line_number,
            start=raw.start,
            end=raw.end,
            phase=raw.phase,
            score=raw.score,
            inverted=getattr(raw, "inverted", False),
        )

        if identifier is None:
            if not profile.permits_anonymous_rows:
                raise NormaliseError(
                    f"line {raw.line_number}: {raw.type} row has no ID, but the "
                    f"{profile.name} profile does not permit anonymous rows. Either the "
                    f"dialect was detected wrongly or the file is truncated."
                )
            ordered.append(_new_feature(raw, identifier, segment))
            continue

        if identifier not in by_id:
            feature = _new_feature(raw, identifier, segment)
            by_id[identifier] = feature
            first_line[identifier] = raw.line_number
            ordered.append(feature)
            continue

        existing = by_id[identifier]
        if raw.type in profile.unique_id_types:
            raise NormaliseError(
                f"line {raw.line_number}: ID {identifier!r} repeats on a {raw.type} row "
                f"(first seen at line {first_line[identifier]}). In the "
                f"{profile.name} dialect every {raw.type} carries its own ID, so a "
                f"repeat is a duplicate rather than a discontinuous feature."
            )
        _assert_same_feature(existing, raw, identifier, first_line[identifier], profile)
        existing.segments.append(segment)

    return ordered


def _new_feature(
    raw: RawFeature, identifier: str | None, segment: Segment
) -> NormalisedFeature:
    parents = raw.parents
    if len(parents) > 1:
        raise NormaliseError(
            f"line {raw.line_number}: {raw.type} row has {len(parents)} parents "
            f"({', '.join(parents)}). cgload stores one parent per feature; a "
            f"multi-parent feature would have to be duplicated or silently truncated."
        )
    return NormalisedFeature(
        source_id=identifier,
        source_type=raw.type,
        feature_type=canonical_type(raw.type),
        seqid=raw.seqid,
        strand=raw.strand,
        source_program=raw.source,
        segments=[segment],
        parent_source_id=parents[0] if parents else None,
        # Set by a tokenizer that computed the edge itself (D-032). None for the
        # GFF3 dialects, where the file states the parent.
        linkage_method=raw.first(LINKAGE_ATTRIBUTE),
        attributes={
            key: list(values)
            for key, values in raw.attributes.items()
            if key != LINKAGE_ATTRIBUTE
        },
    )


def _has_trans_splicing(raw: RawFeature, key: str) -> bool:
    """Does this row carry a trans-splicing exception?

    Matched on the *value* rather than merely on the key, because
    `exception=` covers several unrelated NCBI annotations -- `ribosomal
    slippage`, `rearrangement required for product`, `annotated by transcript or
    proteomic data` -- and only trans-splicing licenses a strand disagreement.
    """
    return any("trans-splicing" in value for value in raw.all(key))


def _assert_same_feature(
    existing: NormalisedFeature,
    raw: RawFeature,
    identifier: str,
    first: int,
    profile: Profile,
) -> None:
    """The D-014 exemption, field by field.

    `phase` and `score` are deliberately excluded: phase varies across CDS
    segments by design, and score varies per segment in AUGUSTUS output.
    """
    # D-068. A trans-spliced gene is assembled from segments on opposite strands,
    # so no single strand describes it. NCBI flags those rows
    # `exception=trans-splicing` and writes `?` on the mRNA -- but not on the gene
    # rows, which is why the gene arrives as two rows differing only in strand.
    # Gated on the qualifier, so an unexplained strand disagreement is still two
    # features sharing an ID and is still refused.
    trans_spliced = (
        profile.strand_exception_qualifier is not None
        and _has_trans_splicing(raw, profile.strand_exception_qualifier)
    )
    if trans_spliced and existing.strand != raw.strand:
        # `?` is GFF3 for "stranded, but not a single strand" -- the honest value
        # for the merged feature, and the one NCBI itself writes on the mRNA.
        existing.strand = "?"

    checks = [
        ("seqid", existing.seqid, raw.seqid),
        ("type", existing.source_type, raw.type),
        ("parent", existing.parent_source_id, raw.parents[0] if raw.parents else None),
    ]
    if not trans_spliced:
        checks.insert(2, ("strand", existing.strand, raw.strand))

    for label, mine, theirs in checks:
        if mine != theirs:
            raise NormaliseError(
                f"line {raw.line_number}: ID {identifier!r} reappears with {label} "
                f"{theirs!r}, but line {first} had {mine!r}. Two different features "
                f"share an ID; cgload will not merge them."
            )

    # D-033. Overlap is a duplicate in GFF3 and a valid record in GenBank: a
    # programmed frameshift reads one base twice, so
    # join(57926..58406,58406..59148) is correct and must load. Gated on the
    # qualifier actually being present, so an unexplained overlap is still an
    # error in every format.
    if profile.overlap_qualifier and raw.first(profile.overlap_qualifier) is not None:
        return

    for segment in existing.segments:
        if raw.start <= segment.end and segment.start <= raw.end:
            raise NormaliseError(
                f"line {raw.line_number}: ID {identifier!r} has segment "
                f"{raw.start}-{raw.end} overlapping {segment.start}-{segment.end} from "
                f"line {segment.line_number}. Segments of one discontinuous feature "
                f"cannot overlap."
            )


def _order_segments(features: list[NormalisedFeature]) -> None:
    """Assign ``segment_index`` in ascending coordinate order.

    Not file order: funannotate emits minus-strand segments 3'->5' (finding 7),
    so file order would make transcript order dialect-dependent. Strand is
    stored separately, so transcript order is coordinate order on `+` and its
    reverse on `-`.
    """
    for feature in features:
        feature.segments.sort(key=lambda segment: (segment.start, segment.end))
        for index, segment in enumerate(feature.segments):
            segment.segment_index = index


def _resolve_parents(
    features: list[NormalisedFeature],
    profile: Profile,
    *,
    id_prefix: str | None,
) -> None:
    """Attach each child to its parent, and record how (D-007's vocabulary).

    Milestone 3 handles explicit `Parent` only. Structural and synthesised
    inference arrive at milestones 6 and 5; a profile declaring either here is a
    programming error rather than a silent fallback.
    """
    if profile.parent_inference is ParentInference.SYNTHESISED:
        raise NormaliseError(
            f"the {profile.name} profile requests synthesised parent inference, which "
            f"is not implemented. Refusing rather than falling back to explicit "
            f"Parent, which would silently produce a wrong gene graph."
        )

    known = {feature.source_id for feature in features if feature.source_id is not None}
    for feature in features:
        if feature.parent_source_id is None:
            continue
        parent = _prefixed(feature.parent_source_id, id_prefix)
        feature.parent_source_id = parent
        if parent not in known:
            line = feature.segments[0].line_number
            raise NormaliseError(
                f"line {line}: {feature.source_type} references Parent {parent!r}, "
                f"which is not defined anywhere in this file. A dangling parent means "
                f"the file is incomplete or was sliced badly; loading it would store an "
                f"orphan that no query would ever find."
            )

        # D-032. Where the format states its own parents, the edge came from the
        # file and `explicit_parent` is the truth. Where it does not -- GenBank
        # has no parent pointers at all -- the tokenizer computed the edge and
        # stamped which rule it used. Overwriting that with `explicit_parent`
        # would be a false claim about provenance, and it is the exact claim the
        # containment algorithm is published on, so it must survive to the
        # database and to `explain`.
        if profile.states_own_parents:
            feature.linkage_method = "explicit_parent"
            continue

        if feature.linkage_method is None:
            raise NormaliseError(
                f"line {feature.segments[0].line_number}: the {profile.name} profile "
                f"does not state its own parents, so every computed edge must carry the "
                f"rule that produced it. {feature.source_id!r} has a parent and no "
                f"linkage method, which means the tokenizer failed to record how it was "
                f"resolved."
            )
        if feature.linkage_method not in LINKAGE_METHODS:
            raise NormaliseError(
                f"line {feature.segments[0].line_number}: unknown linkage method "
                f"{feature.linkage_method!r}. Permitted: {', '.join(LINKAGE_METHODS)}."
            )


def _check_hierarchy(features: list[NormalisedFeature]) -> None:
    """A Tier A feature's parent must also be Tier A, and a permitted one.

    Linkage rules are written per Tier A type and assume Tier A parents (D-017),
    so an `exon` hanging off an inert feature has no defined meaning. Inert
    features may hang off anything, or nothing.
    """
    types = {
        feature.source_id: feature.feature_type
        for feature in features
        if feature.source_id is not None
    }
    #: Whether this file has a transcript layer at all. A prokaryotic GenBank
    #: record has none, and its CDS rows attach straight to their genes.
    #:
    #: **`mRNA` alone, deliberately, and not TRANSCRIPT_TERMS.** This reads like
    #: an oversight next to `genes_with_transcripts` below, which uses all four
    #: transcript terms, and it is not. Every bacterial genome carries tRNA and
    #: rRNA genes, so widening this to TRANSCRIPT_TERMS makes `has_transcripts`
    #: true for every prokaryotic record, and then every `gene`->`CDS` edge --
    #: the normal prokaryotic shape that D-036 exists to permit -- is refused.
    #: Verified by widening it: the whole prokaryotic fixture class fails.
    #:
    #: So this line is not really asking "does the file have transcripts". It is
    #: a proxy asking "is this a eukaryotic annotation", and `mRNA` is the marker
    #: that answers it. The known cost of the proxy is that a eukaryotic file
    #: whose transcript layer is entirely non-coding -- an ncRNA-only or
    #: organellar slice -- reads as prokaryotic, and a genuinely missing
    #: transcript there loads unchecked. That gap is accepted, not fixed, because
    #: no predicate available at this point distinguishes the two cases; it is
    #: recorded in D-059 and pinned by tests in both directions.
    has_transcripts = any(feature.feature_type == "mRNA" for feature in features)

    #: Genes that do have a transcript child. The transcript-bypass check below
    #: is scoped to genes *without* one, because real RefSeq attaches extra
    #: gene-level `exon` rows (`gbkey=exon`, `ID=id-...`) alongside the
    #: transcript-parented ones -- see D-043. Those are redundant output, not a
    #: missing transcript, and a file-wide rule refused a correct RefSeq genome.
    genes_with_transcripts = {
        feature.parent_source_id
        for feature in features
        if feature.feature_type in TRANSCRIPT_TERMS and feature.parent_source_id
    }

    for feature in features:
        if feature.feature_type not in HIERARCHY_TERMS:
            continue
        permitted = ALLOWED_PARENTS[feature.feature_type]
        parent = feature.parent_source_id

        if parent is None:
            # A feature that declares no parent is top-level, and GFF3 permits
            # that for any type -- the specification nowhere requires an exon to
            # be parented. Refusing here conflated two different things: a
            # feature naming a parent it may not have, and a feature naming none
            # at all. Only the first hides something.
            #
            # Nothing is concealed by a parentless row. D-058's case had a
            # parent present with a layer missing between; here the file never
            # claimed a parent, so storing it top-level records exactly what the
            # file says.
            #
            # Found on *Zea mays* B73 GCF_902167145.1 at line 992,090: two
            # chloroplast exons, `ID=id-NC_001666.2:129636..129867-2`, flanking a
            # parentless intron, all three carrying `number=2`. They are the
            # second exon block of a trans-spliced chloroplast gene whose first
            # exon lies elsewhere; NCBI's converter could not express the link
            # and emitted them with coordinate-derived identifiers. Two rows in
            # 1.2 M features, both organellar -- the fourth manifestation of the
            # same converter defect (D-069).
            continue

        parent_type = types[parent]

        # A gene whose transcript is missing, with its parts hanging off the
        # gene instead, hides an incomplete hierarchy -- that is what this rule
        # guards against.
        #
        # **A pseudogene is not that case.** RefSeq attaches exons directly to a
        # `pseudogene` row by design: a pseudogene has no functional transcript
        # to model, so there is no missing layer to hide. Measured on *Zea mays*
        # B73 GCF_902167145.1: 5,195 pseudogenes and 12,757 exons attached to a
        # gene-level parent. Refusing them rejected standard NCBI output at line
        # 247 of chromosome 1 (D-058).
        #
        # This is the second rule in cgload stated about the whole file when it
        # should have been stated about each parent -- the first was the
        # allowlist of types permitted to share an ID across segments, narrowed
        # in milestone 8c. Both were written from what four fixtures happened to
        # contain.
        if (
            feature.feature_type in TRANSCRIPT_LAYER_REQUIRED
            and parent_type == "gene"
            and has_transcripts
            and parent not in genes_with_transcripts
        ):
            line = feature.segments[0].line_number
            raise NormaliseError(
                f"line {line}: {feature.source_type} {feature.source_id!r} attaches "
                f"directly to gene {parent!r}, which has no transcript of its "
                f"own, while the rest of this file does have a transcript layer. "
                f"Attaching past a missing transcript would make an incomplete "
                f"hierarchy look valid. (A pseudogene parent is exempt: it has no "
                f"functional transcript to model.)"
            )

        if parent_type not in permitted:
            line = feature.segments[0].line_number
            raise NormaliseError(
                f"line {line}: {feature.source_type} {feature.source_id!r} has parent "
                f"{parent!r} of type {parent_type!r}, but a {feature.feature_type} may "
                f"only attach to {sorted(permitted) or 'nothing'}."
            )


def counts_by_source_type(normalised: NormalisedFile) -> dict[str, int]:
    """Feature counts grouped by the verbatim column-3 term.

    This is the database-independent half of what ``verify`` will compare. It
    counts *features*, so a four-segment CDS is one, and each anonymous row is
    one (D-024).
    """
    counts: dict[str, int] = defaultdict(int)
    for feature in normalised.features:
        counts[feature.source_type] += 1
    return dict(counts)


def check_declared_proteins(normalised: NormalisedFile) -> list[str]:
    """The AUGUSTUS oracle (D-025). Returns a list of mismatch descriptions.

    ``sum(CDS spans) / 3 - 1`` must equal the declared amino-acid count: the
    minus one is the stop codon, which this dialect places inside the terminal
    CDS (D-022). If a segment were dropped, a span were off by one, or a segment
    were grouped under the wrong transcript, the arithmetic breaks.

    This shares no code with the coordinate handling it checks -- it reads the
    tool's own declared output -- which is the kind of independent check D-008
    asks for and cannot usually get for free.
    """
    if not normalised.proteins:
        return []

    exempt = {
        feature.parent_source_id
        for feature in normalised.features
        if feature.feature_type == "CDS" and is_exempt(feature.attributes)
    }

    coding: dict[str, int] = defaultdict(int)
    transcript_gene: dict[str, str | None] = {}
    for feature in normalised.features:
        if feature.feature_type == "mRNA":
            transcript_gene[str(feature.source_id)] = feature.parent_source_id
        if feature.feature_type == "CDS" and feature.parent_source_id:
            span = sum(segment.end - segment.start + 1 for segment in feature.segments)
            coding[feature.parent_source_id] += span

    by_gene: dict[str, int] = defaultdict(int)
    for transcript, bases in coding.items():
        gene = transcript_gene.get(transcript)
        if gene is not None:
            by_gene[gene] += bases

    exempt_genes = {
        transcript_gene.get(transcript) for transcript in exempt if transcript
    }

    problems: list[str] = []
    for gene, protein in sorted(normalised.proteins.items()):
        if gene in exempt_genes:
            # Pseudogene, partial, artificial location, slippage: the arithmetic
            # is not expected to hold. One shared predicate, so milestone 8's
            # completeness check cannot disagree with this one about which
            # features are exempt.
            continue
        bases = by_gene.get(gene)
        if bases is None:
            problems.append(f"{gene}: a protein is declared but no CDS was parsed")
            continue
        if bases % 3:
            problems.append(f"{gene}: CDS totals {bases} bp, not a whole number of codons")
            continue
        expected = bases // 3 - 1
        if expected != len(protein):
            problems.append(
                f"{gene}: CDS totals {bases} bp implying {expected} aa, but AUGUSTUS "
                f"declares {len(protein)} aa"
            )
    return problems
