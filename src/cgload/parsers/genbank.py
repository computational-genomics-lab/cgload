"""Layer 1 for GenBank flat files: tokenizer plus structural linkage.

GenBank is not GFF3 with different punctuation. Three differences drive
everything in this module:

1. **There is no ``Parent`` attribute.** Hierarchy is implicit. Every parent
   edge here is *computed*, and the computation is the risky part, so it is
   done explicitly, declared per profile, and refuses rather than guesses.
2. **There is no ``ID`` attribute.** Identifiers are recovered from the
   qualifiers the file already carries -- ``locus_tag``, ``transcript_id``,
   ``protein_id``. Nothing is invented; an identifier cgload made up would
   appear in no input file and would break D-008's independent count.
3. **Locations are expressions**, not two integers: ``complement(join(1..96,
   200..300))``. One feature record can therefore be several segments, which
   maps onto the existing ``Segment`` grouping rather than needing new machinery.

The output is a ``TokenizedFile`` of ordinary ``RawFeature`` objects carrying
synthesised ``ID`` and ``Parent`` attributes, so ``normalise()`` and everything
downstream is shared with the GFF3 dialects unchanged.

**No biopython.** D-009 reserved a biopython dependency for this module; it is
not needed. The flat-file grammar exercised here is small and fully covered by
the fixtures, and dropping the dependency keeps the runtime budget at three
packages (P4) and removes a compiled-extension install risk on another
machine. Recorded as D-028.
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

from cgload.parsers.tokenizer import RawFeature, TokenizedFile, TokenizeError

#: See ``cgload.parsers.profiles.LINKAGE_ATTRIBUTE``. Defined here as a literal
#: rather than imported, because ``profiles`` must not depend on this module: the
#: GenBank profile is a declaration, and a tokenizer import would make the
#: profile table depend on one dialect's parser.
LINKAGE_ATTRIBUTE = "cgload_linkage_method"

#: Feature-table records begin at column 5; qualifiers at column 21.
_FEATURE_LINE = re.compile(r"^ {5}(\S+)\s+(\S.*)$")
_QUALIFIER = re.compile(r'^ {21}/(\w+)(?:=(.*))?$')
_LOCUS = re.compile(r"^LOCUS\s+(\S+)\s+(\d+)\s+bp.*?\b(circular|linear)\b", re.I)
_VERSION = re.compile(r"^VERSION\s+(\S+)")
_SPAN = re.compile(r"^<?(\d+)(?:\.\.>?(\d+))?$")

#: ``/note="... from transcript CG12164-RA; ..."``. An exact identifier in a
#: documented position -- not a similarity match. The distinction matters:
#: pairing CDS to mRNA by fuzzy product-name comparison mispairs some and
#: silently drops whatever fails to match.
_FROM_TRANSCRIPT = re.compile(r"from transcript (\S+?)\s*[;\"]")


class GenBankError(TokenizeError):
    """A GenBank record that cannot be read without guessing."""


@dataclass
class _Record:
    """One feature-table record, before any interpretation."""

    line_number: int
    type: str
    location: str
    qualifiers: dict[str, list[str]] = field(default_factory=dict)

    def first(self, key: str) -> str | None:
        values = self.qualifiers.get(key)
        return values[0] if values else None

    @property
    def is_pseudo(self) -> bool:
        return "pseudo" in self.qualifiers or "pseudogene" in self.qualifiers


def _open(path: Path):
    """Open plain or gzipped, detected by magic bytes rather than by suffix.

    NCBI serves ``.gbff.gz``; users rename them. Trusting the extension makes a
    correct file fail with a confusing error.
    """
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_location(
    expression: str, path: Path, line_number: int
) -> tuple[list[tuple[int, int]], str]:
    """Return ``([(start, end), ...], strand)`` from a GenBank location.

    Handles ``complement``, ``join``, ``order`` and the ``<``/``>`` partial
    markers. Partial markers are discarded here and the fact is preserved as an
    attribute by the caller, because a truncated coordinate is still a real
    coordinate -- it is the *completeness* claim that changes, not the position.
    """
    text = expression.strip().replace(" ", "")
    strand = "+"
    if text.startswith("complement("):
        strand = "-"
        text = text[len("complement(") : -1]
    # A nested complement (mixed-strand join) cannot be represented by one
    # strand value, so it is refused rather than flattened.
    if "complement(" in text:
        raise GenBankError(
            path, line_number,
            f"mixed-strand location {expression!r}: cgload stores one strand per "
            "feature, so this cannot be represented without discarding information",
        )
    for wrapper in ("join(", "order("):
        if text.startswith(wrapper):
            text = text[len(wrapper) : -1]
            break

    spans: list[tuple[int, int]] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        match = _SPAN.match(part)
        if not match:
            raise GenBankError(path, line_number, f"cannot parse location {part!r}")
        start = int(match.group(1))
        end = int(match.group(2)) if match.group(2) else start
        spans.append((start, end))
    if not spans:
        raise GenBankError(path, line_number, f"empty location {expression!r}")
    return spans, strand


def _read_records(
    path: Path,
) -> tuple[list[_Record], dict[str, tuple[int, int]], list[tuple[int, str]], str]:
    """Split the file into feature records. One pass, no interpretation."""
    records: list[_Record] = []
    regions: dict[str, tuple[int, int]] = {}
    comments: list[tuple[int, str]] = []
    accession = ""
    in_features = False
    current: _Record | None = None
    qualifier: str | None = None
    locus_name = ""

    with _open(path) as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")

            if line.startswith("LOCUS"):
                match = _LOCUS.match(line)
                if match:
                    locus_name = match.group(1)
                    regions[locus_name] = (1, int(match.group(2)))
                    comments.append((number, f"# topology {match.group(3).lower()}"))
                in_features = False
                continue
            if line.startswith("VERSION"):
                match = _VERSION.match(line)
                if match:
                    accession = match.group(1)
                    # The VERSION accession is the authoritative seqid; LOCUS
                    # drops the version suffix, and the GFF3 for the same
                    # assembly uses the versioned form. Using LOCUS here would
                    # make the two files disagree about the sequence name.
                    if locus_name in regions:
                        regions[accession] = regions.pop(locus_name)
                        locus_name = accession
                continue
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if line.startswith(("ORIGIN", "//", "CONTIG", "BASE COUNT")):
                if current is not None:
                    records.append(current)
                    current = None
                in_features = False
                continue
            if not in_features:
                continue

            feature = _FEATURE_LINE.match(line)
            if feature:
                if current is not None:
                    records.append(current)
                current = _Record(number, feature.group(1), feature.group(2))
                qualifier = None
                continue

            if current is None:
                continue

            qual = _QUALIFIER.match(line)
            if qual:
                qualifier = qual.group(1)
                value = qual.group(2)
                if value is None:
                    # A bare flag such as /pseudo. Recorded as present.
                    current.qualifiers.setdefault(qualifier, []).append("")
                    qualifier = None
                else:
                    current.qualifiers.setdefault(qualifier, []).append(value.strip('"'))
                continue

            # Continuation of a wrapped location or qualifier value.
            continuation = line.strip()
            if qualifier is not None and current.qualifiers.get(qualifier):
                previous = current.qualifiers[qualifier][-1]
                joiner = "" if previous.endswith("-") else " "
                current.qualifiers[qualifier][-1] = (previous + joiner + continuation).rstrip('"')
            elif not current.qualifiers:
                current.location += continuation

    if current is not None:
        records.append(current)
    return records, regions, comments, accession


# --------------------------------------------------------------------------
# Linkage. Every edge below is computed; none is read from the file.
# --------------------------------------------------------------------------

#: Types that are transcripts -- the middle layer, where one exists.
_TRANSCRIPTS = frozenset({"mRNA", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA", "precursor_RNA"})
#: Types that hang off a transcript in a eukaryote, or off a gene in a prokaryote.
_CODING = frozenset({"CDS"})
#: Records that describe the sequence rather than a feature on it.
_NON_FEATURE = frozenset({"source"})


def _identifier(record: _Record) -> str | None:
    """The identifier this record already carries.

    Order is deliberate. ``protein_id`` and ``transcript_id`` are unique per
    isoform; ``locus_tag`` names the *gene*, so it must come last or two
    isoforms would collide onto one identifier.

    When only ``locus_tag`` is available the record type is appended, giving
    ``MAE_RS31305:CDS``. Both halves come from the file -- this is a compound
    key, not an invented identifier (D-030). It is necessary because a
    prokaryotic pseudogene's CDS carries no ``protein_id`` (526 of 5,808 loci in
    *M. aeruginosa* NIES-843), so it would otherwise take its own gene's
    ``locus_tag`` and the two rows would be merged into one feature. Genes keep
    the bare tag, so the parent edge computed below still resolves.
    """
    for key in ("protein_id", "transcript_id"):
        value = record.first(key)
        if value:
            return value
    tag = record.first("locus_tag") or record.first("gene")
    if not tag:
        return None
    if record.type in ("gene", "pseudogene"):
        return tag
    return f"{tag}:{record.type}"


def _gene_identity(record: _Record) -> str | None:
    """Which gene a record belongs to.

    `/gene` where present, `/locus_tag` otherwise. Both are per-gene rather than
    per-transcript, so this narrows containment to one gene's own transcripts --
    it never resolves between isoforms, which is what containment is for.
    """
    return record.first("gene") or record.first("locus_tag")


def _contains(child: list[tuple[int, int]], parent: list[tuple[int, int]]) -> bool:
    """Every child span falls inside some parent span."""
    return all(
        any(pstart <= cstart and cend <= pend for pstart, pend in parent)
        for cstart, cend in child
    )


def _link_by_locus_tag(
    genes: list[tuple[_Record, str]],
    children: list[tuple[_Record, str, list[tuple[int, int]]]],
) -> dict[int, tuple[str, str]]:
    """Prokaryotic linkage. ``locus_tag`` is an explicit, unique key on both rows.

    Measured on *M. aeruginosa* NIES-843: 5,808 loci, every one exactly one gene
    plus exactly one child, zero loci with two CDSs. There is no ambiguity to
    resolve, so containment is not used and cannot introduce an error here.
    """
    by_tag: dict[str, str] = {}
    for record, identifier in genes:
        tag = record.first("locus_tag")
        if tag:
            by_tag[tag] = identifier

    edges: dict[int, tuple[str, str]] = {}
    for record, _identity, _spans in children:
        tag = record.first("locus_tag")
        if tag and tag in by_tag:
            edges[record.line_number] = (by_tag[tag], "locus_tag")
    return edges


def _link_by_containment(
    parents: list[tuple[_Record, str, list[tuple[int, int]]]],
    children: list[tuple[_Record, str, list[tuple[int, int]]]],
    path: Path,
    *,
    use_note: bool,
) -> dict[int, str]:
    """Eukaryotic CDS-to-transcript linkage, in three tiers.

    Measured on *D. melanogaster* Dscam1 (79 transcripts, 75 isoforms on one
    gene): containment alone resolves 75 of 79 uniquely. The 4 failures are
    two-isoform genes whose transcripts differ only in UTR length, so both
    contain the same CDS.

    Positional pairing -- link the nth CDS to the nth transcript -- was
    specified before these files were read, and is wrong: it mispairs 2 of 79,
    silently, because transcripts carry UTRs and therefore sort differently
    from their own coding sequences. It is not used.
    """
    # Transcripts indexed by gene, once, instead of rescanning the whole list per
    # CDS. On a 13,725-gene genome the linear scan made this quadratic -- 173
    # seconds against 10 for the equivalent GFF3 -- which would have been the
    # binding constraint on milestone 10's large-genome demo (D-050).
    by_gene: dict[str | None, list[tuple[_Record, str, list[tuple[int, int]]]]] = {}
    for parent, identifier, pspans in parents:
        by_gene.setdefault(_gene_identity(parent), []).append((parent, identifier, pspans))

    edges: dict[int, tuple[str, str]] = {}
    for record, _identity, spans in children:
        # Gene identity is `/gene` or `/locus_tag`, the same pair used to link a
        # transcript to its gene. They were inconsistent: containment read only
        # `/gene`, and RefSeq's *F. graminearum* GenBank release uses `/locus_tag`
        # exclusively -- so every CDS in the file fell into one group and was
        # compared against all 13,315 transcripts instead of its own gene's. It
        # failed on the first record. See D-050.
        gene = _gene_identity(record)
        candidates = [
            (parent, identifier)
            for parent, identifier, pspans in by_gene.get(gene, ())
            if _contains(spans, pspans)
        ]

        if len(candidates) == 1:
            # Containment resolved it uniquely. Recorded as `structural`, the
            # term D-007 and the paper use for this rule.
            edges[record.line_number] = (candidates[0][1], "structural")
            continue

        if len(candidates) > 1 and use_note:
            note = " ".join(record.qualifiers.get("note", []))
            match = _FROM_TRANSCRIPT.search(note + '"')
            if match:
                wanted = match.group(1)
                named = [
                    identifier
                    for parent, identifier in candidates
                    if wanted in " ".join(parent.qualifiers.get("note", []))
                    or wanted == parent.first("transcript_id")
                ]
                if len(named) == 1:
                    # Containment was ambiguous and /note named the transcript.
                    # A distinct method, not a variant of structural: it is an
                    # exact string in a specific qualifier, and `explain` must
                    # be able to say which of the two rules decided.
                    edges[record.line_number] = (named[0], "note")
                    continue

        if not candidates:
            # An orphan is left unlinked rather than attached to the gene: a
            # CDS whose coordinates escape every transcript of its gene is a
            # fact about the file, and hiding it under the gene would make the
            # hierarchy look complete when it is not.
            continue

        names = ", ".join(identifier for _parent, identifier in candidates)
        raise GenBankError(
            path, record.line_number,
            f"CDS {_identifier(record)!r} of gene {gene!r} is contained by "
            f"{len(candidates)} transcripts ({names}) and the file carries no "
            "qualifier that distinguishes them. cgload will not guess which "
            "transcript this coding sequence belongs to",
        )

    return edges


def tokenize_genbank(path: Path, *, linkage: str = "auto") -> TokenizedFile:
    """Read a GenBank flat file into ``RawFeature`` objects.

    ``linkage`` is ``locus_tag``, ``containment`` or ``auto``. ``auto`` chooses
    ``locus_tag`` when the file contains no transcript records at all, which is
    what a prokaryotic annotation looks like, and reports the choice.
    """
    records, regions, comments, accession = _read_records(path)
    if not records:
        raise GenBankError(path, 1, "no FEATURES table found")

    seqid = accession or next(iter(regions), "unknown")

    parsed: list[tuple[_Record, str | None, list[tuple[int, int]], str]] = []
    for record in records:
        if record.type in _NON_FEATURE:
            continue
        spans, strand = parse_location(record.location, path, record.line_number)
        parsed.append((record, _identifier(record), spans, strand))

    genes = [(r, i) for r, i, _s, _t in parsed if r.type in ("gene", "pseudogene") and i]
    transcripts = [(r, i, s) for r, i, s, _t in parsed if r.type in _TRANSCRIPTS and i]
    coding = [(r, i, s) for r, i, s, _t in parsed if r.type in _CODING and i]

    if linkage == "auto":
        # The discriminator is mRNA specifically, not "any transcript type".
        # A bacterial annotation has tRNA and rRNA records but no mRNA, so
        # testing for transcripts in general selects the eukaryotic strategy
        # for every prokaryote -- which leaves every CDS unlinked, because
        # containment needs a transcript layer that is not there.
        has_mrna = any(r.type == "mRNA" for r, _i, _s, _t in parsed)
        linkage = "containment" if has_mrna else "locus_tag"
    comments.append((1, f"# cgload linkage strategy: {linkage}"))

    edges: dict[int, tuple[str, str]] = {}
    if linkage == "locus_tag":
        edges.update(_link_by_locus_tag(genes, transcripts + coding))
    else:
        # Transcripts hang off genes by shared /gene or /locus_tag; only the
        # CDS-to-transcript edge is ambiguous.
        gene_by_name: dict[str, str] = {}
        for record, identifier in genes:
            for key in ("locus_tag", "gene"):
                value = record.first(key)
                if value:
                    gene_by_name.setdefault(value, identifier)
        for record, _identity, _spans in transcripts:
            for key in ("locus_tag", "gene"):
                value = record.first(key)
                if value and value in gene_by_name:
                    # Transcript-to-gene is never ambiguous: the qualifier is
                    # shared and unique within the record.
                    edges[record.line_number] = (gene_by_name[value], "locus_tag")
                    break
        edges.update(_link_by_containment(transcripts, coding, path, use_note=True))

    features: list[RawFeature] = []
    sources = {"GenBank"}
    for record, identifier, spans, strand in parsed:
        attributes: dict[str, list[str]] = {
            key: list(values) for key, values in record.qualifiers.items()
        }
        if identifier:
            attributes["ID"] = [identifier]
        edge = edges.get(record.line_number)
        if edge:
            parent, method = edge
            attributes["Parent"] = [parent]
            # D-032. GenBank states no parents, so every edge here was computed
            # and must carry the rule that produced it. `explicit_parent` would
            # be a false provenance claim, and it is the claim the containment
            # algorithm is published on.
            attributes[LINKAGE_ATTRIBUTE] = [method]
        if record.is_pseudo:
            # Consumed by the completeness check: a pseudogene CDS has no
            # translation and must not be measured against one.
            attributes["pseudo"] = ["true"]
        if "<" in record.location or ">" in record.location:
            attributes["partial"] = ["true"]

        phase = None
        codon_start = record.first("codon_start")
        if codon_start and codon_start.isdigit():
            phase = int(codon_start) - 1

        for span_start, span_end in spans:
            features.append(
                RawFeature(
                    line_number=record.line_number,
                    seqid=seqid,
                    source="GenBank",
                    type=record.type,
                    start=span_start,
                    end=span_end,
                    score=None,
                    strand=strand,
                    phase=phase,
                    attributes=attributes,
                )
            )

    return TokenizedFile(
        features=features,
        declared_regions=regions,
        comments=comments,
        sources=sources,
        data_line_count=len(features),
    )
