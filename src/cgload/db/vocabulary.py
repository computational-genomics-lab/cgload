"""The controlled vocabulary, its tiers, and its content hash.

Two properties this module exists to guarantee (docs/decisions.md D-017):

1. Normalisation is for linkage, never for counting. ``feature.source_type``
   holds column 3 verbatim and is what ``verify`` counts; ``feature_type``
   holds the canonical term and is what linkage and ``query`` use. The two
   sides of the primary completeness check therefore share no input.
2. The mapping is data, but versioned data. Adding a row changes the meaning
   of linkage, so the seed carries a version and a content hash, the hash is
   recorded on every assembly, and a load against a different hash is refused.
"""

from __future__ import annotations

import hashlib

VOCABULARY_VERSION = 1

#: Tier A -- the canonical terms the loader knows how to link parent to child.
#: Closed, because a parent-child rule cannot be inferred for an unseen term.
#: Every entry here is a *canonical* term; source-side synonyms live in
#: SYNONYMS and must not be repeated here.
HIERARCHY_TERMS: tuple[str, ...] = (
    "gene",
    "pseudogene",
    "mRNA",
    "tRNA",
    "rRNA",
    "ncRNA",
    "exon",
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "repeat_region",
)

#: Terms deliberately kept OUT of Tier A and out of SYNONYMS, with a note, so
#: nobody re-adds them. Tested.
NOT_CANONICALISED: dict[str, str] = {
    "start_codon": "AUGUSTUS and funannotate include the start codon inside "
                   "their terminal CDS rows, so mapping it to CDS would count "
                   "the same bases twice and inflate every CDS count. It loads "
                   "as an inert feature with its explicit Parent.",
    "stop_codon": "Same as start_codon. GTF's stop codon sits OUTSIDE the CDS "
                  "and is absorbed by that profile's coordinate reconciliation "
                  "(D-013), which is a coordinate operation, not a type mapping.",
    "intron": "AUGUSTUS emits introns; RefSeq and funannotate do not. Derived "
              "information, stored verbatim, never synthesised or compared "
              "across dialects.",
    "region": "RefSeq's first row declares the sequence region and carries "
              "taxon and strain. Consumed by the profile as region metadata, "
              "not as a feature in the gene hierarchy.",
}

#: Source term -> canonical term. Only synonyms that resolve *into* Tier A
#: belong here; a synonym pair among inert terms would be tier "B".
SYNONYMS: dict[str, str] = {
    "transcript": "mRNA",
    "primary_transcript": "mRNA",
    # RefSeq annotates a microRNA as `gene -> primary_transcript -> miRNA ->
    # exon`: two transcript tiers, where cgload's model has one. `miRNA` was
    # unmapped, so it loaded inert, and an inert parent makes its Tier A
    # children unloadable -- the open vocabulary is open in one direction only.
    # Measured on *Zea mays* B73 nuclear GFF3: 310 miRNA features under 164
    # `primary_transcript` rows (the 5p/3p arm pattern), every one carrying
    # exons, every one refused (D-062).
    "miRNA": "ncRNA",
    "5'UTR": "five_prime_UTR",
    "3'UTR": "three_prime_UTR",
    "pseudogenic_transcript": "mRNA",
    "pseudogenic_exon": "exon",
    # The Sequence Ontology's pseudogenic transcripts of RNA genes, each a kind of
    # pseudogenic_transcript and mapped as that one is. RefSeq's rice and grape
    # mitochondria write a `pseudogenic_tRNA` with an exon below it, and the exon
    # was refused because its parent loaded inert (D-097).
    "pseudogenic_tRNA": "mRNA",
    "pseudogenic_rRNA": "mRNA",
    # INSDC's `precursor_RNA` feature key, which older NCBI GFF3 writes in
    # column 3 where newer files write `primary_transcript`: a microRNA's
    # precursor, with its exons and its mature miRNA below it. *Xenopus
    # tropicalis* (GCF_000004195.4) has them, and every child was refused (D-097).
    "precursor_RNA": "mRNA",
    "misc_RNA": "ncRNA",
    "snRNA": "ncRNA",
    "snoRNA": "ncRNA",
    "tmRNA": "ncRNA",
    "lnc_RNA": "ncRNA",
    "lncRNA": "ncRNA",
    "SRP_RNA": "ncRNA",
    "RNase_P_RNA": "ncRNA",
    "antisense_RNA": "ncRNA",
    "guide_RNA": "ncRNA",
    "telomerase_RNA": "ncRNA",
    # The rest of the INSDC /ncRNA_class vocabulary, as NCBI's GFF3 writes it in
    # column 3. Each is a transcript that carries exons, like snRNA above. Left
    # unmapped, each loaded inert and its exons were refused, the same failure
    # D-062 fixed for miRNA one class at a time. Vertebrate RefSeq files carry
    # Y_RNA (the RNY genes), vault_RNA (VTRNA) and RNase_MRP_RNA (RMRP) (D-088).
    # `pre_miRNA` is left out: it would put a miRNA under an ncRNA, which the
    # hierarchy does not permit, and RefSeq writes `primary_transcript` instead.
    "Y_RNA": "ncRNA",
    "vault_RNA": "ncRNA",
    "RNase_MRP_RNA": "ncRNA",
    "scRNA": "ncRNA",
    "scaRNA": "ncRNA",
    "piRNA": "ncRNA",
    "siRNA": "ncRNA",
    "rasiRNA": "ncRNA",
    "ribozyme": "ncRNA",
    "hammerhead_ribozyme": "ncRNA",
    # Transcript-level, not gene-level. Measured on RefSeq human GRCh38.p14
    # NC_000014.9, the IGH locus: all 175 segment rows carry a `Parent`, that
    # parent is a `gene`, and 164 CDS plus 162 exon rows hang beneath them. They
    # occupy exactly the position an mRNA occupies.
    #
    # The previous mapping to `gene` was self-contradictory: `ALLOWED_PARENTS`
    # gives `gene` an empty parent set, so every one of these rows was refused
    # with "a gene may only attach to nothing" (D-063).
    "V_gene_segment": "mRNA",
    "C_gene_segment": "mRNA",
    "D_gene_segment": "mRNA",
    "J_gene_segment": "mRNA",
}

#: Which canonical parents each canonical Tier A type may attach to. A Tier A
#: feature whose parent resolves to a term outside Tier A is an error: linkage
#: rules are written per Tier A type and assume Tier A parents.
ALLOWED_PARENTS: dict[str, frozenset[str]] = {
    "gene": frozenset(),
    "pseudogene": frozenset(),
    "repeat_region": frozenset(),
    "mRNA": frozenset({"gene", "pseudogene"}),
    "tRNA": frozenset({"gene", "pseudogene"}),
    "rRNA": frozenset({"gene", "pseudogene"}),
    # `mRNA` is permitted because a precursor transcript genuinely contains a
    # mature one: RefSeq's miRNA sits under `primary_transcript`, which
    # canonicalises to `mRNA`. This is a real nested relationship, not a
    # modelling compromise. Verified uniform on maize: all 310 miRNA parents are
    # `primary_transcript`, none attach directly to a gene (D-062).
    "ncRNA": frozenset({"gene", "pseudogene", "mRNA"}),
    # `pseudogene` is permitted because RefSeq attaches exons directly to a
    # pseudogene row, with no transcript between them: a pseudogene has no
    # functional transcript to model. Measured on *Zea mays* B73: 5,195
    # pseudogenes, 12,757 exons on gene-level parents. `CDS` already listed
    # pseudogene; `exon` did not, and the omission refused standard NCBI output
    # (D-058).
    "exon": frozenset({"mRNA", "tRNA", "rRNA", "ncRNA", "gene", "pseudogene"}),
    # `gene` is permitted because a prokaryotic annotation has no transcript
    # layer at all: in GenBank, a bacterial CDS hangs directly off its gene, and
    # requiring mRNA would refuse every prokaryotic record. It is not a licence
    # to skip the transcript where one exists -- see TRANSCRIPT_LAYER_REQUIRED.
    "CDS": frozenset({"mRNA", "gene", "pseudogene"}),
    # PREDICTED, NOT OBSERVED. No file in evidence attaches a UTR to a gene;
    # this is the same shape D-036 fixed for `CDS` and D-058 for `exon`, and the
    # UTRs were the remaining member of that family. RefSeq human emitted no UTR
    # rows in the IGH region examined, so the failure is inferred from family
    # resemblance rather than measured. Widened anyway because it costs nothing
    # and the forbidden version is easy to state: a UTR may not attach to a CDS,
    # to another UTR, or to a repeat_region (D-064).
    "five_prime_UTR": frozenset({"mRNA", "gene", "pseudogene"}),
    "three_prime_UTR": frozenset({"mRNA", "gene", "pseudogene"}),
}

#: Canonical types that constitute a transcript layer. A gene with any of these
#: as a child has a transcript; a gene with none does not.
TRANSCRIPT_TERMS: frozenset[str] = frozenset({"mRNA", "tRNA", "rRNA", "ncRNA"})

#: Types that must not attach directly to a gene *when the file contains a
#: transcript layer*. Checked per file rather than per profile, because one
#: GenBank profile serves both kingdoms: the bacterial records have no mRNA and
#: the eukaryotic ones do. Without this, widening CDS above would let a
#: eukaryotic CDS silently bypass its transcript and the missing mRNA layer
#: would look like a valid hierarchy.
# `ncRNA` is deliberately NOT here, though D-060's standing rule would suggest
# adding it when ALLOWED_PARENTS["ncRNA"] was widened (D-062). The rule does not
# apply: `ncRNA` is itself a member of TRANSCRIPT_TERMS, so an ncRNA child makes
# its own gene count as transcript-bearing and the check can never fire on it.
# Adding it is a no-op that reads as protection.
#
# More importantly there is nothing to conceal. `gene -> ncRNA` is the normal
# shape of a non-coding gene: the ncRNA *is* the transcript, not a child that
# skipped one. The concealment the rule guards against is a coding gene whose
# mRNA is missing, and that is unaffected by the miRNA widening.
TRANSCRIPT_LAYER_REQUIRED: frozenset[str] = frozenset({"CDS", "exon"})

#: `genome=` values (on RefSeq's `region` rows) that name an organelle. On such a
#: sequence a CDS or exon may attach straight to its gene, because organelle
#: genomes are annotated the prokaryotic way: NCBI writes a vertebrate
#: mitochondrial CDS directly under its gene (AUDIT-5, D-080). A nucleomorph is
#: a reduced nucleus with eukaryotic genes, so it is not listed. Hashed below:
#: it decides whether a `gene`->`CDS` edge is accepted (D-061).
ORGANELLE_GENOMES: tuple[str, ...] = (
    "mitochondrion", "plastid", "chloroplast", "apicoplast", "chromoplast",
    "cyanelle", "kinetoplast", "leucoplast", "proplastid", "chromatophore",
    "hydrogenosome",
)

#: How a feature's parent was determined. NULL means the feature has no
#: parent, which is correct for a top-level or an inert feature.
#:
#: This is the column that makes D-007's rule hierarchy auditable after the
#: fact. `explain` reports it, and it is the difference between "the file said
#: so" and "we worked it out" -- the distinction the paper's claim rests on.
LINKAGE_METHODS: tuple[str, ...] = (
    "explicit_parent",  # GFF3 Parent attribute; the common case
    "locus_tag",  # GenBank: /locus_tag or /gene shared by parent and child
    "qualifier",  # GenBank: /transcript_id present on BOTH mRNA and CDS
    "note",  # GenBank: /note names the parent transcript explicitly
    "structural",  # exon-interval containment (D-007); "containment" in prose
    "positional",  # order of appearance. Never used -- see D-034.
    "synthesised",  # GTF: parent row absent from the file and constructed
)

#: Methods that mean "cgload computed this edge" rather than "the file stated
#: it". A GenBank feature may never be recorded as ``explicit_parent``: the
#: format carries no parent pointers, so that value would be a false claim
#: about provenance.
COMPUTED_LINKAGE: frozenset[str] = frozenset({"structural", "positional"})


def seed_rows() -> list[dict[str, str]]:
    """Every row written to the ``vocabulary`` table at ``init``.

    Terms absent from this list are not rejected: they load verbatim with
    ``source_type == feature_type`` and are counted normally. A type the
    file's dialect profile does not declare is named in the load summary, and
    ``load --strict-vocabulary`` refuses it (D-087).
    """
    rows = [{"term": t, "canonical_term": t, "tier": "A"} for t in HIERARCHY_TERMS]
    rows += [
        {"term": src, "canonical_term": canon, "tier": "A"}
        for src, canon in sorted(SYNONYMS.items())
    ]
    return rows


def content_hash() -> str:
    """SHA-256 over the seed, the hierarchy rules and the linkage vocabulary.

    Recorded on every assembly. Anything that changes how a source term is
    linked must change this hash, so the hierarchy rules are hashed alongside
    the term mapping -- editing ALLOWED_PARENTS without bumping the hash would
    otherwise silently change the meaning of an existing database.
    """
    h = hashlib.sha256()
    h.update(f"version={VOCABULARY_VERSION}\n".encode())
    for row in seed_rows():
        h.update(f"{row['term']}\t{row['canonical_term']}\t{row['tier']}\n".encode())
    for child in sorted(ALLOWED_PARENTS):
        parents = ",".join(sorted(ALLOWED_PARENTS[child]))
        h.update(f"{child}<-{parents}\n".encode())
    h.update(("|".join(LINKAGE_METHODS) + "\n").encode())
    # The transcript-layer constants decide whether a `gene`->`CDS` or
    # `gene`->`exon` edge is accepted at all, so they change how a source term is
    # linked and the docstring above already promises they are covered. They were
    # not: D-058 narrowed the transcript-layer rule and only the ALLOWED_PARENTS
    # half of that change moved the hash. An edit to either set alone would have
    # changed which files load while every existing assembly still reported a
    # matching hash (D-061).
    h.update(("transcript_terms=" + ",".join(sorted(TRANSCRIPT_TERMS)) + "\n").encode())
    h.update(
        ("transcript_layer_required=" + ",".join(sorted(TRANSCRIPT_LAYER_REQUIRED)) + "\n").encode()
    )
    h.update(("organelle_genomes=" + ",".join(sorted(ORGANELLE_GENOMES)) + "\n").encode())
    return h.hexdigest()
