---
title: "cgload: relational loading of genome annotations with dialect reconciliation and completeness verification"
tags:
  - Python
  - genomics
  - genome annotation
  - GFF3
  - GenBank
  - data quality
authors:
  - name: Aditya Upadhyay
    orcid: 0000-0000-0000-0000
    affiliation: "1, 2"
  - name: Arijit Panda
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Sucheta Tripathy
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: "1, 2"
affiliations:
  - name: Structural Biology and Bioinformatics Division, CSIR-Indian Institute of Chemical Biology, Jadavpur, Kolkata 700032, India
    index: 1
  - name: Academy of Scientific and Innovative Research (AcSIR), Ghaziabad 201002, India
    index: 2
  - name: Department of Quantitative Health Sciences, Mayo Clinic, Rochester, MN 55905, USA
    index: 3
date: DD Month 2026
bibliography: paper.bib
---

<!-- ORCIDs are placeholders. JOSS will not proceed without real ones. -->

# Summary

`cgload` loads a genome assembly and its structural and functional annotation
into a queryable relational database, reconciles the incompatible conventions
that different annotation pipelines use, and verifies that every feature in the
input reached the database with its parent–child structure intact.

The verification is the point. A loader that silently drops features leaves a
database that looks populated, answers queries, and is wrong, with nothing
downstream able to detect it. `cgload` treats that as the failure mode to
design against: a load either completes in full or leaves the database and the
filesystem untouched, and every load is followed by a count that must reconcile
against an independent reading of the input.

The guarantee is completeness rather than fidelity: it establishes that no
feature was dropped or invented, not that every field of every feature is
faithful. A count-preserving corruption is outside what a cardinality check can
detect, and `cgload` does not claim otherwise.

# Statement of need

Turning a genome annotation into something queryable is a recurring task for
any group working on non-model organisms, and it is routinely done with ad hoc
scripts. Two properties make that harder than it appears.

**Annotation formats vary by producer, not just by extension.** NCBI RefSeq,
funannotate and AUGUSTUS all emit GFF3, and all three disagree. AUGUSTUS emits
no `exon` rows at all, describing gene structure through `intron` records
instead, and leaves half its child rows without identifiers. funannotate names a
transcript by appending `-T1` to its gene's identifier, so the obvious
normalisation makes a child indistinguishable from its parent. RefSeq
percent-encodes reserved characters in attribute values, and attaches exons
directly to a `pseudogene` with no transcript between them. GenBank flat files
carry no parent pointers whatsoever: which coding sequence belongs to which
transcript must be inferred from coordinates. A parser written against one
producer's output fails quietly on another's.

**Silent partial loading is the dominant failure mode, and per-row validation
cannot see it.** Checks that examine one row at a time — start before end,
strand in the permitted set, coordinates within the sequence — miss defects that
only appear in the relationships between rows. On the *Zea mays* B73 reference
assembly (GCF_902167145.1), `cgload` identifies 36 structural problems in the
annotation as published; **34 of them satisfy every per-row rule** and would
load without complaint, leaving parent extents wrong by up to 281 kb. Only a
check on the parent–child relationship exposes them.

Existing tools address adjacent problems. `gffutils` [@Dale:gffutils] and
`BCBio.GFF` parse GFF3 into local databases but treat the file as authoritative
and do not reconcile across producers. Tripal [@Sanderson:2013] and its Chado
schema [@Mungall:2007] provide a full community database platform, at the cost
of a Drupal deployment. JBrowse 2 [@Diesh:2023] visualises indexed annotation
without a server. None of them counts what it stored and requires that count to
agree with an independent reading.

# Implementation

Input is normalised to a single intermediate representation before any database
contact. A dialect is a declaration — how it identifies itself, which feature
types it provides, how parentage is determined — rather than a branch in the
parser, and adding two dialects after the first required no change to the
reading engine. A file that matches no profile is refused, and the refusal names
the profiles considered.

Feature types use an open vocabulary. Types the profile does not declare are
loaded, counted and reported rather than dropped or refused; the maize load
stores 21 source types of which 13 are undeclared. Each feature retains both the
term the file used and the canonical term, so normalisation is never lossy.

Verification runs in three tiers, ranked by what disagreement can prove:

1. **Source count** — a deliberately naive pass over the input, sharing no code
   with the parser, counted against the database. A verifier built from the
   parser can only detect storage errors; parsing errors are invisible to a
   check that parses the same way.
2. **Cross-parser** — the same assembly read as GFF3 and as GenBank by two
   engines with no shared code. On *Fusarium graminearum* PH-1 the two agree
   exactly on 13,725 genes, 13,315 mRNAs, 13,312 coding sequences, 322 tRNAs and
   88 rRNAs. A reader that loses a feature from one format cannot make the other
   format's count lose it too.
3. **Reference** — the data provider's published per-assembly feature totals.

A tier fails the command only when disagreement implies `cgload` is wrong;
otherwise it reports why it did not run. A check that vanishes when it cannot
run is indistinguishable from one that passed.

Assembly sequence is stored as indexed FASTA beside the database rather than in
it, following the convention of samtools [@Danecek:2021] and every genome
browser. The database holds coordinates and a reference. This keeps the database
portable and makes the output directly consumable by JBrowse 2 without
conversion.

# Performance

The complete *Z. mays* B73 assembly — 687 sequences, 2.18 Gb, 696,928 features
across 21 types, with eggNOG-mapper [@Cantalapiedra:2021] functional annotation
for 55,463 proteins — loads in 2 min 34 s with a peak resident set of 5.4 GB,
and verifies in 0.3 s with an exact source count. *F. graminearum* PH-1, at
36.5 Mb and 78,715 features, loads in 12.9 s with a peak of 0.78 GB.

Maize carries 60 times the bases of *F. graminearum* and 8.9 times the features;
load time rises 11.9-fold and peak memory 6.7-fold. Memory therefore scales with
feature count rather than genome size, since no sequence is held in memory —
approximately 7 GB of resident set per million features, plus a fixed overhead
near 200 MB. Verification is effectively constant: 41 MB on maize against 64 MB
on the smaller genome, because the counting pass streams.

# Availability

<!-- The repository URL and the package name are both provisional: `cgload` is
     taken on PyPI by an unrelated project. See docs/NAME.md. Update this
     section, and the three placeholder ORCIDs above, before submission. -->

`cgload` is available from
`https://github.com/computational-genomics-lab/cgload` under the MIT licence.
Documentation, a tutorial against bundled example data, and a reproducible
protocol for the measurements above are in the repository.

# Acknowledgements

<!-- Funding, computing resources, and anyone who supplied test data. -->

# References
