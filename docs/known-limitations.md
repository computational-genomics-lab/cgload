# Known limitations

Stated here rather than left to be discovered.

## Provider data that cgload refuses or repairs

### Trans-spliced organellar genes

A gene assembled from segments on both strands cannot be represented in GFF3,
which allows one strand per row. Producers work around this in ways that are
individually reasonable and collectively inconsistent, and every one of the
following was found on a single reference assembly (*Zea mays* B73,
GCF_902167145.1), on its organellar sequences:

| symptom | cgload's response |
|---|---|
| parent extent excludes its own children | refuses; `--repair-envelopes` recomputes from children |
| a row where start exceeds end | refuses; repaired from children under the same flag |
| one gene emitted twice, on opposite strands | merged when `exception=trans-splicing` is present |
| exons with no parent at all | loaded as top-level features |

cgload's own count on maize, from the load of 15 October 2026
(`envelope_census.py`): 11 envelope problems in two mitochondrial mRNAs,
rna-ZeamMp017 and rna-ZeamMp186. Ten are children lying outside their mRNA;
the eleventh is rna-ZeamMp017, written backwards. `--repair-envelopes`
recomputes both mRNAs and leaves no problem; the larger change widens
rna-ZeamMp186 by 281,540 bp.

*Arabidopsis thaliana* (GCF_000001735.4) has 13 such problems on its
mitochondrion, in three genes. Repair clears 11. The other two are mRNAs that
NCBI writes past the origin of the circular sequence, with its length added,
under genes written in plain coordinates on several rows. A parent written on
several rows is never recomputed (D-046), so these two are stored as written,
named in a warning and recorded in the manifest (D-091).

`--repair-envelopes` is off by default. Every repair is printed at load time and
recorded in `manifest.json`, so a repaired load is declared rather than merely
successful. Repair proceeds from the innermost features outward; a parent is
never repaired from a child that is itself wrong.

### Published feature counts can disagree with the annotation

The reference tier compares against a provider's summary, which is derived data
and can be inconsistent with the primary files. On GCF_902167145.1 the published
counts include organellar coding sequences but no organellar transcripts, while
the GFF3 contains 274 of each. The tier correctly reports the disagreement; it
is not a cgload defect and cannot be resolved in cgload.

This is why the tiers are ranked as they are. The cross-parser check compares two
*primary* representations of an assembly and is the stronger evidence.

### Rows past the end of a sequence

A row past the end of a linear sequence is refused. On a sequence whose `region`
row says `Is_circular=true`, a row crossing the origin loads normally, and a row
lying wholly past the end is stored as written and named in a warning (D-081).
NCBI writes such rows for trans-spliced organelle genes, for example
*Arabidopsis* *nad5* and rice *rps12*.

### Organelle genes without a transcript layer

A CDS attached directly to its gene is accepted on a sequence whose `region` row
declares an organelle genome (`genome=mitochondrion`, `plastid`, ...), which is
how NCBI annotates vertebrate mitochondria (D-080). A file with no `region` rows
declares no organelles, so the transcript-layer rule then applies to every
sequence in it.

## Format and scope

- **Multi-parent features are refused.** GFF3 permits a feature to name several
  parents, which is how the specification represents an exon shared between
  isoforms. cgload stores one parent per feature, so such a file is refused
  rather than silently attached to one branch. Affects annotations with shared
  exons; not encountered in the fungal, plant or bacterial assemblies tested.
- **GTF is not supported.** It needs stop-codon reconciliation that no supported
  dialect requires. Deferred rather than half-built.
- **Multi-segment parents are never repaired.** Widening one segment of several
  means choosing which, and the file gives no basis for that choice.
- **UTR-to-gene attachment is permitted but unobserved.** Allowed by analogy
  with `CDS` and `exon`, which needed it. No file in evidence exercises it.

## Functional annotation

- **Gene Ontology counts include inherited ancestors.** A protein annotated with
  one specific term also carries every broader term above it, which is why the
  mean is around 90 terms per protein. Counts are always reported as *ontology
  terms including inherited ancestors*; cgload cannot separate direct from
  inherited assignments without the ontology file, which is a large versioned
  download deliberately not bundled.
- **InterProScan support is untested against real output.** The parser follows
  the documented column layout, but no real InterProScan file has been run
  through it.
- **Sequence-similarity search, variants, expression data and orthology are out
  of scope.** No loaders, no tables, no plans.

## Operational

- **Peak memory is roughly 7 GB per million features, plus about 200 MB**
  (docs/scalability.md). Every feature is
  parsed before any row is written, so that cross-file identifier collisions are
  detected before the database is touched. A machine with 8 GB completes a
  700,000-feature genome; one with 4 GB does not.
- **Single-process only.** Concurrent loads against one database are untested.
- **MySQL and PostgreSQL are untested at scale.** PostgreSQL 16 runs a small
  end-to-end test in CI (`tests/test_postgresql.py`); MySQL has no test at all.
  The large runs were SQLite.
- **The tabix index covers positions below 2^29.** A longer sequence gets its
  `annotation.gff3.gz` but no `.tbi`; the load says so, and
  `tabix --csi -p gff` indexes it (D-082).
- **`query` needs recursive common table expressions**: SQLite 3.8.3 or
  newer (released in 2014), PostgreSQL, or MySQL 8.0 or newer.
- **Region queries read whole bins.** Each row is stored with the smallest
  genomic bin that holds it (the scheme of the SAM and CSI indexes, D-094), and
  a region query reads the rows of the few dozen bins that overlap the region.
  A region covering more than about 500 of those bins reads the whole
  sequence's rows instead.
- **An interrupted load leaves its staging directory behind.** A load killed
  part-way changes nothing in the database or the stored assemblies, but the
  copy of the FASTA it was making stays under `<data-dir>/.staging/`. It
  belongs to no assembly and can be deleted when no load is running.
- **A Pfam domain has two spellings.** eggNOG-mapper reports it by name
  (`Pkinase`), InterProScan by accession (`PF00069`). `query --term` matches
  the string given, and cgload holds no table that maps one to the other.
- **Excluding sequences requires filtering the input.** There is no
  `--exclude-sequences` flag yet; to load part of an assembly, filter the GFF3
  and FASTA beforehand.
