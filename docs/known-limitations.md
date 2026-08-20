# Known limitations

Stated here rather than left to be discovered.

## Provider data that cgload refuses or repairs

### Trans-spliced organellar genes

A gene assembled from segments on both strands cannot be represented in GFF3,
which allows one strand per row. Producers work around this in ways that are
individually reasonable and collectively inconsistent, and every one of the
following was found on a single reference assembly (*Zea mays* B73,
GCF_902167145.1), confined to seven genes on the two organellar sequences:

| symptom | cgload's response |
|---|---|
| parent extent excludes its own children (34 rows) | refuses; `--repair-envelopes` recomputes from children |
| a row where start exceeds end | refuses; repaired from children under the same flag |
| one gene emitted twice, on opposite strands | merged when `exception=trans-splicing` is present |
| exons with no parent at all | loaded as top-level features |

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

- **Peak memory is roughly 4.5 GB per million features.** Every feature is
  parsed before any row is written, so that cross-file identifier collisions are
  detected before the database is touched. A machine with 8 GB completes a
  700,000-feature genome; one with 4 GB does not.
- **Single-process only.** Concurrent loads against one database are untested.
- **MySQL and PostgreSQL are untested at scale.** The code path is shared with
  SQLite and covered by tests, but the large runs were SQLite.
- **Excluding sequences requires filtering the input.** There is no
  `--exclude-sequences` flag yet; to load part of an assembly, filter the GFF3
  and FASTA beforehand.
