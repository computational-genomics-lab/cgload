# Changelog

Decisions referred to as D-0xx are in [docs/decisions.md](docs/decisions.md).

## Unreleased

Databases made by 0.1.0 are refused: the schema and the vocabulary changed.
The refusal says so. Load into a new database made with `cgload init`.

### Added

- `cgload query`: features by ID (optionally with everything below), by
  region, or by a GO, KEGG, Pfam, EC or InterPro accession; TSV or GFF3; read
  only (D-085).
- `cgload remove`: an assembly's rows and its data directory go together, or
  neither does; `--dry-run`, `--keep-files`, `--yes` (D-084).
- `load` writes the annotation into the data directory as bgzipped,
  tabix-indexed GFF3 for genome browsers (D-082).
- `verify` checks that the stored FASTA is still the one loaded (D-079), and
  runs the parent-edge tier, which it never did before (D-078). `load` runs the
  parent-edge tier on GFF3 files and rolls back if it fails.
- `load --strict-vocabulary` refuses feature types the dialect profile does
  not declare (D-087).
- The INSDC ncRNA classes (Y_RNA, vault_RNA, RNase_MRP_RNA and others) are
  read as ncRNA, so their exons load (D-088).
- PostgreSQL is tested in CI.
- Region queries use genomic bins, and a PostgreSQL load collects planner
  statistics before it commits; a gene lookup no longer scans the table
  (D-094).
- `benchmark/`: the protocol (`benchmark/PROTOCOL.md`) and the scripts of the
  benchmark the paper reports. Not part of the installed package.

### Fixed

- Six RefSeq genomes of the benchmark's thirty were refused over valid
  content: `precursor_RNA`, `pseudogenic_tRNA` and `pseudogenic_rRNA` now load
  as transcripts (D-097); the rows of a trans-spliced gene may overlap when
  they lie on opposite strands (D-098); a CDS may read a base twice where its
  `exception` says `low-quality sequence region` (D-099); and an exon feature
  NCBI writes on several rows with `gbkey=exon` keeps its shared ID (D-100).
- A load on PostgreSQL failed with a traceback on an attribute value longer
  than about 2,700 bytes, as *Arabidopsis* has; the index on attribute values
  is gone (D-101).
- `verify --against` reported a loss for each feature type NCBI names
  differently in its GFF3 and GenBank files, such as `snoRNA` for an `ncRNA`
  record and `long_terminal_repeat` for a `repeat_region`; the tier now
  compares feature families (D-102).
- `verify --against` reported lost mRNAs on fungal genomes: NCBI's GFF3 gives
  an mRNA to every gene whose CDS its GenBank file states alone, the genes of a
  mitochondrion and pseudogenes such as yeast's FLO8. Those mRNAs are now set
  aside gene by gene, and the report says how many (D-103, D-105). A family
  member written in capitals, as *S. pombe*'s `LONG_TERMINAL_REPEAT`, is
  matched to its family (D-106).
- A feature written on several rows, from the highest coordinate down as NCBI
  writes the minus strand, lost its lowest row's own attributes on load: the
  row came back with the first row's `Target`, `Gap` or `end_range` (D-104).
- A feature written on several rows lost what differed between its rows: each
  row of a `cDNA_match` alignment was exported with the first row's `Target`
  and `Gap`, a partial CDS lost `end_range`, and a trans-spliced gene's rows
  all came back on strand `?`. Each row now keeps its own (D-096).
- A load killed after putting its data directory in place but before
  committing left a directory that made the repeated load fail. The next load
  into the same database now recognises such a directory by the database ID in
  its manifest and replaces it, saying so (D-095).

- Rows without an ID (AUGUSTUS introns and codons) lost their stated parent
  and their attributes (D-077).
- Export wrote a multi-valued attribute such as `Dbxref=a,b` as one value, and
  failed on genomes with more rows than the driver's parameter limit; so did
  loading functional annotation.
- `init` and `stats` failed on PostgreSQL; coordinates above 2,147,483,647
  could not be stored there.
- `load --force` failed on any assembly with features, and could leave rows
  pointing at a directory that had been moved aside (D-083).
- Rows after a `##FASTA` line were read as features.
- Rows past the end of a linear sequence loaded silently, and a child on a
  different sequence from its parent was linked (D-081).
- A vertebrate file including its mitochondrion was refused (D-080).
- `verify` blamed the loader when an input file had been edited after the load
  (D-086).
- Messages named a `--dialect` flag and a `cgload remap` command, neither of
  which exists.
- The CI workflow file was not valid YAML, so CI did not run.
- The parent-edge tier took the file's parents from the parser's own record.
  A parser that dropped a `Parent` attribute stored the feature as
  top-level and the tier passed. The independent pass now counts the
  file's `Parent=` statements per type, and the tier compares those with
  the links stored (D-092).
- A feature written on several rows could lose one of them without any check
  noticing, because the source count counted features only. It now counts
  rows too (D-093).

## 0.1.0

First version under the name cgload. Not released on PyPI.
