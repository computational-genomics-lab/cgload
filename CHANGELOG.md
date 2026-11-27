# Changelog

Decisions referred to as D-0xx are in [docs/decisions.md](docs/decisions.md).

## 0.4.1 (27 November 2026)

What the benchmark's first BLAST run found (PROTOCOL.md section 14). BLAST's
searches, the portal's pages and everything else are those of 0.4.0.

### Fixed

- Two `cgload blastdb` building one assembly at the same time spoiled each
  other: each could clear the other's work, and the one that failed removed
  the databases the other had just built. On the benchmark's server two of 30
  genomes were left without BLAST databases. A build now holds a lock on its
  assembly; a second one is refused at once and changes nothing (D-117).
- When makeblastdb fails, the message includes the end of its own output.

### Changed

- `benchmark/run_blast.py` also records how near a search that did not find
  itself came (PROTOCOL.md section 14, run 2); what counts as found is
  unchanged.

## 0.4.0 (26 November 2026)

The portal searches its genomes with BLAST (D-116). Loading, verification,
export, queries and the portal's other pages are those of 0.3.1, and a
database made by 0.3.1 is served as it is.

### Added

- `cgload blastdb`: two BLAST+ databases per assembly, in its data directory:
  the genome, and every CDS translated as the portal translates it. Each is
  checked by counting: every sequence of the assembly, one protein for every
  CDS feature, and BLAST+'s own count of what it stored must agree. An
  assembly whose counts differ keeps no BLAST database, and the portal does
  not use databases older than the load they were built from.
- A BLAST page (`/blast`) and `POST /api/v1/blast`: blastn (megablast or
  blastn), tblastn, blastp and blastx, against one genome or all of them. Hits
  link to their CDS or to their region, with their alignments.
- A gene's page links its protein and its coding sequence to BLAST here, with
  the sequence filled in, and to NCBI BLAST, by accession where the CDS has an
  NCBI protein_id.
- `cgload serve --blast-jobs`, `--blast-threads` and `--blast-timeout`. A few
  searches run at once and as many again may wait; more are turned away with
  503, so searches cannot take every thread the pages need.
- `benchmark/run_blast.py` (PROTOCOL.md section 14): BLAST on the benchmark's
  30 genomes.

### Changed

- The portal answers POST on its two BLAST URLs; every other URL still
  answers GET and HEAD only. A request body is at most 64 KiB.
- `docs/portal.md` covers BLAST, and its nginx example limits BLAST requests
  per visitor. The examples' addresses are ones reserved for documentation.

BLAST+ is not a Python dependency: without it, the portal says that BLAST is
not set up.

## 0.3.1 (26 November 2026)

What 0.3.0 showed when other people first ran the portal (D-115). Loading,
verification, export, queries and the portal's pages are those of 0.3.0.

### Changed

- The portal installs with cgload: Flask and waitress are runtime
  dependencies. `cgload[portal]` still installs, as an empty extra.
- `cgload serve` listens on every network interface by default, so the other
  computers of a lab can open it, and lists the addresses that open it;
  `--host 127.0.0.1` keeps it to one computer. Behind a web server,
  `docs/portal.md` starts it with `--host 127.0.0.1`.
- When port 8000 is taken, `cgload serve` uses the next free port and says so;
  a port given with `--port` is refused in words if taken. Other start-up
  problems (a port that needs administrator rights, an address the computer
  lacks, a missing package, a database with no genomes) are one sentence, not
  a traceback, and Ctrl-C says the portal stopped.
- The README starts with a quick start and a table of what to do when
  something goes wrong; the tutorial has a step for the portal.

## 0.3.0 (26 November 2026)

The genomes cgload loads can now be put on the web for the people who use
them. Loading, `verify`, `export` and `query` are those of 0.2.2, and a
database made by 0.2.2 is read as it is; `cgload index` adds two tables to it.

### Added

- `cgload serve`: a genome portal over a cgload database
  ([docs/portal.md](docs/portal.md), D-111). Search by name, synonym, locus
  tag, product, cross-reference or functional term, in one genome or all;
  gene pages with the gene model drawn to scale, igv.js at the locus, the
  attributes with cross-references linked, the functional annotation and the
  whole family; regions; genomic, spliced, coding and protein sequences as
  FASTA, translated with the genetic code the file names (D-113); GFF3 and
  the assembly's files with byte ranges; each load's receipt; and a JSON API
  under `/api/v1/`. It needs the `[portal]` extra (Flask and waitress).
- `cgload serve` says who can open the portal: by default only the computer
  it runs on, and how to open it to a network (`--host 0.0.0.0`).
- `cgload index`: builds the portal's search index and checks it by counting:
  for each attribute it indexes, the features that carry it in the database
  against those that reach the index. An assembly whose counts differ keeps
  no index (D-112).
- The portal opens the database read-only, answers only GET and HEAD, binds
  every value it is sent, escapes all text from the files, and sends a
  Content-Security-Policy that allows no inline script and no script from
  another site; igv.js 3.8.10 ships with cgload (D-114).
  `tests/test_portal_security.py` sends injection strings and a seeded fuzz
  of 3,000 requests, and CI runs bandit on the code.
- `benchmark/run_portal.py` (PROTOCOL.md section 13): the portal on the
  benchmark's 30 genomes in one database.

### Changed

- `cgload remove` also removes the assembly's search index.
- `pytest` with no options now runs the default tier: `addopts` leaves out the
  real-genome tier, which fails without its files, and the server tier. CI asks
  for each tier with `-m`, as before.
- `survey/cgsurvey.py` is again byte-identical to tag `survey-v1`, as
  `survey/results/provenance.tsv` records. Two comments edited after the census
  ran had changed its MD5.
- The paper: corrections before submission (`WORKLOG.md`, 23 November 2026),
  and Table S1 now gives each input file's annotation release and MD5.

## 0.2.2 (19 November 2026)

`cgload export` now writes GFF3 columns 2 and 6 exactly as they were loaded.
Loading, `verify` and queries are those of 0.2.1, which the paper's benchmark
ran. Databases made by 0.2.1 are accepted, and features whose rows all name
one program export as before.

### Fixed

- Scores were written with six significant digits (`%g`), so 1234567.5 came
  back as 1.23457e+06 and 0.123456789 as 0.123457. Each score is now written as
  the shortest text that reads back as the stored number (D-109).
- A feature written on several rows kept only its first row's column 2, and
  `export` wrote it on every row, so a later row from another program lost its
  name. Each row now keeps and exports its own (D-110). Rows loaded by 0.2.1
  hold the first row's value; load them again to get each row's own.

Neither change alters a count, so neither could be seen by `verify` or by the
benchmark's round-trip comparison, which did not read columns 2, 6 and 8. A
row-by-row comparison of those columns (`benchmark/columns.py`) found both. In
the benchmark's 30 genomes the first affected 52 scores, in human and mouse,
and the second none.

### Added

- `benchmark/columns.py` and `benchmark/run_columns.py`: the round trip of
  columns 2, 6 and 8, row by row, on the 30 genomes of the benchmark.

## 0.2.1 (11 November 2026)

No change in behaviour: the code that loads, verifies, exports and queries is
that of 0.2.0. This is the version the paper's benchmark ran; the paper
describes 0.2.2.

### Changed

- The README's Scale section gives the benchmark's measurements on 30 genomes
  in place of one early maize run, and `docs/known-limitations.md` and
  `docs/scalability.md` follow it.
- Nine decision numbers were used twice (D-028 to D-036). The later decision of
  each pair is now D-028b to D-036b, in `docs/decisions.md` and in every comment
  that cites it.
- `preflight_gff3.py` moved from the repository root to `tools/`.

### Removed

- `rename.py`, the script that renamed the package once, and `paper/paper.md`,
  an early draft superseded by `paper/manuscript.md`. Both remain in the
  history (`git show v0.2.0:rename.py`).

## 0.2.0 (10 November 2026)

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
- Every RefSeq load but one warned about feature types "not declared by the
  refseq profile", and `--strict-vocabulary` refused 29 of the benchmark's 30
  RefSeq files: the profile declared 9 of the 91 types they use. It now
  declares all 91 (D-107).
- `query --id` on PostgreSQL read the whole feature table twice: 12 s for a
  human gene. It now fetches the features it selects by name, using the index
  (D-108).
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
