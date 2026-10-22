# Benchmark protocol

Version: `benchmark-protocol-v1`. This file and the scripts beside it were
committed and tagged before any benchmark run. Results are reported against
it. A change after the first run is listed in section 12 with its reason and
date, and every measurement it affects is repeated.

## 1. Questions

Each question has one section of the paper's Results.

| | Question | Section |
|---|---|---|
| Q1 | Does cgload load and verify all 30 RefSeq genomes of the survey panel, on SQLite and on PostgreSQL, and at what cost in time, memory and disk? | 7.1 |
| Q2 | Does what cgload stores match the input, feature by feature? | 7.1 |
| Q3 | What do other loaders keep, change or lose from the same files, and do they say so? | 7.2 |
| Q4 | Which damage to a loaded genome, or inside the parser, does cgload notice? | 7.3 |
| Q5 | What does an interrupted load leave behind? | 7.4 |
| Q6 | How long does a user wait for everyday answers? | 7.5 |
| Q7 | Do cgload's structural checks and the survey's checker agree on the same files? | 7.6 |
| Q8 | Does cgload work on one of our own genomes, with its functional annotation? | 7.7 |

## 2. Genomes

The 30 genomes of the survey panel (`survey/panel.tsv`, tag `survey-panel-v1`):
8 vertebrates, 7 insects, 7 plants, 8 fungi. `genomes.tsv` lists them with the
parts of the benchmark each one takes part in:

* all 30: Q1, Q2, Q3, Q7;
* the 8 fungi: also the cross-parser tier of `verify`, with each one's GenBank file;
* *F. graminearum*, *A. thaliana*, human: Q4;
* *F. graminearum*, *A. thaliana*: Q5 (cgload and gffutils); *F. graminearum*: Q5 (Chado);
* *F. graminearum*, *A. thaliana*, maize, human: Q6, chosen to span 0.1 to 4.9
  million rows, with organelles and without.

Q8 uses one genome of ours with funannotate gene models, eggNOG-mapper and
InterProScan output, named in its cgload config when the case study runs.

**Inputs.** The annotation is the exact file the survey read: its MD5 must equal
the one in `survey/results/survey.tsv`. `run.sh` takes it from the survey's
`downloads/` folder, else from the survey's Zenodo record, else from NCBI. A file
whose MD5 differs is used and flagged; Q7 is not computed for it. The assembly
FASTA, and for the fungi the GenBank file, come from NCBI Datasets (REST API v2),
with their MD5 recorded.

## 3. Software

* cgload at the tagged commit of this protocol. `run.sh` refuses to run from a
  checkout with uncommitted changes.
* gffutils 0.14 (PyPI).
* txdbmaker (Bioconductor), through conda-forge and bioconda; the version that
  resolves is recorded. (Before Bioconductor 3.19 the same function,
  `makeTxDbFromGFF`, lived in GenomicFeatures; `run_txdb.R` uses that one only
  if txdbmaker is missing, and records which one ran.)
* python-chado 2.3.9 (PyPI) with SQLAlchemy 1.4.54 and Biopython 1.81, which it
  requires, into the GMOD Chado 1.31 schema (`GMOD/Chado` at commit `ef04588`),
  with the Sequence Ontology (`SO-Ontologies` at commit `4340c14`).
* PostgreSQL 16, htslib (tabix), Python 3.12 (3.11 for python-chado), from
  conda-forge.

`environment.yml`, `environment-r.yml` and `environment-chado.yml` fix the
versions that can be fixed. Every run writes the exact installed builds to
`results/environment/` and the cgload commit to `results/versions.txt`.

Chado has several GFF3 loaders: the GMOD bulk loader (`gmod_bulk_load_gff3.pl`),
Tripal's importer, and python-chado's `load_gff`, which Galaxy Genome Annotation
uses. We use python-chado, the one that installs from PyPI without a Perl or
Drupal stack. The paper says which loader stands for Chado and does not
generalise its speed to the others.

## 4. Hardware

The lab server: Linux 6.5 x86_64, 80 CPUs, 125 GB RAM, 8.7 TB disk. `run.sh`
records the host, CPU count, memory and free disk at the start of each run.
Loads that are timed (Q1, Q6) run one at a time, with no other benchmark job
on the machine. Other loaders run several at a time (8 jobs for gffutils, 6
for Chado, 3 for txdbmaker, whose R processes are large), so their load times
are indicative only and the paper says so.

## 5. Definitions

**Reference count.** `truth.py` reads a GFF3 file with the standard library
only, sharing no code with cgload or any other loader. It follows the GFF3
specification: rows sharing an ID form one feature; a row without an ID is a
feature of its own; `Parent` names the parents; `##sequence-region` and the
`region` rows give each sequence's length, whether it is circular and which
genome it belongs to. For each feature it records the sequence, type, every
row's strand, every row's span, the parents, and a digest of all other
attributes (keys and values sorted, escapes decoded, rows merged).

**A loader's result** is reduced to the same record and compared feature by
feature (`compare.py`). For cgload, the record is read from what `cgload
export` writes (the round trip). For gffutils and Chado it is read from their
databases with SQL; for txdbmaker from what the TxDb holds (`run_txdb.R`).

A feature is

* **kept** when the loader holds it once and every compared field equals the file's;
* **changed** when the loader holds it but a field differs, or holds it as
  several records;
* **lost** when the loader does not hold it.

A change the loader reported is **declared** (cgload lists each envelope it
repairs in its manifest). A loss or change is **announced** when the loader's
own output names the feature, or names its type on a line saying it was
skipped, dropped, ignored, rejected or unknown; otherwise it is **silent**.

Fields compared: sequence, type, strands, spans, parents and attributes for
cgload and gffutils; the same without attributes for Chado (python-chado spreads
them over several tables).

txdbmaker holds transcripts, not features, so it is compared per transcript
(`compare_txdb.py`): sequence, strand, span, gene, exon spans (when the file
gives the transcript exon rows) and CDS spans. What counts as a gene, a
transcript, an exon or a CDS is what txdbmaker itself counts: its own type
lists, which `run_txdb.R` reads from the package, and its rule that a gene with
exons or CDS directly below it is also a transcript. A TxDb does not keep the
file's IDs, so each transcript is paired with the TxDb's by the name txdbmaker
gives it (its `transcript_id`, else its `Name`), and failing that by span.
Features of other types are counted as not modelled, not as lost. txdbmaker
names, by ID, each transcript it drops; such a loss is announced.

**Hard cases.** A feature can be in several.

| Case | Definition |
|---|---|
| multi_row | written on two or more rows sharing one ID |
| overlapping_rows | two of its rows overlap (programmed frameshifts) |
| trans_spliced | a row carries `exception=trans-splicing` |
| strand_unknown | strand `?` |
| crosses_origin | on a circular sequence, ends past its length |
| past_end | otherwise past the end of its sequence |
| inverted | a row with start > end |
| escapes_parent | not inside the span of a parent it names |
| exon_under_pseudogene | an exon whose parent is a pseudogene |
| multi_parent | names two or more parents |
| dangling_parent | names a parent the file does not contain |
| no_id | a row with no ID |
| multi_valued | an attribute with two or more values |
| rows_disagree | its rows carry different attributes |
| organelle | on a mitochondrial or plastid sequence |

## 6. Procedures

All run from `run.sh`, in this order, except that the kill trials of gffutils
and Chado follow their comparisons in step 7, and Chado, the slowest, runs
after the timing. Every step appends to its table in `results/` and skips
what is already there, so a stopped run resumes.

1. **Environments and PostgreSQL.** The three conda environments; a PostgreSQL 16
   server under the run folder, reached by a Unix socket only, with default
   settings. `run.sh` refuses to start while another copy runs in the same
   folder, and stops if PostgreSQL does not start.
2. **Inputs** (section 2) and the **reference count** of every file.
3. **Q7, agreement** (`agreement.py`): every row the survey flagged, with its
   survey class (`survey/classify_hits.py`'s rule) and what cgload does with
   the feature on it, and every feature cgload objects to that the survey did
   not flag. cgload's view is computed with its own code, as `load
   --repair-envelopes` runs it, without a database.
4. **Q1 and Q2, cgload at scale** (`run_cgload.py`), each genome alone, SQLite
   first, then PostgreSQL: `cgload init` on a new database; `cgload load`; if
   refused, the reasons are kept, the database is checked to be still empty, and
   the load is repeated with `--repair-envelopes`; `cgload verify` (with
   `--against` the GenBank file for the fungi); database and data-directory
   sizes; `cgload export` and the round-trip comparison. Wall time, CPU time and
   peak memory come from `measure.py`; the time limit per command is 6 hours.
5. **Q4, faults** (`faults.py`), on copies of the SQLite load. Damage after a
   load: F01 to F19. Damage inside the parser, by wrapping cgload's tokenizer:
   P01 to P08. Each fault's expected outcome is written in `faults.py` before
   any run, from what cgload claims to check: completeness (every feature,
   every row, every stated parent link, the stored assembly unchanged), not
   fidelity. A fault that keeps every count is expected to go unnoticed.
6. **Q5, interrupted loads** (`kills.py`): 10 trials per genome and tool (5 for
   Chado), each killing the whole process group with SIGKILL at a time drawn
   uniformly between 5% and 95% of a full load's duration (measured alone: the
   scale step's load for cgload, the control load below for the others), from
   a seeded generator. cgload: a new load (the database must be empty
   afterwards and the load must succeed when repeated) and a `--force` reload
   (the earlier load must survive and pass `verify`). gffutils and Chado: what the database holds
   afterwards, and, for gffutils, whether a database file is left and whether
   gffutils opens it. What is left is counted in each loader's own units (rows
   of its feature table and of its parent-link table) against the same counts
   from a load that was not killed: the repeated load for cgload, the earlier
   load for `--force`, and one uninterrupted control load per genome for
   gffutils and Chado. Each trial ends empty, partial, complete, or (for
   `--force`) intact or damaged.
7. **Q3, the other loaders** (`run_compare.py`): gffutils with its default
   merge strategy (`error`, to record what a user who changes nothing gets) and
   with `create_unique` (the documented way to load a file whose rows repeat an
   ID); txdbmaker's `makeTxDbFromGFF`, told which sequences are circular;
   python-chado's `load_gff` into a new Chado database. Time limit: 12 hours per
   genome, 24 for Chado. A load that is refused, crashes, runs out of memory or
   time is reported as such; whatever a failed Chado load left in its database
   is still compared.
8. **Q6, query timing** (`query_timing.py`), alone on the machine. Tasks: T1 a
   gene with every feature below it, by ID; T2 every feature overlapping a 1 Mb
   window; T3 feature counts by type; T4 every feature on one sequence, written
   as GFF3; T5, on the case study only, the features carrying a given Pfam, GO,
   KEGG, EC or InterPro accession. Tools: cgload on SQLite and on PostgreSQL,
   gffutils, and tabix on cgload's own bgzipped annotation (T2 and T4; T3 as a
   full scan). gffutils is asked through a few lines of its documented API
   (`gffutils_query.py`), as its own command-line tool cannot do T2 and repeats
   records in T1. 20 genes and 20 windows per genome, three sequences for T4,
   each asked 5 times in a row; targets drawn with a seeded generator. Two views:
   one new process per answer (what a user typing a command waits), and one
   process asking T1 to T3 in turn (what a script waits). The page cache
   is dropped before each first ask where the operating system allows it, and
   the table records whether it was.
9. **Q8, the case study** (`case_study.py`): load with eggNOG-mapper and
   InterProScan output, `verify`, what the functional files added and how many
   of their queries matched no feature, T5 for the most frequent accession of
   each analysis, the round trip, and a JBrowse 2 drawing of cgload's files
   (`jbrowse/render.sh`).
10. **Tables** (`summarise.py`), and the results archive to send back.

## 7. What is reported

7.1 Per genome and backend: outcome (loaded; loaded with repair; refused, with
the reasons), rows, features, load time, peak memory, database and
data-directory size, each `verify` tier, and the round trip (kept, changed,
declared, lost). 7.2 Per loader: genomes loaded, refused, failed; features
(transcripts, for txdbmaker) kept, changed, lost, silent or announced; parent
links kept; per hard case.
7.3 Per fault and genome: detected or not, by which check, against what was
expected. 7.4 Per loader and genome: the state after each kill. 7.5 Per genome,
tool, task and view: the median first ask and the median repeated ask, over
targets. 7.6 Per genome: rows the survey flagged, features cgload objects to,
and how they agree. 7.7 The case study as listed in step 9.

Every number in the paper comes from a file in `results/`; `summarise.py`
writes the tables and `results/tables/headline.json`.

## 8. What counts as a failure of cgload

* `verify` failing a tier after a load that succeeded;
* a round-trip difference that cgload did not declare;
* a crash (a traceback rather than a refusal with reasons);
* a refused load leaving anything in the database;
* an interrupted load leaving part of a genome in the database (a killed load
  must leave nothing, or, if it was killed after committing, the whole genome,
  passing `verify`), or a `--force` reload whose interruption damages the
  earlier load;
* a fault expected to be detected that is not.

A refusal with reasons is not a failure of correctness but is reported, with
the reasons, as a limit of coverage.

## 9. When cgload fails

The failure is reported. If it is a bug, it is fixed with a test, the fix is
tagged (`benchmark-run-2`, ...), and every measurement it can affect is
repeated on every genome. Both runs are reported, as the survey reported both
of its runs.

## 10. Analysis

`summarise.py` computes the tables in section 7 from the result tables, with
no manual steps. Medians, not means, for times. No statistical test is
planned: the comparisons are counts of what each tool did to the same
features.

## 11. Limits known in advance

* One annotation release per genome, all from RefSeq; other producers' files
  are represented only by the case study.
* Load times of the other loaders are measured under parallel load.
* python-chado stands for Chado; the GMOD bulk loader and Tripal are not run.
* txdbmaker models only genes, transcripts, exons and CDS.

## 12. Changes after the first run

None yet.
