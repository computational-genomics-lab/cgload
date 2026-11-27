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
attributes (keys and values sorted, escapes decoded, rows merged). White space
at either end of a tag=value pair is not part of it.

**A loader's result** is reduced to the same record and compared feature by
feature (`compare.py`). For cgload, the record is read from what `cgload
export` writes (the round trip). For gffutils and Chado it is read from their
databases with SQL; for txdbmaker from what the TxDb holds (`run_txdb.R`).
gffutils keeps a space that ends a value, so its values are trimmed as
`truth.py` trims the file's (`run_gffutils.py`).

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
them over several tables). A feature whose rows carry different attributes is
compared row by row: each row's span with that row's own attributes, so a value
given to the wrong row is a change even when no value is lost.

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
   `--against` the GenBank file for the fungi; each report is kept whole in
   `results/verify/`); database and data-directory
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

**First run, 21 October 2026, at `benchmark-protocol-v1`: void.** Two copies of
`run.sh` were started in one folder, six minutes apart. They downloaded,
loaded and wrote the result tables over each other, and the second could not
start PostgreSQL, so it marked Chado "not run" for every genome and Chado
never ran. No number from it is reported as a result. What it showed is
reported, because it is what the benchmark is for, and was fixed with tests
before the second run (tag `benchmark-run-2`):

* cgload refused six of the 30 RefSeq files over valid content: *Xenopus*
  (`precursor_RNA`, D-097), rice and grape (`pseudogenic_tRNA`, D-097),
  *Drosophila* (a trans-spliced gene's rows overlapping across strands,
  D-098), *S. pombe* (CDS rows reading a base twice over a low-quality
  sequence region, D-099) and *Ustilago* (an exon feature written on two rows,
  D-100).
* cgload failed with a traceback on *Arabidopsis* on PostgreSQL only: an
  attribute value too long for a b-tree index (D-101).
* The export round trip changed features silently: a multi-row feature kept
  only its first row's attributes, and a trans-spliced gene's rows came back
  on strand `?` (D-096).
* `verify --against` failed on five fungi over names NCBI writes differently
  in its GFF3 and GenBank files (D-102).
* A killed load could leave a directory that blocked the repeated load
  (D-095; found on the test genomes before the run).
* The scripts: `run.sh` now refuses a second copy in the same folder and
  stops if PostgreSQL does not start; Chado's kill trials run on
  *F. graminearum* as section 2 says (they had taken *A. thaliana*);
  `faults.py` and `query_timing.py` no longer crash on a genome whose load
  failed; every runner reads the assembly ID from the database.

**Preflight, 22 October 2026, at 573d173.** The 14 genomes the first run had
problems with, loaded, verified and exported on both databases, nothing else.
Every one loaded, three with `--repair-envelopes` (*A. thaliana*, rice, *S.
pombe*, each repair declared), and every load passed the source-count,
parent-edge and stored-assembly tiers. It showed, and these were fixed with
tests before the second run:

* The round trip changed 1,938 features silently, every one written on rows
  with different attributes: a minus-strand feature's lowest row came back
  with its first row's attributes (D-104). The comparison had missed part of
  this, since it took the union of a feature's rows; it now compares such a
  feature row by row (section 5).
* `verify --against` still failed on yeast and *S. pombe*: NCBI's GFF3 gives
  an mRNA to pseudogenes that GenBank states as a CDS alone (D-105), and
  writes *S. pombe*'s LTRs in capitals (D-106).
* The survey read a different medaka file from NCBI's: 292,462 rows on 5
  sequences, where NCBI's file, fetched twice on different days with the same
  MD5, has 1,399,050 rows on 25. The other 29 files are the survey's exactly.
  The survey's copy, still on the lab server, ends in the middle of a line: its
  download had been cut short. The survey was run again on the complete file
  and its row replaced; its findings did not change
  (`survey/results/CORRECTIONS.md`).
* The summary reported 0 silent changes for cgload, because the scale table
  had no such column; it now has the announced and silent counts. The load's
  warnings and the whole verify report are kept in `results/`.

**Second preflight, 23 October 2026, at cb651d3.** The same 14 genomes on
SQLite, with the agreement and fault steps. No feature changed in the round
trip except the seven repairs cgload declared, `verify --against` passed on all
eight fungi, and the 54 faults (27 on each of *F. graminearum* and *A.
thaliana*) gave the outcomes `faults.py` predicts. It showed one more thing,
fixed before the second run: every RefSeq load but one warned about feature
types the RefSeq profile did not declare (D-107).

Every measurement is repeated on every genome in the second run.

**Second run, 23 and 24 October 2026, at `benchmark-run-2` (8c999e2).** Every
step finished, in 19.5 hours. Two steps were wrong in the scripts, not in the
tools they measure, and one measurement found a cgload defect. Each was fixed
with tests, and the steps they touch were run again at `benchmark-run-2.1`:

* Chado's database setup. Six genomes were set up at once, and loading the
  Sequence Ontology with one savepoint per term ran PostgreSQL out of lock
  space ("out of shared memory"). The first six genomes, all vertebrates, got
  no ontology, and python-chado stopped at once. Reproduced with eight setups
  at once. The ontology now loads without savepoints, and a setup that fails is
  recorded as such, never as a Chado result. Chado is measured again on all 30.
* txdbmaker's comparison paired two unnamed transcripts that share a span on
  opposite strands crosswise, as it did human MIR103A1 and its antisense
  MIR103B1 precursors, and reported strand changes txdbmaker never made (14 in
  human, 10 in rat). Pairing now tries the same strand first. txdbmaker is
  measured again on all 30.
* Q6: looking up one human gene took 12 s on PostgreSQL (D-108). The timing
  step is measured again. Its library count now counts features as `cgload
  stats` does; it had counted rows, an easier question than the command
  answers.

`forget.py` set aside what the repeated steps replace, in
`results/superseded/`. Every other step stands as measured at
`benchmark-run-2`: the scale, round-trip, fault, interrupted-load, gffutils
and agreement steps do not run the code these fixes changed.

**The repeat, 4 to 6 November 2026, at `benchmark-run-2.1` (9b97178).**
txdbmaker and the timing step finished as expected, and Chado on 29 genomes.
The human Chado load reached its 24-hour limit and left nothing, as python-chado
commits only at the end. Its comparison then ran for 26 hours and was stopped:
`compare.py` searched the load's 15 MB log once for each of 3.2 million lost
features, to tell an announced loss from a silent one. It now answers the same
question from an index of the log, built once (tested to give exactly the old
answer), and `recompare.py` filled in that one comparison without loading
again. The run's other results are as the repeat measured them.

**After the repeat, 6 November 2026.** gffutils' 133 attribute changes in
*S. pombe* were an error in the comparison. NCBI ends some *S. pombe* product
names with a space ("D-glutamate cyclase "). gffutils keeps the space, and so
does cgload, while `truth.py` drops white space at the end of a tag=value pair.
cgload's result is read from its export by `truth.py`, so the space was dropped
on both sides of its comparison. gffutils' result is read from its database, where
the space stayed. `run_gffutils.py` now trims gffutils' values as `truth.py`
trims the file's. It also stops decoding escapes a second time, since gffutils
decodes them itself. Only the gffutils comparison runs this code, and it is
measured again on all 30 genomes at `benchmark-run-2.2`, so that every gffutils
result comes from the same scripts.
Result (10 November): gffutils changes 1,072,125 features, each one a
multi-row feature that `create_unique` stores as several records. It changes
no attribute in any genome, *S. pombe* included. Every other table is the same
as before, byte for byte.

**The results the paper reports.** They come from three tags:
`benchmark-run-2` (`8c999e2`) for cgload, the faults, the interrupted loads,
the agreement with the survey and gffutils with its default settings;
`benchmark-run-2.1` (`9b97178`) for Chado, txdbmaker and the query timing, with
the human Chado comparison filled in at `44c2ca5`; and `benchmark-run-2.2`
(`56f6861`) for gffutils with `create_unique`. Their tables are in
`paper/data/benchmark-run-2.2/`. What is left out there (the PostgreSQL
server log, the reference databases, which `truth.py` rebuilds, and the
superseded results) is kept on the lab server.

**Q8 dropped, 10 November 2026.** No genome of ours with funannotate gene
models, eggNOG-mapper and InterProScan output was available, so the case study
(step 9, section 7.7) was not run and is not reported. `run.sh` still has its
step, which does nothing while `CASE_STUDY` is unset. With it goes the only
file from outside RefSeq (section 11): every result here is on RefSeq GFF3.

**Reporting, 11 November 2026.** No measurement changed. Three things changed
in how the paper reports them, after an outside review of its draft. Table 1
counts each loader's features over the genomes it loaded, and counts a genome
it did not load as not loaded, no longer as so many lost features. Figure 4
shows T4 for the largest sequence of each genome, with all three in a table:
the median over the three hid the cost in human, where two of the three are
small scaffolds. And the paper says how the census's rows map to the reference
reader's features in each hard case.

**Columns 2, 6 and 8, 19 November 2026.** Not in the protocol. A review of the
paper's draft by another language model (ChatGPT) pointed out that the
reference reader leaves out each row's source, score and phase, so the round
trip (section 6) never compared them. `benchmark/columns.py` compares them row
by row between an input file and cgload's export, and
`benchmark/run_columns.py` repeats cgload's part of the round trip on SQLite
to run it. On the test fixtures it found two faults in the export of cgload
0.2.1 (decisions D-109 and D-110), fixed in 0.2.2. It also counts, in each
input file, the values that 0.2.1 would have written differently. Its results
are `columns.tsv`; the other loaders were not checked in this way.

It ran on 20 November 2026 with cgload 0.2.2 (tag `v0.2.2`), results in
`paper/data/columns-v0.2.2/`. In 28,211,625 rows no value of the three columns
differed. The files hold 52 scores with more than six significant digits (49
in human, 3 in mouse), which 0.2.1's export rounded, and no feature with rows
from different programs.

**Clarification, 23 November 2026.** The "outside review of its draft" under
11 November was a review by Claude, the language model that wrote these
scripts, as `paper/llm_use.md` records. Nothing measured changed.

## 13. The portal (added 25 November 2026)

Not part of `benchmark-protocol-v1`, and written after its runs. cgload 0.3.0
adds a web portal (`cgload serve`, `docs/portal.md`), and this section says
how it is measured on the same 30 genomes. It was written, with
`run_portal.py`, before that measurement was made. Nothing in it changes a
result of the earlier sections.

**Questions.**

* P1. Do all 30 genomes load into one database, as a lab's portal would hold
  them, and does each genome's search index pass its count check (D-112)?
* P2. How long does a user wait for the things a user asks for: a search by
  a gene's name, the gene's page, its region, its protein, and the byte ranges
  of the files the genome browser reads? On SQLite with all 30 genomes, and on
  PostgreSQL for human.
* P3. Does a search by a gene's own name list that gene, and on which page?
* P4. How do the times change with 1, 8 and 32 users at once?
* P5. Are the proteins the portal shows right? A CDS translated as the portal
  translates it should have a stop codon only at its end, unless the file
  says why not.
* P6. Does any request make the server fail, return a script it was sent, or
  lack a security header? Does serving change the database?

**Procedure.** `run_portal.py`, whose docstring lists its eight steps, on the
server of section 4, with the inputs `run.sh` downloaded. The portal is
`cgload serve` with 8 threads, one process; the client runs on the same
machine and reaches it over the loopback interface. For each genome, 200
genes with a `Name` are chosen, seeded by the accession; a gene's requests are
the search, its page, its region (the gene with 5 kb on each side, at most
100 kb), the protein of its first CDS, and a 64 KiB range of each of
`annotation.gff3.gz` and `assembly.fa` at a seeded offset. The load test
draws from the same requests, each client with its own seed, for 60 s at each
level. Up to 2,000 CDS per genome, seeded by the accession, are translated. The
security probe sends the injection strings in the script and 3,000 seeded
random requests to the server holding all 30 genomes.

**What counts as a failure of the portal.**

* A genome that loaded in the benchmark does not load here, or its index
  check fails.
* A response with status 500 or above, or no response.
* A search by a gene's name that does not list the gene, where no more than
  1,000 features (20 pages) match that name.
* A protein with a stop codon inside, from a CDS with no `transl_except`,
  `exception` or `ribosomal_slippage`, no `start_range`, `end_range` or
  `partial`, and not part of a pseudogene. Each is listed by ID and looked
  into, and section 9 applies.
* A response without the security headers, a script sent in a request and
  returned in a page, or an answer other than 405 to POST, PUT or DELETE.
* A database file whose SHA-256 after the serving steps differs from before.

**What is reported.** Per genome and in total: the load outcome and the index
check (P1); per kind of request, the median, 95th percentile and maximum time
and the statuses (P2); the share of genes found, and on which page (P3);
requests per second and the 50th, 95th and 99th percentiles at each level
(P4); the translation classes (P5); the probe's statuses and failures and the
database's checksums (P6). The results go to `paper/data/portal-v0.3.0/`.

**Limits known in advance.**

* Server and client share the machine, so the times leave out the network
  and the browser's own drawing of igv.js.
* One server process. Its threads share Python's interpreter lock, so
  requests per second stop rising once that is busy; more needs more
  processes behind the web server (`docs/portal.md`).
* Only genes with a `Name` are sampled.
* The probe is the authors' own, not a test by an outside party.
* No other portal is timed: none serves these 30 files without converting
  them first. Tripal reads Chado, whose loading the earlier sections measure.

## 14. BLAST in the portal (added 26 November 2026)

Not part of `benchmark-protocol-v1`, and written, with `run_blast.py`, before
the measurement it describes. cgload 0.4.0 adds BLAST to the portal
(`cgload blastdb`, `docs/portal.md`, D-116); this section measures it on the
database of section 13, which holds all 30 genomes. Nothing in it changes a
result of the earlier sections.

**Questions.**

* B1. Do BLAST databases build for all 30 genomes, and does each pass its
  count check: every sequence of the assembly in the genome database, one
  protein for every CDS feature in the protein database, and BLAST+'s own
  count of what it stored equal to both? How long does a build take, and how
  much disk do the databases use?
* B2. Does a query find itself? A protein of the portal, searched with blastp
  against every genome, should find its own CDS with the top score at 100%
  identity; a stretch of a genome, searched with blastn (megablast) against
  every genome, its own place, at 100% over its whole length; a protein,
  searched with tblastn against its own genome, an HSP at 100% inside its own
  CDS's span; a CDS, searched with blastx against its own genome's proteins,
  its own protein with the top score. This tests the databases, the
  coordinates and the links, not BLAST's sensitivity, which is NCBI's.
* B3. How long does a visitor wait for a search, by program, against one
  genome and against all 30?
* B4. With 1, 2, 4 and 8 visitors searching at once, how many searches are
  done and how many turned away (503), and do the pages still answer, and how
  fast, against the pages with no search running?
* B5. Does any request to the BLAST URLs make the server fail, return a script
  it was sent, or lack a security header? Does BLAST change the database or
  the BLAST databases?

**Procedure.** `run_blast.py`, on the server of section 4, after
`run_portal.py`, with its work directory. Step 1 runs `cgload blastdb` for
each genome in turn. The portal is then `cgload serve` with 8 threads and its
default BLAST settings (2 searches at once, 4 threads each, 120 s limit), and
the client reaches it over the loopback interface, through `POST
/api/v1/blast`, one search at a time until step 3. Per genome, seeded by the
accession: 20 CDS for blastp, their proteins as the portal serves them
(`protein.fa`), against every genome with 100 hits; 10 stretches of 500 bp
with no N, from sequences of at least 2 kb chosen in proportion to their
length, for blastn against every genome with 100 hits; 3 of the 20 proteins,
those with no stop or X inside first, for tblastn against their own genome,
and their CDS (`cds.fa`) for blastx against their own genome's proteins.
Three proteins, from human, *D. melanogaster* and *S. cerevisiae*, are
searched with tblastn against every genome, the slowest search a visitor can
ask for. The load test runs 60 s with no searching client and then with 1, 2,
4 and 8, each sending blastp searches of the sampled proteins against every
genome, waiting 1 s after a 503, while one more client opens the sampled
CDS's pages. The security probe sends the injection strings of section 13
and strings meant for a shell or for BLAST's own options in each field of
both BLAST URLs, as a form and as JSON, 600 seeded random searches, and
POST, PUT and DELETE where the portal answers only GET.

**What counts as a failure.**

* A genome whose BLAST database does not build, or fails its count check.
* A protein or stretch, with no stop or X inside, that does not find itself
  as B2 states, where nothing else ties with its top score. A query that ties
  (an identical copy elsewhere, a repeat) is counted, not failed. Each failure
  is listed and looked into, and section 9 applies.
* A response with status 500 or above, other than 503 with Retry-After from a
  full portal or a search stopped at its time limit, or no response.
* A script sent in a request and returned in a page, a response without the
  security headers, or an answer other than 405 to a method a URL does not
  take.
* A database file, or a BLAST manifest, whose SHA-256 after the searches
  differs from before.

**What is reported.** Per genome and in total, the build's checks, time,
peak memory and size (B1); per program and scope, the queries that found
themselves, those that tied, those that failed, and those stopped at the
time limit (B2); the median, 95th percentile and maximum time of each
(B3); per load level, searches done and turned away, and the pages' median
and 95th percentile times (B4); the probe's statuses and failures, and the
checksums (B5). The results go to `paper/data/blast-v0.4.0/`.

**Limits known in advance.**

* Server and client share the machine, so the times leave out the network.
* A self-hit shows that the right sequence is in the right database under
  the right name and coordinates; it says nothing about distant homologues.
* The proteins are cgload's translations of the CDS rows (D-113): where the
  file declares a translation exception, the protein in the database is not
  the curated one.
* The E-value of a search against every genome counts all 30 databases;
  against one genome, only that one.
* A search stopped at 120 s is counted as stopped, not timed past it; the
  limit is an option (`--blast-timeout`), and so are the number of searches
  at once and their threads.
* The probe is the authors' own, not a test by an outside party.
* No other BLAST server is timed.

**Run 1 and run 2 (added 27 November 2026).** Run 1, with release v0.4.0 on
26 November, met two of the failures named above, and section 9 applies to
both.

* The BLAST databases of *D. rerio* and *G. gallus* did not build: makeblastdb
  reported that the database it had just written was not where it had
  written it. The messages show the build's working folder renamed and, in
  the second case, removed while makeblastdb ran, which only another build of
  the same genome does, and two builds of one assembly started together fail
  in the same way. cgload had no guard against that; the build that failed
  then removed the databases the other had built, so the portal served 28
  genomes. Release v0.4.1 holds a lock per assembly for the length of a build
  and refuses a second build at once, changing nothing (D-117), with a test.
* Searches whose own sequence was among the hits at 100% identity over its
  whole length, but below another hit's score, counted as not found when
  nothing tied with it; so did tblastn searches with no HSP at 100% identity in
  their CDS's span. What counts as found is unchanged.

Run 2 repeats every measurement of this section on every genome with release
v0.4.1, from new databases, and both runs are reported. It also records, for
each search that did not find itself, how near it came: its own hit's score,
and for tblastn the best HSP in the CDS's span at any identity. Run 1's
results are `paper/data/blast-v0.4.0/`, run 2's `paper/data/blast-v0.4.1/`.
