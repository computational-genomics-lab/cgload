# Work log: toward the cgload paper

What was done, what is open, and the evidence for each decision. Newest phase
last. Decisions that change behaviour are recorded in `docs/decisions.md`; this
file records the work.

## Phase 0: baseline of main at 695fbdd (13 October 2026)

| Check | Result |
|---|---|
| `ruff check .` (ruff 0.17.0) | 17 errors, so CI's lint job failed |
| `.github/workflows/ci.yml` | not valid YAML (line 87), so GitHub could not run CI at all |
| Tests, CI selection (`not slow and not requires_server`) | 437 passed, 1 failed, 16 expected failures |
| The failure | `test_the_bundled_example_loads`: a database committed in `examples/minimal` |
| Expected failures | the 15 audit tests (AUDIT-1 to 10) and D-059's accepted gap |
| Slow tier | needs real genomes that CI downloads; not run here |
| `requires_server` tier | no tests existed |
| README five-minute demo on a fresh clone | failed at `init` and `load` |
| `pip install cgload` (README) | fails: no such package on PyPI |
| `cgload init` on PostgreSQL 16 | failed (AUDIT-9) |
| `cgload stats` on PostgreSQL 16 | failed (GROUP BY) |

## Phase 1: defects fixed (branch `paper-v2`)

One commit per change, each with its tests. After the series (14 October 2026):
lint clean; 474 passed and 1 expected failure (D-059) in the CI selection; the
2 PostgreSQL tests pass against PostgreSQL 16.15; the workflow file parses.

| Commit subject | What it fixes |
|---|---|
| examples: remove the committed demo database | the README demo and the bundled-example test |
| lint: make CI's ruff job pass, and pin ruff | 9 unsorted test imports; frozen survey scripts excluded; ruff 0.17.0 pinned |
| export: keep multi-valued attributes multi-valued | AUDIT-3 |
| export: read attributes with a join | AUDIT-8 |
| schema: quote the end column in its CHECK constraint | AUDIT-9 (END) |
| schema: store coordinates and sequence lengths as 64-bit integers | AUDIT-9 (2^31) |
| stats: group by every selected column; add PostgreSQL tests | `stats` on PostgreSQL; first server tests |
| ci: repair the workflow file, and test the PostgreSQL code path | invalid YAML; PostgreSQL CI job |
| tokenizer: a ##FASTA directive ends the feature section | AUDIT-7 |
| verify: make --help match what --sequences-in-assembly does | AUDIT-10 (help) |
| storage: under --force, keep the old directory in place until the end | AUDIT-4 (crash safety) |
| loader: link and annotate anonymous rows | AUDIT-2; D-077 |
| verify: run the parent-edge tier, and run it inside GFF3 loads | AUDIT-1; D-078 |
| verify: check the stored assembly against its recorded checksum | AUDIT-4 (verify); D-079 |
| init: name a remedy that exists when the vocabulary hash differs | the error pointed at `cgload remap`, which does not exist |
| normalise: accept NCBI's mitochondrial layout | AUDIT-5; D-080 |
| loader: check rows against their sequence's length and their parent's sequence | AUDIT-6; D-081 |
| load: write the annotation as bgzipped, tabix-indexed GFF3 | AUDIT-10; D-082 |
| export: refuse a feature on a sequence the regions do not name | a backstop, so the indexed annotation cannot drop rows silently |
| docs: bring the README and known limitations up to date | install line, the checks table, the memory figure, new limitations |

Measurements behind the decisions:

- The tabix index: 3,000 random region queries through htslib (pysam 0.24.1)
  using cgload's index returned exactly the rows of a brute-force overlap scan,
  on 75,630 rows over three sequences of up to 400 Mb.
- Writing the indexed annotation: 1,000,000 rows in 24.8 s (38 MB compressed).
  The human RefSeq file has 4,932,571 rows (survey), so expect about two
  minutes on top of the load. Real rows carry longer attribute columns than
  the test rows; Phase 4 measures it on real genomes.
- The parent-edge tier on all six fixtures passes; with the AUGUSTUS anonymous
  rows unlinked, as before AUDIT-2, it fails with 175 of 273 edges missing.

## Checkpoint B: the lab server, 14 October 2026

The server: Linux 6.5.0-35-generic, x86_64, 80 CPUs, 125 GB RAM, 8.1 TB free.
Its system Python is 3.10.12; the script used Python 3.11.4 from the PATH.
conda (anaconda3), Docker (snap), `tabix` and `bgzip` are installed; R, a
PostgreSQL client and the NCBI `datasets` CLI are not. The first run stopped at
`pip install` because cgload needs Python 3.11; `checkpoint_b.sh` now finds or
builds one (conda-forge, else uv).

| Check | Result |
|---|---|
| lint | clean |
| tests, CI selection | 474 passed, 1 expected failure (D-059) |
| PostgreSQL tests | not run: no server configured |
| *F. graminearum* GCF_000240135.3 | downloaded through NCBI's REST API v2 |
| load | 78,715 features on 19 sequences (36,458,046 bp); 15.1 s; peak RSS 790 MB |
| verify, against the GenBank file | every tier passes: source count 78,715; parent edges 64,970; stored assembly; cross-parser on 5 types |
| tabix index, read by htslib | 3,800 region queries over 19 sequences and 102,635 rows: 0 mismatches |
| CI's `v2alpha` download URL | HTTP 200, so the nightly job's download works |

The load wrote one warning: `sequence_feature` (1 row) is not declared by the
RefSeq profile.

## Phase 2: the command surface, and defects found on the way

| Commit subject | What it does |
|---|---|
| functional: find annotatable features with a join, not an IN list | AUDIT-8's defect in the functional loader: an eggNOG file for a maize- or human-sized genome passed the driver's parameter limit |
| loader: under --force, delete every row the old assembly owns | `load --force` had never worked on an assembly with features: the foreign key refused it, and the load rolled back. D-083 |
| remove: take an assembly out, rows and data directory together | D-084. It deletes only a directory whose manifest names the assembly and its FASTA checksum |
| query: read-only questions by ID, region or functional term | D-085. Bound parameters only; recursive CTEs walk the parent chain |
| verify: recount an input only if it is still the file that was loaded | D-086. An edited input was reported as features "lost or invented during loading" |
| query: write TSV rows as gene models, parents first | gene, mRNA, then exons and CDS |
| query: complete gene models from any row, and regions across the origin | `--descendants` missed the exons of a trans-spliced mRNA when the region held only its second row; a region near position 1 missed a feature crossing the origin of a circular sequence |
| load: add --strict-vocabulary | D-087. The docs described the flag; no command had it |
| messages: point at config keys that exist, and test that every flag does | two messages named a `--dialect` flag and an unnamed "minimum match rate"; a test now fails on any message naming an undeclared flag |
| docs: rewrite the tutorial from a real run of every step | `--out`, `sequence.fa`, `gene-002` and a break-it step that could not fail |
| docs: drop an unverified claim and a leftover from the rename from the README | "not available elsewhere", and a sentence from the old name |

After the series: lint clean; 515 passed and 1 expected failure (D-059) in the
CI selection; 4 PostgreSQL tests pass against PostgreSQL 16.15.

`checkpoint_b2.sh` loads four real genomes with this code: *F. graminearum*,
*Arabidopsis*, maize and human. For maize and *Arabidopsis* it first runs
`envelope_census.py`, cgload's own count of children outside their parents,
which is what the maize numbers in the docs need. Its timings are a check that
the code works, not benchmark results.

## Checkpoint B2: four real genomes, 15 October 2026

The full log came back. The run used the Phase 1 bundle ("commits on top of
main: 21"), in which `query` and `remove` are stubs ("not implemented in
0.1.0"); that explains every query and remove failure. The script now refuses a
bundle without Phase 2. The loads themselves do not depend on Phase 2.

| Genome | Result |
|---|---|
| *F. graminearum* GCF_000240135.3 | loaded: 78,715 features; 14.8 s; peak RSS 790 MB; every verify tier passes, cross-parser included |
| *Arabidopsis* GCF_000001735.4 | refused, with and without `--repair-envelopes` |
| maize GCF_902167145.1 | refused without the flag (a row written backwards, line 990,621); loaded with it: 696,928 features in 21 types, 138.1 s, peak RSS 5,397 MB; verify passes (source count 696,928; 644,997 parent edges; 2 parentless features, D-069) in 5.5 s; database 697 MB, data directory 2.1 GB |
| human GCF_000001405.40 | refused at line 1,838,264 after 131.4 s, peak RSS 14,997 MB |

`envelope_census.py`, cgload's own count:

- maize: 11 envelope problems on the mitochondrion, in two mRNAs: 10 children
  outside their mRNA (6 under rna-ZeamMp017, 4 under rna-ZeamMp186) and
  rna-ZeamMp017 written backwards. Repair recomputes both mRNAs and leaves
  none; the larger change widens rna-ZeamMp186 by 281,540 bp. The docs said 36
  problems in seven genes; corrected, except in `paper/paper.md`.
- *Arabidopsis*: 13 problems on the mitochondrion in three genes (12 children
  outside, 1 row backwards). Repair clears 11. Two mRNAs written past the origin
  stay outside genes written on several rows, and refused the load. D-091
  stores them as written.

Human: `cds-NP_001165908.1` is written as two rows sharing base 94,664,513, a
programmed frameshift. D-089 permits it where the rows say `ribosomal
slippage`. The load stopped there, so the rest of the file is unchecked; D-090
makes the next run list every problem in a pass at once, and D-088 has already
removed one likely problem further on.

## Since Phase 2

| Commit subject | What it does |
|---|---|
| vocabulary: map the rest of the INSDC ncRNA classes to ncRNA | D-088 |
| loader: after --repair-envelopes, stop suggesting --repair-envelopes | the refusal named the flag the user had already given |
| release: write down how to release, and stop a mismatched tag | `RELEASING.md`, `CHANGELOG.md`, a tag check in `release.yml` |
| refseq: let the CDS rows of a programmed frameshift share a base | D-089, the human refusal |
| normalise: refuse with every problem a pass finds, not only the first | D-090 |
| loader: store an envelope disagreement that exists only across the origin | D-091, the *Arabidopsis* refusal |
| docs: the maize numbers from cgload's own census, and the server run | envelope module, known limitations, README |

Target journal, decided 15 October: *Briefings in Bioinformatics*, as a
Problem Solving Protocol; *Bioinformatics* if it is refused there. Both are on
the ONOS APC list (checked by Aditya Upadhyay on the ONOS site, 15 October).

Pushing from this session is still refused (HTTP 403, 15 October): the org has
not installed the Claude GitHub App. Pushes go from the lab server.

How the work lands, decided 15 October: no other branches and no pull
requests. When the fixes are done and the server run passes, the commits go
straight onto `main` in one push. Then the paper.

## Checkpoint B2 rerun, 15 October 2026

Run on the server with the bundle of `main` (43 commits). Aditya reported that
every step passed and pushed `main`; the log was not sent back, so the numbers
of that run are not recorded here. The commits were then re-authored (author
Aditya Upadhyay, committer Claude, same trees and dates; GitHub's `main` ends
at `119dbad`), because they had been made under the account owner's name.

## Phase 3: the benchmark protocol and its scripts, 16 October 2026

`benchmark/PROTOCOL.md` fixes the questions, genomes, software, hardware,
definitions, procedures, failure criteria and analysis before any run.
`benchmark/run.sh` runs all of it and resumes after a stop. The scripts:

| Script | Does |
|---|---|
| `truth.py`, `gffmodel.py` | the reference count: what a GFF3 file holds, standard library only, no loader code |
| `compare.py`, `compare_txdb.py` | a loader's store against the reference count: kept, changed, lost; declared, announced, silent |
| `run_cgload.py` | Q1 and Q2: load, refuse-then-repair, verify, sizes, export round trip |
| `run_compare.py`, `run_gffutils.py`, `run_txdb.R`, `run_chado.py` | Q3: gffutils (default and create_unique), txdbmaker, python-chado into Chado 1.31 |
| `faults.py` | Q4: 19 database and data-directory faults, 8 parser faults, each with its expected outcome |
| `kills.py` | Q5: loads killed at seeded random moments |
| `query_timing.py`, `gffutils_query.py` | Q6: four tasks, four tools, command and library views |
| `agreement.py` | Q7: cgload's checks against the survey's rows |
| `case_study.py`, `jbrowse/render.sh` | Q8: a lab genome with eggNOG and InterProScan; a JBrowse 2 drawing of cgload's files |
| `summarise.py` | the paper's tables, `headline.json` |

Writing and testing these found five problems in cgload, each fixed with tests
before the protocol was tagged:

| Commit subject | What it fixes |
|---|---|
| verify: count the file's Parent statements independently (D-092) | the parent-edge tier took the file's parents from the parser's own record; a dropped `Parent` passed |
| load: do not call violations unrepairable when no repair was asked for | the plain-load refusal said "could not be repaired" |
| verify: count rows as well as features (D-093) | a dropped row of a multi-row feature passed every tier |
| query: read only the rows a gene or region needs (D-094) | a gene took 3.7 s on PostgreSQL and a 1 Mb window 1.8 s on 120,018 features; now 0.08 s and 0.30 s |
| load: replace a directory an interrupted load left behind (D-095) | a load killed between putting its directory in place and committing (1 kill trial in 10) left a directory that made the repeated load fail |

Measured here, on synthetic genomes made by `benchmark/tests/make_synthetic.py`
(they carry every hard case of the protocol) and on the real-data fixtures,
not on real genomes:

- All 27 faults behave as `faults.py` expects them to.
- A killed cgload load left either nothing, after which the repeated load
  succeeded, or, when the kill came after the commit, the whole genome,
  passing `verify`. Before D-095, one kill in ten landed between the rename
  and the commit and the repeated load was refused. A killed `--force` reload left the earlier load intact and
  passing `verify` in every trial. A killed load leaves its staging copy of the
  FASTA under `<data-dir>/.staging/`; now in known limitations.
- A killed gffutils load left no file, an empty file, or every record; once,
  every record in a file gffutils could not open. python-chado commits once,
  at the end of `load_gff`, so a killed Chado load left nothing in every trial.
- python-chado 2.3.9 loads about 88 rows a second (160,019 rows in 30 min),
  so human (4.9 million rows) should take about 15 hours. It needs SQLAlchemy
  1.4 and Biopython 1.81, hence its own environment. Chado 1.31's schema loads
  on PostgreSQL 16 with 9 errors, all in interval-overlap helper functions that
  use the `box @ box` operator PostgreSQL no longer has; GMOD's Docker image
  gets the same.
- txdbmaker ran here as GenomicFeatures 1.54.3 under R 4.3.3 from Ubuntu
  (`makeTxDbFromGFF` lived there until Bioconductor 3.19; Ubuntu does not
  package txdbmaker). On both synthetic genomes and six real fixtures the
  comparison paired every transcript. It drops a trans-spliced transcript
  whose exons it cannot order and names it; keeps only the first CDS in an
  exon that has two (a frameshift) and names it; and drops exons that have no
  Parent without a word (the maize chloroplast fixture).
- gffutils 0.14's own command, `gffutils-cli`, cannot do the timing tasks: its
  `region` raises NotImplementedError and its `children` prints grandchildren
  twice. The timing calls gffutils' documented API instead.
- A cgload command takes about 0.5 s on a small genome, 0.38 s of it importing
  Python packages (SQLAlchemy alone 0.25 s); the gffutils script about 0.15 s.
  Called from Python, gffutils answered faster than cgload on these small
  genomes (a gene in 1 ms against 4 to 5 ms; a 1 Mb window in 88 ms against
  175 ms on SQLite). Real genomes are Phase 7.
- JBrowse 2's command-line export (`@jbrowse/img` 3.7.0) draws cgload's
  bgzipped annotation through cgload's own .tbi index; two of its dependencies
  needed pinning (`jbrowse/package.json`).
- `run.sh` ran start to finish on the two synthetic genomes, every step
  passing (the case study was not set).

Expected on the lab server: about 1.5 to 2 days, most of it Chado (27.1
million rows over the 30 genomes; human alone about 15 hours), and about
200 GB of disk.

## Phase 4, first attempt: void, 21 October 2026

The run at `benchmark-protocol-v1` came back with four steps failed. Two
copies of `run.sh` had been started six minutes apart in one folder, so every
table, database and download was written twice at once; nothing from it is a
result (PROTOCOL.md, section 12). `diagnose.sh` collected the step logs and
the input lines behind each refusal. What they showed, each fixed with tests:

| Found | Fix |
|---|---|
| a multi-row feature exported with its first row's attributes (3,249 features in human, mostly `cDNA_match`), and trans-spliced rows exported on `?` | D-096 |
| *Xenopus*: 199 `precursor_RNA` with their exons and miRNAs refused; rice, grape: `pseudogenic_tRNA` exons refused | D-097 |
| *Drosophila*: mod(mdg4), five gene rows alternating strands, two overlapping by 149 bases | D-098 |
| *S. pombe*: 18 CDS reading a base twice, `exception=low-quality sequence region` | D-099 |
| *Ustilago*: a tRNA gene's `gbkey=exon` feature on two rows sharing `id-UMAG_16066` | D-100 |
| *Arabidopsis*, PostgreSQL only: an attribute value too long for the (key, value) b-tree index; a traceback | D-101 |
| `verify --against` failing on five fungi over NCBI's two names for one feature | D-102 |
| `verify --against`: 13 to 28 more mRNAs in the GFF3 on four fungi with a mitochondrion; in yeast 19 of 25 on NC_001224.1, whose GenBank record states CDS without mRNA | D-103 |
| `run.sh`: no lock; carried on without PostgreSQL; Chado kills on the wrong genome; faults and timing crashing on a failed load; assembly ID assumed to be 1 | benchmark commit |

Real results it showed in passing, to be measured again: txdbmaker cannot read
maize, *Arabidopsis* or soybean (one row written backwards), and stops with
"exons linked to transcripts not found" on rice, tomato and grape, whose
trans-spliced plastid genes it drops; it drops 54 human exons whose parents it
does not import, and holds no transcript for 10,792 *Xenopus* exon and CDS
rows. gffutils' default merge strategy refuses all 30 files.

## Phase 4, preflight: 22 October 2026

The 14 genomes the void run had trouble with, at 573d173, scale steps only
(1 h 16 min on the lab server). All 14 loaded on SQLite and PostgreSQL, three
with `--repair-envelopes`; every source-count, parent-edge and stored-assembly
tier passed, and the two databases gave the same verify report for every
genome. Nothing was lost in the round trip. What it showed:

| Found | Fix |
|---|---|
| 1,938 features changed silently in the round trip (1,792 `cDNA_match`, 141 CDS, 5 genes): a minus-strand feature's lowest row lost its own attributes | D-104 |
| the comparison took the union of a feature's rows' attributes, so a value moved to the wrong row could pass | per-row digest (`gffmodel.feature_digest`) |
| yeast's six extra mRNAs, and *S. pombe*'s 32: pseudogenes GenBank states as a CDS alone | D-105 |
| *S. pombe*'s 238 `LONG_TERMINAL_REPEAT` rows | D-106 |
| the survey read a medaka file of 292,462 rows on 5 sequences; NCBI's has 1,399,050 on 25 (most likely a download cut short) | survey rerun on the complete file |
| the summary's silent-change count for cgload was always 0 | scale table columns |

*Neurospora*, *Aspergillus* and *Cryptococcus*: their extra mRNAs were all
mitochondrial (28, 20, 13), and their GenBank comparison passed with D-103.

## Phase 4, second preflight: 23 October 2026

The same 14 genomes at cb651d3, SQLite only, with the agreement and fault steps
(42 min). The round trip changed nothing but the 7 declared repairs, so D-104
accounted for all 1,938 changes. `verify --against` passed on all eight fungi.
The faults gave the predicted outcome in 54 of 54. The agreement step compared
13 genomes; medaka was left out because its file was not the survey's.

The survey's copy of medaka, still on the lab server, is 292,477 lines and ends
mid-line: the download was cut short. Surveyed again on the complete file, its
row changed only in the counts that describe the file; no defect before or
after (`survey/results/CORRECTIONS.md`).

The warnings, now kept, showed every RefSeq load but *Aspergillus*'s warning
about types the profile did not declare (D-107).

## Open

- Memory: about 7 GB per million features plus 200 MB, from two loads made
  before the indexed annotation existed. Measure again in Phase 4.
- Three of the 30 RefSeq files (maize, *Arabidopsis*, soybean) each have one
  mitochondrial mRNA row whose start is after its end (survey `hits.tsv`).
  cgload refuses such a file unless `--repair-envelopes` is given (D-056), so
  the Phase 4 loads of those three need it, and the paper must say so.
- `paper/paper.md` still has the old maize numbers (36 problems in seven
  genes). The rest of the docs now use the census.
- The case-study genome (Q8) is not named yet.
- `run_txdb.R` has run with GenomicFeatures 1.54.3 only; the server run uses
  txdbmaker from conda.
- Novelty: check Ensembl's datachecks and InterMine's acceptance tests before
  the paper claims that no other loader recounts against the source.
- `paper/paper.md` says the name cgload is taken on PyPI. PyPI has no project
  of that name (its JSON API returned 404 on 14 October 2026), and the comment
  points at `docs/NAME.md`, which does not exist. Fix when the paper is
  rewritten.
- MySQL has no test. No MySQL server could be started in this session.
- The system Python of Ubuntu 22.04 is 3.10, which cgload does not support.
  The README and tutorial say how to get 3.11 or newer; a bioconda package
  would make it one command (release phase).
- Changing the vocabulary hash and the schema means databases made before this
  branch are refused, with a message saying why. Reload into a new database.
