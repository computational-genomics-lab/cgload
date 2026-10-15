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

## Open

- Memory: about 7 GB per million features plus 200 MB, from two loads made
  before the indexed annotation existed. Measure again in Phase 4.
- Three of the 30 RefSeq files (maize, *Arabidopsis*, soybean) each have one
  mitochondrial mRNA row whose start is after its end (survey `hits.tsv`).
  cgload refuses such a file unless `--repair-envelopes` is given (D-056), so
  the Phase 4 loads of those three need it, and the paper must say so.
- `paper/paper.md` still has the old maize numbers (36 problems in seven
  genes). The rest of the docs now use the census.
- Human and *Arabidopsis* have not yet loaded with the fixes; the next server
  run tests D-088 to D-091 on them.
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
