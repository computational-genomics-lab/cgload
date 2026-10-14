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

## Open, for later phases

- `cgload query` and `cgload remove` are still stubs (Phase 2).
- Memory: the README and known-limitations said 4.5 GB per million features,
  which `docs/scalability.md` had already corrected to about 7 GB plus 200 MB.
  They now agree. All scale numbers predate the indexed annotation, which holds
  one sequence's rows, sorted, while it writes. Measure again in Phase 4.
- Three of the 30 RefSeq files (maize, *Arabidopsis*, soybean) each have one
  mitochondrial mRNA row whose start is after its end (survey `hits.tsv`).
  cgload refuses such a file unless `--repair-envelopes` is given (D-056), so
  the Phase 4 loads of those three need it, and the paper must say so.
- The maize numbers in `src/cgload/envelope.py`, `docs/known-limitations.md`
  and `paper/paper.md` (36 problems, 34 children outside their parents, seven
  genes, 281 kb) came from checkers that read only the last row of a multi-row
  feature. The survey's corrected count is 16 children outside their parents in
  two mitochondrial mRNAs, plus one row whose start is after its end. Confirm
  with a maize load, then correct (Phase 8).
- The README says the three properties are "not available elsewhere". Check
  Ensembl's datachecks and InterMine's acceptance tests before the paper says
  anything like it (Phase 10).
- `paper/paper.md` says the name cgload is taken on PyPI. PyPI has no project
  of that name (its JSON API returned 404 on 14 October 2026). The same comment
  points at `docs/NAME.md`, which does not exist; `docs/RENAMING.md` does.
  Probably left over from the rename (D-076). Fix in Phase 10.
- MySQL has no test at all.
- Nothing on this branch has run on a real genome, because this session cannot
  reach NCBI. `checkpoint_b.sh`, delivered with the patches, runs the slow tier
  and loads *F. graminearum* PH-1 on the server. The slow tier only compares the
  two parsers on that genome; it never loads it. Real loads come in Phase 4.
- CI's nightly job downloads through NCBI's `v2alpha` REST path, with the
  include list comma-separated. NCBI promoted `v2alpha` to `v2` in June 2024.
  Whether the old path still answers is unknown, because the workflow never
  ran. `checkpoint_b.sh` probes it; if it fails, switch `ci.yml` to the `v2`
  path the script uses.
- Changing the vocabulary hash and the schema means databases made before this
  branch are refused, with a message saying why. Reload into a new database.
