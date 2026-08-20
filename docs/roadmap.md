# Roadmap

## Milestones

Each ends with something that runs. Work does not begin on a milestone until
the previous one is green in CI.

| # | Milestone | Done when | Status |
|---|---|---|---|
| 0 | Skeleton: package, `pyproject.toml`, CI, `cgload --help`, paper skeleton | CI green | **done** |
| 1 | Schema + `init` | Tables created in SQLite; schema diagram in `docs/` | **done** |
| 2 | FASTA reading + storage layout + `stats` | Assembly loads via API (D-020), `.fai` generated, regions listed | **done** |
| 3 | Tokenizer + profile engine (D-013), RefSeq profile only, per-file detection (D-014) | Tier-1 tests pass for RefSeq; an unregistered dialect fails detection rather than falling back | **done** |
| 4 | Loader + `load` command | `examples/minimal` loads; a broken file leaves the DB *and the data dir* untouched | **done** |
| 5 | Dialects: funannotate, AUGUSTUS | One Tier-1 test per dialect, no engine change | **done** |
| 5b | GTF profile (stop-codon absorption, parent synthesis) | Bare exon/CDS GTF loads; split stop codon handled | deferred to v2.1 unless a fixture arrives |
| 6 | GenBank parser + isoform containment algorithm + load path | Ambiguous isoform case raises, naming the gene; a `.gbff` loads end to end with its linkage method recorded | **done** |
| 7 | eggNOG + InterProScan loading (D-038..D-045) | Real emapper-2.1.4 output loads against its own assembly; match rate enforced and stored | **done** (7b, 7c) |
| 8 | `verify` + `explain` (D-033, D-034, D-046..D-051) | A parser made to drop one leaf feature fails the load and rolls it back; a matched GFF3/GenBank pair agrees | **done** (8b, 8c) |
| 9 | `export` (D-035, D-052..D-054); `query` and `remove` still pending | Round trip through the database: export, re-import, identical counts, linkage provenance preserved | **export done** |
| 10 | Large-genome demo, Docker, docs, tutorial | A stranger runs one command and sees a loaded genome; maize timings published | in progress — envelope check (D-046, D-055) added, scale run pending |
| 11 | Zenodo archive, v2.0.0 tag, JOSS submission | | |

If schedule pressure appears, cut milestone 5 to v2.1 — deliberately, in
writing here, not by attrition. Milestones 0–4 plus 8–11 is still a complete
submission.

## Out of scope for v2.0

Each exclusion has a reason. Anything argued back in must be argued in
`docs/decisions.md` first.

| Excluded | Why |
|---|---|
| **Web application** | A web layer adds the largest security surface a project like this can acquire, and a deployed demo is a maintenance obligation that outlives the people who set it up. May return as a separate package consuming only the documented query API. |
| **BLAST, EMBOSS, pairwise alignment, weight matrices** | Thin wrappers around existing tools. No novelty, and all of the remote-code-execution and denial-of-service exposure. SequenceServer does the BLAST-frontend job properly. |
| **JBrowse integration** | JBrowse 2 reads indexed FASTA and tabix'd GFF3 directly. Because both are kept on disk, a user points JBrowse at our output with no code from us. Document it; do not build it. |
| **KEGG / GO web links** | An outbound hyperlink is not an integration. We store the identifiers; rendering them is a viewer's job. |
| **User accounts, authentication** | No web layer, so no need. |
| **Chromosome sequence in database columns** | See `docs/decisions.md`, D-004. |

## Deferred to v2.1+

- Additional dialects (Ensembl GFF3, GlimmerHMM, Liftoff output).
- A read-only query server or viewer package.
- Bulk multi-organism load orchestration.
- Cross-assembly ortholog relations.

## After 1.0

- `--exclude-sequences`, so loading part of an assembly does not require
  filtering the input by hand.
- GTF, with the stop-codon reconciliation it requires.
- A `feature_parent` junction table, if multi-parent annotations turn out to
  matter in practice; currently such files are refused.
- Direct-versus-inherited Gene Ontology terms, behind a flag and a user-supplied
  ontology file.
