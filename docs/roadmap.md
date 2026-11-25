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
| 8 | `verify` + `explain` (D-033b, D-034b, D-046..D-051) | A parser made to drop one leaf feature fails the load and rolls it back; a matched GFF3/GenBank pair agrees | **done** (8b, 8c) |
| 9 | `export` (D-035b, D-052..D-054), `remove` (D-083, D-084), `query` (D-085) | Round trip through the database: export, re-import, identical counts, linkage provenance preserved | **done** |
| 10 | Large-genome demo, Docker, docs, tutorial | A stranger runs one command and sees a loaded genome; maize timings published | in progress — envelope check (D-046b, D-055) added, scale run pending |
| 11 | Archived release, paper | Release archived at Software Heritage; paper describes it | in progress: v0.2.2 archived |
| 12 | Read-only web portal: `serve`, `index` (D-111 to D-114) | Search by name, gene pages, genome browser, JSON API; security and load tests in CI; measured on the 30 genomes (PROTOCOL.md section 13) | **done** in 0.3.0; measurement pending |

If schedule pressure appears, cut milestone 5 to v2.1 — deliberately, in
writing here, not by attrition. Milestones 0–4 plus 8–11 is still a complete
release.

## Out of scope for v2.0

Each exclusion has a reason. Anything argued back in must be argued in
`docs/decisions.md` first.

| Excluded | Why |
|---|---|
| **A web application that writes** | Accounts, uploads, annotation editing, a hosted service for others. A web layer that writes is the largest security surface a project like this can acquire, and a hosted demo outlives the people who set it up. The read-only portal of D-111 is in; this stays out. |
| **BLAST, EMBOSS, pairwise alignment, weight matrices** | Thin wrappers around existing tools. No novelty, and all of the remote-code-execution and denial-of-service exposure. SequenceServer does the BLAST-frontend job properly. |
| **A genome browser of our own** | JBrowse 2 and IGV read indexed FASTA and tabix'd GFF3 directly, and the portal embeds igv.js reading the same files (D-114). Nothing is drawn by code of ours beyond the gene model. |
| **KEGG / GO integration** | The portal links an accession to its database's page (D-111); it fetches nothing from them. |
| **User accounts, authentication** | The portal is read-only and has none; access control, where needed, belongs to the web server in front of it (`docs/portal.md`). |
| **Chromosome sequence in database columns** | See `docs/decisions.md`, D-004. |

## Deferred to v2.1+

- Additional dialects (Ensembl GFF3, GlimmerHMM, Liftoff output).
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
