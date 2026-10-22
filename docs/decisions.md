# Decision log

One entry per non-obvious choice: what was decided, what was rejected, and
why. Write the entry when the decision is made, not afterwards — the act of
writing it is how you find out whether you actually understood it. These
entries are also the raw material for the paper's design rationale.

---

## D-001 — Silent partial loading is the failure mode to design against

**Decided:** every subsequent decision in this log is subordinate to one
requirement: a load either stores what the file contains or fails loudly. There
is no third outcome.

Annotation loaders fail silently in characteristic ways, and each of these is a
real pattern rather than a hypothetical: features dropped inside a loop that
returns early; coding records discarded when name-based matching finds no
partner; a parent's children split across two subtrees because the parent was
recorded twice; a bare `except:` followed by a zero exit status. Every one
produces a database that looks complete, and none is visible to the person who
ran the load.

Per-row validation cannot see any of them. Start-before-end, strand in the
permitted set, coordinates within the sequence — all satisfied by a file that has
lost half its exons. What catches them is comparing independently derived counts
(D-008) and checking relationships between features rather than properties of
one (D-046).

**Rejected:** treating completeness as a post-hoc audit the user may run. It is
part of loading, it is not optional, and it can roll a load back (D-046).

---

## D-002 — SQLite is the default backend

**Decided:** SQLAlchemy Core, SQLite by default, MySQL and PostgreSQL through
the same code path via `--db-url`.

No server means no daemon, no credentials, no port and no `sudo`, which is
what makes a one-command demo and a server-free test suite possible. A whole
genome database becomes a single file that can be attached to a paper. Labs
wanting MySQL pass a SQLAlchemy URL and nothing else changes.

**Rejected:** MySQL-only — it makes every test require a running server, and
makes the install instructions long enough that many readers will not finish
them. **Rejected:** the SQLAlchemy ORM — Core is substantially easier to read,
and being able to read your own schema definitions is a project requirement.

---

## D-003 — Credentials never appear in a CLI argument or a tracked file

**Decided:** database credentials come from the environment or from a config
file that is gitignored by default.

A password passed as `--password` lands in shell history and in the process
table. A credential written into a tracked file is one `git add -A` away from
being public, and deployment scripts that emit database passwords into
client-side files are a common way for that to happen.

---

## D-004 — Assembly sequence is stored on disk, not in the database

**Decided:** the assembly is written to the data directory as FASTA with a
`.fai` index; the database stores coordinates plus a path and a SHA-256.

Storing sequence in database columns is a known mistake:
it made large genomes needlessly hard and forced every sequence retrieval
through the database. Indexed FASTA gives O(1) random access to any region and
is what samtools, IGV, JBrowse and Ensembl all do. It also keeps the database
portable, and makes JBrowse compatibility free rather than a feature we build.

**Consequence to hold on to:** the load transaction now spans two systems.
See D-006.

---

## D-005 — No web layer, and this is a security decision

**Decided:** v2 ships no server-side application.

JBrowse is a static application with
effectively zero attack surface, and that adding database queries and
shelled-out tools like BLAST around it creates a significant one — including
denial of service from unbounded heavyweight jobs. Genome portals that add
those layers routinely acquire SQL injection and remote-code-execution defects
along with them.
Writing output in browser-native indexed formats gives users the visualisation
without us operating any of the risk.

**Consequence:** v2 has no screenshots and no user-facing figures. The
contribution must therefore be legible as software and as algorithm, which is
what D-007 and D-008 are for.

---

## D-006 — The transaction boundary covers the filesystem, not just the database

**Decided:** a load stages all files (assembly copy, `.fai`, normalised
annotation, manifest) into a temporary directory, commits the database
transaction, then atomically renames the staging directory into place. On any
failure the staging directory is removed and the database transaction rolled
back.

Without this, "succeeds completely or changes nothing" is false as soon as
sequence moved to disk (D-004): a crash after the file copy but before the
commit leaves orphaned files, and a crash the other way round leaves rows
pointing at nothing. Both are silent-wrongness, which is the one failure class
this project exists to prevent.

---

## D-007 — GenBank isoform linkage is by interval containment, never by string similarity

**Decided:** for a gene with multiple isoforms, link each CDS to its parent
mRNA by exon-structure containment — same strand, every CDS interval contained
within an mRNA exon, splice junctions consistent — accepting only a unique
perfect matching. A unique match is recorded as `linkage_method='structural'`.
Multiple valid matchings or none is an error naming the gene, the counts, and
the candidate assignments.

GenBank has no parent pointers, so some rule is needed. Matching product-name
strings by resemblance is the tempting one, and it silently discards whatever it
cannot match.

**Rejected:** linking by order of appearance. GenBank does not guarantee that
mRNA and CDS blocks for a multi-isoform gene are parallel or interleaved;
features are ordered by coordinate. Positional linkage merely replaces one
silent guess with another. It is retained only as a tiebreak, recorded as
`linkage_method='positional'`, never as the primary rule.

**Rejected:** any form of string similarity, in any guise.

**Validation:** RefSeq publishes both `/protein_id` and `/transcript_id` on
GenBank records, giving ground-truth CDS→mRNA pairs. Accuracy of the
containment algorithm against those pairs, on real multi-isoform genes, is
measurable — and is the results table for the paper.

---

## D-008 — Completeness verification compares against the source file first

**Decided:** the primary, always-on check counts features in the input file
and compares against what the database holds. Comparison against
provider-published counts (NCBI assembly reports) is an opt-in secondary check
with an explicit, documented feature-type mapping.

Source-vs-database is self-consistent, available for every input including
funannotate output that has no published counts, and catches 100% of
silent-drop bugs — which is the stated purpose. NCBI's published counts use
their own biotype and deduplication conventions and will show a systematic
non-zero delta against a naive database count on most assemblies. A headline
feature that reports a spurious mismatch on the first genome anyone loads is
worse than no headline feature.

---

## D-009 — Biopython is the fifth dependency, deliberately

**Decided:** exceed the four-dependency budget for `biopython`, used only for
GenBank flat-file tokenisation.

Multi-line qualifier reassembly and the `join(complement(...))` location
grammar are fiddly, already correct in Biopython, and rewriting them means
rediscovering the same bugs at our own expense. The effort belongs in the
containment algorithm (D-007), which is the actual contribution. Biopython is
the most widely installed package in the field, so the install-failure risk
that motivates the budget is close to nil here.

**Constraint:** the dependency is confined to `parsers/genbank.py`. Biopython
objects must not appear in `NormalisedFeature` or anywhere downstream.

---

## D-010 — TOML configuration, not INI

**Decided:** organism configuration is TOML.

INI is the obvious alternative and is parsed inconsistently in practice — `:`
as the separator in one place and `=` in another, with hand-written parsers that
handle only one.
`tomllib` is in the standard library from Python 3.11, so this costs nothing
and removes an entire bug class.

---

## D-011 — The full command surface is published from the first commit

**Decided:** every v2.0 subcommand is registered immediately; unimplemented
ones exit 2 and name their milestone.

A user cannot otherwise distinguish "this tool does not do that" from "this
version does not do that yet", and a half-built command returning 0 is exactly
the silent wrongness P1 exists to prevent. A test asserts that the registered
commands match the roadmap exactly, so the two cannot drift.

---

## D-012 — The paper is written at milestone 0, not milestone 11

**Decided:** `paper/paper.md` exists from the start and is revised as
milestones land.

The statement of need determines whether the scope is publishable. That is
worth discovering while the scope can still change, rather than after eleven
milestones of work.

---

## Open

- **O-1 — Ground-truth eukaryote for the slow tier.** *S. pombe* is small,
  clean and fast enough for CI. A *Phytophthora* strain is closer to real
  research use, and demonstrates the loader at megabase scale rather than on a
  toy genome. Current intent: both — *S. pombe* in the nightly run, *Phytophthora*
  as the large-genome demo of milestone 10.
- **O-2 — Does a viewer package ever ship?** If yes it is a separate
  repository consuming only the documented query API, and it is not part of
  this submission.

---

## D-017 — The controlled vocabulary is open, and normalisation never touches
the counting column

**Decided:** the vocabulary is open, and the schema carries two type columns.
`feature.source_type` holds column 3 verbatim and is what `verify` counts;
`feature.feature_type` holds the canonical term after synonym mapping and is
what linkage rules and `query` use.

The threat to verification was never unknown types — an unknown term that
loads verbatim and counts verbatim appears identically on both sides of the
comparison. The threat is normalisation that collapses two source types into
one. Storing both terms removes the last shared input between the naive
source-side pass and the database side, which is what makes D-008's
independence claim structural rather than a matter of discipline.

Three tiers, only the first closed:

- **Tier A** — hierarchy-bearing canonical terms (11), the only types the
  loader knows how to link. Closed, because a parent-child rule cannot be
  inferred for an unseen term. A Tier A feature whose parent resolves outside
  Tier A is an error.
- **Tier B** — known synonyms that do not resolve into the hierarchy.
  Reserved; currently empty.
- **Unknown** — anything else. Loaded and counted identically, with
  `source_type == feature_type`, no linkage rules applied, and named in the
  load summary. `--strict-vocabulary` makes them an error; off by default.

**Rejected:** a closed vocabulary. It rejects exotic-but-valid annotation and
buys nothing once counting is done on the verbatim column.

**Rejected:** bundling a pinned Sequence Ontology release to separate "known
inert" from "unrecognised". The distinction buys only summary wording, and
merging them makes `--strict-vocabulary` mean the simpler thing.

**Consequence — the mapping is versioned data, not free-form data.** Adding a
row changes what `feature_type` means, and rows already stored are not
re-mapped, so a bare table would let one insert silently alter the meaning of
features loaded months earlier. Therefore: the seed carries a version and a
SHA-256 over the terms, the parent-child rules and the linkage-method
vocabulary; every assembly records the hash it was loaded under; and `load`,
`verify`, `query` and `export` refuse a mismatch until `cgload remap` has
rewritten `feature_type` across affected rows and appended to the manifest.
`explain` prints the vocabulary version it is reporting under.

**Consequence — `query genes --type mRNA` means the canonical term.** That is
the point of normalising. `stats` reports both columns by default, with
`--by={source,canonical}` to collapse, because `stats` must agree with
`verify` and `query` at the same time and cannot silently pick one.

---

## D-018 — The assembly is copied into the data directory, not referenced

**Decided:** `load` copies the assembly FASTA into the data directory and
indexes the copy. `--link` hard-links instead, for users short on disk who
accept the consequence.

A referenced file that the user later edits, moves or deletes makes every
coordinate in the database silently wrong, with nothing to detect it — the
failure class this project exists to prevent. Copying also makes the database
plus its data directory self-contained, which is what lets an assembly be
archived, moved between machines, or attached to a paper as a unit, and is
what makes the atomic staging rename of D-006 coherent: you cannot atomically
publish a directory whose contents live somewhere else.

**Cost, accepted:** a 500 Mb genome occupies twice the disk. The SHA-256 of
both source and stored copy is recorded in the manifest, so a later integrity
check needs no network and no original.

**Rejected:** referencing by default. It saves disk in exchange for a silent
correctness failure, which is the wrong trade for this project.

---

## D-019 — The rename happens inside the database transaction

**Decided:** amends D-006's ordering. Files are staged, the database
transaction is opened and rows inserted, the staging directory is atomically
renamed into place, and only then is the transaction committed. If the commit
fails, the rename is reversed.

D-006 as originally written committed the database first and renamed second,
which leaves a window where a crash produces rows pointing at a directory that
does not exist — a database that looks loaded and is not. The ordering here
narrows the window to a single rename plus a single commit, and the residual
failure is the harmless direction: a directory in place with no rows
referencing it, which is detectable, removable, and never misreports data as
present. (Since D-095 the next load into the same database replaces it.)

---

## D-020 — Milestone 2 ships an API, not a command

**Decided:** assembly ingestion is `cgload.assembly.register_assembly`, driven
by tests. `load` remains a stub exiting 2 until milestone 4.

D-011 says a command either does its whole job or says plainly that it does
not. A `load` that stored sequence and silently ignored the annotation files
in the same config would be a half-built command returning 0 — the precise
thing P1 forbids. The milestone-2 exit criterion is therefore met by the API
and its tests rather than by a user-facing command.

---

## D-021 — The assembly accession is stored, and it is opaque

**Decided:** `[assembly] accession` in the config, stored on `assembly`,
recorded in the manifest. Optional. Never resolved, never fetched, consistent
with D-016.

`verify`'s primary check is self-contained: a naive pass over the input against
`COUNT(DISTINCT source_id)` (D-008). The *optional secondary* check compares
stored counts against the figures NCBI publishes for an assembly, and that
comparison needs a key. Adding the column at milestone 2 costs one line;
adding it at milestone 8 is a schema migration against loaded databases.

**Not decided here:** whether the NCBI cross-check is offered at all in v2.0.
It requires a network fetch and is therefore opt-in at best. The column is
cheap regardless and does not commit us to the feature.

---

## D-022 — Codon rows are inert features, not CDS

**Decided:** `start_codon` and `stop_codon` are not mapped to canonical `CDS`.
They load as inert features carrying their explicit `Parent`, and are listed in
`vocabulary.NOT_CANONICALISED` with the reason, so the mapping is not
reintroduced.

**Why this reverses an earlier seed.** Real AUGUSTUS output places the start
and stop codons *inside* the terminal CDS rows: for gene `g1`, `stop_codon`
644-646 sits within `CDS` 644-739, and `start_codon` 1838-1840 within `CDS`
1651-1840. Mapping the codon rows to canonical `CDS` would therefore count the
same bases twice, and `query --type CDS` would report more CDS features than
the gene has. This was caught only by reading real AUGUSTUS output; no
synthetic fixture would have exposed it.

GTF is the opposite case — its stop codon sits *outside* the CDS — and that is
handled by the GTF profile's coordinate reconciliation (D-013 field 5), which
is a coordinate operation on the CDS span. Absorption and type mapping are
different mechanisms and must not be conflated.

**Also kept out of the hierarchy, for recorded reasons:** `intron` (AUGUSTUS
emits it, RefSeq and funannotate do not; derived information, never
synthesised) and `region` (RefSeq's first row declares the sequence region and
carries taxon and strain; consumed as region metadata, not as a gene feature).

---

## D-023 — IDs are stored verbatim; there is no ID normalisation

**Decided:** delete field 2 of D-013. No profile strips prefixes or suffixes
from IDs. Cross-references that a user would want to query on — `locus_tag`,
`transcript_id`, `protein_id`, `Name` — are stored as `feature_attribute` rows
and indexed there.

**Why.** The two proposed normalisations are both unsafe against real files.
Stripping funannotate's `-T1` turns mRNA `PHE_000001-T1` into `PHE_000001`,
which is the ID of its own parent gene — a guaranteed collision, and exactly
the injectivity failure the load-time check was added to catch. Stripping
RefSeq's `gene-`/`rna-`/`cds-` prefixes happens to be injective on the files
inspected, but it buys nothing that an attribute lookup does not, while making
`source_id` no longer the string in the file — which weakens D-008, since the
naive counting pass reads the file's own IDs.

Verbatim IDs make `source_id` mean one thing everywhere, keep the injectivity
question from arising, and leave the naive pass with nothing to reproduce.

---

## D-024 — A row with no ID is an anonymous feature; its ID stays NULL

**Decided:** `feature.source_id` is nullable. A GFF3 row carrying no `ID`
attribute is stored with `source_id = NULL`, and the feature-level count is::

    COUNT(DISTINCT source_id) + COUNT(*) WHERE source_id IS NULL

**Why this is blocking, and why it was invisible until real files arrived.**
AUGUSTUS 3.4.0 gives its `intron`, `start_codon` and `stop_codon` rows a
`Parent` and no `ID` — in the 404-line fixture, 175 of 350 child rows are
anonymous. The counting rule the work order specified, "distinct
`(seqid, column-3 type, ID)` tuples", collapses all 49 `stop_codon` rows into a
single tuple. `verify` would have compared 1 against 49 and reported a
mismatch on every AUGUSTUS file, in the feature the paper's central claim rests
on.

**Rejected: synthesising an ID** from position or parent. The naive counting
pass of D-008 would then have to reproduce the synthesis rule, which is shared
logic between the parser and its own check — the circularity D-008 exists to
prevent. NULL requires both sides to agree on one thing only: whether the row
had an ID.

**Consequence for uniqueness.** SQL treats NULLs as distinct, so
`(sequence_region_id, source_id, start, end)` cannot catch two identical
anonymous rows. The duplicate check for anonymous rows is therefore enforced in
the loader, on `(region, source_type, start, end, parent)`, not by a
constraint: MySQL has no partial indexes, and a constraint that works on two
backends of three would void D-002.

**Consequence for the count rule.** It now lives in one place,
`cgload.db.counts`, because `stats`, `verify` and `query` must not drift on
the one number the project's headline claim depends on.

---

## D-025 — AUGUSTUS comment blocks are parsed; the declared protein is a free
oracle, and the translations go to disk

**Decided:** the AUGUSTUS profile reads its `#`-prefixed comment blocks, because
two thirds of the file lives there. Specifically:

1. **Declared protein sequences are extracted.** `# protein sequence = [ ... ]`,
   continuing over `#` lines to a closing `]`, is the *only* place AUGUSTUS
   emits a translation — there is no `protein_id` pointing at a separate FASTA
   as in RefSeq and funannotate.
2. **They are written to the data directory, not the database.** D-004 puts
   sequence on disk; a protein is sequence. They land in
   `annotation/proteins.faa` beside the assembly, indexed the same way, with the
   manifest recording the count. Storing 2 kB strings in `feature_attribute`
   would bloat the table that `explain` reads on every call.
3. **The declared length is checked against the CDS spans at load time.**
   `sum(CDS spans) / 3 - 1` must equal the declared amino-acid count (minus one
   for the stop codon, which AUGUSTUS places inside the terminal CDS, D-022).
   Verified against all 49 genes in the fixture: 49 matches, no mismatches.

**Why point 3 matters more than the other two.** This is a genuine independent
oracle on our coordinate and span arithmetic, produced by AUGUSTUS itself,
requiring no network, no reference database and no second parser. If a CDS span
is off by one base, or a segment is dropped, or the discontinuous-feature
grouping attaches a segment to the wrong transcript, the arithmetic fails and
the load stops. D-008 asks for a check that shares no code with the parser; for
this dialect the file supplies one.

**Consequence for the tokenizer contract.** The skip rule is "first character is
`#`", applied identically by the tokenizer and by D-008's naive counting pass.
No comment line in the fixture contains a tab, so a field-count heuristic is
unnecessary and the simple rule is safe. A pass that counted lines rather than
features would be wrong by a factor of three on this dialect.

**Scope limit.** The hint-support statistics in the same blocks
(`% of transcript supported by hints`, `CDS exons: 0/4`) are not stored in
v2.0. They are prediction diagnostics, not annotation, and no query needs them.

---

## D-026 — `Parent=a,b` is refused, not truncated

**Decided:** a feature naming more than one parent is a load error.

GFF3 permits a comma-separated `Parent` list — an exon shared between two
isoforms is the canonical case. cgload stores one parent per feature, so the
alternatives were to duplicate the feature under each parent, which breaks the
`COUNT(DISTINCT source_id)` count that D-008 rests on, or to take the first
parent and drop the rest, which is a silent wrong answer. Refusing is the only
option that does not misreport.

**Found by a test, not by reading.** The first implementation read `Parent` with
a first-value-only accessor and would have attached such a feature to one
arbitrary branch with no error. None of the three fixtures contains a
multi-parent feature, so only the deliberately hostile test caught it. If a real
multi-parent file arrives, the fix is a `feature_parent` join table and a
revised count rule — a schema change, correctly refused rather than guessed at
now.

---

## D-027 — Recognised-but-unregistered is a distinct outcome from unrecognised

**Decided:** detection has three outcomes, not two. A registered profile loads. A
profile cgload *knows about* but has not yet enabled reports its name and the
milestone it is scheduled for. Anything else is "no confident match".

At milestone 3 only `refseq` is registered, yet `funannotate` and `augustus` are
declared in the profile table. Their files are refused — D-013's rule holds, and
there is no generic fallback — but the refusal says *recognised as funannotate,
scheduled for milestone 5* rather than *unreadable*. The difference matters to a
user holding a perfectly good file, and it lets the fail-detection trap be
exercised from milestone 3 with real files rather than waiting for milestone 5,
by which point a permissive fallback might already have been added to make
something work.

Both files are also loadable now by passing the profile explicitly in code, which
is how their end-to-end shape is already under test — 402 and 1187 lines through
the real engine — before either dialect is user-visible.

---

## D-028 — Parsing completes before any row is inserted

**Decided:** `load` tokenises, detects and normalises *every* annotation file
before it writes a single row. Cross-file ID collisions and duplicate file paths
are checked in the same phase.

The alternative — parse and insert file by file — would leave the second file's
dialect undetected until the first file's features were already in the database.
The transaction would still roll back correctly, but the load would have done
minutes of work on a large genome before discovering something knowable up front.
More importantly, parent resolution across files requires the union of features
(D-014), so the phases cannot interleave without breaking the case
multi-file support exists for.

**Consequence:** peak memory holds every normalised feature for the organism.
Measured at milestone 10 on a genome ≥500 Mb; if it proves excessive, the fix is
to stream the *insertion* while keeping detection and collision-checking up
front, not to interleave the phases.

---

## D-029 — `parent_id` points at the parent's first segment

**Decided:** a child feature's `parent_id` references the parent row with
`segment_index = 0` — its lowest-coordinate segment.

A parent stored as N segment rows has N candidate row ids, so `parent_id` needed
a defined choice or a join table. The first segment makes a walk up the graph
single-valued and deterministic, and the parent's full span remains reachable by
grouping on `source_id`, so the choice costs nothing. Documented because it is
the kind of convention a later refactor breaks silently: pointing at a different
segment would still satisfy the foreign key.

---

## D-030 — Attributes are stored against the feature, not the segment

**Decided:** column-9 attributes attach to the feature's first segment row.

RefSeq repeats every attribute on every segment of a discontinuous feature — a
six-segment CDS carries the same `protein_id`, `locus_tag` and `product` six
times. Storing them per segment would multiply `feature_attribute`, which
`explain` reads on every call, by the segment count for no retrievable
information. `ID` and `Parent` are skipped entirely: both are already modelled as
columns, and a second copy could disagree with the first. (D-096: a later row
whose attributes differ from the first row's keeps its own.)

---

## D-031 — Milestone 5 registers all three profiles; the unregistered mechanism stays

**Decided:** `funannotate` and `augustus` join `refseq` in `REGISTERED`. Each
arrived as a profile plus a real fixture plus an end-to-end test, and no engine
code changed to accommodate either — which is the claim D-013 makes and this is
the evidence for it.

`KNOWN_UNREGISTERED` is now empty and is deliberately kept. It is what lets the
*next* dialect be recognised by name before it is supported (D-027), and its
behaviour stays under test using a synthetic profile rather than by holding a
real dialect back to keep a test passing.

## D-028 — no biopython; the GenBank tokenizer is written directly

D-009 reserved a biopython dependency for `parsers/genbank.py`. Not taken. The
flat-file grammar cgload needs (LOCUS, VERSION, the feature table, and the
location expression language) is small and fully covered by the two fixtures.
Writing it directly keeps the runtime dependency budget at three packages (P4)
and removes a compiled-extension install risk, which matters wherever the
package is installed from scratch rather than from a curated environment.

Reversible: if a later milestone needs GenBank features beyond the feature
table, biopython goes back in behind the same module boundary D-009 defined.

## D-029 — GenBank linkage strategy is declared, not inferred at load time

GenBank carries no `Parent`, so every parent edge is computed. Which computation
is correct depends on the annotation, and choosing wrong is silent, so it is a
declared profile field (`LinkageStrategy`) rather than a runtime guess.

- `LOCUS_TAG` — prokaryotic. Measured on *M. aeruginosa* NIES-843: 5,808 loci,
  each exactly one gene and one child, zero loci with two CDSs. The shared
  `locus_tag` is an explicit unique key; no ambiguity is possible.
- `CONTAINMENT` — eukaryotic. Measured on *D. melanogaster* Dscam1: containment
  resolves 75 of 79 uniquely; the declared `from transcript` note resolves the
  other 4; anything still ambiguous is refused, naming the gene and candidates.
- `AUTO` — selects on the presence of `mRNA` records specifically. Testing for
  "any transcript type" selects the eukaryotic strategy for every prokaryote,
  because tRNA and rRNA are present, leaving every CDS unlinked. Caught by
  running the real bacterial genome, not by any unit test.

**Positional pairing is not used and must not be reintroduced.** The
specification said to link the nth CDS to the nth transcript, "which GenBank
guarantees". It does not: transcripts carry UTRs of differing length and
therefore sort differently from their own coding sequences. Measured wrong for
2 of 79 on Dscam1, silently. `test_positional_pairing_is_not_used` fails if it
returns.

## D-030 — GenBank identifiers may be a compound of two file-provided fields

Identifier order is `protein_id`, `transcript_id`, then `locus_tag`. Where only
`locus_tag` is available and the record is not the gene row, the record type is
appended: `MAE_RS31305:CDS`.

Necessary because a prokaryotic pseudogene's CDS carries no `protein_id` — 526
of 5,808 loci in *M. aeruginosa* — so identified by `locus_tag` alone it takes
its own gene's identifier and `_group_segments` merges the two rows into one
feature. Both halves come from the file; nothing is invented, so D-023 holds and
D-008's independent count can reproduce the key.

---

## D-032 — `linkage_method` records how the edge was found, and GenBank may never say `explicit_parent`

**Decided:** a profile declares `states_own_parents`. Where it is true — the three
GFF3 dialects — the file stated the parent and `linkage_method` is
`explicit_parent`. Where it is false — GenBank, which carries no parent pointers
at all — the tokenizer stamps the rule it used on each raw feature and the
normaliser honours it. A computed edge arriving without a stamped method is an
error, not a default.

**This was a correctness bug, not a gap.** The GenBank tokenizer synthesises
`Parent` attributes so that its output is an ordinary `TokenizedFile`, which is
good design; but `_resolve_parents` then set `explicit_parent` unconditionally, so
every GenBank edge was recorded as *the file told us* when in fact cgload had
computed it. That erases D-007's rule hierarchy from the database, fails its
acceptance check, and leaves `explain` unable to distinguish containment from a
stated pointer — which is the exact claim the paper is published on. It would have
surfaced in the manuscript rather than in a stack trace.

**Method vocabulary added:** `locus_tag` (a qualifier shared by parent and child,
unique within the record) and `note` (`/note` naming the parent transcript
explicitly). Containment is recorded as **`structural`**, the existing term D-007
and the paper use, rather than a new `containment` value — one rule, one name.
`COMPUTED_LINKAGE` marks which values mean *cgload inferred this*.

Measured on the fixtures: bacterial edges are 100% `locus_tag`; the fly's
CDS-to-transcript edges are 67 `structural` and 2 `note`, and its
transcript-to-gene edges are `locus_tag`. `positional` remains unused, per D-034.

---

## D-033 — Overlapping segments are permitted only where a qualifier explains them

**Decided:** `Profile.overlap_qualifier`. GenBank declares `ribosomal_slippage`;
every other profile declares nothing, so an overlap remains an error.

The GFF3 non-overlap rule (D-014) refused a valid record:
`join(57926..58406,58406..59148)` with `/ribosomal_slippage`. Base 58406 is
deliberately translated twice — that is what a programmed frameshift *is*. Two
such records are present in a 28-gene bacterial fixture, so this is common rather
than exotic.

Gated on the qualifier rather than on the format, so an unexplained overlap in a
GenBank file is still a duplicate and still refused. A blanket per-format
exemption would have been the easier fix and would have silently accepted real
duplicates.

---

## D-034 — Positional pairing is not implemented, and a test prevents its return

**Decided:** `positional` stays in the linkage vocabulary as a reserved value and
is never produced. The work order specified it as the tiebreak tier; reading real
records showed it is wrong.

Transcripts carry UTRs of differing length, so a gene's transcripts and its own
coding sequences appear in *different* orders. On *D. melanogaster* CG12164 the
transcripts appear RB then RA while their CDSs appear RA then RB, so nth-to-nth
pairing attaches each protein to the wrong transcript — silently, producing a
database that looks complete. Measured: wrong for 2 of 79.

This is the same class of error as pairing by product-name similarity, and it
fails the same way: plausibly, and without an error. The rule is therefore: match
on something exact, or refuse.

---

## D-035 — File format is settled by content, never by name

**Decided:** `load` reads the first line. A GenBank flat file opens with `LOCUS`;
anything else is treated as tab-delimited. A file *named* like GenBank that does
not begin with `LOCUS` is refused by name rather than silently parsed as GFF3.

The *M. aeruginosa* record downloaded as `.gbff` while actually being gzip. Trusting
the extension produces a confusing error at best and a wrong parse at worst.

---

## D-036 — A prokaryotic CDS attaches to its gene; a eukaryotic one may not bypass its transcript

**Decided:** `ALLOWED_PARENTS["CDS"]` gains `gene` and `pseudogene`, and a
separate check refuses a CDS or exon attaching directly to a gene **when the file
contains a transcript layer**.

A bacterial annotation has no mRNA records at all, so requiring one refused every
prokaryotic file — a eukaryote-shaped assumption that only real data exposed.
Simply widening the rule, though, would let a eukaryotic CDS silently skip its
transcript and make a missing mRNA layer look like a valid hierarchy. The
condition is per *file* rather than per profile, because one GenBank profile
serves both kingdoms.

---

## D-037 — A profile declares which absent types are worth a warning

**Decided:** `Profile.warns_if_missing`, a subset of `provides`. The engine's
hardcoded exclusion list (`region`, `pseudogene`, `ncRNA`) is deleted.

The GenBank profile spans both kingdoms, so `provides` is a union: a bacterial
record legitimately has no mRNA and a eukaryotic slice legitimately has no tRNA.
Warning on either fired on every correct file, and warnings that fire on correct
input train users to ignore warnings — which costs more than the warning gains.
GenBank therefore declares an empty set, which is the honest statement that the
profile cannot predict which record types a given accession carries.

---

## D-038 — Functional annotation is two tables outside the feature-graph engine

**Decided:** `annotation_run` (one row per functional file) and
`protein_annotation` (one row per hit). No dialect profile, no tokenizer, no
entry in the vocabulary.

D-013 set the scope limit and this is it being honoured: eggNOG and InterProScan
produce per-protein tables, not feature graphs. They have no parents, no
segments and no genomic coordinates, so the profile abstraction has nothing to
offer them and forcing them into it would bend it for every real dialect.

**One table rather than a table per tool plus a term table.** A GO term from
eggNOG and a Pfam domain from InterProScan are both *this protein has this
accession*, differing only in whether protein coordinates are known. Uniform rows
mean one index serves every "which proteins have X" query, and adding a third
tool later adds no schema.

**A test asserts that loading functional annotation changes no feature count.**
That is the scope limit made checkable rather than merely stated: `verify` must
be unaffected by whether annotation was attached.

---

## D-039 — The match rate is stored, not just reported

**Decided:** `annotation_run` carries `query_count`, `matched_count` and
`unmatched_count`.

An unmatched query leaves no row behind. Without these columns, the fact that
40% of a functional file matched nothing would be visible only in the load output
and unrecoverable afterwards — so a database could look annotated while half its
annotation had silently gone nowhere. This is the same reasoning as the manifest:
the absence of something must itself be recorded.

---

## D-040 — Every query protein must match a feature, unless a subset is asked for

**Decided:** `minimum_match_rate` defaults to 1.0 at every caller. A file where
any query matches no feature is refused, naming the count, the percentage and the
first few unmatched identifiers. `[[files.functional]] minimum_match_rate = 0.8`
lowers it deliberately.

The join is the whole difficulty of this milestone. eggNOG and InterProScan are
given a protein FASTA, and whatever was in those headers comes back — transcript
ID, protein ID, locus tag, or any of those with a suffix the user's pipeline
added. Nothing in the functional file says which. So a mismatch between the
functional file and the annotation is common, easy to cause, and produces exactly
the failure this project exists to prevent: a load that reports success over a
database that is half-annotated.

**Lookup is by exact key only** — `source_id`, then `protein_id`,
`transcript_id`, `orig_protein_id`, `orig_transcript_id`, `locus_tag`, `Name`,
`gene`; plus the tail after a `|`, because RefSeq writes
`gnl|WGS:AACM|FGSG_11579T0`, and the first whitespace-delimited field, because
that is the identifier by FASTA convention. No prefix matching, no similarity.
Pairing proteins to transcripts by name resemblance mispairs some and drops the
rest, and produces no error either way.

**A query matching two features is refused outright**, not resolved by
preference. Attaching a protein's annotation to the wrong isoform is invisible and
unrecoverable.

**Annotation attaches to CDS or mRNA, never to a gene.** A gene with several
isoforms has several proteins; attaching to the gene would make them
indistinguishable.

---

## D-041 — Functional files load inside the genome's transaction

**Decided:** `[[files.functional]]` entries are attached by `load`, inside the
same transaction as the assembly and the annotation. `load_functional` exists for
adding annotation to an assembly already in the database.

No new command: the eight-command surface is published and fixed (D-011), so
`annotate` would break the promise that the command list does not grow. And a
functional file that fails to attach must roll the genome back with it — a genome
loaded with its annotation half-attached is precisely the partial success that
D-006 exists to prevent. A test loads a genome with a deliberately non-matching
functional file and asserts the feature table, both annotation tables and the data
directory all come out empty.

---

## D-042 — One predicate decides which features are exempt from length arithmetic

**Decided:** `cgload.conformance.is_exempt`. Consulted by the AUGUSTUS protein
oracle now and by milestone 8's completeness check when it arrives.

Three unrelated things in real annotation produce a coding sequence that will not
divide by three or has no translation to check: pseudogenes (526 of 5,808 loci in
*M. aeruginosa*), partial features truncated by the assembly, and
`/artificial_location` where NCBI has adjusted coordinates around a known
sequencing problem (78 records in *F. graminearum*). Ribosomal slippage and
`/transl_except` are two more.

Written as one module now rather than three checks at milestone 8, because three
copies drift and the third is always the one nobody updates. An explicit
`partial=false` is *not* an exemption: reading it as one would silently disable the
check on a feature that should have it.

## D-031 — KEGG `map` pathway identifiers are dropped when the `ko` form is present

eggNOG's `KEGG_Pathway` column emits each pathway under both KEGG identifiers:
`ko00362` (the ortholog reference pathway) and `map00362` (the same pathway's
map). Measured on real `emapper-2.1.4` output for *U. maydis* (6,509 proteins):
2,253 rows carry pathway data, and in all 2,253 the `ko` set and the `map` set
are identical. 8,613 of the 17,226 rows produced — exactly half — were the same
pathway stored twice.

Storing both doubles every pathway count and pushes the convention onto every
downstream query. The `ko` form is kept: it is the identifier eggNOG uses in
`KEGG_ko` and the one KEGG documents as the ortholog pathway.

A `map` entry with **no** `ko` counterpart is kept rather than discarded. That
case does not occur in the measured file, so dropping it would be an assumption
rather than a finding.

Found only by running real tool output; the constructed fixture milestone 7 was
built against contained no `map` identifiers at all.

## D-032 — a KEGG pathway search term is canonicalised before lookup

D-031 stores the `ko` form and drops the redundant `map` form. That creates a
search hole: **`map00362` is the identifier in KEGG's own URLs and in most
published papers**, so it is what a user pastes into a search box. Without
normalisation, a pathway that is in the database returns nothing and the user
concludes it is missing — the same silent-wrong-answer failure the storage fix
was meant to prevent, moved one layer up.

`canonical_kegg_pathway()` rewrites `map<digits>` to `ko<digits>` and leaves
everything else alone, including `ko:K…` orthology terms, which share the prefix
but are a different analysis type.

It lives in `parsers/functional.py`, beside the parser, not in the query layer.
The ko/map equivalence is a fact about the source format; a second copy of it in
the query module would drift from the first.

Every code path that accepts a pathway identifier from a human must call it.
Three tests pin this: the alias resolves, unrelated identifiers are untouched,
and everything the parser stores is already canonical — so normalisation can
never produce a form that is absent from the database.

---

## D-043 — The transcript-bypass check is scoped per gene, not per file

**Decided:** a CDS or exon attaching directly to a gene is an error only when
*that gene* has no transcript child, in a file that otherwise has a transcript
layer. D-036's file-wide condition is withdrawn.

**A real RefSeq genome was refused.** For the tRNA gene `UMAG_16001`,
*M. maydis* GCF_000328475.2 emits the exons twice: once parented to the tRNA
(`ID=exon-UMAG_16001-1`, `gbkey=tRNA`) and again parented straight to the gene
(`ID=id-UMAG_16001`, `gbkey=exon`). The second pair is NCBI's gene-level exon
annotation — redundant, not a bypass. The transcript layer is present for that
gene, so nothing is being hidden.

The per-gene scope is also the stronger rule. The condition worth catching was
never "an exon touches a gene"; it was "a gene has no transcript, yet its parts
attach straight to it, while the rest of the file has transcripts" — a genuinely
missing transcript. The file-wide version caught that case and a large class of
correct RefSeq output with it.

Found only by loading a real matched pair end to end. The rule had passed 272
tests, including three GFF3 fixtures, because none of them contained a tRNA gene
with duplicated gene-level exons.

---

## D-044 — Functional annotation is inserted in bounded batches

**Decided:** `INSERT_BATCH_SIZE = 5000`. Deduplication keeps a key set and flushes
rows as it goes, rather than materialising every row and then a dictionary of all
of them.

Measured on real emapper-2.1.4 output: 57 stored findings per protein. A
13,000-protein fungal genome yields roughly three-quarters of a million rows from
one file — inside a transaction already holding the entire genome. The previous
implementation built the full row list *and* a dictionary keyed on all of it
before the first insert, which made the functional load, not the feature load, the
memory ceiling.

This is D-028's concern arriving in a sharper form. The same fix applies there
when milestone 10 measures a ≥500 Mb genome: batch the insertion, keep detection
and collision-checking up front.

---

## D-045 — GO counts must be reported as ontology terms, not as independent functions

**Decided:** anything that reports a functional-annotation count to a user —
`stats`, `query`, the paper — says *ontology terms including inherited ancestors*.
Recorded as a decision because it constrains wording in three places that do not
exist yet.

eggNOG emits the full inherited GO set, not only the specific assignments: a
protein annotated as one specific hydrolase also carries "hydrolase activity",
"catalytic activity" and every ancestor to the root. Measured on real output:
286,671 of 371,409 findings are GO terms, 89.7 per annotated protein, 570 on the
deepest.

So "the database holds 371,409 functional annotations" is a true row count and a
misleading scientific claim, and showing a user "89 functions" for one protein
misrepresents the data as thoroughly as dropping half of it would. No code change
now; the constraint is written down before the code that would violate it exists.

**Not decided:** whether to store only the most specific terms. That needs the GO
graph, which is a dependency and a staleness problem (D-016's reasoning), and the
inherited terms are what make a query for a broad category work at all.
`canonical_kegg_pathway()` rewrites `map<digits>` to `ko<digits>` and leaves
everything else alone, including `ko:K…` orthology terms, which share the prefix
but are a different analysis type.

It lives in `parsers/functional.py`, beside the parser, not in the query layer.
The ko/map equivalence is a fact about the source format; a second copy of it in
the query module would drift from the first.

Every code path that accepts a pathway identifier from a human must call it.
Three tests pin this: the alias resolves, unrelated identifiers are untouched,
and everything the parser stores is already canonical — so normalisation can
never produce a form that is absent from the database.

## D-033 — a check may fail the command only when disagreement implies cgload is wrong

Three tiers, not equally authoritative:

**Tier 1, source count** — the naive pass versus the database. Disagreement
always means cgload is wrong. Fails, exit 2.

**Tier 2, cross-parser** — the same assembly read as GFF3 and as GenBank by two
engines sharing no code. A GFF3 parser that loses a feature cannot make the
GenBank count lose it too, so agreement cannot be faked and disagreement has no
innocent explanation. Fails, exit 2. Measured on GCF_000240135.3: gene 13,725 =
13,725 (13,714 + 11 pseudogene), mRNA 13,315, CDS 13,312, tRNA 322, rRNA 88 --
exact on every type. Mutation-tested: deleting one gene and dropping three tRNAs
are both caught with the exact delta.

**Tier 3, NCBI published counts** — disagreement has at least four innocent
causes: a partial assembly was loaded; the annotation release differs; NCBI
tallies non-coding and pseudogenes under different classes; re-annotation happens
without the accession changing. Only the first is detectable, so it is the
precondition. When it holds, disagreement *does* imply cgload is wrong and the
tier fails. When it does not, the tier is **SKIPPED with the reason printed** --
not reported as an advisory warning. A check that raises concerns it cannot stand
behind is one people learn to scroll past, which leaves you no better off than
having no check.

Reversal of an earlier assumption: NCBI's counts were thought to be
"different quantities by construction" because the report groups by class rather
than by GFF3 column 3. Summing `Class` within `Feature` gives exactly the
column-3 counts -- measured on GCF_000328475.2, all five types match, with CDS
matching at feature level (6,780) rather than row level (9,739) and thereby
independently validating segment grouping.

## D-034 — three traps in the NCBI feature-count format

All three are present in a real seven-line file (GCF_000328475.2) and each would
silently break the tier-3 comparison:

1. Counts must be **summed across `Class` within a `Feature`**. `gene` appears
   three times: protein_coding 6,762 + rRNA 34 + tRNA 111 = 6,907. Reading one
   row reports 145 missing genes on a correct load.
2. **`Unique Ids` can be the literal string `na`** while `Placements` carries the
   number. The `tRNA` row is exactly this; a reader taking only the first column
   drops the type and reports 111 missing tRNAs.
3. **`Assembly-unit accession` is not the assembly accession** -- GCF_000328485.2
   versus GCF_000328475.2, one digit apart. Matching the wrong column selects no
   rows, and a verifier that silently matched nothing would report every feature
   as missing. `parse_ncbi_feature_counts` raises rather than returning empty.

---

## D-046 — The completeness check runs on every load and can roll it back

**Decided:** `load` runs the primary check itself, after inserting and before
committing. A failure raises and the transaction rolls back — database and data
directory both untouched. The result is written into the manifest either way.

`verify` remains a command for checking an assembly later, when the database may
have been altered by something other than cgload. But the check that catches the
project's defining failure mode is not a step a user has to remember: a load that
silently lost features and reported success is exactly what this tool exists to
prevent, so the check is part of loading rather than adjacent to it.

The failure message names the type and both counts, and says plainly that this is
a bug in cgload rather than in the user's file — because it is. No arrangement of
a valid input file can make the naive count disagree with the stored count.

---

## D-047 — The database half of the check lives outside `verify.py`

**Decided:** `cgload.verify` holds the naive counting pass and the two
file-to-file comparisons, and imports nothing from cgload. The stored-side count
and the file-to-database comparison live in `cgload.db.verification`.

A test enforces the separation by scanning `verify.py` for the strings
`tokenize` and `normalise`, so anything needing the schema had to live elsewhere
regardless. The split is worth having on its own terms: `verify.py` answers *what
does the file say*, `db/verification.py` answers *what did we store*, and neither
can see the other's method. That is the property D-008 rests on, expressed as
module boundaries rather than as discipline.

**Measured:** the naive pass reproduces the parser's per-type counts exactly on
all six real fixtures — RefSeq, funannotate, AUGUSTUS and both GenBank files —
including AUGUSTUS's 175 rows with no `ID`, GenBank's `join()` spans, and the
RefSeq CDS segments that share one identifier. Two independent implementations
agreeing on every real file is what makes the check meaningful rather than merely
present.

---

## D-048 — `explain` reads recorded values and recomputes nothing

**Decided:** every value `explain` prints was stored at load time. It performs no
linkage, no parsing and no inference.

A recomputed answer could differ from the stored one, and neither would be
checkable against the other — so `explain` would become a second opinion rather
than an audit. Reading only what was recorded also makes it a test of D-032: if
the linkage method were not stored, `explain` could not report it, and the fact
that it can is evidence the provenance survived to the database.

It states in words whether a link was **stated by the file** or **computed by
cgload**, which is the distinction the paper is published on made visible to a
user one feature at a time.

---

## D-049 — Shared IDs are governed by a deny-list, not an allow-list

**Decided:** `Profile.unique_id_types` names the types where a repeated ID is
always an error. Everything else may share an ID, subject to the D-014
conditions (same seqid, type, strand and parent, non-overlapping spans).
Replaces `shares_id_across_segments`, which named the types permitted to share.

The allow-list was `{"CDS"}`, derived from three fixtures in which only CDS rows
shared identifiers. The full *F. graminearum* GFF3 refused to load on line 47,781:

```
sequence_feature  5789158..5789507  ID=id-NC_026475.1:5789089..5789507
sequence_feature  5789089..5789107  ID=id-NC_026475.1:5789089..5789507
```

A `misc_feature` with a `join()` location, split across two segments sharing one
ID — a valid discontinuous feature of a type the allow-list had never seen. An
allow-list is wrong here because the set of types that *may* be discontinuous is
open-ended, while the set that must not be is small and observed: exons carry
their own IDs in every dialect examined. Declaring the invariant we measured,
rather than enumerating everything else, keeps the protection and stops refusing
valid records.

---

## D-050 — Gene identity is `/gene` or `/locus_tag`, and transcripts are indexed by it

**Decided:** `_gene_identity` takes `/gene` where present and `/locus_tag`
otherwise, used by both the transcript-to-gene link and the containment grouping.
Transcripts are indexed by gene once per file rather than rescanned per CDS.

**A correctness bug and a performance bug in the same place.** Containment read
only `/gene`; the transcript-to-gene link read `/locus_tag` then `/gene`. RefSeq's
*F. graminearum* GenBank release uses `/locus_tag` exclusively and never `/gene`,
so every CDS in the file had gene identity `None`, fell into one group, and was
tested against all 13,315 transcripts in the genome. It failed on the first
record with an eight-way ambiguity that was an artefact of the grouping, not a
property of the data.

Fixing the identity exposed the second problem: the candidate scan was linear in
the whole transcript list, so linkage was quadratic in genome size — 173 seconds
for the GenBank half against 10 for the equivalent GFF3. Indexing by gene brings
it to 6.4 seconds, unchanged results, and removes what would have been the
binding constraint on milestone 10's large-genome demo.

**Measured on the full genome after both fixes:** 13,725 genes, 13,315
transcripts, 13,312 coding sequences. Every CDS-to-transcript edge resolved by
containment alone — 13,312 `structural`, zero ambiguity failures, zero needing
the `/note` fallback.

---

## D-051 — The matched pair is a fixture, not a footnote

**Decided:** `pair_fgraminearum_NC_026477.gff3` and
`pair_fgraminearum_NC_026477.gbff` ship in the repository — the same 43 genes of
one chromosome region, written in two formats, graph-closed in both.

The cross-format comparison is the strongest evidence this project has, and until
now it could only be demonstrated on files a user supplied. Anyone asked to
believe a claim they cannot reproduce from the repository will discount it.

Both halves are around 60 and 120 kB, which is the price of the project's best
argument being checkable by `pytest` rather than by a download.

**Full-assembly measurement, for the paper:** GCF_000240135.3 read from both
formats agrees exactly on every shared type — mRNA 13,315, CDS 13,312, tRNA 322,
rRNA 88 — and on genes once the naming difference is accounted for: GFF3 reports
13,714 genes plus 11 pseudogenes, GenBank reports 13,725 genes and marks the
pseudogenes with a note. 13,714 + 11 = 13,725. Accessions and commands are in
`PROVENANCE.md`.
## D-035 — export emits a cgload dialect, never the source dialect

Round-tripping back into the original dialect is where export quietly becomes
dishonest. Emitting AUGUSTUS would require synthesising the `exon` rows stage 5
deliberately refused to invent; emitting RefSeq would require fabricating
`gbkey` and `Dbxref` values that were never in the input. A file claiming to be
RefSeq that is not is the same class of failure as a load reporting a success it
cannot back up.

Four rules:

1. **Column 2 is the stored `source_program`, verbatim.** Already per-feature,
   so an assembly loaded from several files exports correctly with no extra
   machinery — an AUGUSTUS gene and a RefSeq gene keep their own column 2 in one
   output file. Tested.
2. **Column 3 is `source_type`, not `feature_type`.** The same reason `verify`
   counts `source_type`: exporting the canonical term rewrites AUGUSTUS's
   `transcript` to `mRNA` and makes the round trip lossy on any dialect using a
   synonym. Tested.
3. **Nothing is synthesised.** AUGUSTUS in, no exons out. A header comment names
   the types the source dialect does not emit, so a downstream operator learns
   it from the file rather than from a support thread.
4. **The header declares `##cgload-export`,** the version and the vocabulary
   hash, so a re-imported file is recognised as cgload output.

**Anonymous features are the one place the round trip cannot be exact by
identifier.** A row with no ID — AUGUSTUS `intron`, `start_codon`, `stop_codon`,
175 of 350 child rows in one fixture — has nothing to write in column 9.
Synthesising one would break D-023. They are emitted ID-less and re-import
correctly, because the parser and the naive counter both treat an ID-less row as
its own feature (D-024).

Consequence, stated so it is not rediscovered the hard way: **the round-trip
assertion is on feature counts per source type, never on identifier sets.** The
stricter assertion fails, and the obvious fix is to synthesise IDs — the exact
thing this decision forbids.

## D-036 — the `cgload` profile must outrank column 2 in detection

Measured: an export of an AUGUSTUS load still carries `AUGUSTUS` in column 2, so
`detect()` returns the `augustus` profile and applies AUGUSTUS's rules to a file
that no longer follows them.

The `cgload` profile's `header_signature` (`##cgload-export`) must therefore be
consulted **before** the column-2 signature, for that profile only. This is the
single edit milestone 9 requires outside its own new files;
`test_reimported_export_detects_as_cgload` is marked `xfail(strict=True)` and
will fail loudly once the profile is registered correctly, at which point the
marker comes off.

---

## D-052 — The export header outranks column 2, and `cgload` is a registered profile

**Decided:** a `cgload` profile joins the registry, and `detect` checks for the
`##cgload-export` pragma *before* consulting column 2.

Column 2 of an export is the original producing program, preserved per feature
(D-014) — an export of an AUGUSTUS load says `AUGUSTUS` there, and correctly so.
But the file no longer follows AUGUSTUS's conventions, so without header
precedence the augustus profile claims it and applies the wrong rules to a file
that states plainly what it is. Milestone 9 identified this and shipped it as a
deliberately failing test rather than a silent gap, which was the right call; the
fix is one profile and four lines in `detect`.

---

## D-053 — An export is not a self-stating format

**Decided:** the `cgload` profile sets `states_own_parents=False`, and the
exporter writes the linkage method under the same attribute name the normaliser
reads.

**This was a live provenance-loss bug, caught by round-tripping through the
database rather than through the parser.** A GenBank load stores 43 CDS features
whose parents were computed by containment. The export writes `Parent` attributes
for them — it must, or the graph is lost — so treating the export as a format that
states its own parentage relabelled every one of them `explicit_parent` on
re-import. Measured before the fix: 148 `structural` rows in the database came
back as 86 `explicit_parent` features.

That is exactly the false-provenance failure D-032 was written to prevent, arriving
by a path D-032 did not cover. The export's truth about parentage is the recorded
stamp, not the presence of a `Parent` attribute, and the profile now says so.

Two contributing details, both worth naming because either alone would have
reintroduced it: the exporter wrote `cgload_linkage` where the normaliser reads
`cgload_linkage_method`, so the stamp was discarded before the profile flag ever
mattered; and the two constants live in modules that deliberately do not import
each other, so a test asserts they agree rather than a comment asking nicely.

---

## D-054 — Export reads the database, and the round trip goes through it

**Decided:** `cgload.db.export_source` reconstructs stored rows into the shape
`normalise` produces, and `cgload export` writes from the database. The
serialiser in `cgload.export` continues to know nothing about the schema.

Milestone 9's round trip went parser → exporter → counter, which tests the
serialiser and not the command: the point of `export` is to write out what is in
the database *after the input files are gone*. The round trip that matters is
parser → database → exporter → parser, and it now runs on four fixtures across
three dialects.

The split is kept because the serialiser stays testable against parser output
directly, and because a duck-typed `StoredFeature` — rather than reusing
`NormalisedFeature` — keeps the exporter free of any import path back into parser
code.

**`--format tsv`** writes functional annotation, one assignment per row, with a
header comment recording that GO rows include inherited ancestors (D-045). Anyone
summing that column needs to know before they do it.
## D-046 — a parent envelope must contain its children; repair is opt-in and bottom-up

Every other structural rule in cgload is per-row. A row can satisfy all of them
and still be wrong, and nothing downstream can tell.

Measured on *Zea mays* B73 GCF_902167145.1, 1.2 M features: **36 structural
problems, of which 1 fails a per-row check and 34 do not.** The 34 load cleanly
with parent envelopes wrong by up to 281 kb. Only the parent-child relationship
exposes them.

**Cause.** All 34 are organellar — mitochondrial `ZeamMp*` and chloroplast
`ZemaCp*` — across seven trans-spliced genes. NCBI's GFF3 converter emits the
parent envelope from **one segment** of the trans-spliced set rather than their
union, and which segment varies: last for `gene-ZeamMp186`, first for
`gene-ZeamMp016`, a middle one for `gene-ZeamMp071`. The GenBank file for the
same assembly has the coordinates right, which is how the children are known to
be correct. To be reported to NCBI.

An early hypothesis that strand `?` marked the affected rows was **wrong**: the
affected parents are six `-`, two `?` and one `+`. A fix scoped to `?` would
have missed two thirds of the cases.

**Policy.** Detect always; refuse by default; `--repair-envelopes` recomputes.

**Repair is bottom-up, and that is not a preference.** For `ZeamMp186` both the
mRNA and the gene envelope are wrong. Repairing the gene from its immediate
child — the still-broken mRNA — gives `50,490..267,232`: plausible, and wrong by
281 kb. Leaves-upward gives `50,490..548,772`, matching the GenBank join exactly.
A top-down pass produces a number that looks fixed and is not, which is worse
than refusing.

Three constraints on repair:

- **Only widens, never narrows.** Narrowing would discard a claim the file
  makes; widening only admits children the parent already claims.
- **Multi-segment parents are never repaired.** Widening one segment of several
  means choosing which, and the file gives no basis for that choice. It remains
  a violation and refuses the load.
- **Re-check after repair.** A repair that does not clear the violation means the
  cause is not a truncated envelope; proceeding would be the silent partial
  success this project exists to prevent.

Every repair is reported at load time with old and new envelope and recorded in
`manifest.json`, so a repaired load is declared rather than merely successful.

---

## D-055 — The envelope check is wired into `load`, and repair is opt-in

**Decided:** `load` runs the containment check on every file. A parent that does
not contain a feature claiming it is a load error. `--repair-envelopes`
recomputes each parent's extent from its contents, leaves-first, and every change
is printed and written to `manifest.json`.

The check itself is the delivered work (D-046); this is the wiring, and three
properties of it are load-bearing.

**Refusing is the default.** A tool that loads coordinates it knows are
contradictory is worse than one that stops. The refusal names every offending
feature and names the flag, so the user is told both what is wrong and what to do.

**Repair is declared, not silent.** The manifest carries an `envelope_repairs`
list on *every* load, empty when nothing was changed. A missing key and an empty
list are different facts, and someone auditing a database later has to be able to
find out that cgload altered coordinates the provider published.

**The repaired span is checked against an independent source.** For maize `nad1`
the repair produces 50,490..548,772, which is exactly what the GenBank release of
the same assembly states. Repairing outside-in instead produces 50,490..267,232 —
plausible, wrong, and indistinguishable without that second file. The test asserts
the GenBank number rather than "some wider span", which is the only version of
that test worth having.

**Measured:** zero violations across all eight real fixtures in the repository —
RefSeq, funannotate, AUGUSTUS and both GenBank files. A check that fired on
correct data would be worse than no check, so that is asserted rather than
assumed.
## D-056 — an inverted row is the same provider bug, and is repairable

*Zea mays* B73 GCF_902167145.1 contains one row where start exceeds end:

    NC_007982.1  RefSeq  mRNA  691776  267232  .  ?  .  ID=rna-ZeamMp017;...

This is not a distinct defect. It is the same NCBI trans-splicing converter bug
as the 34 too-small envelopes (D-046), arriving in the one form a per-row check
happens to catch. Neither number is an extent: 267,232 is the end of one
segment, and 691,776 is not a coordinate on this molecule at all — NC_007982.1
is 569,630 bp. The feature's children span 122,146..548,772, which is correct
and lies inside the sequence.

Before this change cgload **repaired one form of the bug and refused the
other**, though both are recoverable from the same evidence, because the
tokenizer rejected the row long before any envelope logic ran. `--repair-envelopes`
could not help, and the maize load failed with the flag set.

**Coordinates are never swapped.** Reversing 691,776..267,232 gives
267,232..691,776, which runs 122,146 bases past the end of the molecule — so a
swap is not merely a guess, it is a demonstrably wrong one. The written pair is
carried through verbatim with an `inverted` flag, and the extent is recomputed
from the children alone. Widening from the written pair would be meaningless
because neither number is an extent.

- `tokenize(..., tolerate_inverted=True)` records the violation instead of
  raising. Default remains `False`; behaviour without the flag is unchanged.
- `find_violations` reports an inverted feature on its own, before any parent is
  consulted, since it contradicts itself.
- `repair_envelopes` takes the children's span directly for such a feature, and
  clears the flag so the post-repair re-check does not refuse a row it just
  fixed.
- An inverted row **with no children cannot be repaired** — there is nothing to
  recompute from — so it remains a violation and refuses the load.

Repair remains bottom-up: the `ZeamMp017` gene envelope is also wrong and is
repaired from the repaired mRNA, not the broken one.

**Wiring, added after the fact (D-056b).** `--repair-envelopes` has to reach the
*tokenizer*, not only the envelope stage. An inverted row is rejected at read
time, before any envelope logic runs, so a flag that stopped at `apply_policy`
could not help the one maize row it was written for — and the maize load would
have failed with the flag set, which is what prompted this decision in the first
place. `load_organism(..., repair_envelopes=True)` now passes
`tolerate_inverted=True` through to `tokenize`, and a test asserts the whole
chain rather than the envelope stage alone.

**Message shape.** A self-contradicting row is rendered differently from a
containment violation: the parent/child wording describes nothing when the parent
*is* the child, and the first version read
`(self) <anonymous> 691,776..267,232 lies outside gene gene-lonely 691,776..267,232`,
which is noise. It now reads `gene gene-lonely is written 691,776..267,232, which
runs backwards, and has no children to recompute an extent from`, and a test
asserts the containment wording never appears for a self-contradiction.

---

## D-057 — Deliberately malformed fixtures live in their own directory

**Decided:** `tests/fixtures/broken/` holds files that are not valid input.
`tests/fixtures/` holds only files that load cleanly.

Adding `inverted_row_envelope.gff3` alongside the real fixtures broke two tests
immediately — both iterate every fixture and assert it tokenizes and round-trips,
and a file that is *supposed* to be refused fails that by construction.

The tempting fix is a per-name exclusion in each globbing test. That is how those
tests stop asserting anything: every future broken fixture adds another name to
another list, and eventually the "every fixture round-trips" test covers a
shrinking subset nobody has re-read. A directory boundary says the same thing
once, and makes "every fixture in `fixtures/` is valid input" a property rather
than a convention.

It also strengthened a test rather than weakening it: `test_no_correct_fixture_reports_a_violation`
no longer needs to exclude the broken file by name, so it now genuinely means
*every* file it iterates.
# Decisions D-059 – D-061

Append to `docs/decisions.md` after D-058. Written against the four-file drop
`cgload-pseudogene-fix.zip`; see the "Not verified" note at the end of D-060 for
what could not be checked because the rest of the tree was absent.

A correction to D-058's own note first: it describes itself as **the third**
rule in cgload stated about the whole file when it should have been stated
about each parent. It is the **fifth**. The full sequence is D-049 (the
allow-list of types permitted to share an ID across segments, which refused a
valid RefSeq `sequence_feature` with a `join()` location), D-045 (the strand
vocabulary, which refused GFF3's valid `?`), D-036 (`ALLOWED_PARENTS["CDS"]`
lacking `gene`, which refused every prokaryotic GenBank record), D-043 (the
transcript-bypass rule needing per-gene rather than per-file scope), and D-058.
The prose in `normalise.py` says "second" and the change note says "third";
both undercount, and the undercount is itself the problem — each instance has
been recorded as a one-off rather than as a recurring class.

---

## D-059 — `has_transcripts` tests for `mRNA` alone, and that is deliberate

`_check_hierarchy` computes two things that look like they should agree and do
not:

```python
has_transcripts       = any(f.feature_type == "mRNA" for f in features)
genes_with_transcripts = {... if f.feature_type in TRANSCRIPT_TERMS ...}
```

`TRANSCRIPT_TERMS` is `{mRNA, tRNA, rRNA, ncRNA}`. Widening `has_transcripts` to
match is the obvious tidy-up. **It was attempted during this audit and it
refuses every prokaryotic record.** Every bacterial genome carries tRNA and rRNA
genes, so the file reads as transcript-bearing, and then every `gene`→`CDS`
edge — the normal prokaryotic shape that D-036 exists to permit — trips the
transcript-bypass rule. Reverted.

So the line is not asking "does this file have a transcript layer". It is a
proxy asking **"is this a eukaryotic annotation"**, and `mRNA` is the marker
that answers it. Recorded because the code did not say so, and the next reader
to notice the inconsistency will make the same change.

**Accepted cost.** A eukaryotic file whose transcript layer is entirely
non-coding — an ncRNA-only scaffold, an organellar annotation, a filtered slice
containing only lncRNA genes — reads as prokaryotic. The bypass rule never runs,
and a genuinely missing transcript in such a file loads unchecked. This is a
false *negative*, which is the cheaper error: it loads correct data correctly and
fails to catch one class of malformed data. Not fixed, because nothing available
at that point in `_check_hierarchy` distinguishes "eukaryotic ncRNA-only file"
from "prokaryotic record".

**If it needs fixing later**, the predicate does not belong in `normalise.py` at
all — the profile knows the kingdom, or the `region` row's `taxon` does, and
either is a fact about the source rather than a guess from the contents.

Pinned by `tests/test_transcript_layer_scope.py`: a prokaryotic fixture that must
load, and a strict `xfail` on the ncRNA-only gap that will start failing the day
someone improves the predicate.

---

## D-060 — permitted-set audit register

The recurring defect, stated once so it does not need rediscovering:

> **A permitted-set is a guess about the world; a forbidden-set is a claim about
> a specific observed failure.** A rule that enumerates what is allowed, derived
> from whichever fixtures were on hand, refuses correct data the first time a
> genome exceeds them — and it refuses it as a hard error, at load, on a file
> that is not wrong.

Every remaining rule of that shape in `normalise.py` and `db/vocabulary.py`,
with a recommendation. Classification: **(a)** invert to a forbidden-set,
**(b)** keep as a refusal because the forbidden version is genuinely statable,
**(c)** downgrade to a warning because the forbidden version cannot be stated
and refusing correct data is the worse error.

### 1. `SYNONYMS` missing `miRNA` — (a), highest priority

*Permits:* the 21 source terms listed. Unmapped terms are not rejected — they
load inert — but **an inert parent makes its Tier A children unloadable**, so
the open vocabulary is open in one direction only.

*Derived from:* four fungal fixtures plus whatever RefSeq rows the maize work
touched. Fungi have almost no miRNA.

*Exceeded by:* RefSeq's miRNA shape, `gene → primary_transcript → miRNA → exon`.
`primary_transcript` maps to `mRNA`, but `miRNA` maps to nothing, so it is inert,
and its `exon` child is refused with *"a exon may only attach to ['gene', 'mRNA',
'ncRNA', 'pseudogene', 'rRNA', 'tRNA']"*. Reproduced.

*Why it matters now:* **Zea mays B73 — the genome in D-058 — contains miRNA
genes in this shape.** The load that stopped at `gene-LOC103630212` will very
likely stop again a few thousand lines later on a miRNA exon. This is the sixth
instance of the class and it is already queued.

*Not fixed here, because the one-line fix does not work.* Mapping
`miRNA → ncRNA` leaves `ALLOWED_PARENTS["ncRNA"] == {gene, pseudogene}`, and this
miRNA's parent is a `primary_transcript` (canonically `mRNA`), so it is still
refused. RefSeq's miRNA is genuinely **two transcript tiers deep** and cgload's
model has one. The real options are (i) map `miRNA → ncRNA` *and* add `mRNA` to
`ALLOWED_PARENTS["ncRNA"]`, admitting nested transcripts; or (ii) collapse
`primary_transcript`/`miRNA` to a single tier in the profile. Both are model
changes, neither is low-risk, and choosing needs a real RefSeq file.

### 2. `SYNONYMS` mapping the four `*_gene_segment` terms to `gene` — (a)

*Permits:* `V_gene_segment`, `C_gene_segment`, `D_gene_segment`,
`J_gene_segment` → canonical `gene`.

*Derived from:* not from any fixture visible here — it looks like a reading of
the SO term names rather than of a file.

*Exceeded by:* its own vocabulary. `ALLOWED_PARENTS["gene"]` is `frozenset()`,
i.e. a `gene` may have no parent at all. So the moment one of these rows carries
a `Parent` — which it does in RefSeq human and mouse, where the segment sits
*under* a gene and *above* its exons — normalisation fails with *"a gene may
only attach to nothing"*. Reproduced.

This one is unambiguous as a defect regardless of what RefSeq actually emits:
the two tables contradict each other. Four terms are declared canonicalisable
into a type that cannot appear in a parented position.

*Not fixed here* because the correct target depends on the real file and the two
branches fail in opposite directions. If the segment is transcript-level
(parented under a gene, exons beneath it) the mapping should be `mRNA`. If it is
gene-level (top-level, no parent) the current mapping is right and the bug is
elsewhere. Mapping to `mRNA` when it is actually top-level converts a working
load into *"an mRNA must attach to one of ['gene', 'pseudogene']"*. **One row of
real RefSeq human GFF3 settles it** — see "Not verified".

### 3. `ALLOWED_PARENTS["five_prime_UTR"]` and `["three_prime_UTR"]` — (a)

*Permits:* `{mRNA}` only.

*Derived from:* the fungal GFF3 fixtures, where UTRs always hang off an mRNA.

*Exceeded by:* the same shape D-036 fixed for `CDS`. When a dialect has no
transcript layer, a UTR attaches to the gene. D-036 widened `CDS` and D-058
widened `exon`; the UTRs were left behind, so they are the next member of the
same family to fail. Reproduced: a `five_prime_UTR` on a `gene` is refused.

*Recommendation:* invert. The forbidden version is easy to state and does not
depend on future fixtures — **a UTR may not attach to a `CDS`, to another UTR,
or to a `repeat_region`.** Everything else is a plausible parent. Add both UTR
types to `TRANSCRIPT_LAYER_REQUIRED` at the same time, so widening does not
reopen the concealed-transcript hole D-043 closed.

*Not implemented here* only because no fixture or real file in evidence actually
attaches a UTR to a gene; the failure is predicted from the family resemblance,
not observed. It is the cheapest of the three to do once someone has one.

### 4. `ALLOWED_PARENTS` as a whole — (b), keep

*Permits:* a fixed parent set per Tier A child type.

*Assessment:* keep it as a refusal. The forbidden version *is* statable here —
"a CDS may not attach to a `repeat_region`" — and the check earns its place: it
is what catches a genuinely scrambled hierarchy. The defect has never been the
table's existence, only that individual rows were written narrow. The durable
protection is the fixture sweep, not inversion.

*Standing rule for edits:* when widening a row, ask whether the widening reopens
a concealment the transcript-layer rule was closing, and if so add the child type
to `TRANSCRIPT_LAYER_REQUIRED` in the same commit. D-036 and D-058 both did this
correctly; it is the pattern to copy.

### 5. `HIERARCHY_TERMS` under `--strict-vocabulary` — (c)

*Permits:* 11 canonical terms.

*Assessment:* harmless by default — an unlisted term loads verbatim, is counted
normally, and is named in the load summary, which is the correct open-vocabulary
behaviour. `--strict-vocabulary` converts that open set into a permitted-set that
refuses. Since the whole point of the flag is to be strict, keep it, but it
should refuse with the *list of unrecognised terms in this file* rather than one
term at a time, so a user adding a new dialect learns the whole gap in one run
instead of by repeated failed loads. Cheap, and it is the difference between the
flag being usable and being switched off permanently after the second failure.

### 6. `_new_feature`, multi-parent refusal — (b), keep, improve the message

*Permits:* exactly one `Parent`.

*Exceeded by:* GFF3 explicitly permits multiple `Parent` values, and the natural
use is one exon shared by two transcripts. Reproduced: refused.

*Assessment:* keep. This is not a guess about the world, it is a true statement
about the storage model — one parent column per feature — and silently picking
the first parent or duplicating the feature would corrupt counts. It is the
*right* kind of refusal: a claim about a specific known limitation. NCBI and
Ensembl both emit one `Parent` per row, so the practical exposure is low. Worth
adding to the message which transcripts were named, so a user who hits it can
decide whether to split the file.

### 7. `_assert_same_feature`, overlap gated on `profile.overlap_qualifier` — (b)

*Permits:* overlapping segments of one feature only when the profile names a
qualifier and the row carries it (D-033, programmed frameshift).

*Exceeded by:* a GFF3 file describing a frameshifted feature. GFF3 has no
equivalent qualifier, so the escape hatch is unreachable in that dialect.

*Assessment:* keep. D-033 reasoned this deliberately and an unexplained overlap
really is an error in GFF3. Noted rather than actioned.

### 8. `_resolve_parents`, dangling-parent refusal — (b), keep

*Permits:* only parents defined in the same file.

*Exceeded by:* a legitimately sharded annotation — per-chromosome GFF3 where a
`Parent` lives in a sibling file.

*Assessment:* keep as a refusal; the message already says the file is incomplete
or was sliced badly, which is exactly right. Worth naming multi-file loads in the
message as the other explanation, since that is the case where the user's file is
fine and their invocation is not.

### 9. `LINKAGE_METHODS` — (b), keep

*Permits:* seven method names. Unlike everything else here, the values are
produced by cgload's own tokenizer rather than read from external data, so it
is a closed internal vocabulary, not a guess about the world. Correctly a
refusal.

### 10. `profile.unique_id_types`, `permits_anonymous_rows`, `overlap_qualifier` — **not audited**

Per-profile permitted-sets, and `profiles.py` was not in the drop. `unique_id_types`
in particular is the D-049 allow-list, i.e. a known member of this class that was
already narrowed once. It should be re-read against this register.

---

## D-061 — `content_hash()` now covers the transcript-layer constants

The docstring on `content_hash()` promises that *anything that changes how a
source term is linked must change this hash*, because editing the hierarchy
rules without bumping the hash would silently change the meaning of an existing
database. The hash covered `seed_rows()`, `ALLOWED_PARENTS` and
`LINKAGE_METHODS`. It did **not** cover `TRANSCRIPT_TERMS` or
`TRANSCRIPT_LAYER_REQUIRED`, and both decide whether a `gene`→`CDS` or
`gene`→`exon` edge is accepted at all.

D-058 is the demonstration. It made two changes: widening
`ALLOWED_PARENTS["exon"]`, which moved the hash, and narrowing the
transcript-layer rule, which did not. Had the change been only the second half —
which is the half that actually unblocked the maize load — every existing
assembly would have kept reporting a matching hash while a different set of files
became loadable.

Both sets are now hashed. Covered by `tests/test_vocabulary_hash.py`, which
asserts movement for `ALLOWED_PARENTS`, `SYNONYMS`, `TRANSCRIPT_TERMS` and
`TRANSCRIPT_LAYER_REQUIRED` by monkeypatching each in turn, so the property is
tested rather than a literal digest being pinned.

**Migration consequence, which D-058 did not record:** the hash changed. Every
existing assembly now mismatches and is refused until `cgload remap` runs, even
though `seed_rows()` is byte-identical and the `vocabulary` table content has not
moved. That is the designed behaviour and it is correct — the *meaning* of
linkage did change — but it is a breaking change for existing databases and
belongs in the release note rather than being discovered at the next load.
`VOCABULARY_VERSION` is still `1`; it arguably should move too, since the version
is the human-readable half of the same signal.
## D-058 — a pseudogene is a legitimate terminal parent

Two rules refused standard NCBI output, both written from what four fungal
fixtures happened to contain:

1. `TRANSCRIPT_LAYER_REQUIRED` fired when a transcript-bearing type attached to
   a `gene` **or `pseudogene`** with no transcript of its own.
2. `ALLOWED_PARENTS["exon"]` omitted `pseudogene`, although
   `ALLOWED_PARENTS["CDS"]` already listed it.

Measured on *Zea mays* B73 GCF_902167145.1: **5,195 pseudogenes and 12,757
exons on a gene-level parent.** The load stopped at line 247 of chromosome 1 --
`gene-LOC103630212`, an ordinary nuclear pseudogene with two exons and no
transcript.

RefSeq attaches exons directly to a `pseudogene` row by design. A pseudogene has
no functional transcript to model, so there is no missing layer being concealed
-- which is the only thing rule 1 exists to catch. The rule now applies to
`gene` alone, and the refusal message names the exemption so nobody has to read
the source to discover it.

**This is the third rule in cgload stated about the whole file when it should
have been stated about each parent.** The first was the allowlist of types
permitted to share an ID across segments (milestone 8c, AUGUSTUS
`sequence_feature`); the second was the strand vocabulary (D-045). Each was a
list of what is *permitted*, and each was wrong the first time real data
exceeded the fixtures. Worth auditing the remaining rules in `normalise.py` for
the same shape before the next unfamiliar genome finds a fourth.

## D-062 — RefSeq microRNAs are two transcript tiers deep

RefSeq annotates a microRNA as `gene -> primary_transcript -> miRNA -> exon`.
cgload's model has one transcript tier; `miRNA` was unmapped, so it loaded
inert, and **an inert parent makes its Tier A children unloadable** — the open
vocabulary is open in one direction only.

Measured on *Zea mays* B73 nuclear GFF3: **310 miRNA features under 164
`primary_transcript` rows** (1.9 each, the 5p/3p arm pattern), every one
carrying exons, every one refused. Predicted by audit before it was hit, and
confirmed by counting rather than by another failed load.

Two changes, one commit:

- `SYNONYMS["miRNA"] = "ncRNA"`
- `ALLOWED_PARENTS["ncRNA"]` gains `"mRNA"`

The second admits a nested transcript, which is a real relationship rather than
a modelling compromise: a precursor transcript genuinely contains a mature one.
Verified uniform — all 310 miRNA parents are `primary_transcript`; none attach
directly to a gene.

**D-060's standing rule was applied and then reverted, deliberately.** The rule
says a widened parent set should be matched by adding the child type to
`TRANSCRIPT_LAYER_REQUIRED`. Here that is wrong twice over: `ncRNA` is itself in
`TRANSCRIPT_TERMS`, so an ncRNA child makes its own gene count as
transcript-bearing and the check can never fire — a no-op that reads as
protection. And there is nothing to conceal: `gene -> ncRNA` is the normal shape
of a non-coding gene, where the ncRNA *is* the transcript rather than a child
that skipped one. Caught by a test asserting the guard fires, which it did not.

## D-063 — immunoglobulin gene segments are transcript-level

`SYNONYMS` mapped `V_gene_segment`, `C_gene_segment`, `D_gene_segment` and
`J_gene_segment` to canonical `gene`, while `ALLOWED_PARENTS["gene"]` is
`frozenset()`. The two tables contradicted each other: four terms declared
canonicalisable into a type that cannot appear in a parented position.

Settled from real data — RefSeq human GRCh38.p14, NC_000014.9, the IGH locus at
105,586,437–106,879,844:

- **175 segment rows, all carrying a `Parent`** (126 V, 27 D, 13 C, 9 J)
- every parent is a `gene`
- **164 CDS and 162 exon rows hang beneath them**
- zero segments without children

They occupy exactly the position an mRNA occupies, so they canonicalise to
`mRNA`. No other table changes: `mRNA` already permits a `gene` parent, and
`exon`/`CDS` already permit an `mRNA` parent.

## D-064 — UTR parents widened; predicted, not observed

`ALLOWED_PARENTS["five_prime_UTR"]` and `["three_prime_UTR"]` were `{mRNA}`.
This is the same shape D-036 fixed for `CDS` and D-058 for `exon`: when a
dialect has no transcript layer, the child attaches to the gene. The UTRs were
the remaining member of that family.

**No file in evidence attaches a UTR to a gene.** RefSeq human emitted no UTR
rows at all in the IGH region examined, so this is inferred from family
resemblance rather than measured — recorded as such so nobody later mistakes it
for an observed case. Widened anyway because it costs nothing and the forbidden
version is easy to state: a UTR may not attach to a `CDS`, to another UTR, or to
a `repeat_region`.

Both sets become `{mRNA, gene, pseudogene}`.

## D-065 — NCBI feature counts must not be summed across assembly units

A fourth trap in the feature-count format, and the one that produced a false
failure on a correct load.

A file may carry several groupings under one `Full Assembly`: an `all` rollup
plus one row-set per assembly unit. *Zea mays* GCF_902167145.1 has **24 `all`
rows, 23 `Primary Assembly` rows and 8 `non-nuclear` rows**, and the first is
the sum of the other two.

`parse_ncbi_feature_counts` summed every row matching `Full Assembly`, so the
rollup was added to its own components. Measured: mRNA is 57,071 in both `all`
and `Primary Assembly`, and the naive sum reported **114,142 against a correct
load of 57,071** — a 100% overstatement, presented as cgload having lost half
the genome.

Missed because the *U. maydis* file the parser was written against carries one
unit grouping per feature, where summing is safe. The same shape as every other
defect in this project: a rule derived from the fixtures on hand.

**Fix.** Group by assembly unit and never sum across groups. The default takes
the `all` rollup when present; `assembly_unit` selects a specific unit, which is
what a partial load needs. With no rollup and several units, refuse and name
them rather than guessing.

After the fix, against the `Primary Assembly` rows — the unit actually loaded,
the two organellar sequences having been excluded — mRNA, CDS, rRNA and tRNA all
match the loaded counts exactly.

## D-066 — a pseudogene is reported under `gene` by every external reference

Both references cgload compares against fold pseudogenes into `gene`:

- **GenBank** writes a pseudogene as a `gene` record carrying `/pseudo`.
- **NCBI's feature-count file** reports `Feature=gene, Class=pseudogene`.

cgload keeps `pseudogene` as a distinct `source_type`, because that is what the
GFF3 said and D-035 requires the source term to survive.

The cross-parser tier already applied this fold. The reference tier did not, so a
correct maize load reported a 5,195-gene shortfall — exactly its pseudogene
count. Every other type matched exactly.

Measured twice, exact both times:

| | gene | pseudogene | sum | published |
|---|---|---|---|---|
| *F. graminearum* GFF3 vs GenBank | 13,714 | 11 | 13,725 | 13,725 |
| *Z. mays* Primary Assembly vs NCBI | 44,321 | 5,195 | 49,516 | 49,516 |

`GENBANK_TYPE_EQUIVALENCE` is renamed `EXTERNAL_TYPE_EQUIVALENCE` — the mapping
was never GenBank-specific — with the old name kept as an alias. Both tiers now
use it.

**The fold cannot mask a real loss.** It moves counts between type names within
one side of the comparison; it does not change the total. Tested: removing 100
genes still fails with the exact delta.

---

## D-067 — The reference tier's precondition may not default

**Decided:** `--reference-counts` requires `--sequences-in-assembly`. It is not
defaulted, and the tier refuses to run without it.

The code read
`sequences_in_assembly=sequences_in_assembly or sequences_loaded(...)`, so
omitting the flag made the precondition `loaded == loaded` — always true. The tier
could therefore never skip, and a partial load was compared against totals that
may describe the whole assembly: it would pass by luck, or fail for a reason
cgload already classifies as innocent.

That is D-033 inverted. D-033 says a check that vanishes when it cannot run is
indistinguishable from one that passed, which is why every tier is always
reported. The same argument applies to a precondition: one that cannot fail is
indistinguishable from one that is not being evaluated.

**The legitimate case has to be asserted, not inferred.** Comparing against counts
that describe exactly the unit loaded is valid — that is precisely what
`--assembly-unit` selection exists for (D-065). But cgload cannot distinguish
"these totals describe my two sequences" from "these totals describe all
twenty-three and I loaded two", so the user states it by passing the count
explicitly. Consistent with `--accession`, required on the same tier for the same
reason: guessing wrong there selects nothing and reports every feature missing.

**Found by reading a passing maize run, not a failing one.** The tier reported
PASS, which was the right answer for that run — 23 Primary Assembly sequences
against Primary Assembly totals — arrived at by logic that would have reported
PASS regardless of what was loaded.

---

## Numbering note

D-059 is used twice in this log: once for the `has_transcripts` predicate scope
and once, in an earlier draft, for this precondition rule. The precondition is
**D-067**; the `has_transcripts` entry keeps D-059. Recorded because a decision
log whose numbers are ambiguous is worse than one with a gap.

---

## D-068 — A trans-spliced gene spans both strands, and says so with `?`

**Decided:** `Profile.strand_exception_qualifier`. RefSeq declares `exception`;
every other profile declares nothing. When two rows share an ID, differ in
strand, and carry `exception=trans-splicing`, they merge into one feature whose
strand becomes `?`. Without the qualifier a strand disagreement is still two
features sharing an ID, and is still refused.

**Third form of the same provider defect.** After 34 too-small envelopes (D-046)
and one inverted row (D-056), NCBI's GFF3 converter writes a trans-spliced gene as
several rows sharing one ID with differing strands — `gene-ZeamMp186` (*nad1*) at
lines 990,466 and 990,469 of the maize annotation, `+` then `-`.

The rows are correct. A trans-spliced gene *is* assembled from segments on
opposite strands, so no single strand describes the feature and GFF3's answer is
`?` — which NCBI itself writes on the mRNA row of that same gene, but not on the
gene rows. Taking `?` on merge makes cgload's stored value the honest one and
agrees with the provider's own annotation of the transcript.

**Matched on the value, not the key.** `exception=` also covers `ribosomal
slippage`, `rearrangement required for product` and `annotated by transcript or
proteomic data`. Only trans-splicing licenses a strand disagreement, so a
key-only match would silently merge genes flagged for something unrelated — the
permitted-set failure of D-060's register, arriving inside a fix for it.

**Found only by including the organelles.** The earlier scale run loaded the 23
Primary Assembly sequences and excluded the mitochondrion and chloroplast, so this
row was never read. Every trans-splicing defect in this project lives on those two
sequences, which is an argument for running the full assembly at least once
rather than the nuclear unit alone.

## D-069 — a feature that declares no parent is top-level, whatever its type

GFF3 nowhere requires an `exon` to be parented; a feature with no `Parent`
attribute is top-level, and that is legal for any type. cgload refused
parentless `exon` rows because `ALLOWED_PARENTS["exon"]` is non-empty.

The rule conflated two different things:

- **a feature naming a parent it may not have** — an error, because the file
  asserts a relationship the model cannot hold;
- **a feature naming no parent at all** — not an error, because the file asserts
  nothing.

Only the first hides something. D-058's case had a parent present with a layer
missing between them; here the file never claimed a parent, so storing the row
top-level records exactly what the file says.

**Found on *Zea mays* B73 GCF_902167145.1, line 992,090.** Two chloroplast exons
flanking a parentless intron, all three carrying `number=2` and
coordinate-derived identifiers (`ID=id-NC_001666.2:129636..129867-2`). They are
the second exon block of a trans-spliced chloroplast gene whose first exon lies
elsewhere; NCBI's converter could not express the link and emitted them
unparented. **Two exons in 1.2 M features, both organellar** — the fourth
manifestation of the same converter defect, after 34 too-small envelopes
(D-046), one inverted row (D-056) and the both-strand gene rows (D-068).

Naming a disallowed parent is still refused, unchanged.

**Eighth instance of D-060's class.** Every one has been a rule stated over a
permitted set, derived from the fixtures then on hand. This one differs in shape:
the permitted set was correct, and the defect was applying it to a case the set
was never about.

## D-070 — sum the component units; the rollup can be incomplete

D-065 established that rows must not be summed across assembly units, because
an `all` rollup already contains its components. The corollary, found on the
same file: **when components exist, sum them rather than reading the rollup.**

Measured on GCF_902167145.1, rRNA:

| unit | Unique Ids | Placements |
|---|---|---|
| all | 2427 | 2439 |
| Primary Assembly | 2427 | 2427 |
| non-nuclear | **na** | 12 |

The rollup's `Unique Ids` total silently omits the 12 organellar rRNAs, because
NCBI computed no unique count for that component. 2,427 + 12 = 2,439, which is
what the GFF3 contains and what cgload loaded. Summing the components recovers
the right answer with no per-column special-casing.

This is not D-065's double-count: components are summed with each other, never
with the rollup that already contains them.

## D-071 — NCBI's published counts omit organellar mRNA (their defect, not ours)

After D-070, one disagreement survives on *Zea mays* B73 and it cannot be fixed
in cgload, because the reference is wrong.

The feature-count file for GCF_902167145.1 carries a `CDS non-nuclear 274` row
and an `rRNA non-nuclear 12` row, but **no `mRNA non-nuclear` row at all**. It
therefore reports 274 organellar coding sequences and zero organellar
transcripts.

NCBI's own GFF3 for the same assembly contains **274 organellar mRNA rows** —
counted directly:

```console
$ awk -F'\t' '!/^#/ && $3=="mRNA" && ($1=="NC_007982.1" || $1=="NC_001666.2")' \
      genomic.gff | wc -l
274
```

57,071 nuclear + 274 organellar = 57,345, which is what cgload loaded and what
the source-count tier independently confirms is in the file.

**Consequence for the tiers.** The reference tier is correctly reporting a
disagreement; the disagreement is between two NCBI products. This is why D-033
ranks the tiers as it does — the cross-parser check compares two *primary*
representations of an assembly, while the reference check compares against a
derived summary that can be wrong in ways the primary data is not.

The fifth defect found in one provider's description of one assembly, after 34
too-small envelopes, an inverted row, both-strand gene rows, and unparented
exons. Reported to NCBI alongside the others.

## D-072 — the comparison tiers compare the union, not the intersection

`compare_across_formats` and `compare_with_reference` both narrowed to
`set(expected) & set(actual)` before comparing. A type present on one side and
absent from the other was silently dropped.

That made **the most alarming failure the one case that could not fail**: if
every feature of a type was gone, the type vanished from the comparison. Losing
321 of 322 tRNAs failed; losing all 322 passed. Tier 1 always used a union;
these two inverted it.

Both now compare the union with a missing type counted as zero. The
`structural_only` allowlist already carried the types that legitimately exist in
one format and not the other, so the intersection was doing nothing the
allowlist was not already doing — while suppressing total loss.

## D-073 — `--force` moves the superseded assembly aside, and deletes its rows

Two defects that combined into data corruption produced by the tool's own flag.

`staged_load` ran `shutil.rmtree(target)` **before** opening the transaction, and
no exception path undid it. `_reject_duplicate_assembly` returned silently under
`force`, leaving the previous assembly's rows in place.

So a `--force` load that failed for any reason left an `assembly` row whose
`fasta_path` pointed at a directory that no longer existed; a `--force` load that
succeeded left two rows sharing one directory, one of them describing
coordinates against a FASTA that had been replaced. `fasta_sha256` is recorded
but never re-checked, so nothing downstream could detect either state.

Now: the existing directory is renamed aside within the same filesystem, restored
if anything raises, and removed only after the transaction commits. The
superseded `assembly` rows are deleted inside that same transaction, so a
rollback restores both the rows and the directory together.

## D-074 — parent edges are counted, and an absent hierarchy fails

Every tier counted features; none counted edges. The GenBank containment linker
is the most novel and most fallible component here, and both linkers fail
quietly: `_link_by_containment` skips an orphan CDS, `_link_by_locus_tag` emits
no edge for a child without a `locus_tag`.

**A GenBank load that produced zero parent edges passed every check and reported
PASS.** Every feature was present and counted; only the structure connecting
them was missing, and cardinality cannot see that.

`check_parent_edges` counts edges by the rule that produced them and fails when
any feature of a parent-taking type is unlinked. For GenBank it additionally
fails on any edge recorded as `explicit_parent`: the format carries no parent
pointers, so such an edge means the tokenizer's synthesised attribute was read
back as the file's own — the linkage-provenance bug arriving from the other
direction.

The `linkage_method` column has been stored since that fix. Nothing asserted on
it until now.

## D-075 — the claim is completeness, not fidelity

Tier 1 compares distinct `(type, ID)` counts from an independent pass against
`COUNT(DISTINCT source_id)` in the database. The implementations share no code,
but they encode the same rule, so what they establish is that **no feature was
dropped or invented** — cardinality.

They cannot see a count-preserving corruption: a wrong strand, a swapped
coordinate, a mangled identifier, a dropped segment of a multi-segment feature.
(The last is caught since D-093, which counts rows as well as features.)
`check_declared_proteins` catches the last of these only where the source
declares a protein length, which is AUGUSTUS and funannotate, not RefSeq.

"Prove that nothing was lost" reads as the stronger claim. The wording in the
README and the paper is narrowed to what is actually verified rather than the
code being extended to meet it, because fidelity cannot be checked cheaply
without a second independent derivation of the same values.

Tiers 3 and 4 are also reclassified in the documentation. They compare two
descriptions of an assembly to each other rather than the database to its input,
so they are consistency checks on the provider's outputs — genuinely useful, and
how four defects in a published reference annotation were found, but a failure
there does not imply a cgload defect.

## D-076 — the package name is `cgload`, and a rename can corrupt its own documentation

`cgload` — comparative genomics load. Two earlier names were rejected after
collision checks found them taken late: `tessera` (a tissue-segmentation
algorithm on bioRxiv) and `sextant` (PyPI 0.5.0, a Helm CLI, whose console
script also collides).

**A mechanical rename rewrote the document describing the collision.**
`docs/NAME.md` recorded that `sextant` was taken by a Helm CLI. `rename.py`
substituted every occurrence of the old name, including inside that prose,
producing a document asserting that the *new* name was taken by a Helm CLI —
directly contradicting `pyproject.toml` and the README two directories away.

The script's own verification could not catch it: no occurrence of the old name
survived, and the cross-module constants agreed. Both checks passed on a tree
containing a false claim.

`NAME.md` is merged into `docs/RENAMING.md`, which now carries the collision
history with the rejected names written as data rather than as prose that a
substitution can rewrite, plus a warning to read that file by hand after any
rename rather than trusting the script's output.

Related: the substitution also produced `Cgload` wherever the old name had been
capitalised mid-sentence. Harmless, but it is the same class — a rename knows
about strings, not about the language around them.

## D-077 — every row records its input line, so anonymous rows keep their parent

**Found by audit (AUDIT-2).** `_resolve_parent_ids` and `_insert_attributes`
found the rows they had just bulk-inserted by `source_id`, and skipped any
feature that had none. An anonymous row (D-024) therefore lost the parent its
`Parent=` attribute stated, while still being stored with
`linkage_method = 'explicit_parent'`: 175 of 322 features in the AUGUSTUS
fixture. Its attributes were dropped the same way. No tier counted edges, so the
load reported success (D-074 covered the counting; nothing wired it in).

**Decided:** `feature.source_line` records the input line each row came from.
The loader maps rows back by `(source_line, segment_index, start, end)`, the
span included because the segments of an anonymous GenBank `join()` share one
line, and links and annotates named and anonymous features the same way.

**Rejected: relying on insertion order** to pair anonymous rows with their
generated keys. It holds on SQLite and in practice on PostgreSQL and MySQL, but
no backend promises it, and a wrong pairing would attach rows to the wrong
parent silently.

**Cost:** one nullable integer per row, and a schema change. A database created
before it is refused by every command with a message naming the missing column
(`REQUIRED_COLUMNS` in `cgload.db.init`), rather than failing on the first insert
with a driver error that names no cause.

## D-078 — the parent-edge tier runs, in `verify` and inside every GFF3 load

**Found by audit (AUDIT-1).** `check_parent_edges` (D-074) was written and unit
tested, and nothing called it. `verify` printed three tiers where the README
listed four, and a database whose `parent_id` column had been set to NULL
throughout still verified PASS.

**Decided:** `verify` reports a parent-edge tier for every annotation file,
computed from the database alone, so it runs even when the input file is gone.

- **GFF3 states its own parents.** Every feature the file gave a parent must be
  linked, whatever its type: the loader records the statement as
  `linkage_method`, so a stated parent missing from the database is a loss by
  cgload. Counting only `LINKED_TYPES` would have missed AUDIT-2, whose unlinked
  rows were introns and codons. A feature the file left parentless is top-level
  by its own statement (D-069) and is reported, not failed.
- **GenBank edges are computed.** Every feature of a `LINKED_TYPES` type is
  expected to have one, as D-074 decided.

The tier also runs inside every GFF3 load, next to the source count (D-046),
and a failure rolls the load back: a test that drops the anonymous rows'
parents, the AUDIT-2 bug, now fails the load instead of storing an unlinked
hierarchy. For GenBank the load leaves the tier to `verify`, because an
unlinked GenBank feature need not be cgload's fault (D-033).

## D-079 — `verify` checks the stored assembly; the assembly row records the data root

**Found by audit (AUDIT-4).** `fasta_sha256` was recorded at load time (D-018)
and never read again. Deleting or replacing the stored FASTA left a database
whose every coordinate referred to a file that was gone or different, and
`verify` still reported PASS.

**Decided:** `verify` has a "stored assembly" tier. It fails when the stored
FASTA or its `.fai` is missing, or when the FASTA's SHA-256 differs from the one
recorded at load. To find the file, `assembly.data_root` records where the data
directory was at load time, absolute. `--data-dir` overrides it after the
directory moves, and `fasta_path` stays relative, so moving the directory still
works.

**This bends the comment on `fasta_path`** that the database should hold no
absolute paths. `source_file.path` already held absolute input paths, which the
source-count tier reads; the alternative, making the user always name the data
directory, would have meant the tier silently skipped in every ordinary run,
which is D-033's failure.

**Cost:** hashing the stored FASTA. About as long as reading it once: a few
seconds for a 2 Gb assembly, where the other tiers take under a second.

## D-080 — on an organelle sequence a CDS may attach straight to its gene

**Found by audit (AUDIT-5).** NCBI writes a vertebrate mitochondrial CDS
directly under its gene, with no mRNA between: chimpanzee GCF_000001515.7, ND1,
and every other vertebrate RefSeq assembly that includes its mitochondrion. The
transcript-bypass rule (D-036, scoped per gene by D-043) refuses a gene with no
transcript whose CDS attaches to it, in a file that has a transcript layer
elsewhere. A vertebrate file has one in its nuclear genes, so the whole file was
refused.

**Decided:** the rule does not apply on a sequence whose `region` row declares an
organelle genome (`genome=mitochondrion`, `plastid`, `chloroplast` and the other
values in `ORGANELLE_GENOMES`). Organelle genomes are annotated the prokaryotic
way. The exemption is per sequence: a nuclear gene missing its transcript in the
same file is still refused, and a test pins that.

**Not used: the `Is_circular` flag.** Some nuclear-like replicons are circular and
some organelle assemblies are linear; `genome=` says what the sequence is.

**Known gap:** a file with no `region` rows (funannotate, AUGUSTUS) declares no
organelles, so the rule applies to all of its sequences, as before.

`ORGANELLE_GENOMES` decides whether a `gene`->`CDS` edge is accepted, so it is
part of the vocabulary hash (D-061), and the hash changed. Databases created
before this change are already refused for lacking `feature.source_line`
(D-077), so no further migration step is added.

## D-081 — a row must lie within its sequence, and a child on its parent's sequence

**Found by audit (AUDIT-6).** Two coordinate checks were never made. A row
ending at 20,800 on a 20,000 bp sequence loaded without a word, so the database
held a coordinate that pointed past the end of the stored FASTA. And the
envelope check compared numbers only, so an exon on one sequence was linked to
an mRNA on another.

**Decided:**

- The FASTA is authoritative for length (D-004). A row past the end of a linear
  sequence is refused, naming the rows.
- On a sequence whose `region` row says `Is_circular=true`, a row that starts
  inside the sequence and ends past it, by less than one length, crosses the
  origin. GFF3 allows that and it loads silently.
- A row lying wholly past the end of a circular sequence is not GFF3-legal
  either. NCBI writes them for trans-spliced organelle genes (Arabidopsis
  *nad5*, rice *rps12*; the RefSeq survey's class C rows), and on a circle the
  position is still readable, so they are stored as written and named in a
  warning. Refusing them would refuse those two reference genomes outright.
- A child on a different sequence from its parent is refused. The RefSeq survey
  found none in 25.5 million parent links, so this costs no real file.

**Test fixtures changed with it.** Several tests built FASTAs shorter than their
annotation said the sequence was: the *F. graminearum* fixture's `region` row
declares the chromosome's full 8,033,942 bp. Those FASTAs now match the
annotation, which is what a real download has.

## D-082 — the data directory gets the annotation as BGZF plus a tabix index, written in Python

**Found by audit (AUDIT-10).** The README promised `annotation.gff3.gz` and its
`.tbi` beside the assembly, so that JBrowse 2 could read cgload's output with no
code from us, and nothing wrote either file.

**Decided:** every load writes the normalised annotation as cgload-dialect GFF3
(D-035), its rows sorted for tabix, compressed as BGZF and indexed as `.tbi`,
inside the staging directory so it lands with the rows or not at all. Both
formats are written by `cgload.bgzf` with `zlib` and `struct`.

**Rejected:**

- *Calling `bgzip` and `tabix`.* Two more programs to install, and a
  subprocess per load, which is where D-005 says genome tools acquire their
  attack surface.
- *pysam.* A compiled dependency with no Windows wheels, against the
  dependency budget, for two short published formats.

**Measured:** on 75,630 rows over three sequences of up to 400 Mb, 3,000 random
region queries through htslib (pysam 0.24.1) using cgload's index returned
exactly the rows a brute-force overlap scan returned, as did htslib's own index
over the same file. The test suite repeats this with a reader written from the
specification, and with htslib where pysam is installed.

**Limit:** a `.tbi` index addresses positions below 2^29. Past that the file is
still written, the index is not, and the load says so; `tabix --csi -p gff`
indexes it. No sequence in the 30-genome RefSeq panel comes near that size.

## D-083 — deleting an assembly deletes every row it owns, children first

**Found while building `remove` (Phase 2).** `load --force` replaced an
assembly by deleting its `assembly` row and nothing else. Foreign keys are
enforced, on SQLite too (`build_engine` turns them on), so with any features
stored that delete was refused with "FOREIGN KEY constraint failed" and the load
stopped with a traceback. The transaction rolled back, so nothing was damaged,
but `--force` had never replaced an assembly that had features. No test loaded
the same assembly twice through the loader. The older `register_assembly` path
did not delete the old row at all, so with no version set it left two assembly
rows sharing one data directory.

**Decided:** one function, `cgload.db.removal.delete_assembly_rows`, deletes
everything an assembly owns, inside the caller's transaction, in dependency
order: functional annotation and its runs, feature attributes, the parent links
among the assembly's own features, the features, sequence regions and source
files, then the assembly row. `load --force`, `register_assembly(force=True)`
and `remove` all use it.

- Parent links are cleared before the features are deleted. SQLite and
  PostgreSQL check a foreign key at the end of a statement, so deleting parents
  and children in one statement would pass there, but MySQL's InnoDB checks row
  by row and refuses to delete a parent before its children.
- Only rows the assembly owns are touched. A row of another assembly pointing
  at this one should not exist; if one did, the foreign key error is the right
  outcome, not a cascade into someone else's data.
- Every statement selects its rows with a subquery, not an IN list built in
  Python (AUDIT-8).

## D-084 — `remove` takes the rows and the data directory together, and only a directory that proves it is the assembly's

**The stub said** "transactional removal of an assembly". Rows alone would
leave the copied FASTA, its index and the annotation behind, and files alone
would leave rows pointing at nothing, which is the state AUDIT-4 was about.

**Decided:**

- One transaction for both, built like `staged_load` in reverse
  (`storage.staged_removal`). The rows are deleted (D-083) and counted again;
  anything left stops the removal. Then, inside the still-open transaction and
  immediately before the commit, the data directory is renamed into `.staging`,
  and it is renamed back if the commit fails. It is deleted only after the
  commit. A process killed between the rename and the commit leaves the rows in
  place and the files under `.staging`; `verify` then reports the stored
  assembly as missing, so the failure is visible.
- The directory is found from what the database recorded at load: the data
  root (D-079) and the FASTA's relative path. `--data-dir` overrides the root
  after a move, as it does for `verify`.
- Nothing is deleted on a guess. The directory must lie inside the data root,
  and its `manifest.json` must name the same organism, strain and assembly and
  record the same FASTA checksum as the database. Otherwise `remove` stops and
  changes nothing. A directory that is not where the database says also stops
  it; `--keep-files` removes the rows alone, on purpose.
- The strain and then the organism are deleted when no other assembly or
  strain uses them, so removing an organism's only assembly leaves the database
  as it was before the load. Empty organism and strain directories go too.
- It asks before deleting unless `--yes` is given, and `--dry-run` prints the
  counts and the directory without changing anything.

## D-085 — `query` answers three questions with whole features, through bound parameters only

**The stub promised** "region, annotation and KO queries". The paper's timing
comparison (Phase 7) needs the everyday questions answered by cgload itself: a
gene with its children, the features in a region, the genes carrying a domain.

**Decided:**

- Three questions, exactly one per call: `--id` (with `--descendants` for
  everything below it), `--region` in samtools notation, and `--term` for an
  accession from eggNOG-mapper or InterProScan. Answers are whole features,
  every segment included.
- A region matches rows, not envelopes. A position inside an intron returns
  the gene and the mRNA, not the CDS whose rows all lie elsewhere; a CDS with
  any row in the region comes back whole. A sequence name containing a colon
  is tried whole first, as samtools does. An end past the sequence is cut to
  its length; a start past it is an error.
- On a circular sequence a region also matches rows that wrap round the
  origin. Such a row is stored as GFF3 writes it, ending past the sequence's
  length (D-081), so it covers position p at p + length, and the region is
  tried a second time one length further on. On a linear sequence that second
  test never holds, because the loader refuses rows past the end.
- `--descendants` starts from every row of each matched feature, not only from
  the rows that matched. A child's `parent_id` names its parent's first row,
  and a region can overlap only a later row of a trans-spliced mRNA; starting
  from the matched rows alone returned that gene model without the exons that
  hang from the first row.
- Accessions match exactly. The only rewrites are spellings of one identifier
  that the source formats are known to differ on: `K00799` and `ko:K00799`,
  KEGG `map` and `ko` pathways (D-032), and GO and EC prefixes. Nothing is
  matched by prefix or similarity. GO matches include the ancestor terms
  eggNOG-mapper adds (D-045), and the TSV says so when GO rows are present.
- `--type` keeps only features of the given canonical types. It is refused with
  `--term`, where the annotated feature is the CDS or mRNA the functional file
  named; the `top_level` column names its gene instead.
- TSV has one row per feature (one per assignment with `--term`), with the
  parent and the identifier at the top of the parent chain. GFF3 adds every
  ancestor, so each `Parent` resolves inside the file, and it carries the
  `##cgload-export` header (D-035) plus a comment naming the question.
- Every value the user types is a bound parameter. Sets of rows pass between
  steps as subqueries, never as IN lists built in Python (AUDIT-8); a test runs
  a whole-sequence query under a parameter limit of 10. The parent chain is
  walked with recursive CTEs using UNION, not UNION ALL, so a cycle in the
  parent links cannot make a query run forever.
- Read-only: nothing is committed, and a test checks that the database file is
  byte-identical after queries.
- Exit codes: 0 for an answer, including an empty one; 1 when the ID asked for
  does not exist, since that is usually a typo; 2 for a malformed question.
- `protein_annotation` gains an index on `accession` alone. The existing index
  leads with `analysis`, so it cannot serve a term asked without `--analysis`.

**Not done:** no server and no web layer (roadmap). Query speed on large
genomes is measured in Phase 7, against gffutils and against tabix on cgload's
own files, before anything is claimed about it.

## D-086 — `verify` recounts the input only if it is still the file that was loaded

**Found while checking the tutorial.** The source-count tier recounts the
original input file at the path recorded at load. It never checked that the
file there was still the one loaded, although its SHA-256 is stored. A gene
appended to the file after the load was reported as a feature "lost or invented
during loading", which blames the loader for an edit, and an edit that kept
every count the same, such as a moved coordinate, passed.

**Decided:** before counting, `verify` compares the file's SHA-256 with the one
recorded at load. If they differ, the tier is skipped with a reason that names
the file and says how to recover: restore the original, or reload with
`--force`. A skipped tier is reported, never hidden (D-033). The cost is one
extra read of each input file, next to the read the count already makes.

## D-087 — `load --strict-vocabulary` refuses the types the load summary names

**Found while checking messages against the commands.** `docs/schema.md`, the
vocabulary module and D-017 all describe a `--strict-vocabulary` flag that turns
unknown feature types into an error. No command had it. A test now reads every
string in the package and fails if one names a flag that no command declares.

**Decided:** `load --strict-vocabulary` refuses a file containing feature types
that its dialect profile does not declare: the types the load summary already
names in its "not declared by the ... profile" warning. Every such type in every
file is listed in one refusal, as item 5 of the D-060 audit asked, and
the database is not touched. Without the flag nothing changes: such types load
verbatim with a warning.

**Rejected:** refusing every term absent from the vocabulary table, as D-017
first described. `region`, `intron`, `start_codon` and `stop_codon` are kept
out of the vocabulary on purpose (`NOT_CANONICALISED`), so that rule would
refuse every RefSeq and every AUGUSTUS file, and the flag would be switched
off after the first try. The profile already records which types each producer
is expected to write, so a type outside it is the one worth stopping for.

## D-088 — the INSDC ncRNA classes map to `ncRNA`

**Found while reading why the human load could fail.** A Tier A feature whose
parent loads inert is refused (D-062), and an `exon` may attach only to a
transcript, a gene or a pseudogene. NCBI writes each ncRNA class from the
INSDC `/ncRNA_class` vocabulary as its own column-3 term. cgload mapped some of
them (snRNA, snoRNA, lnc_RNA, SRP_RNA, RNase_P_RNA, antisense_RNA, guide_RNA,
telomerase_RNA) but not the rest, so a file containing them was refused at
the first such exon: "a exon may only attach to ['gene', 'mRNA', 'ncRNA',
...]". The human RNY, VTRNA and RMRP genes belong to the Y_RNA, vault_RNA and
RNase_MRP_RNA classes. Whether this is what stopped the human load on the lab
server on 15 October is not yet confirmed. D-062 fixed miRNA the same way, one
class at a time; this adds the rest at once.

**Decided:** `Y_RNA`, `vault_RNA`, `RNase_MRP_RNA`, `scRNA`, `scaRNA`, `piRNA`,
`siRNA`, `rasiRNA`, `ribozyme` and `hammerhead_ribozyme` are synonyms of
`ncRNA`. The vocabulary hash changes, so databases made before this are refused
with a message saying why.

**Left out:** `pre_miRNA`, which would put a miRNA under an ncRNA, which the
hierarchy does not permit (RefSeq writes `primary_transcript` there instead),
and `autocatalytically_spliced_intron`, an intron rather than a transcript.

## D-089 — RefSeq GFF3 may overlap a CDS where it says `ribosomal slippage`

**Found on the lab server, 15 October 2026.** The human RefSeq file
(GCF_000001405.40) was refused at line 1,838,264: `cds-NP_001165908.1` is
written as two rows, 94,663,334-94,664,513 and 94,664,513-94,665,682, which
share one base. D-033 had permitted overlapping segments only for GenBank's
`/ribosomal_slippage`, and declared overlap "always an error, which is correct
for GFF3". It is not: NCBI writes the same frameshift in GFF3, the base read
twice appearing in both rows.

**Decided:** `Profile.overlap_exception_value`. The RefSeq and cgload profiles
permit the segments of one feature to overlap when the row's `exception`
attribute says `ribosomal slippage`. Matched on the value, as D-068 matches
`trans-splicing`, because `exception` also carries RNA editing, trans-splicing
and other annotations that explain nothing about an overlap. An overlap without
it is still refused.

**To confirm on the next server run:** that the human rows carry this
`exception` value. If they do not, the load is still refused, and the run
prints the refused rows so the next change can be made from them.

## D-090 — a refusal lists every problem its pass found

**Found on the lab server, 15 October 2026.** Normalisation stopped at the
first problem. On the human RefSeq file that meant two minutes of parsing per
problem, and a server run, and a message back, for each one.

**Decided:** the three normalisation passes (grouping rows into features,
resolving parents, checking the hierarchy) each finish and then refuse with
every problem they found: the first ten in full and a count of the rest. A
single problem reads exactly as before. The passes still run in order, because
each assumes the one before it succeeded, so a file with problems in two passes
needs two runs.

## D-091 — an envelope disagreement that exists only across a circle's origin is stored as written

**Found on the lab server, 15 October 2026.** *Arabidopsis thaliana*
(GCF_000001735.4) was refused even with `--repair-envelopes`. On its
mitochondrion (NC_037304.1, 367,808 bp, `Is_circular=true`), NCBI writes two
trans-spliced genes in plain coordinates, each on several rows, and their
mRNAs past the origin with the sequence length added: rna-DA397_mgp37 is
234,352..396,144 under a gene spanning 26,817..360,128. Both describe the same
bases, but compared as numbers the mRNA lies outside its gene, and a parent
written on several rows is never recomputed (D-046).

**Decided:** with `--repair-envelopes`, a disagreement that remains after
repair is stored as written when the sequence's region row says
`Is_circular=true`, the parent or the child runs past the sequence's length,
and one of them carries `exception=trans-splicing`. Each one is named in a
warning and recorded in `manifest.json` under `envelope_exceptions`. Anything
else still refuses, and without the flag nothing changes.

**Rejected:** recomputing the gene, which means choosing which of its rows to
widen, the guess D-046 forbids; and comparing positions modulo the sequence
length. A parent and a child written in different conventions can describe
different arcs of the circle (rna-DA397_mgp31 is written the long way round,
161,832..460,814), so arithmetic alone does not settle containment, and
choosing the arc is another guess.

## D-092 — the parent-edge tier counts the file's `Parent=` statements itself

**Found while writing the fault-injection protocol, 15 October 2026.** The
tier read only the database. Whether a feature had a parent "in the file" was
taken from `linkage_method`, which the parser writes. A parser that dropped a
`Parent` attribute therefore stored the feature as top-level, and the tier
reported it as top-level by the file's own statement (D-069). Checked on the
code before this change: removing one exon's link and its `linkage_method` from
a loaded database left `verify` at PASS, with "1 more named no parent and are
top-level". The source count could not see it either, because the exon was
still there. D-008's rule, that a check must not share the parser's view of the
file, held for the source count and not for this tier.

**Decided:** the naive pass, in the same single read, also counts per column-3
type the features whose rows say `Parent=`, by the same rule as the feature
count: distinct (type, ID), rows with no ID one by one. It reads no parent
names and resolves nothing. For a GFF3 file the tier compares those counts with
the features stored with an explicit parent link, per type, and fails on any
difference. `load` runs it inside the transaction, as before. `verify` runs it
when the file is still the one loaded; when the file is missing or has changed,
the tier checks the database alone and says so. GenBank states no parents, so
nothing changes there.

**Still not detected, by design (D-075):** a link moved to the wrong parent of
the same type, which keeps every count. The tier shows that every stated
parent is stored as a link, not that each link points where the file says.

## D-093 — the source count counts rows as well as features

**Found while testing the fault-injection scripts, 16 October 2026.** The
naive pass counts features: rows sharing an ID are one. A CDS written on six
rows that lost one row during parsing is still one CDS, so the source count
passed, and so did every other tier. On the synthetic test genome a parser
that dropped one row of a multi-row CDS loaded with a PASS. For a protein-coding
gene that is a missing exon's worth of coding sequence, and nothing said so.

**Decided:** the same single read also counts data rows per column-3 type,
every row counted, and the source-count tier compares them with the rows
stored for that type, after the feature counts agree. cgload stores one row
per segment and refuses a file in which rows repeat a span under one ID, so
for any file it accepts the two numbers must be equal; the full suite, every
real fixture included, agrees. A difference fails the tier and, inside a
load, rolls it back. GenBank is unchanged: a GenBank record is one feature
however many segments its location joins, and its rows are not lines of the
file.

## D-094 — region and gene queries read only the rows they need

**Found while testing the query-timing script, 16 October 2026.** On a
synthetic genome of 120,018 features, one gene and its descendants took 0.28 s
on SQLite and 3.7 s on PostgreSQL; a 1 Mb window took 0.5 s and 1.8 s. The
plans showed four causes:

1. "Every row of the matched features" was an OR of two IN tests. No planner
   turns that into index lookups; PostgreSQL scanned the feature table.
2. The top of each feature's parent chain was found with a recursive CTE that
   PostgreSQL planned as a hash join over the whole table at every step.
3. A region on a linear sequence was also tested one sequence length further
   on, for rows crossing a circular origin. That second test,
   `start <= end + length`, is true of every row, so the query read the whole
   sequence. The index on (sequence, start, end) bounded a region only from
   above in any case, and PostgreSQL, unable to estimate how few rows overlap,
   chose table scans.
4. Right after a load, PostgreSQL had no statistics for the new rows.

**Decided:** whole features are selected by one IN over a UNION. The top of
each chain is walked in Python over the rows already fetched, looking up only
ancestors outside them, 500 at a time (AUDIT-8's bound). `sequence_region` records
whether the file marks a sequence circular, and only those get the second test.
Each row stores the smallest genomic bin that holds it, in the hierarchical
scheme of the SAM and CSI indexes (16 kb leaves, eight-fold steps, coordinates
up to 2^38), indexed with its sequence; a region query matches the few dozen
bins that overlap it, or reads the whole sequence when it would need more than
500. A PostgreSQL load runs ANALYZE on its tables before it commits.

On the same genome: a gene with its descendants 0.006 s on SQLite and 0.08 s
on PostgreSQL; the 1 Mb window, 3,005 features, 0.19 s and 0.30 s, most of it
now spent building the answer in Python.

**Cost:** two columns and an index, so databases made before this are refused
with the message REQUIRED_COLUMNS gives.

## D-095 — a directory left by an interrupted load does not block the next load

**Found by the benchmark's kill trials, 16 October 2026.** A load renames its
directory into place inside the open transaction and then commits (D-019's
ordering). One trial in ten on a small genome was killed between the two. The
database was empty, as designed, but the directory stayed, and the repeated
load was refused: "already exists. Use a different assembly name, or --force
to replace it." D-019 called this the harmless direction, and it is harmless
to the data, but a user who reruns a killed load gets a refusal that names
neither the cause nor the fact that nothing uses the directory.

**Decided:** `init` gives each database an identity, a UUID in
`vocabulary_meta.database_id`, and every load records it in the directory's
`manifest.json`. When a load finds its directory already there, it replaces
it, as `--force` would, only if both hold: the manifest names this database,
and no assembly row of this database keeps its files in that directory. Then
no committed load of this database put it there, and no other database claims
it. The load says what it replaced. A directory whose manifest names another
database, or none, still stops the load, and the refusal says when it is
another database's: two databases may share a data directory, and replacing
another database's assembly would break it.

**Cost:** one column, so databases made before this are refused with the
message REQUIRED_COLUMNS gives.

## D-096 — a feature written on several rows keeps each row's attributes and strand

**Found by the benchmark's round trip, 21 October 2026.** Exporting each
RefSeq genome and comparing it with its input, feature by feature, changed
features no check had reported: 3,249 in human, 7,729 in zebrafish, 976 in
maize. Two causes. The rows of one feature were stored with the first row's
attributes only (D-030), but some rows say something of their own: each row of
a `cDNA_match` alignment has its own `Target` and `Gap`, and a partial CDS
has `start_range` on its first row and `end_range` on its last. And the rows
of a trans-spliced gene that lies on both strands were all written back as `?`
(D-068), losing which strand each piece is on.

**Decided:** a row whose attributes differ from the first row's keeps its own,
stored against its own segment; a row that repeats them, as RefSeq's CDS rows
usually do, stores nothing more, so D-030's reason still holds. Each row keeps
the strand it was written with. The feature as a whole still reads `?` when its
rows lie on both strands, in `query` and `explain`; export, and the indexed
annotation in the data directory, write each row as it was. Readers of a
feature's attributes (functional matching, `explain`) read the first segment's,
as before.

**Cost:** none to the schema; databases loaded before this export the old way.

## D-097 — `precursor_RNA`, `pseudogenic_tRNA` and `pseudogenic_rRNA` are transcripts

**Found by the benchmark's first run, 21 October 2026.** Three RefSeq files of
the 30 were refused because a transcript type was missing from the synonym
table, so it loaded inert and its Tier A children were refused. *Xenopus
tropicalis* (GCF_000004195.4) writes a microRNA precursor with INSDC's feature
key, `precursor_RNA`, where newer NCBI files write `primary_transcript`; 398
exons and miRNAs below such precursors were refused. Rice and grape
mitochondria (GCF_034140825.1, GCF_030704535.1) write a `pseudogenic_tRNA` with
an exon below it.

**Decided:** `precursor_RNA` maps to `mRNA`, as `primary_transcript` does
(D-062). `pseudogenic_tRNA` and `pseudogenic_rRNA` map to `mRNA`, as
`pseudogenic_transcript` does: the Sequence Ontology makes each a kind of
pseudogenic transcript. `source_type` keeps the file's term.

**Cost:** the vocabulary hash changes, so databases made before this are refused
with the message `assert_vocabulary_current` gives.

## D-098 — the rows of a trans-spliced feature may overlap across strands

**Found by the benchmark's first run, 21 October 2026.** The *Drosophila
melanogaster* RefSeq file (GCF_000001215.4) was refused at line 370,106:
`gene-Dmel_CG32491`, the trans-spliced gene mod(mdg4), is written as two rows,
21,365,609-21,370,249 and 21,370,101-21,377,401, which overlap by 149 bases.
D-068 already merged the rows of a trans-spliced gene that lie on opposite
strands; the overlap check then refused them, because segments of one feature
may not overlap (D-033).

**Decided:** when a row carries the trans-splicing qualifier, a segment on the
opposite strand may overlap it. The pieces come from different pre-mRNAs, so
neither duplicates the other. An overlap on one strand, or across strands
without the qualifier, is still refused.

## D-099 — a CDS may read a base twice over a low-quality sequence region

**Found by the benchmark's first run, 21 October 2026.** The *S. pombe*
RefSeq file (GCF_000002945.2) was refused: 18 CDS are written with rows that
overlap by one or two bases, for example `cds-NP_592796.4`, 101,872-103,544
and 101,836-101,872. Every one carries `exception=low-quality sequence
region`: the model reads a base twice to step over an error in the genome
sequence, the way a programmed frameshift does (D-089).

**Decided:** `Profile.overlap_exception_values` replaces the single value. The
RefSeq and cgload profiles accept `ribosomal slippage` and `low-quality
sequence region`; the GenBank profile accepts the second as
`/exception="low-quality sequence region"`. Matched on the value, as before;
any other overlap is still refused.

## D-100 — a GenBank exon feature written on several rows shares its ID

**Found by the benchmark's first run, 21 October 2026.** The *Ustilago maydis*
RefSeq file (GCF_000328475.2) was refused at five tRNA genes with an intron.
Beside the tRNA's own exons, NCBI writes the gene's GenBank `exon` feature,
which has a join() location, as two rows sharing one ID, `ID=id-UMAG_16066;
Parent=gene-UMAG_16066;gbkey=exon`, exactly as it writes a misc_feature
(D-049). The RefSeq profile refuses any exon ID that repeats.

**Decided:** `Profile.shared_id_marker`. In the RefSeq and cgload profiles,
exon rows whose `gbkey` is `exon` may share an ID; they are then one feature
held to D-014's conditions (same sequence, strand and parent, no overlap). An
exon ID repeated on any other rows is still a duplicate and still refused.

## D-101 — no index on attribute values

**Found by the benchmark's first run, 21 October 2026.** The *Arabidopsis*
RefSeq file loaded on SQLite and failed on PostgreSQL with a traceback: "index
row size 3736 exceeds btree version 4 maximum 2704 for index
ix_feature_attribute_key_value". One attribute value is longer than
PostgreSQL can hold in a b-tree index entry. SQLite has no such limit, and
MySQL cannot index a TEXT column without a prefix length at all.

**Decided:** the index on (key, value) goes. No query looks attributes up by
value; functional matching reads them through the (feature, key) index.

**Cost:** none to the data. A database made before this keeps the index and can
still fail on such a value; load into a new one.

## D-102 — the cross-parser tier compares feature families, not names

**Found by the benchmark's first run, 21 October 2026.** `verify --against`
failed on five fungal RefSeq genomes. NCBI's GFF3 writes Sequence Ontology
terms where its GenBank file writes INSDC keys, and for `ncRNA` and
`repeat_region` the term follows a qualifier: yeast's 115 GenBank `ncRNA`
records are 76 `snoRNA`, 22 `ncRNA`, 7 `antisense_RNA`, 6 `snRNA` and four
single others in the GFF3, and its 384 `repeat_region` records are 384
`long_terminal_repeat` rows. `transcript` is GenBank's `misc_RNA`,
`origin_of_replication` its `rep_origin`, `mobile_genetic_element` its
`mobile_element`, `regulatory_region` its `regulatory`. GenBank also has
`assembly_gap` records, which the GFF3 leaves out. The tier compared names and
reported each as a loss.

**Decided:** `CROSS_FORMAT_FAMILIES` names each family by its GenBank key and
counts it once on each side, and `assembly_gap` joins the structural types
that only one format has. The measured families are listed as measured; the
others follow the same INSDC-to-SO naming. A feature missing from one side
still fails the tier.
