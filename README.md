# cgload

**Load genome assemblies and their annotations into a queryable relational
database, reconcile the incompatible formats different pipelines produce, and
check that every feature in the file reached the database.**

A cgload fixes your position by measuring against an independent fixed
reference. That is what this tool does with a genome annotation: it counts what
it stored, and checks that count against sources it does not control.

## Why

A parser that silently drops a feature leaves a database that looks populated,
answers queries, and is wrong. Nothing downstream can tell. cgload is built
around one property — **a load either completes in full or changes nothing** —
and completeness is checked rather than assumed.

Three things, together, that are not available elsewhere:

1. **Dialect reconciliation.** RefSeq, funannotate, AUGUSTUS/BRAKER and GenBank
   disagree about feature types, parent–child linkage and identifier
   conventions, and GenBank carries no parent pointers at all. Every input is
   normalised to one representation before it reaches the database, with a
   documented resolution rule and a test per dialect.
2. **Completeness verification.** After every load, cgload counts what it
   stored and compares it against the source and against the provider's own
   published totals. A mismatch fails the command, whether or not an exception
   was raised.
3. **Functional annotation as queryable relations.** eggNOG output becomes one
   row per GO term, KO, PFAM domain or EC number — not a text blob to grep.

## The three checks

| tier | compares | what it can prove |
|---|---|---|
| source count | an independent pass over your file against the database | every feature in the file reached the database, and none was invented |
| parent edges | the stored hierarchy against the features that require one | the parent–child structure exists, not merely the features |
| cross-parser | the same assembly as GFF3 and as GenBank | the two files describe the same feature counts |
| reference | the provider's published totals against the file | the file matches what the provider says it published |

**What this does and does not establish.** Tiers 1 and 2 compare the database
against the input, and a failure there means cgload is wrong. Tiers 3 and 4
compare two descriptions of an assembly to each other — they are consistency
checks on the *provider's* outputs, which is useful and is how four defects in a
published reference annotation were found, but a failure there is not
necessarily a cgload defect.

The guarantee is about **completeness, not fidelity**: every feature is present
and counted. A corruption that preserves counts — a wrong strand, a swapped
coordinate — is not something a counting check can see, and cgload does not
claim otherwise.

The naive counting pass shares no code with the parser. A verifier built from
the parser can only catch storage mistakes; the mistakes that matter happen
during parsing and are invisible to a check that parses the same way.

When a tier cannot run it says so, with the reason. A tier that vanishes when it
cannot run looks identical to one that passed.

## What it reads

| Source | Notes |
|---|---|
| NCBI RefSeq GFF3 | the reference case |
| funannotate GFF3 | `-T1` transcript suffixes, `DBxref` spelling |
| AUGUSTUS / BRAKER | no `exon` rows, `intron` instead, unnamed child rows |
| GenBank flat file | no parent pointers; hierarchy computed by containment |
| eggNOG-mapper | functional annotation |

Each is a declared profile — what the dialect provides and how it identifies
itself — not a pile of special cases. Adding funannotate and AUGUSTUS required
no change to the reading engine.

**An unrecognised file is refused, not guessed at.** The refusal names the
profiles that were considered.

## Install

Python 3.11 or newer. Three runtime dependencies: `sqlalchemy`, `pyfaidx`,
`click`.

```console
pip install cgload
```

From source:

```console
git clone https://github.com/computational-genomics-lab/cgload
cd cgload
pip install -e ".[dev]"
pytest
```

## Five minutes

```console
cd examples/minimal
cgload init
cgload load --config organism.toml
cgload stats
cgload verify --assembly-id 1
```

```
Created 11 tables in sqlite:///cgload.db
Seeded vocabulary version 1 (33 terms, c8c88be997de)
Loaded Demonstratus minimalis / DEMO-1 / Demo_v1 (assembly_id=1)
  1 sequence regions, 5,000 bp
  annotation.gff3: funannotate (high) -- 12 features
    CDS 2, exon 4, gene 3, mRNA 2, tRNA 1

  source count      PASS  12 features in the file, 12 in the database, across 5 type(s)
```

The database is a single SQLite file — no server, no credentials, nothing
running in the background. MySQL and PostgreSQL work through the same code path
via `--db-url`.

Full walkthrough: [docs/tutorial.md](docs/tutorial.md).

## Your own data

```toml
[organism]
genus = "Fusarium"
species = "graminearum"
strain = "PH-1"
ncbi_taxon_id = 229533

[assembly]
name = "ASM24013v3"
accession = "GCF_000240135.3"
fasta = "genomic.fna"

[[files.annotation]]
path = "genomic.gff"

[[files.functional]]
path = "eggnog.emapper.annotations"
```

```console
cgload load --config organism.toml
cgload verify --assembly-id 1 --against genomic.gbff
```

## Sequence lives on disk

The assembly is stored as indexed FASTA beside the database, not inside it. The
database holds coordinates and a reference.

```
data/<organism>/<strain>/<assembly>/
    assembly.fa          your file, copied
    assembly.fa.fai      faidx index
    annotation.gff3.gz   normalised, bgzip'd
    annotation.gff3.gz.tbi
    manifest.json        checksums, provenance, any repairs made
```

This is what samtools, IGV, JBrowse and Ensembl all do. It gives O(1) access to
any region, keeps the database small enough to attach to a paper, and makes the
output directly usable by JBrowse 2 without any code from us.

## Scale

*Zea mays* B73 (GCF_902167145.1), complete assembly including organelles, with
eggNOG functional annotation:

| | |
|---|---|
| sequences | 687, 2.18 Gb |
| features | 696,928 across 21 source types |
| load | 2 min 34 s, peak 5.4 GB RSS |
| verify | 0.3 s |
| database | 1.03 GB; sequence 2.1 GB on disk |
| source count | **PASS**, exact |

Memory scales with feature count, not genome size — sequence is never held in
memory. Budget roughly 4.5 GB per million features. Protocol:
[docs/scalability.md](docs/scalability.md).

## Commands

```
cgload init      create the schema
cgload load      load an organism; transactional
cgload verify    completeness checks
cgload explain   why one feature was linked the way it was
cgload stats     counts per organism, assembly, feature type
cgload export    write back out as GFF3

cgload query     NOT YET IMPLEMENTED - genes by identifier, region, annotation
cgload remove    NOT YET IMPLEMENTED - transactional removal
```

The two unimplemented commands exist in the surface, exit non-zero, and name
what is missing. They never return success.

## Deliberately out of scope

No web layer. Server-side query execution and shelled-out command-line tools are
where genome portals acquire their attack surface. Because the assembly and its
annotation are written in browser-native indexed formats, JBrowse 2 — a static
application with effectively no server-side attack surface — can be pointed
directly at cgload's output. BLAST is served better by SequenceServer.

See [docs/roadmap.md](docs/roadmap.md).

## Design notes

Every non-obvious decision is recorded in [docs/decisions.md](docs/decisions.md)
with the evidence that produced it, usually a real file that broke an
assumption. It is long, and it is the honest record of why the code looks the
way it does.

## Licence

MIT. See [LICENSE](LICENSE).
