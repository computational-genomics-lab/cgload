# cgload

**Load genome assemblies and their annotations into a queryable relational
database, reconcile the incompatible formats different pipelines produce,
check that every feature in the file reached the database, and put the result
on the web for people to search, browse and download.**

cgload counts what it stored and checks that count against the file it read,
using code that shares nothing with its parser. `cgload portal` then makes the
database a genome portal, in one command: search by gene name, product or GO
term, a gene's model drawn to scale, a genome browser, sequences and proteins,
BLAST, a JSON API.

## Quick start

You need Python 3.11 or newer. Install cgload from GitHub, best into a fresh
virtual or conda environment:

```console
pip install git+https://github.com/computational-genomics-lab/cgload
```

That installs the loader and the portal. Then one command takes a genome to a
running portal with BLAST. Try it on the small genome in the repository:

```console
git clone https://github.com/computational-genomics-lab/cgload
cd cgload/examples/minimal
cgload portal organism.toml
```

```
cgload portal: 1 genome to load, then the search index, then BLAST, then the portal.
  database: created, 11 tables, at sqlite:///cgload.db
  load: Demonstratus minimalis / DEMO-1 / Demo_v1 ...
    assembly 1: 1 sequence, 12 features; every check passed (parent edges [annotation.gff3], source count [annotation.gff3])
  search index: assembly 1, 17 terms; every count matches
  BLAST+: not on this computer; cgload installs it for itself, once (about 900 MB)
  downloading micromamba 2.10.0-0, which installs BLAST+ from Bioconda
  installing BLAST+ 2.17.0 into .../share/cgload/blast (once; a few minutes)
  BLAST databases: building for 1 genome
    assembly 1 (Demonstratus minimalis · DEMO-1 · Demo_v1): 1 sequence, 2 proteins: PASS
cgload portal: the portal is running, read-only, with 1 genome.
  Open it on this computer:   http://127.0.0.1:8000/
  or from another computer:   http://192.0.2.10:8000/
  Anyone who can reach this computer can open it. To keep it to this computer: --host 127.0.0.1;
  to serve it beyond your network, put a web server with HTTPS in front (docs/portal.md).
  Stop it with Ctrl-C.
  BLAST: 1 of 1 genome, 2 searches at a time.
```

Open one of the addresses it prints in a web browser. Each step is done once:
run `cgload portal` again and it starts the portal at once; give it another
organism file and it loads that genome first. BLAST+ goes into cgload's own
folder (inside the environment), not onto the system; `--no-blast` leaves
BLAST out. On Windows, install BLAST+ yourself, or run cgload in WSL.

**Your own genomes.** An organism file names the assembly's FASTA and its
GFF3 annotation; copy `examples/minimal/organism.toml` and change the names
and paths. Then `cgload portal mygenome.toml`, or several at once:
`cgload portal genomes/*.toml`. The database is a single SQLite file
(`cgload.db`), with no server and no credentials; PostgreSQL works through
the same code with `--db-url` and `pip install "cgload[postgres] @ git+https://github.com/computational-genomics-lab/cgload"`.

**Step by step.** `cgload portal` runs commands you can also run yourself:

```console
cgload init                          # an empty database, ./cgload.db
cgload load --config organism.toml   # load, count, keep it only if the counts agree
cgload verify --assembly-id 1        # check it again, at any time
cgload index                         # the portal's search index
cgload blastdb                       # the BLAST databases (BLAST+ needed)
cgload serve                         # the portal, read-only
```

### If something goes wrong

| What you see | What to do |
|---|---|
| `cgload: command not found` | The environment you installed cgload into is not active. Activate it, or run `python -m cgload`. |
| `no database at .../cgload.db` | Run the command in the folder that holds `cgload.db`, or give `--db-url`; `cgload portal organism.toml` makes one. |
| `Port 8000 was in use, so it is on port 8001` | Nothing to do: open the address it prints. To choose the port, give `--port 9000`. |
| `port 9000 is already in use` | The port you asked for is taken. Choose another, or stop the program that holds it. |
| Another computer cannot open the portal | The server's firewall closes the port. On Ubuntu: `sudo ufw allow from 192.168.1.0/24 to any port 8000 proto tcp` (your network's range). And do not start it with `--host 127.0.0.1`. |
| A search for a gene's name finds nothing | Run `cgload portal` (or `cgload index`), and again after loading more genomes. |
| Your data are not public | Anyone who can reach the computer can open the portal. Keep it to this computer with `--host 127.0.0.1`, or put it behind a web server with a password ([docs/portal.md](docs/portal.md)). |
| `BLAST+: cgload cannot install BLAST+ on this computer` | On Windows, or a processor Bioconda does not serve: install BLAST+ by hand (the message says how), then run `cgload portal` again. |
| `BLAST+: could not download ...` | No internet, or a proxy in the way. Run it again where there is, or install BLAST+ by hand; `--no-blast` starts the portal without it. |
| `BLAST is off: no genome has BLAST databases yet` | Run `cgload portal` (or `cgload blastdb`), and again after loading or reloading a genome. |
| `the portal needs Flask and waitress` | The installation is incomplete. Reinstall: `pip install --force-reinstall git+https://github.com/computational-genomics-lab/cgload`. |

Full walkthrough: [docs/tutorial.md](docs/tutorial.md). The portal, its JSON
API, deployment behind nginx and its security model:
[docs/portal.md](docs/portal.md).

## Why

A parser that silently drops a feature leaves a database that looks populated,
answers queries, and is wrong. Nothing downstream can tell. cgload is built
around one property — **a load either completes in full or changes nothing** —
and completeness is checked rather than assumed.

What it does:

1. **Dialect reconciliation.** RefSeq, funannotate, AUGUSTUS/BRAKER and GenBank
   disagree about feature types, parent–child linkage and identifier
   conventions, and GenBank carries no parent pointers at all. Every input is
   normalised to one representation before it reaches the database, with a
   documented resolution rule and a test per dialect.
2. **Completeness verification.** After every load, cgload recounts the
   source file with code that shares nothing with the parser, compares the
   counts with what it stored, and checks that every parent the file names is
   stored as a link. A mismatch rolls the load back, whether or not an
   exception was raised. `verify` runs these checks again later, and can also
   compare the file with the provider's own published totals.
3. **Functional annotation as queryable relations.** eggNOG-mapper and
   InterProScan output becomes one row per GO term, KO, Pfam domain or EC
   number, not a text blob to grep.
4. **A portal for the people who use the genome.** `cgload serve` puts the
   loaded genomes on the web: search, gene pages, igv.js, BLAST against the
   genomes and their proteins, FASTA and GFF3 downloads, and the load's own
   receipt for every genome. It opens the
   database read-only and was built against the attacks genome portals
   suffer, injection first ([docs/portal.md](docs/portal.md)).

## The checks

| tier | compares | what it can prove |
|---|---|---|
| source count | an independent pass over your file against the database | every feature, and every row of a feature written on several rows, reached the database, and none was invented |
| parent edges | the stored links against the `Parent=` statements the same independent pass counts in your file | every parent the file names is stored as a link, not merely the features |
| stored assembly | the stored FASTA and its index against the checksum taken at load | the sequence every coordinate refers to is still there, unchanged |
| cross-parser | the same assembly as GFF3 and as GenBank | the two files describe the same feature counts |
| reference | the provider's published totals against the file | the file matches what the provider says it published |

`load` runs the source count, and for GFF3 the parent-edge tier, itself, and
rolls back if either fails. `verify` runs every tier again later.

**What this does and does not establish.** The first three tiers compare the
database against the input, and a failure there means cgload is wrong, or that
something changed the database or the data directory afterwards. The last two
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
    assembly.fa             your file, copied
    assembly.fa.fai         faidx index
    annotation.gff3.gz      the normalised annotation, bgzip'd, sorted for tabix
    annotation.gff3.gz.tbi  its tabix index
    manifest.json           checksums, provenance, any repairs made
```

This is what samtools, IGV, JBrowse and Ensembl all do. It gives O(1) access to
any region, keeps the database small enough to attach to a paper, and makes the
output directly usable by JBrowse 2 without any code from us. The BGZF file and
its index are written by cgload itself, with no external tools (D-082).

## Scale

Measured on 30 RefSeq genomes, from *Candida albicans* (25,889 GFF3 rows) to
human (4.9 million), on SQLite and on PostgreSQL 16, on one server with 80
cores and 125 GB of memory. Three of them:

| | rows | features | load, SQLite | load, PostgreSQL | peak memory | `verify` | database |
|---|--:|--:|--:|--:|--:|--:|--:|
| *Fusarium graminearum* | 102,635 | 78,715 | 14 s | 43 s | 0.8 / 0.9 GB | 1.7 / 2.9 s | 0.07 / 0.15 GB |
| *Zea mays* | 990,782 | 696,928 | 132 s | 366 s | 5.3 / 6.5 GB | 7.3 / 10.3 s | 0.52 / 1.10 GB |
| *Homo sapiens* | 4,932,571 | 3,203,237 | 641 s | 1,978 s | 29.5 / 36.4 GB | 29.1 / 47.8 s | 3.00 / 5.86 GB |

Where two values are given, the first is SQLite and the second PostgreSQL.

Time and memory grow with the number of GFF3 rows, not with the size of the
genome: the sequence is never held in memory. A load holds the whole
annotation, about 5.7 KB per row on SQLite and 6.9 KB on PostgreSQL (medians
over the 30 genomes), so a human-sized annotation needs 30 to 37 GB. SQLite
loads about 8,000 rows a second and PostgreSQL about 2,700.

The portal was measured on the same 30 genomes in one 20.6 GiB SQLite
database by [benchmark/run_portal.py](benchmark/run_portal.py) (PROTOCOL.md
section 13): a search by a gene's name took a median of 23 ms and a gene's
page 20 ms, all 6,000 genes searched for by name were found, every one of
60,000 translated proteins was either free of internal stop codons or
explained by its file, and 3,063 injection attempts and random requests
caused no server error and left the database unchanged.

All 30 genomes, the protocol and the comparison with other loaders are in
[benchmark/PROTOCOL.md](benchmark/PROTOCOL.md) and
[paper/data/benchmark-run-2.2/](paper/data/benchmark-run-2.2/).
[docs/scalability.md](docs/scalability.md) records two earlier runs, step by
step, with eggNOG-mapper output loaded as well.

## Commands

```
cgload init      create the schema
cgload load      load an organism; transactional
cgload verify    completeness checks
cgload explain   why one feature was linked the way it was
cgload stats     counts per organism, assembly, feature type
cgload export    write back out as GFF3
cgload query     features by identifier, region or functional term; read-only
cgload remove    remove an assembly, rows and files together; transactional
cgload index     build and check the portal's search index
cgload serve     the web portal; read-only
```

`query` answers with whole features, as TSV or as GFF3:

```console
cgload query --assembly-id 1 --id gene-FGSG_11614 --descendants        # a gene model
cgload query --assembly-id 1 --region NC_026477.1:1-20000 --type gene   # genes in a region
cgload query --assembly-id 1 --term PF03198                             # proteins with a Pfam hit
cgload query --assembly-id 1 --term K00799 --format gff3                # GFF3, with ancestors
```

The TSV has one row per feature, with its parent and the top of its parent
chain (`top_level`, usually the gene). With `--term` there is one row per
functional assignment. The GFF3 adds every ancestor, so each `Parent` resolves
inside the file. Accessions match exactly, apart from a few listed spellings
(`K00799` finds eggNOG's `ko:K00799`; D-085).

`remove` asks before it deletes anything; `--dry-run` shows what would go. It
deletes the data directory only when that directory's `manifest.json` names
the assembly and its FASTA checksum (D-084). `--keep-files` removes only the
rows.

## Deliberately out of scope

The portal reads; it does not edit, upload or keep accounts. Writing to the
database from the web is where genome portals acquire their attack surface, so
loading stays on the command line, behind the checks. BLAST is the one program
the portal runs, and it runs as [D-116](docs/decisions.md) sets out: an
argument list of fixed choices, never a shell, the query as its input and
never an argument, a few searches at a time. A lab that wants all of BLAST's
options, or databases of its own, is served better by SequenceServer. The
assembly and annotation are also written in browser-native indexed formats, so
JBrowse 2 can be pointed at cgload's output directly as well.

See [docs/roadmap.md](docs/roadmap.md).

## Developing cgload

```console
git clone https://github.com/computational-genomics-lab/cgload
cd cgload
pip install -e ".[dev]"
pytest
```

## Design notes

Every non-obvious decision is recorded in [docs/decisions.md](docs/decisions.md)
with the evidence that produced it, usually a real file that broke an
assumption. It is long, and it is the honest record of why the code looks the
way it does.

## Licence

MIT. See [LICENSE](LICENSE).
