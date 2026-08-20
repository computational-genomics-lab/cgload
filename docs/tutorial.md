# Tutorial

Fifteen minutes, start to finish, using data bundled with the package. No
downloads, no database server, no credentials.

By the end you will have loaded a genome, proved the load was complete, broken
it deliberately to watch the proof fail, and exported the data back out.

## Before you start

```console
pip install cgload
cgload --version
```

Python 3.11 or newer. Three dependencies. Nothing else.

---

## 1. Create a database

```console
cd examples/minimal
cgload init
```

```
Created 11 tables in sqlite:///cgload.db
Seeded vocabulary version 1 (33 terms, c8c88be997de)
```

`cgload.db` is a single file in the current directory. That file *is* the
database — there is no server to start and no password to set. You can copy it,
email it, or attach it to a paper.

The vocabulary line matters later: it is a fingerprint of the feature-type rules
this database was built with. If those rules change, cgload will refuse the
database rather than reinterpret it silently.

## 2. Look at the inputs

Three files:

```console
ls
```

```
annotation.gff3   organism.toml   sequence.fa
```

`sequence.fa` is 5,000 bases of DNA. `annotation.gff3` says where three genes
are. `organism.toml` tells cgload what they are:

```toml
[organism]
genus = "Demonstratus"
species = "minimalis"
strain = "DEMO-1"

[assembly]
name = "Demo_v1"
fasta = "sequence.fa"

[[files.annotation]]
path = "annotation.gff3"
```

## 3. Load

```console
cgload load --config organism.toml
```

```
Loaded Demonstratus minimalis / DEMO-1 / Demo_v1 (assembly_id=1)
  1 sequence regions, 5,000 bp
  annotation.gff3: funannotate (high) -- 12 features
    CDS 2, exon 4, gene 3, mRNA 2, tRNA 1
  data directory: .../data/demonstratus_minimalis/demo-1/demo_v1
```

Three things happened worth noticing.

**The format was identified, not assumed.** `funannotate (high)` means cgload
recognised the dialect with high confidence. It did not guess: an unrecognised
file stops the load and names the profiles it considered.

**The DNA was copied, not referenced.** Look in the data directory:

```console
ls data/demonstratus_minimalis/demo-1/demo_v1/
```

```
assembly.fa   assembly.fa.fai   manifest.json
```

The copy costs disk space and buys something: a file you later move, edit or
delete cannot silently invalidate every coordinate in your database. The `.fai`
index means any region can be fetched instantly without reading the whole file.

**A receipt was written.** `manifest.json` records which version of cgload
loaded this, when, from which files, what each was recognised as, and a
checksum of every input. Later you can prove nothing was swapped.

## 4. See what is in there

```console
cgload stats
```

```
Demonstratus minimalis / DEMO-1 / Demo_v1 1.0
  1 sequence regions, 5,000 bp
  by source type: exon 4, gene 3, CDS 2, mRNA 2, tRNA 1
  by canonical type: exon 4, gene 3, CDS 2, mRNA 2, tRNA 1
```

Two counts, and the distinction is deliberate.

**Source type** is the word your file used. **Canonical type** is cgload's
normalised term. Here they match. They would not for an AUGUSTUS file, which
writes `transcript` where everyone else writes `mRNA` — cgload knows they mean
the same thing, and keeps both so nothing is lost and nothing is invented.

## 5. Prove the load was complete

```console
cgload verify --assembly-id 1
```

```
source count      PASS  12 features in the file, 12 in the database, across 5 type(s)
cross-parser      ----  no second format supplied; pass --against with a GenBank or GFF3 file describing the same assembly
reference (NCBI)  ----  no published counts supplied; pass --reference-counts and --accession
```

The first line is the guarantee. cgload read your file a second time, with a
deliberately simple counter that shares no code with the parser, and got the
same number.

That independence is the point. A checker built from the parser reads the file
the same wrong way and reports everything fine.

The other two tiers say **why** they did not run. A check that disappears when
it cannot run looks exactly like one that passed.

## 6. Break it on purpose

A check that never fires proves nothing.

```console
cp annotation.gff3 /tmp/backup.gff3
grep -v "ID=gene-002" annotation.gff3 > /tmp/broken.gff3
cp /tmp/broken.gff3 annotation.gff3

cgload init --force
cgload load --config organism.toml
cgload verify --assembly-id 1
```

The load succeeds — the file is still valid, just missing a gene — and
verification reports fewer features than before. Nothing about the load looked
wrong; only the count revealed it.

Restore:

```console
cp /tmp/backup.gff3 annotation.gff3
```

## 7. All or nothing

```console
printf 'chr1\tX\tCDS\t500\t100\t.\t+\t.\tID=broken\n' >> annotation.gff3
cgload init --force
cgload load --config organism.toml
```

```
cgload load: annotation.gff3: line 13: end 100 precedes start 500
```

Now check what survived:

```console
cgload stats
ls data/ 2>/dev/null
```

Nothing. No half-loaded genome, no orphaned files, no rows pointing at things
that are not there. The database and the data directory are exactly as they were.

This is harder than it sounds, because two separate things change during a load:
the database and a folder of files. cgload stages everything in a temporary
directory and moves it into place at the last possible moment, inside the
database transaction.

Restore and reload:

```console
cp /tmp/backup.gff3 annotation.gff3
cgload init --force && cgload load --config organism.toml
```

## 8. Ask why

```console
cgload explain gene-002
```

`explain` shows how one feature was linked and which rule produced the link.
That matters most for GenBank files, where nothing states which coding sequence
belongs to which transcript and cgload has to work it out from coordinates —
`explain` reports whether a link came from the file or was computed, and by
which method.

## 9. Get it back out

```console
cgload export --assembly-id 1 --out roundtrip.gff3
head -5 roundtrip.gff3
```

```
##gff-version 3
##cgload-export 2.0.0
```

cgload writes **its own dialect**, and declares it in the header. It does not
write back out as the format the data came in as. Producing an AUGUSTUS file
would mean inventing exon rows that AUGUSTUS does not emit; producing a RefSeq
file would mean fabricating cross-reference identifiers that were never there.

What it does preserve: the producing program per feature, and the feature type
exactly as the source wrote it. Load a mixed set from three pipelines and each
feature still says where it came from.

Round trip:

```console
cgload init --force
sed -i 's|annotation.gff3|roundtrip.gff3|' organism.toml
cgload load --config organism.toml && cgload stats
```

Same counts.

## Where to go next

- **Your own genome** — see the config example in the [README](../README.md).
- **Functional annotation** — add a `[[files.functional]]` block pointing at
  eggNOG-mapper output. cgload matches proteins to genes by exact identifier
  and refuses the load if too many fail to match, since a functional file
  produced from a different annotation is the commonest way to end up with a
  database that looks annotated and is not.
- **Large genomes** — [docs/scalability.md](scalability.md) has a full protocol
  for a 2.18 Gb, 696,928-feature assembly.
- **Browsing** — point JBrowse 2 at `assembly.fa` and the bgzip'd annotation in
  the data directory. cgload writes browser-native indexed formats so no
  conversion is needed.

## When something goes wrong

Every error names the file, the line where possible, and what to do:

```
cgload load: mystery.gff3: no confident match: column 2 is ['MysteryTool'], and
no registered profile claims it (registered: refseq, funannotate, augustus,
genbank, cgload). Pass --dialect to force one.
```

```
cgload load: eggnog.annotations: 847 of 1,204 query proteins (70.3%) match no
feature in assembly 1. First unmatched: FUN_000001-T1, FUN_000002-T1, ...
Either the file was produced from a different annotation, or the FASTA headers
were renamed.
```

If a message is unclear, that is a bug worth reporting.
