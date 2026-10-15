# Tutorial

Fifteen minutes, start to finish, using data bundled with the package. No
downloads, no database server, no credentials.

By the end you will have loaded a genome, proved the load was complete, broken
the database on purpose to watch the proof fail, asked it questions, exported
the data back out, and removed it again.

Every output below was produced by the commands shown, on the bundled example.

## Before you start

```console
pip install git+https://github.com/computational-genomics-lab/cgload
cgload --version
```

cgload is not on PyPI yet. It needs Python 3.11 or newer and three packages,
which pip installs. If your system Python is older (Ubuntu 22.04 ships 3.10),
make a newer one first, for example with
`conda create -n cgload -c conda-forge python=3.12`.

---

## 1. Create a database

```console
cd examples/minimal
cgload init
```

```
Created 11 tables in sqlite:///cgload.db
Seeded vocabulary version 1 (43 terms, d5fa82cadbce)
```

`cgload.db` is a single file in the current directory. That file *is* the
database. There is no server to start and no password to set. You can copy it,
email it, or attach it to a paper.

The vocabulary line matters later: it is a fingerprint of the feature-type rules
this database was built with. If those rules change, cgload will refuse the
database rather than reinterpret it silently.

## 2. Look at the inputs

```console
ls
```

```
annotation.gff3
assembly.fa
cgload.db
organism.toml
```

`assembly.fa` is 5,000 bases of DNA. `annotation.gff3` says where three genes
are. `organism.toml` tells cgload what they are:

```toml
[organism]
genus = "Demonstratus"
species = "minimalis"
strain = "DEMO-1"
# An opaque identifier. cgload stores it and resolves nothing (D-016).
ncbi_taxon_id = 999999

[assembly]
name = "Demo_v1"
version = "1.0"
fasta = "assembly.fa"

[[files.annotation]]
path = "annotation.gff3"
# dialect = "funannotate"   # optional; detection reports its verdict either way
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
  data directory: ./data/demonstratus_minimalis/demo-1/demo_v1
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
annotation.gff3.gz
annotation.gff3.gz.tbi
assembly.fa
assembly.fa.fai
manifest.json
```

The copy costs disk space and buys something: a file you later move, edit or
delete cannot silently invalidate every coordinate in your database. The `.fai`
index lets any region be fetched without reading the whole file. The
annotation is there too, compressed and indexed the way genome browsers expect,
so JBrowse 2, IGV and `tabix` can read it directly.

**A receipt was written.** `manifest.json` records which version of cgload
loaded this, when, from which files, what each was recognised as, and a
checksum of every input.

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
writes `transcript` where others write `mRNA`. cgload knows they mean the same
thing, and keeps both, so nothing is lost and nothing is invented.

## 5. Prove the load was complete

```console
cgload verify --assembly-id 1
```

```
  source count      PASS  12 features in the file, 12 in the database, across 5 type(s)
  parent edges      PASS  9 edge(s), all accounted for (explicit_parent 9)
  stored assembly   PASS  ./data/demonstratus_minimalis/demo-1/demo_v1/assembly.fa and its index are present; SHA-256 matches the load
  cross-parser      ----  no second format supplied; pass --against with a GenBank or GFF3 file describing the same assembly
  reference (NCBI)  ----  no published counts supplied; pass --reference-counts and --accession
```

The first line is the guarantee. cgload read your file a second time, with a
deliberately simple counter that shares no code with the parser, and got the
same number. A checker built from the parser would read the file the same wrong
way and report everything fine.

The second line says every parent the file names is stored as a link. The
third says the DNA in the data directory is the DNA that was loaded.

The last two tiers say why they did not run. A check that disappears when it
cannot run looks exactly like one that passed.

## 6. Break it on purpose

A check that never fires proves nothing. Delete one exon from the database
behind cgload's back:

```console
python3 - <<'EOF'
import sqlite3

db = sqlite3.connect("cgload.db")
(row,) = db.execute(
    "SELECT feature_id FROM feature WHERE source_id = 'DEMO_000003-T1.exon1'"
).fetchone()
db.execute("DELETE FROM feature_attribute WHERE feature_id = ?", (row,))
db.execute("DELETE FROM feature WHERE feature_id = ?", (row,))
db.commit()
EOF
cgload verify --assembly-id 1
```

```
  source count      FAIL  1 feature type(s) differ between the file and the database; features were lost or invented during loading
      exon                 expected        4   found        3   (-1)
  parent edges      PASS  8 edge(s), all accounted for (explicit_parent 8)
  stored assembly   PASS  ./data/demonstratus_minimalis/demo-1/demo_v1/assembly.fa and its index are present; SHA-256 matches the load
  cross-parser      ----  no second format supplied; pass --against with a GenBank or GFF3 file describing the same assembly
  reference (NCBI)  ----  no published counts supplied; pass --reference-counts and --accession
```

`verify` exits with status 2. The parent-edge tier still passes, because the
exon took its own link with it. Only the count against the file shows the loss.

Replace the damaged assembly with a fresh load:

```console
cgload load --config organism.toml --force
cgload verify --assembly-id 1
```

The source count passes again.

## 7. All or nothing

Add a gene that runs past the end of the 5,000 bp sequence, and try to load the
file again:

```console
cp annotation.gff3 backup.gff3
printf 'demo_scaffold_1\tfunannotate\tgene\t4900\t5200\t.\t+\t.\tID=DEMO_000004;\n' >> annotation.gff3
cgload load --config organism.toml --force
```

```
cgload load: annotation.gff3: 1 row(s) extend past the end of their sequence, so their coordinates are wrong:
  gene DEMO_000004 at line 15: 4,900-5,200 on demo_scaffold_1 (5,000 bp)
Only a feature that crosses the origin of a circular sequence may end past its length. If the sequence is circular, its region row should carry Is_circular=true.
```

Now check what survived:

```console
cgload stats
```

```
Demonstratus minimalis / DEMO-1 / Demo_v1 1.0
  1 sequence regions, 5,000 bp
  by source type: exon 4, gene 3, CDS 2, mRNA 2, tRNA 1
  by canonical type: exon 4, gene 3, CDS 2, mRNA 2, tRNA 1
```

The previous load is untouched, and so is its data directory. A load changes
two things, the database and a folder of files, and cgload changes both or
neither: it stages the files in a temporary directory and moves them into
place inside the database transaction, at the last moment.

`verify` notices that the file on disk is no longer the one that was loaded,
and says so instead of blaming the loader:

```console
cgload verify --assembly-id 1
```

```
  source count      ----  ./annotation.gff3 has changed since it was loaded (its SHA-256 differs), so counting it would compare the database with a different file. Restore the original, or reload with --force
```

The other tiers are as before. Restore the file:

```console
cp backup.gff3 annotation.gff3
```

## 8. Ask why

```console
cgload explain DEMO_000001-T1.exon1
```

```
exon DEMO_000001-T1.exon1
  organism    Demonstratus minimalis / DEMO-1 / Demo_v1
  location    demo_scaffold_1:1,000-1,300 (+)
  type        exon in the file, stored canonically as exon
  from file   ./annotation.gff3
              recognised as funannotate (high confidence)
              sha256 74baa173424a84cd...
  parent      mRNA DEMO_000001-T1
  linked by   explicit_parent — the file stated it (GFF3 Parent attribute)
              this link was stated by the file
  recorded by cgload 0.1.0 under vocabulary version 1
```

`explain` shows how one feature was linked and which rule produced the link.
That matters most for GenBank files, where nothing states which coding sequence
belongs to which transcript and cgload works it out from coordinates. `explain`
reports whether a link came from the file or was computed, and by which method.

## 9. Ask questions

`query` finds features by ID, by region, or by a functional term from
eggNOG-mapper or InterProScan. It changes nothing. A gene with everything below
it:

```console
cgload query --assembly-id 1 --id DEMO_000001 --descendants
```

```
# cgload query: ID DEMO_000001 and its descendants
id	type	feature_type	seqid	start	end	strand	segments	parent	top_level
DEMO_000001	gene	gene	demo_scaffold_1	1000	2400	+	1		DEMO_000001
DEMO_000001-T1	mRNA	mRNA	demo_scaffold_1	1000	2400	+	1	DEMO_000001	DEMO_000001
DEMO_000001-T1.exon1	exon	exon	demo_scaffold_1	1000	1300	+	1	DEMO_000001-T1	DEMO_000001
DEMO_000001-T1.cds	CDS	CDS	demo_scaffold_1	1000	2400	+	2	DEMO_000001-T1	DEMO_000001
DEMO_000001-T1.exon2	exon	exon	demo_scaffold_1	1500	2400	+	1	DEMO_000001-T1	DEMO_000001
```

The CDS is written on two rows in the file; `segments` says so, and its
`start` and `end` span both. A count of the features found goes to standard
error, so standard output holds only the table. Everything in a region:

```console
cgload query --assembly-id 1 --region demo_scaffold_1:2900-3700
```

```
# cgload query: features overlapping demo_scaffold_1:2,900-3,700
id	type	feature_type	seqid	start	end	strand	segments	parent	top_level
DEMO_000002	gene	gene	demo_scaffold_1	3000	3600	-	1		DEMO_000002
DEMO_000002-T1	mRNA	mRNA	demo_scaffold_1	3000	3600	-	1	DEMO_000002	DEMO_000002
DEMO_000002-T1.exon1	exon	exon	demo_scaffold_1	3000	3600	-	1	DEMO_000002-T1	DEMO_000002
DEMO_000002-T1.cds	CDS	CDS	demo_scaffold_1	3000	3600	-	1	DEMO_000002-T1	DEMO_000002
```

A region matches the rows of a feature, not its outer span, so a position inside
an intron returns the gene and the mRNA but not the CDS. `--format gff3` writes
GFF3 instead, with every ancestor of each feature found, so each `Parent` in the
output resolves within it. `--term` takes an accession such as `GO:0005524`,
`K00799` or `PF00069`; the example has no functional annotation to search.

## 10. Get it back out

```console
cgload export --assembly-id 1 --output roundtrip.gff3
head -8 roundtrip.gff3
```

```
##gff-version 3
##cgload-export 0.1.0
##cgload-vocabulary d5fa82cadbce07c9f15895aad21658c1e83705c0b902f49e867f9c39818edea5
##sequence-region demo_scaffold_1 1 5000
# The source dialect does not emit these types, so their absence is a property of the input and not a loss during loading:
#   intron, start_codon, stop_codon
# source programs present: funannotate
demo_scaffold_1	funannotate	exon	1000	1300	.	+	.	ID=DEMO_000001-T1.exon1;Parent=DEMO_000001-T1;cgload_linkage_method=explicit_parent
```

cgload writes **its own dialect**, and declares it in the header. It does not
write the data back out in the format it came in. Producing an AUGUSTUS file
would mean inventing exon rows that AUGUSTUS does not emit; producing a RefSeq
file would mean fabricating cross-reference identifiers that were never there.

What it does keep: the producing program of each feature, in column 2, and the
feature type exactly as the source wrote it. Load a mixed set from three
pipelines and each feature still says where it came from.

Round trip:

```console
sed 's|annotation.gff3|roundtrip.gff3|' organism.toml > roundtrip.toml
cgload load --config roundtrip.toml --force
cgload stats
```

```
Loaded Demonstratus minimalis / DEMO-1 / Demo_v1 (assembly_id=1)
  1 sequence regions, 5,000 bp
  roundtrip.gff3: cgload (high) -- 12 features
    CDS 2, exon 4, gene 3, mRNA 2, tRNA 1
  data directory: ./data/demonstratus_minimalis/demo-1/demo_v1
Demonstratus minimalis / DEMO-1 / Demo_v1 1.0
  1 sequence regions, 5,000 bp
  by source type: exon 4, gene 3, CDS 2, mRNA 2, tRNA 1
  by canonical type: exon 4, gene 3, CDS 2, mRNA 2, tRNA 1
```

Same counts, and the file was recognised as cgload's own dialect.

## 11. Take it out again

`remove` deletes an assembly's rows and its data directory together, or
neither. First see what it would delete:

```console
cgload remove --assembly-id 1 --dry-run
```

```
assembly 1: Demonstratus minimalis / DEMO-1 / Demo_v1 1.0
  12 features in 13 rows, 3 attributes
  0 functional assignments from 0 file(s)
  1 sequence regions, 2 source file record(s)
  data directory ./data/demonstratus_minimalis/demo-1/demo_v1
Dry run: nothing was removed.
```

Then remove it. Without `--yes` it asks first:

```console
cgload remove --assembly-id 1 --yes
cgload stats
```

```
...
Removed assembly 1.
No assemblies loaded.
```

It deletes only a directory whose `manifest.json` names this assembly and
records the same FASTA checksum as the database. `--keep-files` removes the
rows alone.

## Where to go next

- **Your own genome**: see the config example in the [README](../README.md).
- **Functional annotation**: add a `[[files.functional]]` block pointing at
  eggNOG-mapper or InterProScan output. cgload matches proteins to features by
  exact identifier and refuses the load if any fail to match, since a
  functional file produced from a different annotation is the commonest way to
  end up with a database that looks annotated and is not.
- **Large genomes**: [docs/scalability.md](scalability.md) has a full protocol
  for a 2.18 Gb, 696,928-feature assembly.
- **Browsing**: point JBrowse 2 at `assembly.fa` and `annotation.gff3.gz` in
  the data directory. Both are indexed, so no conversion is needed.

## When something goes wrong

Every error names the file, the line where possible, and what to do. A file
whose producer cgload does not recognise:

```
cgload load: mystery.gff3: no confident match: column 2 is ['MysteryTool'], and no registered profile claims it (registered: refseq, funannotate, augustus, genbank, cgload). To force one, set dialect = "<name>" under this file's [[files.annotation]] entry in the config.
```

An eggNOG-mapper file made from a different genome than the one being loaded:

```
cgload load: eggnog.annotations: 20 of 20 query proteins (100.0%) match no feature in assembly 1. First unmatched: XP_011326024.1, XP_011326026.1, XP_011326027.1, XP_011326029.1, XP_011326036.1. Either the file was produced from a different annotation, or the FASTA headers were renamed. To load only the matching subset, on purpose, set minimum_match_rate below 1.0 under this file's [[files.functional]] entry in the config.
```

Neither load changed anything. If a message is unclear, that is a bug worth
reporting.
