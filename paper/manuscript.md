---
title: "cgload: genome annotations loaded with a completeness check and served through a read-only web portal"
---

<!--
Source of the manuscript, laid out as a research article for BMC
Bioinformatics. paper/build.py makes the Word file: it puts in the tables that
paper/figures.py made from paper/data/benchmark-run-2.2, numbers the citations
(written here by key, as [@key]) in order of first citation, puts in the
portal's numbers (written here as {{portal:KEY}}) from the portal's results,
and shows the notes to the authors. Do not type a result number into a table
by hand; change the data or the script.
-->

**Aditya Upadhyay**^1,2^, **Arijit Panda**^3^ and **Sucheta Tripathy**^1,2\*^

^1^ Structural Biology and Bioinformatics Division, CSIR-Indian Institute of Chemical Biology, Jadavpur, Kolkata 700032, India.
^2^ Academy of Scientific and Innovative Research (AcSIR), Ghaziabad 201002, India.
^3^ Department of Quantitative Health Sciences, Mayo Clinic, Rochester, MN 55905, USA.

\* Corresponding author.

Corresponding author e-mail address: tsucheta@gmail.com

# Abstract

**Background:** Genome annotations in GFF3 format are loaded into databases,
and most people who use them reach them through a web portal built on such a
database. A loader can drop or alter features as it reads, and the database
still answers queries, so the loss reaches the portal unseen. We know of no
study that measures how often this happens with real files.

**Results:** We developed cgload, a Python tool that loads an assembly and its
annotation into SQLite or PostgreSQL in one transaction, compares what it
stored with a count of the input made by separate code, and rolls back a load
that falls short. Two more commands serve the database as a read-only web
portal: search by name, product or functional term, gene pages with a genome
browser, and sequences. Under a protocol fixed
beforehand, we compared cgload with three other loaders on 30 RefSeq genomes
(28.2 million GFF3 rows, 18.0 million features), checking each result against
a reference reading that shares no code with any loader. cgload stored
18,017,606 features unchanged in every field compared, declared the 12 it
repaired and lost none. gffutils with `create_unique` stored 1,072,125
multi-row features as separate records without a message. txdbmaker stopped on six plant genomes. python-chado did not
finish nine genomes and, in the other 21, lost 5.0% of the features, four in
five without a message. cgload reported all 57 planted faults it was expected
to catch and, as predicted, none of the 24 that kept every count. One portal
then served all 30 genomes from one database: searches by name found all
{{portal:genes}} genes sampled, a gene page took a median of
{{portal:gene_page_median}} ms, and none of {{portal:probe_requests}} injection
attempts and random requests caused a server error or changed the database.

**Conclusions:** Each of the other three loaders changed or lost features
without a message. A count of the input by independent code catches such loss
at load time, though not a changed value or name, and a read-only portal puts
the checked database in front of its users. cgload is open source (MIT
license):
https://github.com/computational-genomics-lab/cgload.

**Keywords:** genome annotation; GFF3; relational database; data integrity;
benchmarking; RefSeq; web portal; genome browser

# Background

Genome projects publish their annotations as GFF3 files [@gff3], such as those
of NCBI RefSeq [@refseq]. To ask questions of an annotation, many groups first
load it into a database. gffutils builds an SQLite database from a GFF3 file
[@gffutils]. Chado is a relational schema for genome data, used by many
organism databases [@chado] and filled by loaders such as the one in Tripal
[@tripal]. Bioconductor's txdbmaker builds a TxDb, a database of transcripts
for analysis in R [@granges; @txdbmaker], and data warehouses such as InterMine
combine annotations from many sources [@intermine]. Most people who use an
annotation never open its file. They meet it through a web portal built on
such a database, with a search box, a page for each gene, a genome browser and
sequences to download, as in Ensembl [@ensembl] and the organism databases
built with Tripal [@tripal].

A loader can drop or change features while it reads a file. GFF3 lets one
feature be written on several rows that share an ID, which is how RefSeq writes
coding sequences, so a reader that expects one row per feature, or one ID per
row, gets these features wrong. Some features cross the origin of a circular
sequence, some genes are trans-spliced from distant parts of an organelle
genome, and some strands are unknown. We counted these in the 30 RefSeq
annotations used here (see "What the files contain"). Rows that break a GFF3
rule outright are rare. Rows that follow a documented convention, or carry an
attribute that documents an exception, are more common, and a loader has to
take them as written. A loader can finish without a message and still leave
such features out, or store them in another form. The database then answers
queries, and nothing downstream can see what is missing. A portal built on it
shows its users the reduced annotation as if it were whole.

Large resources test their databases before each release. InterMine's
acceptance tests are SQL queries with expected answers, such as a minimum
number of GO annotations [@intermine-tests]. Ensembl's datachecks test row
counts and metadata in each database before it is handed over [@ensembl;
@datacheck]. Curators write these tests for a resource they know, and the tests
compare the built database with their expectations. As documented, neither
compares a database with an independent reading of the files it was built
from. Recent tools aim at speed: gffbase stores GFF3 and GTF files in DuckDB
and, after a load, checks the database against its own structural rules, such
as unique IDs and parents that resolve [@gffbase], and GFFx indexes GFF files
for fast extraction and assumes well-formed input [@gffx]. Toolkits such as
GenomeTools [@genometools] and AGAT [@agat] check or repair a file before it is
loaded, and do not see what a loader then stores.

A group that has annotated a genome of its own needs a portal as well, and the
options differ in what they need and what they offer. Tripal builds a
community database on the Drupal content management system, with Chado in
PostgreSQL, and needs that stack installed and kept up to date [@tripal;
@tripal-docs]. GenomeHubs sets up an Ensembl database and website for a new
species in software containers [@genomehubs]. Genome browsers such as JBrowse 2
[@jbrowse2] and igv.js [@igvjs] show indexed files directly, without a
database, and find features by the names in those files. Each such site is also a web application on
the open internet. Injection, in which text that a visitor sends reaches a
database or a shell as code, is among the commonest flaws of web applications
[@owasp]. A survey of more than 20,000 bioinformatics websites found that the
most common server version among them was outdated and vulnerable [@tao], and
a fifth
of 2,396 web tools published over a decade never answered during 133 days of
tests [@kern].

We wrote cgload to check every load against an independent reading of its
file, and to put the checked database in front of the people who use it. It
reads the input file twice,
once to load it and once, with separate code, to count it, and it keeps the
load only if the database holds what the count says. Two more commands serve
the database as a portal that only reads: it cannot change the database, runs
no other program, and passes every value a visitor sends to the database as
data, never as code. To find out whether such a check is needed, and whether it
works, we ran a benchmark whose protocol we fixed before the first run. On 30
RefSeq genomes we asked what the files contain that a loader has to handle,
whether cgload stores all of it, what three widely used loaders keep, change
and lose from the same files, which kinds of damage cgload's check notices, and
what the check costs in time and memory. When the portal was added, we asked,
under a section of the protocol written before the measurement, how quickly
and how safely one portal serves all 30 genomes.

# Methods

## The cgload software

cgload is a Python package (3.11 or later) with three runtime dependencies:
SQLAlchemy, pyfaidx and click. The portal adds two more, Flask [@flask] and
waitress [@waitress], installed as an optional extra. cgload is used from the
command line, and a genome passes through four commands. `cgload init`
creates an empty database. `cgload load` reads a short configuration file that
names an assembly and its annotation files, and writes them into the database
in one transaction. `cgload verify` checks a loaded genome against its input
files at any later time. `cgload query`, `export`, `explain` and `stats` read
the database, and `cgload remove` deletes a genome from it. `cgload index` and
`cgload serve` put the database on the web (see "The portal"). Figure 1a shows
every step of a load,
what stops it and what it writes, and Figure 1b the checks of `cgload verify`;
the sections below follow the figure.

The database is reached through SQLAlchemy, and the same code runs on SQLite
and PostgreSQL. SQLAlchemy also reaches MySQL, but cgload has no test on
MySQL, and we do not claim it works there. Beside the database, cgload keeps a
data folder that holds, for each assembly, a copy of the FASTA file with its
index, the annotation as cgload stored it, and a manifest of the load. A
genome is therefore a set of database rows and a folder, and cgload treats the
two as one unit.

### The database

The schema has 11 tables (Figure 2). It is written once, in
`src/cgload/db/schema.py`, and SQLAlchemy creates it on either database system.
Five tables describe the genome and where it came from. `organism`, `strain`
and `assembly` name the genome; the assembly row also records the path and
SHA-256 checksum of the stored FASTA file, the cgload version that loaded it,
and the version of the type vocabulary in force at the time. `source_file`
records every input file, with its checksum and the dialect cgload detected in
it. `sequence_region` holds one row per sequence of the assembly, with its
length and whether it is circular.

The annotation itself is in two tables. `feature` holds one row per row of a
GFF3 file (for a GenBank file, one row per segment of a feature's location).
A feature written on several rows, such as a CDS split by introns,
keeps one stored row per segment: the rows share the feature's ID, and
`segment_index` numbers them in coordinate order. Each row records its file,
its sequence, its coordinates (1-based and inclusive, as in GFF3), strand,
phase and score, and the input line it came from. Each row also has two type
columns. `source_type` is column 3 of the file exactly as written, and is what
`verify` counts. `feature_type` is the canonical term after synonyms are
mapped (for example, an AUGUSTUS `transcript` is an `mRNA`), and is what
queries use. The hierarchy is stored in `parent_id`, which points at the first
segment of the feature's parent. Because `parent_id` holds one value, a
feature that names more than one parent cannot be stored, and cgload refuses a
file that has one. `linkage_method` records how the parent was found: stated
by `Parent=`, or computed, as it must be for GenBank files. The other
attributes of column 9 are rows of `feature_attribute`, one per key and value,
held against the feature's first segment; a later segment whose attributes
differ keeps its own. `ID` and `Parent` are not repeated there, because they
are columns.

Constraints in the schema refuse some impossible rows outright: a start below
1, an end before its start, a strand other than `+`, `-`, `.` or `?`, or a
phase other than 0, 1 or 2. On SQLite, foreign keys are switched on for every connection, so a row
cannot name a parent, file or sequence that does not exist. The indexes are
plain B-trees, the only kind that both database systems share. To find the
rows in a region quickly, each row also stores a `bin`: the smallest bin that
holds it whole, in the hierarchical binning scheme of the SAM and CSI indexes
[@htslib]. A region query reads only the rows whose bins overlap the region.

Two tables hold functional annotation. `annotation_run` records each output
file of eggNOG-mapper [@eggnog] or InterProScan [@interproscan], with how many
of its proteins matched a stored feature. `protein_annotation` holds one row
per assignment, such as a GO term, an EC number, a KEGG ortholog or pathway,
or a Pfam domain with its position in the protein, linked to the feature it
describes. The last two tables, `vocabulary` and
`vocabulary_meta`, are written when the database is created and hold the map
of type synonyms with its version and checksum.

### Creating a database and describing a genome

`cgload init` creates the 11 tables, writes the 46 terms of the type
vocabulary, and gives the database a random identifier. It refuses a database
that already holds cgload's tables, unless told to replace them. `load`,
`query`, `stats` and `remove` first check that the vocabulary stored in the
database is the one the installed cgload expects, and stop if it is not,
because a changed synonym map would change what earlier loads mean. `verify`
needs no such check, because it counts the types as the file wrote them.

A genome is described in a TOML file. It names the organism (genus, species,
strain and, optionally, the NCBI taxon), the assembly (name, version,
accession and FASTA file), one or more annotation files, and any functional
annotation files. The repository's `examples/minimal/` holds such a file with
a small genome; Table 1 shows how cgload stores its first gene.

### Reading and checking the files

`cgload load` first reads and checks every annotation file named in the
configuration, and writes nothing while it does so (Figure 1a, top). A stop at
this stage names the file and the line, and the option that would help where
one exists. There are five steps.

*Tokenise.* Each line is split into its nine columns and its attributes. A line without nine columns, a coordinate, strand or phase that is not valid, or
a start below 1 stops the load. So does a start after its end, unless
`--repair-envelopes` is given; with that option the row is passed on to the
check of parent spans below. GenBank flat files are read
by a parser of their own into the same rows.

*Detect the dialect.* Annotation pipelines write GFF3 differently. cgload reads
RefSeq GFF3, the GFF3 of funannotate [@funannotate] and AUGUSTUS [@augustus],
GenBank flat files, and its own export. Each dialect is a declared profile that
says which feature types it provides, how a child finds its parent, and how a
file identifies itself: by the program named in column 2 of a GFF3 file, by the
LOCUS line of a GenBank file, or by the header line of a cgload export. A file
that no profile claims, or that two profiles claim, is refused with the names
of the profiles tried; the configuration can name the dialect instead. The
dialect, the confidence of the match and its evidence are recorded in
`source_file` and in the manifest.

*Normalise.* Every dialect is brought to one representation. Rows that share an
ID become one feature, with its segments in coordinate order, types are mapped
to canonical terms, and each feature's parent is found. In GFF3 the parent is
the one named by `Parent=`. GenBank has no such attribute, so the parent is
computed: from a `/locus_tag` or `/gene` that parent and child share, from a
`/transcript_id` that an mRNA and its CDS both carry, from a note that names
the transcript, or by exon containment. `linkage_method` records which, and a
CDS that fits no transcript is refused rather than guessed. Normalisation also
refuses a feature with two parents, a parent that is not in the file, a parent
of a type that may not hold the child, a parent on another sequence, and one
ID used for two different features.

*Check parent spans.* A feature must lie inside the span of the parent that
claims it, and no row may run backwards. By default a file that breaks these
rules is refused, and the message names each offending feature and the option
that repairs it. With `--repair-envelopes`, cgload recomputes each parent's
span from its children, leaves first, and writes every change to the load's
manifest. A disagreement that exists only across the origin of a circular
sequence is stored as written and named in a warning.

*Check across files.* An ID defined in two files of one genome, or a file
listed twice, stops the load; an `id_prefix` in the configuration keeps the IDs
of two files apart. In an AUGUSTUS file, each protein written in the comments
must have the length that its CDS implies. Types that the dialect's profile
does not declare are loaded and named in a warning, or refused with
`--strict-vocabulary`.

### Populating the database

`cgload load` then writes the genome in one transaction, in nine steps (Figure
1a, bottom):

1. *Stage the assembly.* A staging directory is made inside the data folder.
   The FASTA file is copied into it, or hard-linked with `--link`, and indexed
   with pyfaidx, which writes a samtools-compatible `.fai` index. A duplicate
   or empty sequence stops the load.
2. *Write the genome rows.* The organism and strain rows are found or created,
   and the assembly row is written with the checksum of the FASTA file, the
   cgload version and the vocabulary checksum. An assembly already loaded under
   the same name stops the load, unless `--force` is given. A row for each input
   file follows, then a row for each sequence of the FASTA index.
3. *Write the features.* The features of each file are written in bulk, one row
   per segment, each with its bin. A sequence name that the FASTA file does not
   contain stops the load, because the annotation then describes another
   assembly, and so does a row that runs past the end of a linear sequence. On
   a circular sequence, a feature that crosses the origin may end past the
   sequence length, as GFF3 allows.
4. *Link the parents.* The parent links are set in a second pass. A child may
   come before its parent in the file, so the parent's row may not exist when
   the child is written. Each stored row is found again by its input line,
   segment and span, and every row of a child is pointed at the first row of
   its parent.
5. *Write the attributes* to `feature_attribute`. The proteins that an
   AUGUSTUS file writes in its comments go to the staging directory as FASTA.
6. *Attach functional annotation*, if any. eggNOG-mapper and InterProScan
   report the protein names they were given, which may be a transcript ID, a
   protein ID or a locus tag, so each name is looked up among the stored IDs
   and the attributes `protein_id`, `transcript_id`, `locus_tag` and a few
   others. A name that matches two features stops the load. By default every
   protein must match, and the configuration can lower that threshold;
   `annotation_run` records how many matched.
7. *End-of-load checks.* The stored counts are compared with the independent
   count of each file (see "Verification"). A short count stops the load.
8. *Write out.* On PostgreSQL, the planner's statistics are collected for the
   new rows. The annotation is written back out as GFF3 in the staging directory,
   compressed as BGZF and indexed in the tabix format [@tabix; @htslib], for
   genome browsers; cgload writes both formats itself, without bgzip or tabix. A manifest records the input files and their checksums, the detected
   dialects, the feature counts, every repair and the results of the checks.
9. *Commit.* The staging directory is renamed into place, and the transaction
   is committed. If the commit fails, the directory is moved back.

An error at any step rolls the transaction back and removes the staging
directory, so the database and the data folder are left as they were. A killed
process cannot clean up after itself. Its transaction still rolls back, but
the staging directory stays in the data folder, and so does the folder that a
`--force` reload had set aside. Neither belongs to a stored assembly, a later
load does not remove them, and they have to be deleted by hand (see
"A killed load leaves the database unchanged"). With `--force`, an assembly already loaded under the same
name is replaced; its folder is moved aside only at step 9 and deleted only
after the commit. A successful load prints the assembly's identifier, the
number and total length of its sequences, and for each file the dialect, the
confidence of the match and the features by type, followed by any warnings. A
stopped load prints the reason and exits with status 2.

**Table 1. How cgload stores the first gene of the example genome.**

<!-- include tables/table0_example.md -->

The rows of table `feature` for gene `DEMO_000001` of `examples/minimal/`, made
by loading that example into a new SQLite database. The CDS is written on two
rows of the file, so it is stored as two rows that share its ID, numbered by
`segment_index`. Each row's `parent_id` is the `feature_id` of its parent's
first row. The mRNA's attribute `product=hypothetical protein` is stored as a
row of `feature_attribute` with `feature_id` 2. All rows are on the + strand;
columns not shown are listed in Figure 2.

### Verification

At the end of the load, cgload counts what it stored and compares the counts
with an independent count of the file. If any count differs, it rolls the
transaction back. A load therefore either stores as many features, rows and
parent links of each type as the file holds, or leaves the database as it was.

The command `cgload verify` repeats these checks at any later time, in five
tiers (Figure 1b). The source-count tier counts the file again with a module
that imports no other part of cgload, so it shares no code with the parser:
features of each type, rows of each feature written on several rows, and
parent statements. The end-of-load check of step 7 uses the same count. These counts must equal the database's,
counted by `source_type`. The parent-edge tier requires every feature that the
file gives a parent to have a stored link, and as many such features of each
type as the file has. The stored-assembly tier compares
the FASTA and its index with the checksums taken at load. Two optional tiers
compare the file with other descriptions of the same assembly: the GenBank
file (cross-parser tier) and the totals NCBI publishes for the annotation
release (reference tier). A tier that cannot run says so and gives the reason.
`verify` exits with status 2 when a tier fails and 0 otherwise.

These checks compare counts. They establish that the database holds as many
features and rows of each type as the file, and as many features with a
parent, and that no feature the file gives a parent was left without one. They
cannot tell which feature is which, or whether a value is right. A renamed
feature, a child linked to the wrong parent, two swapped IDs, a removed
attribute, a wrong strand or a shifted coordinate leaves every count as it was
(see "cgload reported every planted fault that changes a count"), and so does
a feature lost while another of its type is invented.
Nor do these checks read columns 2, 6 and 8 of the file. In the benchmark,
agreement feature by feature was measured separately, by comparing cgload's
export with a reference reading of the file (see "Reference reading and
comparison"); that comparison is not part of `cgload verify`.

### Reading the database

`cgload query` returns features by ID, by region or by functional annotation
term, as a table or as GFF3, and with `--descendants` also everything below
each match, such as a gene's transcripts, exons and CDS. A region query reads
the rows whose bin overlaps the region and keeps those that overlap it; on a
circular sequence it also finds the rows that cross the origin. `cgload export`
writes an assembly's annotation back out as GFF3. It writes column 3 as the
file had it, from `source_type`, and adds no row that the input did not have.
`cgload explain` shows where a feature came from: its segments, the file and
the dialect detected in it, its parent, and how the parent was found. `cgload stats` counts features by their original and
canonical types, and `cgload remove` deletes an assembly's rows and folder in
one transaction. Because the database is plain SQL, it can also be read with
any SQL client.

### The portal

`cgload serve` turns a cgload database into a website for the people who use
its genomes (Figure 3). It is written with Flask [@flask] and served by
waitress [@waitress], a production web server written in Python, and it reads
the same database and data folder as the other commands, on SQLite or
PostgreSQL. One portal serves every genome in the database, and it only reads.

*What a visitor can do.* A feature's page draws its gene with every
transcript to scale and opens the igv.js genome browser [@igvjs] at its locus
(Figure 3a). It then lists the feature's attributes, with cross-references
linked to the databases they name, the functional annotation of its proteins,
and its whole family from gene to exons. A visitor finds features by searching
one genome or all of them, by a gene's name or synonym, a locus tag, a
transcript or protein ID, the words of a product or description, a
cross-reference such as `GeneID:672`, or a term that eggNOG-mapper or
InterProScan assigned: a GO term, a KEGG ortholog, a Pfam or InterPro
accession or an EC number (Figure 3b). Exact IDs are listed first, then exact
names in any case, then names that begin with the search, then features whose
text holds every word of it, then functional matches; within each group, genes
come before transcripts and transcripts before coding sequences. A region page
shows the browser and the genes in a window of up to 5 Mb. The home page lists
the genomes, and a genome's page gives its counts by feature type and the
record of its load, read from the manifest: the checks that `cgload load` ran
and their outcome, any repairs, and the input files with their SHA-256
checksums (Figure 3c).

The portal gives as FASTA a feature's genomic sequence, a transcript spliced,
a coding sequence from its first whole codon, as the phase of its first row
places it, and the protein. The protein is translated with the NCBI genetic
code [@gencode] that the file names in `transl_table`, or with the standard
code when it names none. Where the file declares something that a plain
translation does not apply, such as a selenocysteine or a ribosomal slippage,
the page names it. A JSON interface under `/api/v1/` answers the same
questions for programs, and the four files of each assembly that a genome
browser reads (the FASTA file, its index, the compressed annotation and its
tabix index) are served in byte ranges.

*The search index.* Names and words are looked up in an index that
`cgload index` builds from the features' IDs and attributes: a row for each
lower-cased term and the feature it names, in two tables kept apart from the
eleven. The index is checked by counting, as a load is. For each attribute key
it covers, the number of features that carry the key must equal the number of
features that reach the index under that key, and an assembly whose counts
differ is left without an index. Without one, the portal still finds features
by ID and by functional term, and says that names cannot be searched.

*Safety.* The portal's defences are limits on what it can do. It opens the
database read-only, SQLite in read-only mode
and PostgreSQL in read-only transactions, so no request can change it, and it
answers only GET and HEAD requests. Every value that a request sends reaches
SQL as a bound parameter, and its length is limited. The portal runs no other
program, keeps no session and has no accounts. Text from the annotation files
is escaped wherever a page shows it, and every response carries a Content
Security Policy [@csp] that lets a page run scripts from the portal only, and
no inline script. igv.js is served from the cgload package itself, so a page
loads nothing from another site. Of each assembly's folder, only the four
browser files can be read, by name. The test suite sends injection strings and
a seeded fuzz of 3,000 random requests to the portal's routes and checks that
none fails the server or comes back as code, and continuous integration runs
the Bandit security linter [@bandit] on the code. `docs/portal.md` in the
repository describes how to run the portal behind a web server that adds
HTTPS.

## Benchmark design

The benchmark follows a protocol (`benchmark/PROTOCOL.md` in the repository)
that we committed and tagged `benchmark-protocol-v1` on 16 October 2026, with
every script it names, five days before the first run. The protocol states the
questions, the genomes, the software, the hardware, the definitions and, for
every planted fault, the outcome expected from what cgload claims to check.
Benchmarks written by a tool's own authors tend to favour that tool [@weber;
@buchka], and fixing these things first was our way to limit that. A change
made after the first run is listed in the protocol with its date and reason,
and every measurement it affects was repeated (see "Changes after the first run"). The protocol asks seven questions:

- Q1. Does cgload load and verify all 30 genomes, on SQLite and on PostgreSQL,
  and at what cost in time, memory and disk?
- Q2. Does what cgload stores match the input, feature by feature?
- Q3. What do other loaders keep, change or lose from the same files, and do
  they say so?
- Q4. Which damage to a loaded genome, or inside the parser, does cgload notice?
- Q5. What does an interrupted load leave behind?
- Q6. How long does a user wait for everyday answers?
- Q7. Do cgload's structural checks agree with an independent census of the
  same files?

It also listed an eighth question, a case study on one of our own genomes with
its functional annotation. No such genome was available, and the question was
dropped before any of it was run. A ninth question came with the portal, after
these runs:

- Q9. Can one portal serve all 30 genomes from one database, how long do its
  users wait, and does anything they send make it fail?

We wrote its procedure into the protocol as section 13, with what counts as a
failure of the portal, before the portal was measured.

### Genomes and input files

The panel is 30 RefSeq genomes: 8 vertebrates, 7 insects, 7 plants and 8 fungi
(Additional file 1). They were chosen to span clades, assembly sizes and
annotation ages, and the list was fixed in a tagged commit (`survey-panel-v1`,
6 October 2026) before any of them was read. Together the annotations have
28,211,625 GFF3 rows, 18,017,618 features and 16,566,107 parent links. A
parent link here is a feature and the parent it names, however many of the
feature's rows name it; counted per row, as the census counts them, the files
hold 26,599,375 parent statements.
The smallest is *Candida albicans* (25,889 rows) and the largest is human
(4,932,571 rows, 3,203,237 features). Every file is identified by its MD5,
which the benchmark checks before each run. One download, medaka, had been cut
short, and this check found it; the file was fetched again and every count of
it redone (`survey/results/CORRECTIONS.md`). The assembly FASTA files, and the
GenBank files of the eight fungi, came from NCBI Datasets [@datasets].

### Census of the files

Before the benchmark we counted what the 30 files contain, with a script
(`survey/cgsurvey.py`) that reads a GFF3 file as text, line by line, with the
Python standard library only. It shares no code with cgload and is not
installed with it. For each file it counts the rows, the parent statements, the
sequences and the rows that repeat an ID already seen, and it flags four kinds of
row: a feature outside the parent that claims it, a start after its end, a row
past the end of its sequence, and a strand written `?`. It counts separately
the rows that follow two documented conventions, an exon placed directly under
a pseudogene and a row that crosses the origin of a circular sequence, which
GFF3 and NCBI both write by letting the end run past the sequence length. For
every flagged row it records the raw GFF3 line and any attribute that documents
an exception, such as `exception=trans-splicing` or `partial=true`, so that
each can be read by hand.

We then classed every flagged row into one of four: (A) valid GFF3 that a
loader must accept, (B) a documented convention or an exception the file itself
states, (C) an error in the file, and (D) a finding no convention explains.
Classes A and B are what a loader has to handle; C and D are what a file ought
not to contain.

### Reference reading and comparison

A reference reader (`benchmark/truth.py`) reads each GFF3 file with the Python
standard library only and shares no code with cgload or with any other loader.
It follows the GFF3 specification [@gff3]: rows that share an ID form one
feature, a row without an ID is a feature of its own, and `Parent` names the
parents. For each feature it records the sequence, the type, the strand and
span of every row, the parents, and a digest of all other attributes (keys and
values sorted, escapes decoded). White space at either end of a tag=value pair
is not part of it. It does not record columns 2, 6 and 8 of a row (the program
that wrote it, its score and its phase); for cgload these were compared
separately, row by row, as described below.

A loader's result is reduced to the same record and compared with the
reference feature by feature. For cgload, the record is read from the GFF3 file
that `cgload export` writes from the database, so the comparison covers the
round trip. For gffutils and Chado it is read from their databases with SQL. A
feature is kept when the loader holds it once and every compared field equals
the file's. It is changed when the loader holds it with a field that differs,
or holds it as several records. It is lost when the loader does not hold it. A
change that the loader reported in its own record is declared (cgload lists
each repair in its manifest). A change or loss is announced when the loader's
output names the feature, or names its type on a line saying it was skipped,
dropped, ignored, rejected or unknown. Otherwise it is silent.

The fields compared are sequence, type, strands, spans, parents and attributes
for cgload and gffutils. Columns 2, 6 and 8 are not among them, so for the
other loaders the counts of changed features are lower bounds. For Chado the attributes are left out, because
python-chado spreads them over several tables. txdbmaker holds transcripts, not
features, so it is compared per transcript: sequence, strand, span, gene, exon
spans and CDS spans. Transcripts are paired by the name txdbmaker gives them,
or failing that by span, and features of types txdbmaker does not model are
counted apart, not as lost. A feature whose rows carry different attributes is
compared row by row. Each feature is also placed in the hard cases the
protocol defines, such as a feature written on several rows, a trans-spliced
feature, a feature on an organelle or one that crosses the origin. Twelve of
these cases occur in the 30 files.

For cgload we also compared columns 2, 6 and 8, row by row. A second program
(`benchmark/columns.py`), which likewise uses the Python standard library only
and shares no code with any loader, reads an input file and cgload's export of
it, pairs the rows of each feature in coordinate order, and compares column 2
as text after percent-decoding, column 6 as a number and column 8 as text.
`benchmark/run_columns.py` repeats cgload's part of the round trip on SQLite
for each of the 30 genomes. This check was added after all other runs, outside
the protocol (see "Changes after the first run"), and ran with version 0.2.2 on
20 November 2026. It also counts, in each input file, the values that version
0.2.1 writes differently (see "Loaders and versions"): scores with more than
six significant digits, and rows whose column 2 differs from the first row of
their feature.

### Loaders and versions

Table 2 lists the versions used. gffutils was run twice: with its default merge
strategy, which is what a user who changes nothing gets, and with
`create_unique`, the documented way to load a file whose rows repeat an ID
[@gffutils]. txdbmaker's `makeTxDbFromGFF` was told which sequences are
circular [@txdbmaker]. For Chado we used python-chado [@pychado], the
`load_gff` command used by Galaxy Genome Annotation, which reads GFF3 with
Biopython [@biopython], into a new Chado database with the Sequence Ontology
[@so]. Chado has other GFF3 loaders, such as the GMOD bulk loader and Tripal's
importer, and our results stand for python-chado only. Each load had a time
limit of 12 hours, or 24 hours for Chado. A load that was refused, crashed, or
ran out of memory or time is reported as not loaded; none of those loads left
anything in its database.

cgload was benchmarked at two commits of its repository, `8c999e2` and
`9b97178`, and released as version 0.2.1. The code that loads, verifies and
exports is the same in all three, apart from the version string and comments.
Version 0.2.2 differs from 0.2.1 only in how it keeps and writes back columns
2 and 6. The export of 0.2.1 wrote scores with
six significant digits, so 1234567.5 came back as 1.23457e+06, and it wrote a
feature's first-row column 2 on all its rows. Neither changes a count. Version
0.2.2 keeps both as loaded (decisions D-109 and D-110 in the repository).
Version 0.3.0, the one this paper describes, adds the portal and its search
index (D-111 to D-114). It loads, verifies and exports as 0.2.2 does, and its
`remove` also deletes a genome's index.
The query code changed between the two commits (see "Changes after the first
run"), and the timing was measured after that change.

**Table 2. Software versions used in the benchmark.**

<!-- include tables/table6_versions.md -->

Versions are those recorded in the server's environments at the time of the
runs; the Chado schema and the Sequence Ontology are the files that
`benchmark/run.sh` fetches, at the commits shown. pyfaidx is described in
[@pyfaidx] and R in [@rcore]. GenomicFeatures is RRID:SCR_016960 and HTSlib is
part of SAMtools, RRID:SCR_002105.

### Faults, interrupted loads, timing and agreement

For Q4 we planted 27 faults, each in three genomes (*F. graminearum*, *A.
thaliana* and human). Faults F01 to F19 damage a copy of the SQLite database or
its data folder after the load: deleting a feature or a row, removing a parent
link, editing a base of the stored FASTA, or editing the input file. Faults P01
to P08 damage the load itself, through a wrapper around cgload's tokenizer that
drops, adds or alters rows. The protocol states each fault's expected outcome.
A fault that changes a count, the stored assembly or the input file should be
detected, and one that changes only a value, an ID or a parent, such as a strand
or a coordinate, should not (Additional file 2).

For Q5 we killed loads with SIGKILL, sent to the whole process group at a time
drawn uniformly between 5% and 95% of the duration of an uninterrupted load,
from a seeded generator. cgload and gffutils had 10 trials per genome on *F. graminearum* and *A.
thaliana*, and Chado 5 trials on *F. graminearum*; txdbmaker was not part of
this test. For
cgload we killed new loads and `--force` reloads over an earlier load. After
each kill we recorded what the database held, counted in the loader's own
tables and compared with an uninterrupted load.

For Q6 we timed four tasks on four genomes spanning 0.1 to 4.9 million rows
(*F. graminearum*, *A. thaliana*, maize and human): T1, one gene and every
feature below it; T2, every feature overlapping a 1 Mb window; T3, feature
counts by type; and T4, every feature on one sequence, written as GFF3. Targets
were drawn with a seeded generator: 20 genes that have children, 20 windows on
sequences of at least 1 Mb, and three sequences, the one with the most features
and two drawn from those with at least 1,000. Each target was asked 5 times, in
two ways: one new process per answer, as a user typing a command waits, and one
process asking T1 to T3 in turn, as a script waits. As commands, cgload was
asked with `cgload query --id GENE --descendants`, `cgload query --region
SEQUENCE:START-END`, `cgload stats --by source` and `cgload query --region
SEQUENCE --format gff3`, on SQLite and on PostgreSQL. gffutils was asked
through its Python interface (`children`, `region`, `count_features_of_type`),
and tabix on cgload's compressed annotation for T2 and T4, with T3 as a full
scan. From one process, cgload was asked through the functions behind those
commands. We report the median of the repeated asks; the first asks are in
Additional file 3, and the sequences drawn for T4 in Additional file 4. The server did not let us drop the page cache, so no
ask started from a cold disk.

For Q7 we took every row the census flagged, with the class the census gave
it, and ran cgload's structural checks on the same files, without a database,
as `cgload load --repair-envelopes` runs them. The two find problems by
different routes: the census reads each file as text and checks rows, while
cgload checks features it has already parsed and assembled.

### Measuring the portal

For Q9, cgload 0.3.0 loaded all 30 genomes into one SQLite database, as a
group's portal would hold them, with `--repair-envelopes` where a file was
refused, as for Q1, and `cgload index` built and checked the index of each.
The portal ran as `cgload serve` with 8 threads in one process, and a client
on the same server asked it, over the loopback interface, what a user asks.
For 200 genes with a `Name` in each genome, drawn with a seed, the client
searched by the gene's name and opened the gene's page, its region (the gene
with 5 kb on either side, at most 100 kb), the protein of its first coding
sequence, and a 64 KiB range of each of the two files that the genome browser
reads, at a seeded offset. It recorded the status, size and time of every
response, and on which page of results, if any, the search listed the gene.
The same requests were then sent by 1, 8 and 32 clients at once, for 60 s at
each level. Up to 2,000 coding sequences per genome, drawn with a seed, were
translated as the portal translates them, and every protein with a stop codon
before its end was checked against what the file declares: a translation
exception, a partial coding sequence or a pseudogene. A security probe sent
injection strings and 3,000 seeded random requests to the running server,
recorded every status and checked every response for the security headers.
The SHA-256 checksum of the database file was taken before the portal started
and after its last read. Human was also loaded, indexed and asked the same
requests on PostgreSQL. The protocol counts as a failure of the portal any
response with a status of 500 or above, a search by a gene's name that does not
list the gene while no more than 1,000 features match the name, a stop codon
inside a protein that the file does not explain, a response without the
security headers or with a script it was sent, and any change to the database.

### Hardware

All runs used one lab server with Linux 6.5, 80 CPUs, 125 GiB of memory and 8.7
TiB of disk, and PostgreSQL with default settings. Loads that were timed for Q1
and Q6 ran alone on the machine. The other loaders ran several at a time (8
jobs for gffutils, 6 for Chado and 3 for txdbmaker), so their times are
indicative only. Q9 ran on the same server. Wall time, CPU time and peak memory were measured for every
command. Memory and sizes are given in binary units: 1 GiB is 2^30^ bytes, 1
MiB 2^20^ and 1 KiB 2^10^.

### Changes after the first run

The first run, on 21 October 2026, is void: two copies of the run script were
started in one folder and wrote over each other's results. It still showed
problems, and we fixed each with tests before the second run. cgload refused
six RefSeq files over valid content, failed on one attribute value that was too
long for a PostgreSQL index, and changed some features on the way back out
through `cgload export`. Its cross-parser tier failed on five fungi over names
that NCBI writes differently in GFF3 and GenBank. A killed load could leave a
directory that blocked the next load. Two preflights on the 14 genomes involved
then found more of the same kind, and also showed that the medaka file we had
counted was incomplete.

The second run (23 and 24 October) finished every step. It found two errors in
our scripts, in the Chado database setup and in txdbmaker's pairing of
transcripts, and one cgload defect: a gene lookup on PostgreSQL took 12 s for
human. We fixed these, and repeated Chado, txdbmaker and the timing on all 30
genomes. The human Chado load then reached its 24-hour limit with nothing
stored, and we finished its comparison without loading again. Last, we found
that the comparison had charged gffutils with 133 attribute changes in *S.
pombe* that were ours: gffutils keeps a space at the end of a value, and the
reference reader dropped it. We fixed the comparison and measured gffutils
again on all 30 genomes. The protocol lists each change with its date, its
commit and the tests added.

cgload was therefore corrected on the same 30 genomes it is measured on
(see Limitations). The results reported here come from three tagged commits
(`benchmark-run-2`, `benchmark-run-2.1` and `benchmark-run-2.2`), and the
protocol says which result comes from which.

One check was added after all runs, outside the protocol: the row-level
comparison of columns 2, 6 and 8 described above. On the test files of the
repository it found the two faults in the export of version 0.2.1 that version
0.2.2 corrects. On the 30 genomes it was run once, with version 0.2.2.

The portal and Q9 came last. Section 13 of the protocol, which sets out how
the portal is measured and what counts as its failure, was committed with the
code of version 0.3.0, before the portal was run on the 30 genomes.

## Use of large language models

We used Claude, a large language model made by Anthropic [@claude], mainly the
model Claude Opus 5.5, from 20 August to 25 November 2026. Under the authors'
direction it wrote and changed code, tests and documentation of cgload,
including its portal, and of the census; drafted the benchmark protocol with
us and wrote the benchmark scripts; analysed the results of each run; wrote
the scripts that make every figure and table of this paper from the results;
searched the literature for related tools; drafted the text of this
manuscript; and checked drafts and the submission files for errors. No figure
or number was made by asking the model for it directly. Figures and tables are
the output of `paper/figures.py`, and Figure 3 of `paper/screenshots.py`,
which photographs pages of the portal in a browser and crops them as its
legend says; no image was edited by hand.

We also asked ChatGPT (OpenAI) [@chatgpt] to review drafts of this paper, with
its data and code, for errors. Its review led to corrections in the text and to
the row-level comparison of columns 2, 6 and 8, which found the two faults that
version 0.2.2 corrects. We checked every point that either model raised against
the code and the results before changing anything.

The authors ran every benchmark step on their own server, from the scripts in
the repository, and every result in this paper comes from those runs. The
authors reviewed the code, the results and this text, and are fully responsible
for all of them. Claude's work is public in the repository's history: from 14
October 2026 each commit it wrote names it in the commit message (87 of the
100 commits up to the version submitted), and its draft of this manuscript is
`paper/manuscript.md`. Additional file 5 gives the details.

# Results

## What the files contain

A loader has to cope with whatever real files contain, so we first counted what
is in the 30 files. The annotations hold 28,211,625 GFF3 rows on 23,793 sequences, of which 33
are organellar. Writing one feature on several rows is not an edge case:
1,072,125 features are written on more than one row, and the rows beyond the
first of each number 10,194,007, over a third of all rows. Exons placed
directly under a pseudogene, which some readers do not expect, occur 91,009 times in 20 of the
30 files. No file had a dangling parent statement, an ID used for two different
features, or a non-numeric coordinate.

The census flagged 110 findings on 102 distinct rows, about 4 per million. Most
are things a loader must accept: 67 rows lie outside the parent that claims
them, 63 of them in genes whose own attributes state trans-splicing or RNA
editing, and 31 rows carry the strand `?`, which GFF3 allows. Of the rest, 9
rows run past the end of their sequence and 3 have a start after their end.
Sixteen findings in 5 files carry no documented exception: 12 on plant
organelles, and 4 in *Schizosaccharomyces pombe*, the only ones outside an
organelle, where a gene's span leaves out part of a child. A further 9 rows
cross the origin of a circular sequence, which the census counts as a
convention. Nearly everything that looks wrong at first reading is therefore
valid GFF3 or a documented convention that a loader has to store as written,
but not everything: the 16 findings without a documented exception remain.

## cgload loaded every genome and kept every feature

cgload loaded all 30 genomes on both SQLite and PostgreSQL (Q1; Additional
file 1). Twenty-three loaded as they are. Seven were refused at first, each over
features lying outside the parent that claims them, or a start after its end,
and loaded with `--repair-envelopes`, which repaired 12 features in all and
listed each one. After every one of the 60 loads, `verify` passed the
source-count, parent-edge and stored-assembly tiers, and for the eight fungi the cross-parser tier as well. The reference tier did
not run: it needs NCBI's published counts, which the benchmark did not supply.

The round trip (Q2) kept 18,017,606 of the 18,017,618 features unchanged in every
field compared. The other 12 were the repaired ones, changed in their spans as
the manifest declared. No feature was lost, none changed without a declaration,
and all 16,566,107 parent links came back, each to the parent the file names.
SQLite and PostgreSQL gave the same result for every genome.

Columns 2, 6 and 8, which the round trip does not read, were compared
afterwards with version 0.2.2 (see "Reference reading and comparison"). It
wrote back the source, score and phase of every row in all 30 genomes as the files had them:
28,211,625 rows, with no difference in any of the three columns. The 12 rows
whose span changed are the 12 declared repairs. No feature in the files has
rows from different programs, but 52 scores have more than six significant
digits, 49 in human and 3 in mouse. In the 13 that the comparison kept as
examples, the score lies between 1.1 and 3.2 million and within 41 of the
row's length. The export of version 0.2.1, which the benchmark ran, writes six significant
digits, so it changed these 52 scores without a message: 1,231,658 became
1.23166e+06. It changed no other value in these columns.

## cgload's structural checks agree with the census

Seven of the 30 files were refused at first. To see whether cgload objects to
the right rows, we ran its structural checks on the rows the census flagged
(Q7). cgload accepts the 9 rows that cross the origin of a circular sequence,
as GFF3 allows. Of the other 102 rows (Table 3), both checkers accepted the 23 rows
whose only finding is the strand `?`, and both objected to the 13 the census
classed as errors or unexplained, though cgload loaded 6 of these with a
warning. For the other 66 rows, nearly all in trans-spliced organelle genes,
cgload refused the file unless `--repair-envelopes` was given. With it, cgload
repaired 12 features. Sixty-three of the flagged rows are children of those 12
features and were accepted once their parent was repaired, up to 10 rows under
one transcript; 8 flagged rows are themselves among the repaired features. No
row that the census classed as an error was accepted by cgload without comment.
cgload objected to 8 features that the census did not flag. Each is a gene or
pseudogene whose span leaves out a child that the census did flag, so the two
found the same 11 parent and child pairs from opposite ends.

**Table 3. cgload's structural checks and the census on the same rows.**

<!-- include tables/table5_agreement.md -->

Rows are those the census listed in the 30 files, by the census's classes.

## Other loaders changed or lost features, mostly without a message

Given the same files, the other three loaders changed or lost features, and
mostly did not say so (Q3). Table 4 sums up the results and Figure 4 shows them
per genome.

**Table 4. What each loader kept, changed and lost.**

<!-- include tables/table1_loaders.md -->

Counts are features, except for txdbmaker, which holds transcripts, and cover
the genomes each loader loaded. cgload gave the same result on SQLite and on
PostgreSQL, and loaded 7 of its 30 genomes with `--repair-envelopes`. A change
is declared when the loader recorded it, announced when its output names the
feature or its type, and silent otherwise. The nine genomes python-chado did
not load hold 7,244,078 features, and the six txdbmaker did not load are all
plants; neither is counted in its column.

gffutils with its default settings stops at the first ID written on two rows
(`ValueError: Duplicate ID`), the documented behaviour of its default merge
strategy, and every one of the 30 files has such IDs. With `create_unique` it
loaded all 30 and lost no feature, but stored each of the 1,072,125 features
written on several rows (6.0% of all features) as separate records, the second
and later ones under new IDs. A query for such a feature by its ID returns one
part. The share ranged from 1.2% in *Saccharomyces cerevisiae* to 13.1% in
*Aspergillus nidulans* and *Neurospora crassa*. gffutils printed nothing about
it. It changed no attribute and kept every parent link.

txdbmaker stopped on six of the seven plant genomes. On maize, *Arabidopsis
thaliana* and soybean it stopped at a row whose end precedes its start, an
inverted mitochondrial *nad1* transcript in each. The census classes these
three rows as errors in the files, and cgload refuses the same three files at
the same rows unless `--repair-envelopes` is given. On rice, tomato and
grapevine txdbmaker stopped while building the TxDb, reporting 2 exons in each
file linked to transcripts it could not find. On the 24 genomes it loaded, it
kept 1,169,686 of 1,169,826 transcripts. It changed 104 and dropped 36, and its
messages named all but 3 of the changes. Another 647,583 features have types
that a TxDb does not model.

python-chado stopped with an error about a feature's location on all seven
plant genomes and on *Drosophila melanogaster*, each time at a trans-spliced
gene: organelle genes in the plants and CG32491 in the fly. The human load
reached its 24-hour limit. None of these nine loads stored anything. In the 21
loads that finished, python-chado lost 540,255 of 10,773,540 features (5.0%),
of three kinds. It dropped all 109,181 features of type lnc_RNA, a type the
Sequence Ontology now names lncRNA (SO:0001877) [@so], and its log named the
type. It dropped their 430,836 exons without a message. It also dropped 238
long terminal repeats that the *S. pombe* file writes in capitals, and named
that type too. So four in five of its losses were silent, from none in 9
genomes to 14.6% of the features in *S. pombe* (Figure 4). These loads kept
9,627,024 of 10,167,041 parent links, and python-chado added 670,102
polypeptide features that the files do not contain.

The hard cases show where the alterations come from (Table 5). gffutils split
every feature written on several rows, and these include the features of most
other hard cases. txdbmaker's changes and losses fell on rows that overlap, on
unknown strands and on trans-spliced and organelle transcripts, and its
messages named them. In the loads that finished, python-chado kept every
hard-case feature; the genomes with trans-spliced genes or unknown strands were
the ones on which it stopped. cgload kept every hard-case feature or declared
its repair.

**Table 5. Hard cases and what each loader altered.**

<!-- include tables/table2_cases.md -->

For each hard case the protocol defines: the number of features in the 30
files, what each loader altered (changed or lost) and how many of those
alterations were silent. A feature can be in several cases, and nearly every
feature has an attribute with several values. Chado's counts cover its 21
finished loads and txdbmaker's the transcripts of the 24 genomes it loaded,
each with its own total; a dash means the case does not occur there.

Table 5 counts features, and the census above counts rows, so their numbers
differ for the same cases. The census reads rows as they are written. The
reference reader behind Table 5 first joins the rows of each feature and puts
each row's coordinates in order. The 31 rows with strand `?` belong to 30
features, because one gene is written on two rows. The 67 rows outside their
parent belong to 47 features, since coding sequences written on several rows
are flagged row by row. The reference reader counts 46. The two counts need not
agree, because the census compares each row with its parent as written, and the
reference reader compares whole features after putting coordinates in order:
it finds 8 such features in maize where the census finds 10, and 13 in *A.
thaliana* where the census finds 12. All 9 rows that the census counts as past the end
of a sequence lie on circular organelle sequences. The census keeps them apart
from its 9 origin-crossing rows because they do not have the form of that
convention. The reference reader has one rule for both: on a circular
sequence, a feature that begins before the end of the sequence and runs past
it crosses the origin, and one written wholly beyond the end is past the end.
Its 14 origin-crossing features are therefore the 9 that the census counts as
such, the 3 transcripts written with start after end, and 2 coding sequences
with one or two rows written beyond the end. Its 3 features past the end are 3
exons written wholly beyond it.

## cgload reported every planted fault that changes a count

We then damaged loads and loaded genomes on purpose, to see what cgload's
checks notice (Q4). All 81 fault trials (27 faults in 3 genomes) ended as the protocol predicted
(Table 6, Additional file 2). `verify` detected every fault that deleted,
duplicated or retyped a feature, removed a parent link, or damaged or moved
the stored assembly. Every fault planted in the parser that dropped, added or
retyped a row stopped the load: one by its structural check (a child whose
parent was dropped), the others by the end-of-load counts, which rolled the
load back. As stated in advance, cgload did not notice the eight faults that
change a value, an ID or a parent but no count: a flipped strand, a shifted start, a renamed feature, a
parent pointed elsewhere, a deleted attribute, and in the parser a flipped
strand, a shifted start and two swapped IDs.

One fault needs a caveat. When the input file was edited after the load
(F16), `verify` named the file, said it had changed since the load, and skipped
the source count. It still exited with status 0, because a skipped tier is not
a failure. A script that reads only the exit status would miss this edit, and
Table 6 counts it as reported but not as a failing exit.

**Table 6. Planted faults, grouped by what they do.**

<!-- include tables/table3_faults.md -->

Each fault was planted in *Fusarium graminearum*, *A. thaliana* and human.
"Reported" counts the trials in which cgload reported the fault, and "Exit 2"
those in which its command also exited with a failing status. "Caught by"
names the check that reported it.

## A killed load leaves the database unchanged

A load can also be cut short, by a crash, a full disk or a user (Q5). Every
killed cgload load left the database as it was (Table 7). After the 20
killed new loads the database was empty, and the same load then succeeded each
time. After the 20 killed `--force` reloads, the earlier load was intact and
passed `verify`. In 28 of the 40 trials, cgload left files on disk outside any
stored assembly, from its staging folder or from the folder a `--force` reload
sets aside (up to 124 MiB). After a killed new load, the next load did not
remove them. We list this among cgload's known limitations.

python-chado commits only at the end, so its five killed loads left an empty
database. gffutils left a database file after all 20 kills, and gffutils could
open none of them. Seventeen were empty. Three held every record of a complete
load, because the kill came after the last insert but before gffutils wrote the
metadata it needs to open the file.

**Table 7. What a killed load left behind.**

<!-- include tables/table4_kills.md -->

Each kill came at a random time between 5% and 95% of an uninterrupted load.

## Time, memory and query waits

cgload's checks and its storage cost memory and time (Q1), and it answers
everyday questions more slowly than gffutils (Q6).

Load time, memory and size grew in proportion to the number of rows (Figure 5).
On SQLite, cgload loaded a median 7,961 rows per second; the quickest load, *S.
cerevisiae*, took 4.8 s and human took 641 s. PostgreSQL was about three times
slower (median 2,700 rows per second; human 1,978 s). `verify` took 0.8 to 29 s
on SQLite and 1.5 to 48 s on PostgreSQL. The load held about 6 KiB of memory per
row (median 5.7 KiB on SQLite and 6.9 KiB on PostgreSQL), so human needed
29.5 GiB and 36.4 GiB at its peak. Its database took 3.0 GiB on SQLite and
5.9 GiB on PostgreSQL, and its data folder 3.2 GiB.

One gene with its parts took 0.26 to 0.30 s as a command on SQLite and 0.37
to 0.41 s on PostgreSQL, against 0.13 to 0.24 s for gffutils (Figure 6,
Additional file 3). Most of that time is
the start of a new Python process: asked from one running process, the same
lookup took 4 to 5 ms on SQLite, 8 to 25 ms on PostgreSQL and 1 ms or less in
gffutils. For a 1 Mb window, tabix answered in 8 to 36 ms, at least six times
faster than any of the databases asked as a command; from a running process
gffutils answered about as fast as tabix on maize (17 ms). Counting features
by type took cgload longest: 3.8 s for human on SQLite and 19.5 s on
PostgreSQL, against 0.3 s for gffutils. The two answers differ, since gffutils
counts records and cgload counts features.

For T4, writing out every feature on one sequence, the protocol takes the
sequence with the most features in each genome and two others drawn at random
from those with at least 1,000 features. In human the two drawn were small
scaffolds, so a median over its three sequences would hide the cost. Figure 6
therefore shows the largest sequence of each genome, and Additional file 4
all twelve. Writing human chromosome 1 (437,536 rows) took 25.0 s on SQLite,
88.6 s on PostgreSQL, 14.4 s with gffutils and 0.5 s with tabix.

## One portal served all 30 genomes

cgload 0.3.0 loaded all {{portal:genomes}} genomes into one SQLite database in
{{portal:load_hours}} h (Q9), {{portal:repaired}} of them with
`--repair-envelopes`, as for Q1. `cgload index` then built
{{portal:index_terms_million}} million index terms for
{{portal:index_features_million}} million features in
{{portal:index_minutes}} min, and its check passed for all
{{portal:index_pass}} genomes. With its index, the database took
{{portal:database_gib}} GiB.

The client sent {{portal:requests}} requests about {{portal:genes}} genes
(Table 8). Every one was answered without an error and within a second. A
search by a gene's name took a median of {{portal:search_median}} ms, a gene's
page {{portal:gene_page_median}} ms and a region's page
{{portal:region_page_median}} ms; 95% of all requests took
{{portal:all_p95}} ms or less, and the slowest took {{portal:gene_page_max}} ms.
On PostgreSQL, a human gene's page took a median of
{{portal:gene_page_pg_median}} ms and a search {{portal:search_pg_median}} ms.
The search by name listed every one of the {{portal:genes}} genes,
{{portal:found_first_page}} of them on the first page of results. The other
{{portal:found_later}} were tRNA genes of zebrafish and *X. tropicalis*, whose
names are each shared by {{portal:found_later_matches_min}} to
{{portal:found_later_matches_max}} features, and they appeared by page
{{portal:found_later_max_page}}.

Users who come at the same time share the one process. One client asking
alone received {{portal:c1_per_second}} answers a second, with a median wait of
{{portal:c1_p50}} ms. More clients did not raise that much:
{{portal:c8_per_second}} answers a second with 8 clients and
{{portal:c32_per_second}} with 32, when the median wait was
{{portal:c32_p50}} ms and the 95th percentile {{portal:c32_p95}} ms. No request
failed.

Of {{portal:cds_checked}} coding sequences translated,
{{portal:cds_no_internal_stop}} ({{portal:cds_clean_percent}}%) gave a protein
without a stop codon before its end. The file explained each of the other
{{portal:cds_explained}}: {{portal:cds_stop_declared_exception}} declare an
exception to plain translation, and {{portal:cds_stop_pseudo}} are marked as
pseudogenes. None was unexplained.

The security probe sent {{portal:probe_requests}} requests to the running
portal. It answered {{portal:probe_200}} normally and redirected
{{portal:probe_redirects}}, refused {{portal:probe_400}} as malformed, answered
{{portal:probe_404}} with "not found", and refused the {{portal:probe_405}}
that used a method other than GET or HEAD. None caused a server error, every
response carried the security headers, and none returned a script it was sent.
The database file was unchanged: its SHA-256 checksum was the same before the
portal started and after its last read.

**Table 8. How long the portal took to answer.**

<!-- include tables/table7_portal.md -->

Times run from the client's request to the end of the response, measured on
the server itself, for 200 genes in each of the 30 genomes held in one SQLite
database (`cgload serve`, 8 threads, one process), and for 200 human genes on
PostgreSQL. A gene with no coding sequence has no protein to ask for.
Percentiles are nearest-rank. Additional file 6 gives the medians per genome.

# Discussion

On the same 30 RefSeq files, each of the other three loaders stored some
features in another form, or not at all, without a message: gffutils 1,072,125
features, python-chado 430,836 in its finished loads and txdbmaker 3
transcripts. cgload kept every feature in every field compared and declared its
12 repairs. The causes of the silent changes differ. RefSeq writes coding
sequences on several rows that share an ID, as GFF3 allows.
gffutils refuses such files by default, which is a safe choice, and with the
setting its documentation offers for them it stores each row as a feature of
its own. A user who then asks for a coding sequence by its ID gets one segment
of it, and nothing in the database says that the rest exists under other
names. python-chado stopped at trans-spliced genes, and where it finished, it
dropped a whole class of transcripts over the spelling of a type name and
their exons along with them. txdbmaker was the most careful of the three: it
stopped on the files it could not model, three of them files that cgload also
refuses by default, and it named nearly every transcript it changed or dropped.

The other tools do some things better. gffutils loaded the human annotation
with 140 MiB of memory, where cgload needed 29.5 GiB, and it answered lookups
faster. Each tool makes a choice about rows that share an ID, and the trouble
is that the choice is silent. A benchmark of this kind can only be fair to each
tool's documented use, and we ran each as its documentation describes; other
settings, or other Chado loaders, may do better.

cgload adds a check that does not trust its own loader. The source count
is made by code that shares nothing with cgload's parser, so an error in the
parser cannot hide itself in the check. The faults show both sides of this.
Every fault that changed a count was caught, including those planted inside the
parser, which a check built from the same parser could not see. Every fault
that changed only a value, an ID or a parent went unnoticed, as the protocol
predicted. A user who needs those checked has to compare the database with the
file feature by feature, as our round trip does, and that comparison is not
part of `cgload verify`. Checking identity as well, for instance by comparing
the set of IDs and of child and parent pairs, would close much of this gap at
little cost; it is not in version 0.3.0. The same blind spot hid a fault in
cgload itself. Its export rounded scores to six significant digits, which
changed 52 scores in human and mouse and no count; only the later comparison
of columns 2, 6 and 8 showed it.

For anyone who builds a database from an annotation file, the lesson does not
depend on cgload. The database should be compared with the file it came from,
by code that does not share the loader's reading of the file. A count is cheap
and catches loss; a comparison feature by feature also catches changes. cgload
makes the first part of every load. The reference reader and comparison scripts
of this benchmark, which already read gffutils, txdbmaker and Chado databases,
can be adapted for the second.

The portal is where this check reaches the people who use a genome. A
biologist who searches for a gene does not read a load log, and a portal on a
database that lost features gives no sign of it. Had the Chado databases that
python-chado filled been put behind a portal, it would have offered the 21
genomes they hold without their 109,181 long non-coding RNAs and the 430,836
exons of those RNAs, and no page would have said so. cgload's portal shows on each
genome's page how the genome was loaded and what the checks found, and every
feature it lists is one that the load stored and counted.

The portal also takes a position among the ways a group can put a genome on
the web. Tripal runs a whole community website, with user accounts and pages
that curators edit [@tripal; @tripal-docs], and GenomeHubs brings the Ensembl
website with a BLAST server beside it [@genomehubs]; both do much more than
cgload's portal. They also bring a large software stack that has to be
installed and kept up to date, with security updates of its own; Drupal, on which Tripal is built, issued a highly critical fix for remote
code execution in 2018 [@drupal]. JBrowse 2 can read cgload's indexed files
with no code on the server at all [@jbrowse2] and searches the names it
indexes from them, but the functional annotation and the record of each load
stay in the database, out of its reach. cgload's portal lies between these:
two commands on a database that the load has already checked, one Python
process, and no request that can change the database or start a program.
Read-only access, bound parameters, escaping and a strict Content Security
Policy are each standard defences against injection [@owasp]. What the portal
leaves out matters as much: it has no form that writes and no shell that a
visitor's text could reach.

Refusing by default has a cost. Seven of the 30 RefSeq files were refused
until `--repair-envelopes` was given, over rows that carry a documented
exception or follow a documented convention. We chose to stop and name the
features, and to repair only on request, because a parent span that does not
contain its children is contradictory whatever its cause. A user who wants
every RefSeq file to load without a flag will find this strict.

## Limitations

cgload was corrected on the same 30 genomes it is measured on. Fourteen of its
recorded decisions (D-095 to D-108 in the repository) fixed problems that these
files exposed in the first run, the preflights and the second run. We did not
run a second panel of genomes that cgload had never seen, so its results on
other files, or on later RefSeq releases, may be worse. The planted faults and
the other loaders were not tuned in this way.

All 30 genomes are RefSeq. cgload's profiles for funannotate and AUGUSTUS
files, and its loading of functional annotation, are covered by its tests but were not part of this benchmark, and neither was
the reference tier of `verify`. MySQL was not tested at all. cgload refuses a
feature that names more than one parent, which GFF3 allows and RefSeq does not
use. Our results therefore show what cgload does with RefSeq files of the kind
tested, not with every valid GFF3 file.

We did not benchmark gffbase or GFFx. Neither was in the protocol: we found
them in a literature search after the runs were complete. gffbase in particular treats a multi-row feature as one feature, and belongs in
any repeat of this benchmark. Nor did we test BioPerl's
Bio::DB::SeqFeature::Store [@bioperl], an older loader that stores GFF3 in
MySQL, PostgreSQL or SQLite for genome browsers such as GBrowse.

The comparison is not symmetric in one respect. cgload's result was read from
the GFF3 file that `cgload export` writes, while gffutils' and Chado's were
read from their databases with SQL. Our numbers therefore test cgload's load
and export together and cannot tell the two apart.

`cgload verify` and the end-of-load check compare counts, not identities: they
cannot see a renamed feature, a wrong parent or a changed value (see
"Verification"). The round-trip comparison can, but it belongs to this
benchmark, not to cgload. It also left out columns 2, 6 and 8 of every row.
We compared those afterwards, for cgload 0.2.2 only, so for the other loaders
the counts of changed features are lower bounds. The reference reader shares no
code with any loader, but we wrote it, as we wrote cgload; a second,
independent reading of some of the files would make the reference stronger.

cgload holds the whole file in memory during a load, about 6 KiB per row. A
human-sized annotation needs 30 to 37 GiB, which rules out a laptop. Counting
features on PostgreSQL took 19.5 s for human, and gene lookups as commands
took about a quarter of a second, mostly to start Python; query waits were
measured on four genomes, not 30. A killed load leaves its staging files on
disk, to be deleted by hand. `verify` exits 0 when it skips the source count
because the input file has changed, so a script that reads only the exit
status misses the change; this is unchanged in version 0.3.0.

The portal and Q9 came after the benchmark. Q9's procedure was fixed before it
was run, but on the same 30 genomes that cgload had been corrected on, and
Figure 3 shows the repository's own test record. The portal's times were
measured with the client on the server, so they leave out the network and the
time a browser takes to draw igv.js. They are also the times of one process,
which answered about 40 requests a second however many clients asked at once;
its threads share Python's interpreter lock, so serving more needs more
processes behind a web server, which the portal allows because it keeps no
state. The security probe and the tests are our own, not an audit by an
outside party. We did not time another portal, since none serves these files
without converting them first, and installing Tripal or GenomeHubs to do so was
outside the protocol. The portal translates coding sequences itself, as the
file describes them. It names the exceptions that a file declares but does not
apply them, so for such a coding sequence the protein it shows is not the
curated one. It has no BLAST, no editing and no accounts.

The benchmark was written by cgload's authors, and its protocol was registered
as a tagged commit in our own public repository, not in an independent
registry. We tried to limit the bias that follows [@weber; @buchka] by fixing
the protocol first, by using a reference reader that shares no code with any
loader, and by reporting our own errors. Three errors in our scripts made other
loaders look worse until we found and fixed them
(see "Changes after the first run"). Chado stands here
for python-chado alone. The other loaders ran several at a time, so their load
times are indicative only, and we draw no conclusion from them.

# Conclusions

Loaders can change or lose annotation features without saying so: on 30
RefSeq genomes, each of the three other loaders we tested did. A database built
from a GFF3 file deserves the same check as any other copy of data: does it hold
what the file holds? cgload answers part of that question after every load. It
compares counts of what it stored with counts of the file, made by code that
does not share the loader's assumptions, and stops when they differ. On the 30
genomes it kept every feature in every field we compared, declared its 12
repairs and reported every planted fault that changed a count; version 0.2.2
also returns every row's source, score and phase as written. Counts cannot show
that a name, a link or a value is right; that takes a comparison feature by
feature, like the one this benchmark made. With `cgload index` and
`cgload serve` the same database becomes a portal that its users can search
and browse, that shows how each genome was loaded, and that cannot change what
it serves; one such portal served all 30 genomes from one database, with a
median of {{portal:gene_page_median}} ms for a gene's page. The
protocol, scripts and results are public, so the same questions can be put to
other loaders and other files.

# Availability and requirements

- Project name: cgload
- Project home page: https://github.com/computational-genomics-lab/cgload
- Operating system(s): Linux, macOS and Windows (the test suite passes on all
  three in continuous integration; the benchmark and the portal's measurement
  ran on Linux); the portal is used from any current web browser
- Programming language: Python 3.11 or later
- Other requirements: SQLAlchemy 2.0 or later, pyfaidx 0.7 or later, click 8.1
  or later; SQLite 3.8.3 or later (part of Python), or PostgreSQL (tested with
  version 16) with psycopg 3; for the portal, Flask 3.0 or later and waitress
  3.0 or later
- License: MIT
- Any restrictions to use by non-academics: none
- Version described: 0.3.0 (tag `v0.3.0`; source archive in Additional file 7)
- RRID:SCR_029282
- bio.tools ID: biotools:cgload
- Archived version: Software Heritage, swh:1:rel:26af2fd022ceebaff1a602bf15c036c107713aad (release v0.3.0),
  https://archive.softwareheritage.org/swh:1:rel:26af2fd022ceebaff1a602bf15c036c107713aad;origin=https://github.com/computational-genomics-lab/cgload.
  The loading results rest on release v0.2.2,
  swh:1:rel:e0db1e5b682047ae65d4e9e21b25749e886a8962,
  https://archive.softwareheritage.org/swh:1:rel:e0db1e5b682047ae65d4e9e21b25749e886a8962;origin=https://github.com/computational-genomics-lab/cgload;visit=swh:1:snp:ce34fb036eb8063458ec04aa5f64c3c8231d176a

cgload is installed with pip from its repository, at the version described
here, with the portal:

    pip install "cgload[portal] @ git+https://github.com/computational-genomics-lab/cgload@v0.3.0"

For PostgreSQL, add its driver as well: `cgload[portal,postgres]`.

`examples/minimal/` in the repository loads and verifies a small genome in
under a minute, as a check that the installation works; `cgload index` and
`cgload serve` then show it in a web browser.

The repository holds test data (`examples/minimal/`, `tests/fixtures/`) with
expected outputs, a test suite of 655 tests, the benchmark protocol and
scripts (`benchmark/`), and the portal's documentation (`docs/portal.md`).
`paper/figures.py` and `paper/screenshots.py` make every figure and table of
this paper from the results.

# List of abbreviations

BGZF: blocked GNU zip format; BLAST: Basic Local Alignment Search Tool; CDS:
coding sequence; CPU: central processing unit; CSI: coordinate-sorted index;
EC: Enzyme Commission; GFF3: Generic Feature Format version 3; GO: Gene
Ontology; GTF: Gene Transfer Format; HTTPS: Hypertext Transfer Protocol Secure;
ID: identifier; JSON: JavaScript Object Notation; KEGG: Kyoto Encyclopedia of
Genes and Genomes; MD5: message digest 5; NCBI: National Center for
Biotechnology Information; RRID: Research Resource Identifier; SAM: Sequence
Alignment/Map; SHA: Secure Hash Algorithm; SO: Sequence Ontology; SQL:
Structured Query Language; TOML: Tom's Obvious, Minimal Language; TxDb:
transcript database of Bioconductor.

# Declarations

## Ethics approval and consent to participate

Not applicable.

## Consent for publication

Not applicable.

## Availability of data and materials

The datasets supporting the conclusions of this article are included within
the article and its additional files (Additional files 1 to 4 and 6). They are
available in full in Additional file 7, a copy of the cgload repository as
submitted, and in the repository itself,
https://github.com/computational-genomics-lab/cgload: the results folder of
the benchmark, as its scripts left it, is `paper/data/benchmark-run-2.2/`, the
check of columns 2, 6 and 8 is `paper/data/columns-v0.2.2/`, and the portal's
measurement is `paper/data/portal-v0.3.0/`. The census scripts and their
results are in `survey/`. The reference readings are not
included, because `benchmark/truth.py` rebuilds them from the input files. All
input files are public RefSeq annotations, identified by accession, annotation
release and MD5 in Additional file 1. Release v0.2.2 of the repository, which
holds the software, the benchmark results and the census, is archived at
Software Heritage as swh:1:rel:e0db1e5b682047ae65d4e9e21b25749e886a8962
(https://archive.softwareheritage.org/swh:1:rel:e0db1e5b682047ae65d4e9e21b25749e886a8962;origin=https://github.com/computational-genomics-lab/cgload;visit=swh:1:snp:ce34fb036eb8063458ec04aa5f64c3c8231d176a);
the check of columns 2, 6 and 8 was added to the repository after that release.
Release v0.3.0, which adds the portal and its measurement, is archived as
swh:1:rel:26af2fd022ceebaff1a602bf15c036c107713aad
(https://archive.softwareheritage.org/swh:1:rel:26af2fd022ceebaff1a602bf15c036c107713aad;origin=https://github.com/computational-genomics-lab/cgload).

## Competing interests

The authors declare no competing financial or personal interests that could
have influenced the work reported in this manuscript.

## Funding

Funding for this work was provided by Department of Biotechnology (DBT) BIC and
NNP projects with ids BT/PR40137/BTIS/137/35/2022 and BT/PR40243/BTIS/137/75/2023
to ST. Fellowship for AU was provided by Council of Scientific and Industrial
Research, Government of India. The authors acknowledge the financial support for
the Article Processing Charges (APC) provided by the One Nation One Subscription
(ONOS) scheme, an initiative of the Ministry of Education, Government of India.

## Authors' contributions

AU: Conceptualization, Methodology, Software, Investigation, Validation, Formal
analysis, Data curation, Visualization, Writing – original draft.

AP: Conceptualization, Methodology, Writing – review & editing.

ST: Conceptualization, Methodology, Supervision, Project administration,
Funding acquisition, Writing – review & editing.

All authors read and approved the final manuscript.

## Acknowledgements

Not applicable.

# References

{{references}}

[@gff3]: Stein LD. Generic Feature Format Version 3 (GFF3), version 1.26. The Sequence Ontology Project. 2020. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md. Accessed 10 Nov 2026.
[@refseq]: O'Leary NA, Wright MW, Brister JR, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 2016;44:D733–45. https://doi.org/10.1093/nar/gkv1189.
[@gffutils]: Dale R. gffutils: GFF and GTF file manipulation and interconversion. Version 0.14. https://github.com/daler/gffutils. Accessed 10 Nov 2026.
[@chado]: Mungall CJ, Emmert DB, The FlyBase Consortium. A Chado case study: an ontology-based modular schema for representing genome-associated biological information. Bioinformatics. 2007;23:i337–46. https://doi.org/10.1093/bioinformatics/btm189.
[@tripal]: Sanderson LA, Ficklin SP, Cheng CH, et al. Tripal v1.1: a standards-based toolkit for construction of online genetic and genomic databases. Database. 2013;2013:bat075. https://doi.org/10.1093/database/bat075.
[@granges]: Lawrence M, Huber W, Pagès H, et al. Software for computing and annotating genomic ranges. PLoS Comput Biol. 2013;9:e1003118. https://doi.org/10.1371/journal.pcbi.1003118.
[@txdbmaker]: Pagès H, Carlson M, Aboyoun P, et al. txdbmaker: tools for making TxDb objects from genomic annotations. R package version 1.6.2. Bioconductor. https://doi.org/10.18129/B9.bioc.txdbmaker.
[@intermine]: Smith RN, Aleksic J, Butano D, et al. InterMine: a flexible data warehouse system for the integration and analysis of heterogeneous biological data. Bioinformatics. 2012;28:3163–5. https://doi.org/10.1093/bioinformatics/bts577.
[@intermine-tests]: InterMine. Acceptance tests. InterMine documentation. https://intermine.readthedocs.io/en/latest/database/data-integrity-checks/acceptance-tests/. Accessed 10 Nov 2026.
[@ensembl]: Dyer SC, Austine-Orimoloye O, Azov AG, et al. Ensembl 2025. Nucleic Acids Res. 2025;53(D1):D948–57. https://doi.org/10.1093/nar/gkae1071.
[@datacheck]: Ensembl. ensembl-datacheck: code for checking Ensembl databases. https://github.com/Ensembl/ensembl-datacheck. Accessed 10 Nov 2026.
[@gffbase]: Chao K-H. gffbase: a genomic-annotation engine for GFF3 and GTF. Version 0.3.0. 2026. https://github.com/Kuanhao-Chao/gffbase. Accessed 10 Nov 2026.
[@gffx]: Chen B, Wu D, Zhang G. GFFx: a Rust-based suite of utilities for ultra-fast genomic feature extraction. GigaScience. 2025;14:giaf124. https://doi.org/10.1093/gigascience/giaf124.
[@genometools]: Gremme G, Steinbiss S, Kurtz S. GenomeTools: a comprehensive software library for efficient processing of structured genome annotations. IEEE/ACM Trans Comput Biol Bioinform. 2013;10:645–56. https://doi.org/10.1109/TCBB.2013.68.
[@agat]: Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. https://doi.org/10.5281/zenodo.3552717.
[@weber]: Weber LM, Saelens W, Cannoodt R, et al. Essential guidelines for computational method benchmarking. Genome Biol. 2019;20:125. https://doi.org/10.1186/s13059-019-1738-8.
[@buchka]: Buchka S, Hapfelmeier A, Gardner PP, et al. On the optimistic performance evaluation of newly introduced bioinformatic methods. Genome Biol. 2021;22:152. https://doi.org/10.1186/s13059-021-02365-4.
[@tabix]: Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. Bioinformatics. 2011;27:718–9. https://doi.org/10.1093/bioinformatics/btq671.
[@htslib]: Bonfield JK, Marshall J, Danecek P, et al. HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience. 2021;10:giab007. https://doi.org/10.1093/gigascience/giab007.
[@eggnog]: Cantalapiedra CP, Hernández-Plaza A, Letunic I, et al. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Mol Biol Evol. 2021;38:5825–9. https://doi.org/10.1093/molbev/msab293.
[@interproscan]: Jones P, Binns D, Chang H-Y, et al. InterProScan 5: genome-scale protein function classification. Bioinformatics. 2014;30:1236–40. https://doi.org/10.1093/bioinformatics/btu031.
[@funannotate]: Palmer JM, Stajich J. Funannotate v1.8.1: eukaryotic genome annotation. Zenodo. 2020. https://doi.org/10.5281/zenodo.1134477.
[@augustus]: Stanke M, Morgenstern B. AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints. Nucleic Acids Res. 2005;33:W465–7. https://doi.org/10.1093/nar/gki458.
[@datasets]: O'Leary NA, Cox E, Holmes JB, et al. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. Sci Data. 2024;11:732. https://doi.org/10.1038/s41597-024-03571-y.
[@pychado]: Galaxy Genome Annotation. python-chado: a Python library for interacting with Chado. Version 2.3.9. https://github.com/galaxy-genome-annotation/python-chado. Accessed 10 Nov 2026.
[@so]: Eilbeck K, Lewis SE, Mungall CJ, et al. The Sequence Ontology: a tool for the unification of genome annotations. Genome Biol. 2005;6:R44. https://doi.org/10.1186/gb-2005-6-5-r44.
[@biopython]: Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics. 2009;25:1422–3. https://doi.org/10.1093/bioinformatics/btp163.
[@claude]: Anthropic. Claude [large language model]. https://www.anthropic.com/claude. Used 20 Aug to 25 Nov 2026.
[@pyfaidx]: Shirley MD, Ma Z, Pedersen BS, Wheelan SJ. Efficient "pythonic" access to FASTA files using pyfaidx. PeerJ PrePrints. 2015;3:e970v1. https://doi.org/10.7287/peerj.preprints.970v1.
[@rcore]: R Core Team. R: a language and environment for statistical computing. Version 4.5.3. Vienna: R Foundation for Statistical Computing; 2026. https://www.R-project.org/.
[@bioperl]: Stajich JE, Block D, Boulez K, Brenner SE, Chervitz SA, Dagdigian C, et al. The Bioperl toolkit: Perl modules for the life sciences. Genome Res. 2002;12:1611–8. https://doi.org/10.1101/gr.361602.
[@chatgpt]: OpenAI. ChatGPT [large language model]. https://chatgpt.com. Accessed 19 Nov 2026.
[@tripal-docs]: Tripal. Tripal 4 documentation: requirements. https://tripaldoc.readthedocs.io/en/latest/install/requirements.html. Accessed 25 Nov 2026.
[@genomehubs]: Challis RJ, Kumar S, Stevens L, Blaxter M. GenomeHubs: simple containerized setup of a custom Ensembl database and web server for any species. Database. 2017;2017:bax039. https://doi.org/10.1093/database/bax039.
[@jbrowse2]: Diesh C, Stevens GJ, Xie P, et al. JBrowse 2: a modular genome browser with views of synteny and structural variation. Genome Biol. 2023;24:74. https://doi.org/10.1186/s13059-023-02914-z.
[@igvjs]: Robinson JT, Thorvaldsdottir H, Turner D, Mesirov JP. igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). Bioinformatics. 2023;39:btac830. https://doi.org/10.1093/bioinformatics/btac830.
[@owasp]: OWASP Foundation. OWASP Top 10:2025. https://owasp.org/Top10/2025/. Accessed 25 Nov 2026.
[@tao]: Tao T, Chen Y, Liu B, Jin X, Yan M, Ji S. Security analysis of bioinformatics web application. In: Yang CN, Peng SL, Jain LC, editors. Security with Intelligent Computing and Big-data Services. SICBS 2018. Advances in Intelligent Systems and Computing, vol 895. Cham: Springer; 2020. https://doi.org/10.1007/978-3-030-16946-6_30.
[@kern]: Kern F, Fehlmann T, Keller A. On the lifetime of bioinformatics web services. Nucleic Acids Res. 2020;48:12523–33. https://doi.org/10.1093/nar/gkaa1125.
[@flask]: Pallets. Flask. Version 3.1. https://flask.palletsprojects.com. Accessed 25 Nov 2026.
[@waitress]: Pylons Project. Waitress. Version 3.0. https://docs.pylonsproject.org/projects/waitress/. Accessed 25 Nov 2026.
[@gencode]: Elzanowski A, Ostell J. The genetic codes. National Center for Biotechnology Information. Updated 23 Sep 2024. https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. Accessed 25 Nov 2026.
[@csp]: West M, Sartori A, editors. Content Security Policy Level 3. W3C Working Draft. https://www.w3.org/TR/CSP3/. Accessed 25 Nov 2026.
[@bandit]: PyCQA. Bandit: a tool to find common security issues in Python code. Version 1.8.6. https://github.com/PyCQA/bandit. Accessed 25 Nov 2026.
[@drupal]: Drupal Security Team. Drupal core – highly critical – remote code execution – SA-CORE-2018-002. 28 Mar 2018. https://www.drupal.org/sa-core-2018-002. Accessed 25 Nov 2026.

# Figure legends

**Figure 1. The steps of a cgload load and of cgload verify.** (a) `cgload
load` first reads and checks every annotation file and writes nothing (top,
blue). It then writes the genome in one transaction, in the nine steps of the
text (bottom): grey steps write and the blue step checks. Left of each step,
in orange, is what stops the load at that step; right of it is what the step
writes. The independent count is made from each annotation file by separate
code and gives the numbers that step 7 compares with the database. Any stop in
the transaction rolls the database back and deletes the staging folder; a
killed process rolls the database back but leaves the staging folder on disk.
(b) The five tiers of `cgload verify`, which can be run at any later time. They
compare counts, so they cannot tell which feature is which.

**Figure 2. The database schema of cgload.** The 11 tables, with selected
columns; `src/cgload/db/schema.py` defines every column, constraint and index.
The italic line under each name says what one row of the table stands for.
Grey tables describe the genome and its input files, dark blue tables hold the
annotation, and light blue tables hold functional annotation. The two outlined
tables are written by `cgload init`. A line joins each foreign key to the key
it names, marked with a dot; `feature.parent_id` names a row of the same table.
The figure is drawn from the schema itself by `paper/figures.py`.

**Figure 3. Three pages of the cgload portal.** The genome is the repository's
test record of the Dscam1 region of *Drosophila melanogaster* chromosome 2R
(NCBI NT_033778.4, bases 7,312,924 to 7,386,899, numbered from the start of
the region), loaded with the record's own sequence. (a) A gene's page:
its three transcripts drawn to scale, with exons as boxes, coding sequence as
taller boxes and introns as lines, and igv.js reading the assembly and the
annotation that cgload wrote, served by the portal itself. The page goes on
with the gene's attributes, sequences and family. (b) A search by the words of
a product, with how each feature matched. (c) The genome's page records how it
was loaded: the checks that `cgload load` ran and their outcome, any repairs,
and the input files with their checksums. `paper/screenshots.py` photographs
the pages in Chromium and crops them; it removes the empty rows below the
browser's tracks and a one-line key under the gene model.

**Figure 4. Features changed or lost without a message, per genome and
loader.** Each cell gives the percentage of the genome's features (for
txdbmaker, transcripts) that the loader changed or lost without a message.
Dark grey cells mark genomes a loader did not load. cgload's results were
identical on SQLite and PostgreSQL. † cgload loaded these genomes with
`--repair-envelopes` and declared every repair.

**Figure 5. The cost of loading the 30 genomes with cgload.** Against the
number of GFF3 rows: (a) load time, (b) time of `cgload verify`, (c) peak
memory of the load and (d) size of the database. Each point is one genome on
one database system. Both axes are logarithmic. The values for human, the
largest genome, are labelled.

**Figure 6. Time to answer four everyday tasks on four genomes.** Median time
of repeated asks, on genomes from 0.1 to 4.9 million rows. Top: one new process
per answer, as a user typing a command waits. Bottom: one process asking T1 to
T3 in turn, as a script waits. T4 is shown for the largest sequence of each
genome. tabix answers T2 and T4 from cgload's compressed annotation, and T3 by
a full scan. Medians at or under 1 ms are drawn at 1 ms. The page cache was not dropped.

# Additional files

**Additional file 1 (CSV).** Table S1. The 30 genomes, with the annotation
release and MD5 of each input file, and cgload's load and verify times, peak
memory and sizes on each database system.

**Additional file 2 (CSV).** Table S2. Every planted fault, its expected
outcome and what happened in the three genomes.

**Additional file 3 (CSV).** Table S3. Query timing: median of the first asks
and of the repeated asks, per genome, tool, task and view.

**Additional file 4 (CSV).** Table S4. The sequences written for T4, and the
time each tool took to write them.

**Additional file 5 (Markdown).** Use of large language models
(`paper/llm_use.md` in the repository): which model, what it did, how the
authors checked it, and where its output is kept.

**Additional file 6 (CSV).** Table S5. The portal on each of the 30 genomes:
the load, the search index and its check, the median time of each kind of
request, how often a search by a gene's name found it, and the translation of
its coding sequences.

**Additional file 7 (ZIP).** Source code of cgload, version 0.3.0: the
repository at the time of submission, with its tests, the benchmark protocol
and scripts, and the results of the benchmark and of the portal's
measurement.
