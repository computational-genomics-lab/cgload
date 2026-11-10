---
title: "cgload: loading genome annotations into a database and checking that every feature arrived"
subtitle: "Problem Solving Protocol, for Briefings in Bioinformatics. Draft of 10 November 2026."
---

<!--
Built by paper/build.py: each "include" line below is replaced by the table
that paper/figures.py made from paper/data/benchmark-run-2.2. Do not type a
result number into a table by hand; change the data or the script.
-->

**Aditya Upadhyay**^1,2^, **Arijit Panda**^3^ and **Sucheta Tripathy**^1,2,\*^

^1^ Structural Biology and Bioinformatics Division, CSIR-Indian Institute of Chemical Biology, Jadavpur, Kolkata 700032, India.
^2^ Academy of Scientific and Innovative Research (AcSIR), Ghaziabad 201002, India.
^3^ Department of Quantitative Health Sciences, Mayo Clinic, Rochester, MN 55905, USA.

\* Corresponding author: Sucheta Tripathy, [institutional e-mail].

Running title: Loading genome annotations without losing features

# Abstract

A loader that reads a GFF3 annotation into a database can drop or alter
features, and the database it leaves still looks complete. We present cgload,
which loads an assembly and its annotation into SQLite or PostgreSQL in one
transaction, compares the stored features and parent links with a count of the
input made by separate code, and rolls back a load that falls short. We tested
cgload and three other loaders on 30 RefSeq genomes (28.2 million GFF3 rows,
18.0 million features) under a protocol committed before the first run,
comparing each result feature by feature with a reference reading that shares
no code with any loader. cgload stored 18,017,606 features unchanged on both
database systems, declared the 12 it repaired, and lost none. gffutils refused
all 30 files with its default settings, and with create_unique it split
1,072,125 multi-row features into separate records without a message.
txdbmaker failed on six plant genomes. python-chado failed or timed out on nine
genomes and, in the 21 loads that finished, lost 5.0% of the features, four in
five of them without a message. cgload noticed all 57 planted faults that its
checks cover and, as documented, none of the 24 that changed values but kept
every count. Every killed cgload load left the database as it was. Loading the
human annotation took 641 s and 29.5 GB of memory on SQLite. cgload is free
software under the MIT license.

**Keywords:** genome annotation; GFF3; database; data integrity; benchmarking; RefSeq

# Key Points

- A database built from a GFF3 file can miss features, or hold them in altered
  form, without any message to the user.
- On 30 RefSeq genomes, gffutils (create_unique) split 1.07 million multi-row
  features into separate records, and python-chado lost 5% of the features in
  the loads that finished, mostly without a message.
- cgload compares every load with an independent count of the input file and
  rolls back a load that falls short. It stored all 18.0 million features and
  declared the 12 it repaired.
- cgload checks that every feature and stated parent link is present. It does
  not check their values, and it missed every planted fault that changed a
  feature's values while keeping the counts.
- The protocol, scripts, results and figures are public, and every change made
  after the first run is recorded with its date and reason.

# Introduction

Genome projects publish their annotations as GFF3 files [1], such as those of
NCBI RefSeq [2]. To ask questions of an annotation, many groups first load it
into a database. gffutils builds an SQLite database from a GFF3 file [3]. Chado
is a relational schema for genome data, used by many organism databases [4]
and filled by loaders such as the one in Tripal [5]. Bioconductor's txdbmaker
builds a TxDb of transcripts for analysis in R [6, 7], and data warehouses such
as InterMine combine annotations from many sources [8].

A loader can drop or change features while it reads a file. GFF3 lets one
feature be written on several rows that share an ID, which is how RefSeq writes
coding sequences, so a reader that expects one row per feature, or one ID per
row, gets these features wrong. Some features cross the origin of a circular
sequence, some genes are trans-spliced from distant parts of an organelle
genome, and some strands are unknown. In a survey of the 30 RefSeq genomes used
here [9], rows that break a GFF3 rule were rare (9 of 28.2 million, and 4 more
that no convention explains). Rows that follow a documented convention or carry
a documented exception were more common (98), and a loader has to take them as
written. A loader can finish without a message and still leave such features
out, or store them in another form. The database then answers queries, and
nothing downstream can see what is missing.

Large resources test their databases before each release. InterMine's
acceptance tests are SQL queries with expected answers, such as a minimum
number of GO annotations [10]. Ensembl's datachecks test row counts and
metadata in each database before it is handed over [11, 12]. Curators write
these tests for a resource they know, and the tests compare the built database
with their expectations. As documented, neither compares a database with an
independent reading of the files it was built from. Recent tools aim at speed:
gffbase stores GFF3 and GTF files in DuckDB and validates a sample of records
[13], and GFFx indexes GFF files for fast extraction and assumes well-formed
input [14]. Toolkits such as GenomeTools [15] and AGAT [16] check or repair a
file before it is loaded, and do not see what a loader then stores.

We wrote cgload to make that check part of every load. It reads the input file
twice, once to load it and once, with separate code, to count it, and it keeps
the load only if the database holds what the count says. We then asked how
cgload and other loaders behave on real files. Benchmarks written by a tool's
own authors tend to favour that tool [17, 18]. To limit this, we fixed the
questions, the genomes and the scripts in a tagged protocol before the first
run, compared every loader with the same reference reading, and report every
change made after the first run. This paper describes cgload and the benchmark,
and reports its answers to seven questions:

- Q1. Does cgload load and verify all 30 genomes, on SQLite and on PostgreSQL,
  and at what cost in time, memory and disk?
- Q2. Does what cgload stores match the input, feature by feature?
- Q3. What do other loaders keep, change or lose from the same files, and do
  they say so?
- Q4. Which damage to a loaded genome, or inside the parser, does cgload notice?
- Q5. What does an interrupted load leave behind?
- Q6. How long does a user wait for everyday answers?
- Q7. Do cgload's structural checks agree with the survey's checker on the same
  files?

<!-- AUTHORS: the protocol also lists Q8, a case study on one of our own
genomes with eggNOG-mapper and InterProScan output. It is not run yet. Either
run it (it adds a section 3.7 and the JBrowse figure) or record in PROTOCOL.md
section 12 that it was dropped, and why. -->

# Materials and methods

## cgload

cgload is a Python package (3.11 or later) with three runtime dependencies:
SQLAlchemy, pyfaidx and click. It stores a genome in a relational database
through SQLAlchemy, on SQLite, PostgreSQL or MySQL. The schema has 11 tables. A
feature written on several rows keeps one stored row per segment, each with
its own attributes; features are identified by their ID, and parent links are
kept in a table of their own. The assembly FASTA is compressed with bgzip and
indexed [19, 20] in a data folder beside the database, and the database records
its checksum. Functional annotation from eggNOG-mapper [21] and InterProScan
[22] becomes one row per GO term, KEGG ortholog, Pfam domain or EC number.

Annotation pipelines write GFF3 differently. cgload reads RefSeq GFF3, the GFF3
of funannotate [23] and AUGUSTUS [24], and GenBank flat files. Each dialect is
a declared profile that says which feature types it provides and how the file
identifies itself, and a file that matches no profile is refused with the names
of the profiles that were tried. Every dialect is normalised to one
representation before anything reaches the database.

Before it writes anything, cgload checks the structure of the file. A feature
must lie inside the span of the parent that claims it, a start must not exceed
its end, and a row must not run past the end of its sequence unless the
sequence is circular. By default a file that breaks these rules is refused, and
the message names each offending feature and the option that repairs it. With
`--repair-envelopes`, cgload recomputes each parent's span from its children,
leaves first, and writes every change to the load's manifest. A disagreement
that exists only across the origin of a circular sequence is stored as written.

The load is one database transaction. The data folder is prepared in a staging
directory and renamed into place just before the commit. At the end of the
load, cgload counts what it stored and compares the counts with an independent
count of the file. If any count is short, it rolls the transaction back. A load
therefore stores the whole file or leaves the database as it was.

The command `cgload verify` repeats these checks at any later time, in tiers
(Figure 1). The source-count tier counts the file again with code that imports
nothing from cgload: features of each type, rows of each feature written on
several rows, and parent statements. These counts must equal the database's.
The parent-edge tier requires every `Parent=` statement in the file to be a
stored link. The stored-assembly tier compares the FASTA and its index with the
checksums taken at load. Two optional tiers compare the file with other
descriptions of the same assembly: the GenBank file (cross-parser tier) and
the totals NCBI publishes for the annotation release (reference tier). A tier
that cannot run says so and gives the reason. `verify` exits with status 2 when
a tier fails and 0 otherwise.

These checks establish completeness: every feature, every row and every stated
parent link is present, and nothing was invented. They do not establish that a
value is right. A wrong strand or a shifted coordinate keeps every count the
same, and the documentation says so.

![**Figure 1.** How cgload loads and verifies an annotation. Grey boxes load and blue boxes check. A separate reader counts the input file with code that imports nothing from cgload. At the end of the load its counts decide whether the transaction is committed or rolled back, and `cgload verify` compares them with the database again at any later time.](figures/fig1_design.png)

## Benchmark design

The protocol (`benchmark/PROTOCOL.md` in the repository) was committed and
tagged `benchmark-protocol-v1` before the first run, with every script it
names. It states the questions, the genomes, the software, the hardware, the
definitions and, for every planted fault, the outcome expected from what cgload
claims to check. A change made after the first run is listed in the protocol
with its date and reason, and every measurement it affects was repeated
(section "Changes after the first run").

## Genomes and input files

The panel is the 30 RefSeq genomes of the survey [9]: 8 vertebrates, 7
insects, 7 plants and 8 fungi (Supplementary Table S1). The list was fixed in a
tagged commit before the survey ran. Together the annotations have 28,211,625
GFF3 rows, 18,017,618 features and 16,566,107 stated parent links. The smallest
is *Candida albicans* (25,889 rows) and the largest is human (4,932,571 rows,
3,203,237 features). Each annotation is the exact file the survey read,
checked by its MD5. One survey download, medaka, had been cut short; the survey
was repeated on the complete file and its findings did not change [9]. The
assembly FASTA files, and the GenBank files of the eight fungi, came from NCBI
Datasets [25].

## Reference reading and comparison

A reference reader (`truth.py`) reads each GFF3 file with the Python standard
library only and shares no code with cgload or with any other loader. It
follows the GFF3 specification [1]: rows that share an ID form one feature, a
row without an ID is a feature of its own, and `Parent` names the parents. For
each feature it records the sequence, the type, the strand and span of every
row, the parents, and a digest of all other attributes (keys and values sorted,
escapes decoded). White space at either end of a tag=value pair is not part of
it.

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
for cgload and gffutils. For Chado the attributes are left out, because
python-chado spreads them over several tables. txdbmaker holds transcripts, not
features, so it is compared per transcript: sequence, strand, span, gene, exon
spans and CDS spans. Transcripts are paired by the name txdbmaker gives them,
or failing that by span, and features of types txdbmaker does not model are
counted apart, not as lost. A feature whose rows carry different attributes is
compared row by row. Each feature is also placed in the hard cases the
protocol defines, such as a feature written on several rows, a trans-spliced
feature, a feature on an organelle or one that crosses the origin. Twelve of
these cases occur in the 30 files (Table 2).

## Other loaders

gffutils 0.14 was run twice: with its default merge strategy, which is what a
user who changes nothing gets, and with `create_unique`, the documented way to
load a file whose rows repeat an ID [3]. txdbmaker's `makeTxDbFromGFF` was told
which sequences are circular [7]. For Chado we used python-chado 2.3.9 [26],
the `load_gff` command used by Galaxy Genome Annotation, which reads GFF3 with
Biopython 1.81 [28], into a new Chado 1.31 database with the Sequence Ontology
[27]. Chado has other GFF3 loaders, such as
the GMOD bulk loader and Tripal's importer, and our results stand for
python-chado only. Each load had a time limit of 12 hours, or 24 hours for
Chado. A load that was refused, crashed, or ran out of memory or time is
reported as such, and whatever a failed Chado load left in its database was
still compared.

## Faults, interrupted loads, timing and agreement

For Q4 we planted 27 faults, each in three genomes (*Fusarium graminearum*,
*Arabidopsis thaliana* and human). Faults F01 to F19 damage a copy of the
SQLite database or its data folder after the load: deleting a feature or a
row, removing a parent link, editing a base of the stored FASTA, or editing the
input file. Faults P01 to P08 damage the load itself, through a wrapper around
cgload's tokenizer that drops, adds or alters rows. The protocol states each
fault's expected outcome. A fault that changes a count should be detected, and
one that changes only a value, such as a strand or a coordinate, should not
(Supplementary Table S2).

For Q5 we killed loads with SIGKILL, sent to the whole process group at a time
drawn uniformly between 5% and 95% of the duration of an uninterrupted load,
from a seeded generator. Each loader had 10 trials per genome on *F.
graminearum* and *A. thaliana*, and Chado 5 trials on *F. graminearum*. For
cgload we killed new loads and `--force` reloads over an earlier load. After
each kill we recorded what the database held, counted in the loader's own
tables and compared with an uninterrupted load.

For Q6 we timed four tasks on four genomes spanning 0.1 to 4.9 million rows
(*F. graminearum*, *A. thaliana*, maize and human): T1, one gene and every
feature below it; T2, every feature overlapping a 1 Mb window; T3, feature
counts by type; and T4, every feature on one sequence, written as GFF3. Each
task was asked of 20 genes or windows, or 3 sequences, 5 times each, in two
ways: one new process per answer, as a user typing a command waits, and one
process asking T1 to T3 in turn, as a script waits. The tools were cgload on
SQLite and on PostgreSQL, gffutils through its Python interface, and tabix on
cgload's compressed annotation for T2 and T4, with T3 as a full scan. We report
the median of the repeated asks; the first asks are in Supplementary Table S3.
The server did not let us drop the page cache, so no ask started from a cold
disk.

For Q7 we took every row the survey flagged and the survey's class for it, and
ran cgload's checks on the same files with its own code, without a database, as
`cgload load --repair-envelopes` runs them.

## Hardware and software

All runs used one lab server with Linux 6.5, 80 CPUs, 125 GB of memory and 8.7
TB of disk, PostgreSQL 16 with default settings, and Python 3.12 (3.11 for
python-chado). Loads that were timed for Q1 and Q6 ran alone on the machine.
The other loaders ran several at a time (8 jobs for gffutils, 6 for Chado and 3
for txdbmaker), so their times are indicative only. Wall time, CPU time and peak
memory were measured for every command. Memory and sizes are given in GB of
2^30^ bytes. The exact versions installed are kept with the results.

## Changes after the first run

The first run, on 21 October 2026, is void: two copies of the run script were
started in one folder and wrote over each other's results. It still showed
problems, and we fixed each with tests before the second run. cgload refused
six RefSeq files over valid content, failed on one attribute value that was too
long for a PostgreSQL index, and changed some features on the way back out
through `cgload export`. Its cross-parser tier failed on five fungi over names
that NCBI writes differently in GFF3 and GenBank. A killed load could leave a
directory that blocked the next load. Two preflights on the 14 genomes involved
then found more of the same kind, and also showed that the survey had read an
incomplete medaka file.

The second run (23 and 24 October) finished every step. It found two errors in
our scripts, in the Chado database setup and in txdbmaker's pairing of
transcripts, and one cgload defect: a gene lookup on PostgreSQL took 12 s for
human. We fixed these, and repeated Chado, txdbmaker and the timing on all 30
genomes. The human Chado load then reached its 24-hour limit with nothing
stored, and we finished its comparison without loading again. Last, we found
that the comparison had charged gffutils with 133 attribute changes in
*Schizosaccharomyces pombe* that were ours: gffutils keeps a space at the end
of a value, and the reference reader dropped it. We fixed the comparison and
measured gffutils again on all 30 genomes. The protocol lists each change with
its date, its commit and the tests added.

cgload was therefore corrected on the same 30 genomes it is measured on, and we
return to this under Limitations. The results reported here come from three
tagged commits, and the protocol says which result comes from which.

# Results

## cgload loads and verifies all 30 genomes (Q1, Q2)

cgload loaded all 30 genomes on both SQLite and PostgreSQL (Supplementary Table
S1). Twenty-three loaded as they are. Seven were refused at first, each over
features lying outside the parent that claims them, or a start after its end,
and loaded with `--repair-envelopes`, which made 12 repairs in all and listed
each one. After every one of the 60 loads, `verify` passed the source-count,
parent-edge and stored-assembly tiers, and for the eight fungi the cross-parser
tier as well.

The round trip kept 18,017,606 of the 18,017,618 features. The other 12 were
the repaired ones, changed in their spans as the manifest declared. No feature
was lost, none changed without a declaration, and all 16,566,107 parent links
came back. SQLite and PostgreSQL gave the same result for every genome.

Load time, memory and size grew in proportion to the number of rows (Figure 2).
On SQLite, cgload loaded a median 7,961 rows per second; the quickest load, *S.
cerevisiae*, took 4.8 s and human took 641 s. PostgreSQL was about three times slower
(median 2,700 rows per second; human 1,978 s). `verify` took 0.8 to 29 s on
SQLite and 1.5 to 48 s on PostgreSQL. The load held about 6 KB of memory per row
(median 5.7 KB on SQLite and 6.9 KB on PostgreSQL), so human needed 29.5 GB and
36.4 GB at its peak. Its database took 3.0 GB on SQLite and 5.9 GB on
PostgreSQL, and its data folder 3.2 GB.

![**Figure 2.** The cost of loading the 30 genomes with cgload, against the number of GFF3 rows: (a) load time, (b) time of `cgload verify`, (c) peak memory of the load and (d) size of the database. Each point is one genome on one database system. Both axes are logarithmic. The values for human, the largest genome, are labelled.](figures/fig2_scale.png)

## What other loaders keep, change and lose (Q3)

Table 1 sums up the results and Figure 3 shows them per genome.

<!-- include tables/table1_loaders.md -->

**Table 1.** What each loader kept, changed and lost from the 30 annotations.
Counts are features, except for txdbmaker, which holds transcripts. A change is
declared when the loader recorded it, announced when its output names the
feature or its type, and silent otherwise. cgload's seven repaired genomes
needed `--repair-envelopes`. Chado's 9 failed loads stored nothing, so every
feature of those genomes counts as lost.

gffutils with its default settings stopped on every file at the first ID
written on two rows (`ValueError: Duplicate ID`), so a user who changes nothing
gets no database. With `create_unique` it loaded all 30 and lost no feature,
but stored each of the 1,072,125 features written on several rows (6.0% of all
features) as separate records, the second and later ones under new IDs. A query
for such a feature by its ID returns one part. The share ranged from 1.2% in
*Saccharomyces cerevisiae* to 13.1% in *Aspergillus nidulans* and *Neurospora
crassa*. gffutils printed nothing about it. It changed no attribute and kept
every parent link.

txdbmaker stopped on six of the seven plant genomes. On three it rejected a row
whose end precedes its start, the inverted *nad1* transcripts the survey
reported [9], and on the other three it stopped while building the TxDb. On the
24 genomes it loaded, it kept 1,169,686 of 1,169,826 transcripts. It changed 104
and dropped 36, and its messages named all but 3 of the changes. Another 647,583
features have types that a TxDb does not model.
<!-- AUTHORS: the three makeTxDbFromGRanges errors (rice, tomato, grapevine)
are in compare/txdbmaker/<accession>/steps/load.err on the server; add the
reason in a few words. -->

python-chado stopped with an error about a feature's location on all seven
plant genomes and on *Drosophila melanogaster*, each time at a trans-spliced
gene: organelle genes in the plants and CG32491 in the fly. The human load
reached its 24-hour limit. None of these nine loads stored anything. In the 21
loads that finished, python-chado lost 540,255 of 10,773,540 features (5.0%),
of three kinds. It dropped all 109,181 features of type lnc_RNA, a type the
current Sequence Ontology names lncRNA [27], and its log named the type. It
dropped their 430,836 exons without a message. It also dropped 238 long
terminal repeats that the *S. pombe* file writes in capitals, and named that
type too. So four in five of its losses were silent, from none in 9 genomes
to 14.6% of the features in *S. pombe* (Figure 3). These loads kept 9,627,024
of 10,167,041 parent links, and python-chado added 670,102 polypeptide features
that the files do not contain.

![**Figure 3.** Features changed or lost without a word, as a percentage of each genome's features (txdbmaker: transcripts). Grey cells mark loads that failed or ran out of time. cgload's results were identical on SQLite and PostgreSQL. † cgload loaded these genomes with `--repair-envelopes` and declared every repair.](figures/fig3_loaders.png)

The hard cases show where the alterations come from (Table 2). gffutils split
every feature written on several rows, and these include the features of most
other hard cases. txdbmaker's changes and losses fell on rows that overlap, on
unknown strands and on trans-spliced and organelle transcripts, and its
messages named them. In the loads that finished, python-chado kept every
hard-case feature; the genomes with trans-spliced genes or unknown strands were
the ones on which it stopped. cgload kept every hard-case feature or declared
its repair.

<!-- include tables/table2_cases.md -->

**Table 2.** Hard cases, as defined in the protocol: for each, the number of
features in the 30 files, what each loader altered (changed or lost) and how
many of those alterations were silent. A feature can be in several cases, and
nearly every feature has an attribute with several values. Chado's counts cover
its 21 finished loads and txdbmaker's the transcripts of the 24 genomes it
loaded, each with its own total; a dash means the case does not occur there.

## Faults (Q4)

All 81 fault trials (27 faults in 3 genomes) ended as the protocol predicted
(Table 3, Supplementary Table S2). `verify` detected every fault that deleted,
duplicated or retyped a feature, removed a parent link, or damaged or moved
the stored assembly. Every fault planted in the parser that dropped, added or
retyped a row stopped the load: one by its structural check (a child whose
parent was dropped), the others by the end-of-load counts, which rolled the
load back. As stated in advance, cgload did not notice the eight faults that
change only a value: a flipped strand, a shifted start, a renamed feature, a
parent pointed elsewhere, a deleted attribute, and in the parser a flipped
strand, a shifted start and two swapped IDs.

One fault needs a caveat. When the input file was edited after the load
(F16), `verify` named the file, said it had changed since the load, and skipped
the source count. It still exited with status 0, because a skipped tier is not
a failure. A script that reads only the exit status would miss this edit.

<!-- include tables/table3_faults.md -->

**Table 3.** Planted faults, grouped by what they do. Each fault was planted in
*F. graminearum*, *A. thaliana* and human. "Caught by" names the check that
reported it.

## Interrupted loads (Q5)

Every killed cgload load left the database as it was (Table 4). After the 20
killed new loads the database was empty, and the same load then succeeded each
time. After the 20 killed `--force` reloads, the earlier load was intact and
passed `verify`. In 28 of the 40 trials, cgload left files on disk outside any
stored assembly, from its staging folder or from the folder a `--force` reload
sets aside (up to 124 MB). After a killed new load, the next load did not
remove them. The documentation lists this.

python-chado commits only at the end, so its five killed loads left an empty
database. gffutils left a database file after all 20 kills, and gffutils could
open none of them. Seventeen were empty. Three held every record of a complete
load, because the kill came after the last insert but before gffutils wrote the
metadata it needs to open the file.

<!-- include tables/table4_kills.md -->

**Table 4.** What a killed load left behind. Each kill came at a random time
between 5% and 95% of an uninterrupted load.

## Query waits (Q6)

cgload answered everyday questions more slowly than gffutils (Figure 4). One
gene with its parts took 0.26 to 0.30 s as a command on SQLite and 0.37 to 0.41
s on PostgreSQL, against 0.13 to 0.24 s for gffutils. Most of that time is
the start of a new Python process: asked from one running process, the same
lookup took 4 to 5 ms on SQLite, 8 to 25 ms on PostgreSQL and 1 ms or less in
gffutils. For a 1 Mb window, tabix answered in 8 to 36 ms, at least six times
faster than any of the databases. Counting features by type took cgload
longest: 3.8 s for human on SQLite and 19.5 s on PostgreSQL, against 0.3 s for
gffutils. The two answers differ, since gffutils counts records and cgload
counts features. Writing out every feature on one sequence took up to 9.4 s on
SQLite and 30.6 s on PostgreSQL (*A. thaliana*, median of three sequences),
against 6.6 s for gffutils.

![**Figure 4.** Median time of repeated asks for four everyday tasks, on four genomes from 0.1 to 4.9 million rows. Top: one new process per answer, as a user typing a command waits. Bottom: one process asking T1 to T3 in turn, as a script waits. tabix answers T2 and T4 from cgload's compressed annotation, and T3 by a full scan. Medians at or under 1 ms, the resolution of the clock, are drawn at 1 ms. The page cache was not dropped.](figures/fig4_queries.png)

## Agreement with the survey (Q7)

The survey flagged 111 rows in the 30 files: 102 listed one by one, and 9
crossing the origin of a circular sequence, which it only counted [9]. cgload
accepts rows of the second kind, which GFF3 allows. Of the 102 listed rows
(Table 5), both checkers accepted 23 rows with an unknown strand, and both
objected to the 13 the survey classed as errors or unexplained, though cgload
loaded 6 of these with a warning. For 66 rows that the survey classed as
documented conventions or exceptions, nearly all in trans-spliced organelle
genes, cgload refused the file unless `--repair-envelopes` was given. With it,
cgload repaired the parents of 59 of these rows, repaired 5 rows itself and
stored 2 as written. No row that the survey classed as an error was accepted by
cgload without comment. cgload objected to 8 features that the survey did not
list. Each is a gene or pseudogene whose span leaves out a child that the
survey did list, so the two checkers found the same 11 parent and child pairs
from opposite ends.

<!-- include tables/table5_agreement.md -->

**Table 5.** cgload's checks and the survey's checker on the same 30 files.
The survey's classes are its own [9].

# Discussion

Each of the other three loaders stored some features in another form, or not
at all, without a message: gffutils 1,072,125 features, python-chado 430,836
in its finished loads and txdbmaker 3 transcripts. The causes differ. RefSeq
writes coding sequences on several rows that share an ID, as GFF3 allows.
gffutils refuses such files by default and, with the setting its documentation
suggests, stores each row as a feature of its own. A user who then asks for a
CDS by its ID gets one segment of it, and nothing in the database says that
the rest exists under other names. python-chado stopped at trans-spliced genes,
and where it finished, it dropped a whole class of transcripts over the
spelling of a type name and their exons along with them. txdbmaker was the
most careful of the three: it stopped on the files it could not model and
named nearly every transcript it changed or dropped.

The other tools do some things better. gffutils loaded the human annotation in
7 minutes with 140 MB of memory, where cgload needed 11 minutes and 29.5 GB,
and it answered lookups faster. Each tool makes a choice about rows that share
an ID, and the trouble is that the choice is silent. A benchmark of this kind
can only be fair to each tool's documented use, and we ran each as its
documentation describes; other settings, or other Chado loaders, may do better.

cgload adds a check that does not trust its own loader. The source count
is made by code that shares nothing with cgload's parser, so an error in the
parser cannot hide itself in the check. The faults show both sides of this.
Every fault that changed a count was caught, including those planted inside the
parser, which a check built from the same parser could not see. Every fault
that changed only a value went unnoticed, as the protocol predicted. A user
who needs values checked has to compare them with the file, as our round trip
does, and that comparison is not part of `cgload verify`.

Refusing by default has a cost. Seven of the 30 RefSeq files were refused
until `--repair-envelopes` was given, over rows that the survey classed as
documented conventions or exceptions. We chose to stop and name the features,
and to repair only on request, because a parent span that does not contain its
children is contradictory whatever its cause. A user who wants every RefSeq
file to load without a flag will find this strict.

## Limitations

cgload was corrected on the same 30 genomes it is measured on. Fourteen of its
recorded decisions (D-095 to D-108 in the repository) fixed problems that these
files exposed in the first run, the preflights and the second run. Its results
on files from other sources, or on later RefSeq releases, may be worse. The planted faults and the other loaders were not
tuned in this way. All genomes are RefSeq; GFF3 from other annotation pipelines
was not part of this benchmark.

cgload holds the whole file in memory during a load, about 6 KB per row. A
human-sized annotation needs 30 to 37 GB, which rules out a laptop. Counting
features on PostgreSQL took 19.5 s for human, and gene lookups as commands
took about a quarter of a second, mostly to start Python. A killed load leaves
its staging files, and `verify` exits 0 when it skips the source count because
the input file has changed.

The benchmark was written by cgload's authors. We tried to limit the bias that
follows [17, 18] by committing the protocol first, by using a reference reader
that shares no code with any loader, and by reporting our own harness errors.
Three of those errors made other loaders look worse, until we found and fixed
them. Chado stands here for python-chado alone, and the other loaders ran in
parallel, so their load times are indicative only.

# Conclusion

A database built from a GFF3 file deserves the same check as any other copy of
data: does it hold what the file holds? cgload answers that question after every
load, with code that does not share the loader's assumptions, and stops when the
answer is no. On 30 RefSeq genomes it kept every feature, declared its 12
repairs and noticed every fault that changed a count. The protocol, scripts,
results and figures of this benchmark are public, so the same questions can be
put to other loaders and other files.

# Data availability

The benchmark results reported here are in `paper/data/benchmark-run-2.2/` of
the repository. The full results, including the superseded ones, the PostgreSQL
server log and the reference databases, are deposited at Zenodo
[DOI to be added]. The genome panel and the survey results are in `survey/` of
the same repository [9]. All input files are public RefSeq releases, identified
by accession and MD5 in `paper/data/benchmark-run-2.2/inputs.tsv`.

# Code availability

cgload is free software under the MIT license at
https://github.com/computational-genomics-lab/cgload. The benchmark protocol and
scripts are in `benchmark/`, and `paper/figures.py` makes every figure and
table from the deposited results. The version described is tagged
[release tag to be added] and archived at Zenodo [DOI to be added].

# Use of AI tools

Claude (Anthropic) wrote code, ran the analyses and drafted text for cgload,
the benchmark and this manuscript, under the authors' direction. The authors
checked the code, the results and the text, and take responsibility for them.
<!-- AUTHORS: confirm this wording against the journal's current policy. -->

# Funding

[To be completed by the authors.]

# Conflict of interest

None declared.

# Author contributions

[To be completed by the authors, for example with CRediT roles.]

# Biographical notes

[One or two sentences per author, as the journal asks: position, institution
and research interest. To be completed by the authors.]

# References

1. Stein LD. Generic Feature Format Version 3 (GFF3), version 1.26. The Sequence Ontology Project, 2020. https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
2. O'Leary NA, Wright MW, Brister JR, et al. Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Res* 2016;**44**:D733–45. doi:10.1093/nar/gkv1189
3. Dale R. gffutils: GFF and GTF file manipulation and interconversion, version 0.14. https://github.com/daler/gffutils
4. Mungall CJ, Emmert DB, FlyBase Consortium. A Chado case study: an ontology-based modular schema for representing genome-associated biological information. *Bioinformatics* 2007;**23**:i337–46. doi:10.1093/bioinformatics/btm189
5. Sanderson LA, Ficklin SP, Cheng CH, et al. Tripal v1.1: a standards-based toolkit for construction of online genetic and genomic databases. *Database* 2013;**2013**:bat075. doi:10.1093/database/bat075
6. Lawrence M, Huber W, Pagès H, et al. Software for computing and annotating genomic ranges. *PLoS Comput Biol* 2013;**9**:e1003118. doi:10.1371/journal.pcbi.1003118
7. Pagès H, Carlson M, Aboyoun P, et al. txdbmaker: tools for making TxDb objects from genomic annotations. Bioconductor. doi:10.18129/B9.bioc.txdbmaker
8. Smith RN, Aleksic J, Butano D, et al. InterMine: a flexible data warehouse system for the integration and analysis of heterogeneous biological data. *Bioinformatics* 2012;**28**:3163–5. doi:10.1093/bioinformatics/bts577
9. Upadhyay A, Panda A, Tripathy S. Few structural errors in RefSeq genome annotations: a survey of 30 eukaryotes. *Database*, submitted.
10. InterMine documentation. Acceptance tests. https://intermine.readthedocs.io/en/latest/database/data-integrity-checks/acceptance-tests/ (10 November 2026, date last accessed).
11. Dyer SC, Austine-Orimoloye O, Azov AG, et al. Ensembl 2025. *Nucleic Acids Res* 2025;**53**:D948–57. doi:10.1093/nar/gkae1071
12. Ensembl. ensembl-datacheck: code for checking Ensembl databases during release production. https://github.com/Ensembl/ensembl-datacheck
13. Chao K-H. gffbase: a genomic-annotation engine for GFF3 and GTF, version 0.3.0. 2026. https://github.com/Kuanhao-Chao/gffbase
14. Chen B, Wu D, Zhang G. GFFx: a Rust-based suite of utilities for ultra-fast genomic feature extraction. *GigaScience* 2025;**14**:giaf124. doi:10.1093/gigascience/giaf124
15. Gremme G, Steinbiss S, Kurtz S. GenomeTools: a comprehensive software library for efficient processing of structured genome annotations. *IEEE/ACM Trans Comput Biol Bioinform* 2013;**10**:645–56. doi:10.1109/TCBB.2013.68
16. Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format. Zenodo. doi:10.5281/zenodo.3552717
17. Weber LM, Saelens W, Cannoodt R, et al. Essential guidelines for computational method benchmarking. *Genome Biol* 2019;**20**:125. doi:10.1186/s13059-019-1738-8
18. Buchka S, Hapfelmeier A, Gardner PP, et al. On the optimistic performance evaluation of newly introduced bioinformatic methods. *Genome Biol* 2021;**22**:152. doi:10.1186/s13059-021-02365-4
19. Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files. *Bioinformatics* 2011;**27**:718–9. doi:10.1093/bioinformatics/btq671
20. Bonfield JK, Marshall J, Danecek P, et al. HTSlib: C library for reading/writing high-throughput sequencing data. *GigaScience* 2021;**10**:giab007. doi:10.1093/gigascience/giab007
21. Cantalapiedra CP, Hernández-Plaza A, Letunic I, et al. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. *Mol Biol Evol* 2021;**38**:5825–9. doi:10.1093/molbev/msab293
22. Jones P, Binns D, Chang H-Y, et al. InterProScan 5: genome-scale protein function classification. *Bioinformatics* 2014;**30**:1236–40. doi:10.1093/bioinformatics/btu031
23. Palmer JM, Stajich J. Funannotate v1.8.1: eukaryotic genome annotation. Zenodo 2020. doi:10.5281/zenodo.1134477
24. Stanke M, Morgenstern B. AUGUSTUS: a web server for gene prediction in eukaryotes that allows user-defined constraints. *Nucleic Acids Res* 2005;**33**:W465–7. doi:10.1093/nar/gki458
25. O'Leary NA, Cox E, Holmes JB, et al. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. *Sci Data* 2024;**11**:732. doi:10.1038/s41597-024-03571-y
26. Galaxy Genome Annotation. python-chado, version 2.3.9. https://github.com/galaxy-genome-annotation/python-chado
27. Eilbeck K, Lewis SE, Mungall CJ, et al. The Sequence Ontology: a tool for the unification of genome annotations. *Genome Biol* 2005;**6**:R44. doi:10.1186/gb-2005-6-5-r44
28. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 2009;**25**:1422–3. doi:10.1093/bioinformatics/btp163

<!-- AUTHORS: check every reference against the journal's style before
submission. Reference 18's article number and reference 11's pages were taken
from memory; their DOIs were checked. -->

# Supplementary material

**Supplementary Table S1.** The 30 genomes, with cgload's load and verify times,
peak memory and sizes on each database system.

<!-- include tables/tableS1_genomes.md -->

**Supplementary Table S2.** Every planted fault, its expected outcome and what
happened in each of the three genomes.

<!-- include tables/tableS2_faults.md -->

**Supplementary Table S3.** Query timing: median of the first asks and of the
repeated asks, per genome, tool, task and view.

<!-- include tables/tableS3_timing.md -->
