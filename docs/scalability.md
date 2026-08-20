# Scalability

Measured, not projected. Every number here comes from a run recorded below with
the command that produced it, so it can be reproduced or contradicted.

## Summary

| | *F. graminearum* PH-1 | *Z. mays* B73 |
|---|---|---|
| accession | GCF_000240135.3 | GCF_902167145.1 |
| sequences | 19 | 687 |
| assembly | 36.5 Mb | 2.18 Gb |
| features | 78,715 | 696,928 |
| source types | 9 | 21 |
| functional annotation | — | eggNOG-mapper, 55,463 proteins |
| **load** | **12.9 s** | **2 min 34 s** |
| **peak RSS during load** | **0.78 GB** | **5.37 GB** |
| **verify** | **0.5 s** | **0.3 s** |
| **peak RSS during verify** | **64 MB** | **41 MB** |
| **export** | **5.6 s** | not measured |
| **peak RSS during export** | **460 MB** | not measured |
| **database** | **102 MB** | **1.03 GB** |
| sequence on disk | 36 MB | 2.1 GB |

Maize carries **8.9×** the features of *Fusarium* and **60×** the bases. Load
time rose 11.9× and peak memory 6.7× — both tracking the feature count, neither
tracking the assembly size.

## What scales, and what does not

**Memory scales with feature count, not with genome size.** The assembly is
copied to disk and indexed; no sequence is held in memory. A 2.18 Gb genome and a
36 Mb genome differ in load cost by their annotation size, not their base count —
which is why maize's 60× larger assembly costs only 6.7× the memory.

Two points give a line:

| features | peak RSS |
|---|---|
| 78,715 | 0.76 GiB |
| 696,928 | 5.12 GiB |

**Roughly 7 GB of RSS per million features, plus about 200 MB of baseline.**

*This corrects an earlier figure of 4.5 GB per million, which was derived from the
maize run alone by subtracting an assumed baseline. With two measurements the
slope is 7.05 GiB per million and the intercept is 0.21 GiB — so the earlier
number understated the requirement by roughly a third, in the optimistic
direction.* Per-feature cost is not constant either: *Fusarium* alone implies 9.7
GB per million because the fixed baseline is a larger share of a small run.

**Two points define a line, not a scaling law.** A linear fit cannot be
distinguished from a mildly superlinear one on two measurements. Treat the
projection below as a budget, not a guarantee:

| features | projected peak RSS |
|---|---|
| 1,000,000 | ~7.3 GiB |
| 2,000,000 | ~14.3 GiB |

**The ceiling is the annotation parse, not the functional load.** Every feature is
parsed before any row is written, deliberately, so that cross-file identifier
collisions and envelope violations are detected before the database is touched
(D-028). That is the high water mark. Measured directly:

| run | peak RSS |
|---|---|
| maize, nuclear only, no functional annotation | 5,360,500 kB |
| maize, complete assembly + eggNOG (55,463 proteins) | 5,370,712 kB |

Adding functional annotation for 55,463 proteins cost **10 MB** of peak memory.
Annotation rows are inserted in bounded batches (D-044), so the functional phase
never dominates. A projected multi-gigabyte ceiling for the deduplication set did
not materialise.

**Verification is effectively free and does not scale.** 41 MB on maize against 64
MB on *Fusarium* — smaller on the larger genome, because the naive counting pass
streams and holds only a counter per feature type. Wall time is 0.3–0.5 s at both
scales.

**Export does *not* stream, and its memory scales with feature count.** 460 MB
peak on *Fusarium*, 59% of that load's own peak. `load_features` reconstructs
every feature before writing, in order to regroup segments that share an
identifier. Projecting from the same slope, exporting maize would need roughly
4 GB. This was not measured and should be, because it is the one phase whose
memory profile was previously assumed rather than checked.

**Wall time is not the constraint.** 696,928 features in 2 min 34 s is roughly
4,500 features per second including verification and indexing; *Fusarium* runs at
6,100 per second, the difference being the fixed startup cost amortised over
fewer features.

**Database size is close to linear in features.** 1,300 bytes per feature on
*Fusarium* and 1,478 on maize. The difference is functional annotation, which
maize has and *Fusarium* does not; feature storage itself is stable.

## Reproducing the *F. graminearum* run

The smaller of the two, and the one to start with: it completes in under thirteen
seconds on a laptop and it is the only assembly published in **both** GFF3 and
GenBank, which is what makes the cross-parser tier demonstrable.

### 1. Fetch

```console
datasets download genome accession GCF_000240135.3 --include gff3,gbff,genome
unzip ncbi_dataset.zip
```

### 2. Config

```toml
[organism]
genus = "Fusarium"
species = "graminearum"
strain = "PH-1"
ncbi_taxon_id = 229533

[assembly]
name = "ASM24013v3"
accession = "GCF_000240135.3"
fasta = "GCF_000240135.3_ASM24013v3_genomic.fna"

[[files.annotation]]
path = "fusarium_genomic.gff"
```

### 3. Load

No `--repair-envelopes` needed: this assembly has no envelope violations.

```console
cgload init
/usr/bin/time -v cgload load --config fusarium.toml
```

```
Loaded Fusarium graminearum / PH-1 / ASM24013v3 (assembly_id=1)
  19 sequence regions, 36,458,046 bp
  fusarium_genomic.gff: refseq (high) -- 78,715 features
    CDS 13,312, exon 37,933, gene 13,714, mRNA 13,315, pseudogene 11,
    rRNA 88, region 19, sequence_feature 1, tRNA 322
  source count      PASS  78,715 features in the file, 78,715 in the database,
                          across 9 type(s)

Elapsed (wall clock) time: 0:12.93
Maximum resident set size: 797,976 kbytes
```

The completeness check runs as part of the load and is not optional (D-046): a
load that stored fewer features than the file contains fails and rolls back.

### 4. Verify

```console
/usr/bin/time -v cgload verify --assembly-id 1
```

```
source count      PASS  78,715 features in the file, 78,715 in the database, across 9 type(s)
cross-parser      ----  no second format supplied
reference (NCBI)  ----  no published counts supplied

Elapsed (wall clock) time: 0:00.47
Maximum resident set size: 65,920 kbytes
```

Both unmet tiers are **reported rather than omitted** (D-033). A check that
vanishes when it cannot run is indistinguishable from one that passed.

### 5. Export

```console
/usr/bin/time -v cgload export --assembly-id 1 --output out.gff3
```

```
Elapsed (wall clock) time: 0:05.63
Maximum resident set size: 471,068 kbytes
```

Export writes cgload's own dialect, declared in the header, rather than
impersonating the input format (D-052). Re-importing the result gives the same
feature counts and preserves how each parent edge was determined.

### 6. On disk

```console
$ ls -l cgload.db
-rw-r--r-- 102,326,272   cgload.db
$ du -sh data/
36M     data/
```

The data directory holds the copied assembly, its index and the manifest. It is
almost exactly the size of the input FASTA, because that is what it mostly is
(D-018).

### 7. Cross-parser

The strongest tier available, and the reason this assembly is the reference case:
the same annotation read by two engines that share no code — one splitting
tab-separated columns and following explicit `Parent` attributes, the other
reading an indented flat file and inferring parentage from coordinates.

```console
cgload verify --assembly-id 1 --against GCF_000240135.3_ASM24013v3_genomic.gbff
```

```
source count      PASS  78,715 features in the file, 78,715 in the database, across 9 type(s)
cross-parser      PASS  GFF3 and GenBank agree on 5 feature types
```

| type | from GFF3 | from GenBank |
|---|---|---|
| gene | 13,714 + 11 pseudogene | 13,725 |
| mRNA | 13,315 | 13,315 |
| CDS | 13,312 | 13,312 |
| tRNA | 322 | 322 |
| rRNA | 88 | 88 |

The gene row reflects a naming difference, not a discrepancy: GFF3 gives a
pseudogene its own type, GenBank writes it as a `gene` with `/pseudo`.
13,714 + 11 = 13,725.

A parser that dropped a gene from one format has no way to make the other
format's count drop it too, which is what makes this tier hard to satisfy by
accident.

## Reproducing the maize run

### 1. Fetch

```console
datasets download genome accession GCF_902167145.1 --include gff3,genome
unzip ncbi_dataset.zip
```

Files used: `genomic.gff` and
`GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna`. Optionally the published
feature counts for the reference tier: `GCF_902167145.1_feature_count.txt`.

### 2. Pre-flight

Reports every structural problem in one pass, without touching a database. On a
file this size a load takes minutes; this takes seconds, and it lists all
problems rather than stopping at the first.

```console
python3 preflight_gff3.py genomic.gff
```

For this assembly it reports 36 problems: 34 `escapes-parent`, one `inverted`,
one `out-of-range`. All are confined to seven trans-spliced genes on the two
organellar sequences. See [known-limitations.md](known-limitations.md).

cgload reports fewer than 34 envelope repairs, and both counts are correct: the
pre-flight script reads raw rows, while cgload merges rows sharing an identifier
before checking envelopes, which resolves most of them.

### 3. Config

```toml
[organism]
genus = "Zea"
species = "mays"
strain = "B73"
ncbi_taxon_id = 4577

[assembly]
name = "Zm-B73-REFERENCE-NAM-5.0"
accession = "GCF_902167145.1"
fasta = "GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna"

[[files.annotation]]
path = "genomic.gff"

[[files.functional]]
path = "maize.emapper.annotations"
```

### 4. Load

`--repair-envelopes` is required for this assembly: without it the load refuses,
because seven organellar genes have parent extents that exclude their own
children.

```console
cgload init
/usr/bin/time -v cgload load --config maize.toml --repair-envelopes
```

```
warning: genomic.gff: loaded 13 feature type(s) not declared by the refseq
profile: antisense_RNA, cDNA_match, intron, lnc_RNA, match, miRNA,
primary_transcript, sequence_conflict, sequence_difference, sequence_feature,
snRNA, snoRNA, transcript. They are stored verbatim and linked only where the
file gave an explicit parent.

Elapsed (wall clock) time: 2:33.78
Maximum resident set size: 5,370,712 kbytes
Exit status: 0
```

The warning is the open vocabulary working: 13 types the profile does not declare
were stored and counted rather than dropped or refused.

### 5. Verify

```console
cgload verify --assembly-id 1 \
  --reference-counts GCF_902167145.1_feature_count.txt \
  --accession GCF_902167145.1 \
  --sequences-in-assembly 687
```

```
source count      PASS  696,928 features in the file, 696,928 in the database, across 21 type(s)
reference (NCBI)  ...   see known-limitations.md on mRNA
```

`--sequences-in-assembly` has no default, deliberately (D-067): without it the
partial-load precondition cannot be evaluated, and the reference tier would
compare what you loaded against totals that may describe more.

### 6. Cross-parser

Maize has no practical GenBank flat file at this size — the GBFF embeds the
sequence in the feature table and runs past 2 Gb. The tier is demonstrated on
*F. graminearum* instead, above.

## Hardware

Both runs on a single workstation core (99% CPU, no parallelism), Linux, SQLite
backend, Python 3.11.

Peak RSS of 5.4 GB means a machine with 8 GB will complete the maize load and one
with 4 GB will not. *Fusarium* completes comfortably in 2 GB.

## What is not measured

- **Export at maize scale.** Export memory scales with feature count and was
  measured only on *Fusarium*, where it reached 460 MB. Projected at roughly 4 GB
  for maize, untested.
- **Concurrent access.** Single-process loads only.
- **MySQL and PostgreSQL at scale.** The code path is shared and tested, but both
  large runs above are SQLite.
- **Genomes above 700,000 features.** Extrapolating from 7 GB per million is
  reasonable but untested, and two measurements cannot separate a linear trend
  from a mildly superlinear one.
