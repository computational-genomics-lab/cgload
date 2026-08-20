# Fixture provenance

Real annotation output, sliced to one sequence and to a graph-closed set of
records (every row's full ancestor chain is present). Kept small deliberately;
each is a few hundred lines.

| File | Source | Sequence | Notes |
|---|---|---|---|
| `refseq_fgraminearum_NC_026477.gff3` | NCBI RefSeq, *Fusarium graminearum* PH-1, assembly ASM24013v3 (`GCF_000240135.3`) | `NC_026477.1` | protein-coding, tRNA, rRNA and pseudogene loci; URL-encoded attribute values; partial features |
| `funannotate_aphaeospermum_pheo_arth1.gff3` | funannotate, *Apiospora phaeospermum* | `pheo_arth1` | `-T1` transcript suffix; `DBxref` spelling; trailing `;`; minus-strand segments in descending file order |
| `augustus_aphaeospermum_pheo_arth1.gff3` | AUGUSTUS 3.4.0, *Apiospora phaeospermum* | `pheo_arth1` | no `exon` rows; `intron`/`start_codon`/`stop_codon` rows carry no `ID`; real score values; **per-gene comment blocks retained verbatim, including the protein-sequence blocks** — 788 comment lines to 399 data lines |

Provided by the project authors from their own analyses. Retained in the
repository, and therefore in the Zenodo archive, as test data.

## genbank_maeruginosa_NC_010296.gbff

*Microcystis aeruginosa* NIES-843, RefSeq NC_010296.1 (GCF_000010625.1),
downloaded as gzipped GBFF. Whole locus blocks only: no gene is missing its
child and no child is missing its gene.

28 loci: 15 plain protein-coding, 4 pseudogenes (no `/translation`, truncated
coordinates), 2 with `join()` CDS spans, and one block each of tRNA, rRNA,
ncRNA and tmRNA. Carries the prokaryotic shape — gene → CDS with no transcript
layer — which is what selects the `locus_tag` linkage strategy.

## genbank_dscam1_CG12164.gbff

*Drosophila melanogaster* Release 6, NT_033778.4 region 7,312,924..7,386,899.
Five genes: three Dscam1 isoforms plus the whole of CG12164.

CG12164 is the reason this fixture exists. Its two isoforms are listed with
transcripts in the order RB, RA and coding sequences in the order RA, RB, so
pairing by file order attaches each protein to the wrong transcript. Both
transcripts also contain both coding sequences, so containment alone cannot
separate them either — only the `/note` "from transcript" text can.

Measured on the full region before trimming: 79 transcripts, 75 isoforms on
Dscam1 alone; containment resolves 75 of 79, the note resolves the remaining 4,
positional pairing is wrong for 2.

---

## Reproducing the measurements

Every figure quoted as a *full-record* number in `MILESTONE-6.md` comes from the
complete NCBI record, not from the trimmed fixture beside it. The fixtures are
smaller on purpose. Commands below regenerate each number from a fresh download,
so nothing rests on trust.

### Downloads

```console
# Microcystis aeruginosa NIES-843
datasets download genome accession GCF_000010625.1 --include gbff
# Fusarium graminearum PH-1
datasets download genome accession GCF_000240135.3 --include gbff,gff3
# Drosophila melanogaster, Dscam1 region
efetch -db nuccore -id NT_033778.4 -format gbwithparts \
       -seq_start 7312924 -seq_stop 7386899 > dscam1_region.gbff
```

### Full-record counts

```console
# coding records, and how many use join() -- the slippage rate
grep -c '^     CDS  ' NC_010296.gbff                     # 5756
grep -A1 '^     CDS  ' NC_010296.gbff | grep -c 'join('   # 37   (0.64%)
grep -c '/ribosomal_slippage' NC_010296.gbff              # 37

# genes, and pseudogenes among them
grep -c '^     gene ' NC_010296.gbff                      # 5808
grep -c '/pseudo' NC_010296.gbff                          # 526

# the qualifier question: /transcript_id by feature type
awk '/^     (mRNA|CDS) /{f=$1} /\/transcript_id/{print f}' dscam1_region.gbff \
  | sort | uniq -c        # mRNA only; never CDS

# legitimately non-conforming CDSs in F. graminearum
grep -c '/artificial_location' GCF_000240135.3.gbff       # 78
```

### How each fixture was sliced

Whole locus blocks only, so that no feature is missing an ancestor and no
ancestor is missing a child. Specifically:

- **`genbank_maeruginosa_NC_010296.gbff`** — first 28 `gene` blocks that
  together cover every shape present in the genome, plus a deliberate selection
  rule of `join()` span count ≥ 2 to capture two `/ribosomal_slippage` records.
  **That rule is why the fixture holds two: it does not reflect their frequency
  in the genome, which is 0.64%.** Any claim about how common a feature is must
  come from the full record, never from this file.
- **`genbank_dscam1_CG12164.gbff`** — region 7,312,924..7,386,899 of
  `NT_033778.4`, chosen because it contains the whole of CG12164, whose two
  isoforms defeat both positional pairing and containment.
- **The three GFF3 fixtures** — one sequence each, cut at locus boundaries.

### Known gap in the fixture set

No fixture currently pairs a GBFF and a GFF3 **for the same assembly**. The
GenBank fixtures are *M. aeruginosa* and *D. melanogaster*; the RefSeq GFF3 is
*F. graminearum*. A matched pair would give a free cross-parser check — load both
representations of one assembly and compare feature counts per type, where the
two parsers share no code. `GCF_000240135.3` is available in both formats and one
download would close this. Worth doing at milestone 8, when the completeness
check needs every independent oracle it can get.

---

## Functional annotation fixtures (milestone 7) — SYNTHETIC

`eggnog_fgraminearum.emapper.annotations` and `interproscan_fgraminearum.tsv`
are **constructed, not downloaded.** This is the first milestone whose fixtures
are not real tool output, and that is a weaker footing than milestones 3–6 —
stated here rather than left to be discovered.

What they are: the real column layouts of eggNOG-mapper 2.1.12 and
InterProScan 5 TSV, filled with plausible values, and **keyed to the 20 real
`protein_id` values in `refseq_fgraminearum_NC_026477.gff3`** so that the join
against real features is genuinely exercised.

Traps deliberately included, each of which a naive reader gets wrong:

- `COG_category = EGP` — three categories in a separator-free letter run. A
  comma split yields one bogus accession named `EGP`.
- `-` placeholders in most columns, which must not become accessions.
- Multi-valued `GOs`, `KEGG_ko`, `eggNOG_OGs` cells.
- An InterProScan score column holding an e-value (`2.1E-9`) for HMM analyses
  and a bit score (`43.7`) for others.
- The same Pfam signature (`PF00096`) matching one protein at two positions,
  which must survive as two rows.
- Optional InterProScan columns 12–14 present on some rows and absent on others.

**What still needs real files.** Real eggNOG and InterProScan output for one of
these assemblies would confirm the column layouts against the actual tool
versions, and would exercise the match-rate path against genuine FASTA headers —
which is where the join is most likely to break, because whatever header the user
fed the tool is what comes back. Until then, milestone 7 is validated against the
format as documented rather than as emitted.

## eggnog_mycmay_emapper214.annotations

Real `emapper-2.1.4` output, run by the cgload authors on the RefSeq protein
set of *Mycosarcoma maydis* (*Ustilago maydis*) GCF_000328475.2. The complete
file is 6,509 proteins and 5.2 MB; the shipped fixture is an 11-row slice with
the original `##` header block and trailer preserved.

Rows were selected to carry one trap each, not at random:

| query | trap |
|---|---|
| XP_011385997.1, XP_011386076.1 | `COG_category = IJT` — three categories, no separator |
| XP_011386956.1 | 570 GO terms on one protein |
| XP_011386301.1, XP_011386313.1 | CAZy present (rare: 95 of 6,509) |
| XP_011385992.1 | `KEGG_Pathway` carries `ko00362` and `map00362` — the same pathway twice |
| XP_011385994.1, XP_011386008.1 | comma inside `Description` |
| XP_011385990.1 and two others | `EC` and `PFAMs` both populated |

Measured on the **full** file, not the fixture:

- 6,509 queries, all with `XP_` RefSeq protein accessions as the query ID
- 380,022 hits before the D-031 pathway fix, 371,409 after
- GO is 286,671 of those — mean 89.7 terms per annotated protein, max 570
- 61,425 cells contain `-`; zero contain an empty string
- 97 of 121 distinct `COG_category` values are multi-letter
- 2,253 rows carry `KEGG_Pathway`; in **all 2,253** the `ko` and `map` sets are
  identical
- zero duplicate (query, analysis, accession) triples

Reproduce with:

```console
python - <<'PY'
from cgload.parsers.functional import parse_eggnog
from collections import Counter
ff = parse_eggnog("Mycmay_chromosome_emapper.annotations")
print(len(ff.queries), len(ff.hits))
print(Counter(h.analysis for h in ff.hits).most_common())
PY
```

**Known gap:** this file is *U. maydis* and no genome fixture in this repository
is that assembly, so it exercises parsing but not the join to real features. A
matching `.faa` or RefSeq GFF3 would close that.

## refseq_mycmay_GCF_000328475.gff3

RefSeq GFF3 for *Mycosarcoma maydis* (*Ustilago maydis*) GCF_000328475.2 —
**the same assembly as `eggnog_mycmay_emapper214.annotations`**. This pair is the
only one in the repository where a genome and a functional-annotation file
describe the same proteins, so it is the only fixture that tests the join rather
than just the parsing.

Measured on the full pair before trimming:

- 6,780 `protein_id` values on CDS rows in the GFF3
- 6,509 query IDs in the eggNOG file
- **all 6,509 eggNOG queries appear as a CDS `protein_id`** — zero unmatched
- 271 proteins in the GFF3 have no eggNOG hit

The 271 matter: they are the opposite direction from the unmatched-query guard
and must not trigger it. A genome where some genes carry no functional
annotation is the normal case.

The shipped slice is 15 whole gene blocks, graph-closed: 11 whose proteins are in
the eggNOG fixture, 3 whose proteins eggNOG returned no hit for, and one tRNA gene
with no protein at all.

---

## The matched pair — the only real join in the repository

`refseq_mycmay_GCF_000328475.gff3` and `eggnog_mycmay_emapper214.annotations`
describe the **same assembly** (*Mycosarcoma maydis* / *Ustilago maydis*,
GCF_000328475.2). Every other functional fixture is parse-tested only: no genome
in the repository describes its proteins.

Loading this pair end to end through `cgload load` is what exposed D-043 — the
transcript-bypass rule refusing a correct RefSeq genome. That rule had passed 272
tests, including three GFF3 fixtures, because none of them contained a tRNA gene
whose exons NCBI emits twice.

The remaining paired-fixture gap is a **GBFF and a GFF3 for one assembly**, which
would give a cross-parser count comparison between two parsers sharing no code —
the strongest independent oracle available to milestone 8. `GCF_000240135.3` is
published in both formats and one download closes it.
and must not trigger it. A genome where some genes carry no functional
annotation is the normal case.

The shipped slice is 15 whole gene blocks, graph-closed: 11 whose proteins are in
the eggNOG fixture, 3 whose proteins eggNOG returned no hit for, and one tRNA gene
with no protein at all.

## ncbi_feature_count_GCF_000328475.txt

NCBI's published per-assembly feature totals for *Mycosarcoma maydis*
GCF_000328475.2, verbatim -- seven data rows, unmodified. Downloaded from the
assembly's FTP directory alongside the GFF3.

Carries all three format traps documented in D-034: `gene` split across three
`Class` rows, `na` in the `Unique Ids` column of the `tRNA` row, and an
`Assembly-unit accession` (GCF_000328485.2) that differs from the assembly
accession by one digit.

Measured against the full GFF3 for the same assembly -- all five shared types
match exactly:

| type | GFF3 rows | GFF3 features | NCBI |
|---|---|---|---|
| gene | 6,907 | 6,907 | 6,907 |
| mRNA | 6,780 | 6,780 | 6,780 |
| CDS | 9,739 | **6,780** | 6,780 |
| tRNA | 111 | 111 | 111 |
| rRNA | 34 | 34 | 34 |

CDS matching at feature level rather than row level is what makes this an
independent check on segment grouping.

## The cross-parser pair (not shipped)

`GCF_000240135_3_ASM24013v3_genomic.gff` and `.gbff` -- *F. graminearum* PH-1,
the same assembly in both formats. Too large to ship; download both from the
assembly's FTP directory to reproduce.

| type | GFF3 features | GBFF records |
|---|---|---|
| gene | 13,714 + 11 pseudogene = 13,725 | 13,725 |
| mRNA | 13,315 | 13,315 |
| CDS | 13,312 | 13,312 |
| tRNA | 322 | 322 |
| rRNA | 88 | 88 |

---

## The matched pair — same assembly, two formats (D-051)

`pair_fgraminearum_NC_026477.gff3` and `pair_fgraminearum_NC_026477.gbff` are the
same 43 genes of *Fusarium graminearum* PH-1 chromosome 4 (`NC_026477.1`,
1–120,000), sliced from the full GCF_000240135.3 release in both formats and
graph-closed in each: no feature is missing an ancestor, no ancestor a child.

They exist so the cross-format check — the strongest evidence in the project —
can be reproduced by running `pytest`, not only by downloading 100 MB.

```console
datasets download genome accession GCF_000240135.3 --include gff3,gbff
```

### Measured on the full release, before slicing

Both formats read by parsers sharing no code:

| | GFF3 | GenBank |
|---|---|---|
| gene | 13,714 (+11 pseudogene) | 13,725 |
| mRNA | 13,315 | 13,315 |
| CDS | 13,312 | 13,312 |
| tRNA | 322 | 322 |
| rRNA | 88 | 88 |
| exon | 37,933 | — (implicit in join spans) |
| gap | — | 414 |

Exact on every shared type. The gene row differs only in naming: GFF3 gives
pseudogenes their own type, GenBank calls them genes and marks them with a note,
and 13,714 + 11 = 13,725.

Every GenBank CDS-to-transcript edge — all 13,312 — was resolved by containment
alone, with no ambiguity failures and no use of the `/note` fallback.

```console
cgload load --config fg.toml
cgload verify --assembly-id 1 --against fg.gbff
```

```
  source count      PASS  78,715 features in the file, 78,715 in the database
  cross-parser      PASS  GFF3 and GenBank agree on 5 feature types
```

Load: 30 s for 78,715 features across 19 sequences. Verification: 1.3 s.

### Two bugs this pair exposed

Both were found by loading it, not by reading it, and both had passed 311 tests.

1. A `misc_feature` with a `join()` location shares one ID across its segments,
   exactly as a CDS does. The shared-ID rule was an allow-list of `{"CDS"}` and
   refused the record (D-049).
2. This release uses `/locus_tag` and never `/gene`, so containment's gene
   grouping — which read only `/gene` — put every CDS in the genome into a single
   group and failed on the first record (D-050).

**Note on `eggnog_mycmay_emapper214.annotations`:** the third `##` header line
records the emapper command line, and in the original file that line contained
local filesystem paths including a username. It has been replaced with an
equivalent generic path. Nothing else in the file is altered; the version line,
the column header and every data row are verbatim. The substitution preserves
the property the fixture tests -- that the command-line row contains the string
`emapper` and must not be mistaken for the version declaration.
