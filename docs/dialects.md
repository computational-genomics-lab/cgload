# Dialect profiles, from the real files

Written from the three fixtures in `tests/fixtures/`, and verified by the
assertions in `tests/test_fixtures_real.py` — every claim below is checked by a
test, so a fixture swapped for a tidier one fails loudly rather than silently
removing the trap it was carrying.

Scale: 1187 / 402 / 406 lines, one sequence each, graph-closed (every row's
full ancestor chain is present). Findings that **changed the design** are
marked as such, with the decision they produced.

## Detection

All three are identifiable from column 2 alone, corroborated by a header line.
Detection is cheap and high-confidence, which is what makes D-013's
error-on-ambiguity rule affordable rather than obstructive.

| Profile | Column 2 | Header corroboration |
|---|---|---|
| `refseq` | `RefSeq` | `#!processor NCBI annotwriter`, `#!genome-build-accession` |
| `funannotate` | `funannotate` | none (only `##gff-version 3`) |
| `augustus` | `AUGUSTUS` | `# This output was generated with AUGUSTUS (version …)` |

Column 2 is `source_program` and is stored verbatim per feature (D-014), so
detection and provenance read the same field.

## What each dialect actually provides

| | RefSeq | funannotate | AUGUSTUS |
|---|---|---|---|
| gene row | yes | yes | yes |
| transcript row | `mRNA` | `mRNA` | `transcript` |
| `exon` rows | yes | yes | **no** |
| `CDS` rows | yes | yes | yes |
| `intron` rows | no | no | yes |
| codon rows | no | no | `start_codon`, `stop_codon` |
| explicit `Parent` | yes | yes | yes |
| CDS ID shared across segments | yes | yes | yes |
| exon ID shared across segments | no | no | n/a |
| rows with **no** `ID` | none | none | **175 of 350 children** |
| non-coding genes | tRNA, rRNA, pseudogene | tRNA | none |
| minus-strand segments in file order | mixed | **descending** | ascending |
| comment lines : data lines | 6 : 400 | 1 : 401 | **788 : 399** |
| protein sequence | via `protein_id`, external FASTA | via `protein_id` | **inside comment blocks only** |
| URL-encoded attribute values | yes | no | no |
| score column | `.` | `.` | real values |

## The findings

### 1. 175 of 350 AUGUSTUS child rows carry no `ID` at all — **changed, blocking**

`intron`, `start_codon` and `stop_codon` have only a `Parent`. The counting rule
originally specified — distinct `(seqid, type, ID)` tuples — collapses all 49
`stop_codon` rows into a single tuple, so `verify` would compare 1 against 49
and fail on every AUGUSTUS file, in the one feature the paper's central claim
rests on.

`source_id` is now nullable, and the feature count is
`COUNT(DISTINCT source_id) + COUNT(*) WHERE source_id IS NULL`. **D-024**,
including why synthesising an ID would have broken D-008's independence.

### 2. Two thirds of an AUGUSTUS file is comments, and the proteins are in there — **changed**

788 comment lines to 399 data lines. `# protein sequence = [ ... ]`, continuing
over `#` lines to a closing `]`, is the only place AUGUSTUS emits a
translation — there is no `protein_id` pointing at an external FASTA.

And the free oracle: `sum(CDS spans) / 3 - 1` equals the declared amino-acid
count for **all 49 genes**, no mismatches. AUGUSTUS audits our coordinate
arithmetic using its own output — no network, no reference database, no second
parser. A dropped segment or an off-by-one span breaks the equality and stops
the load. **D-025**; translations go to `annotation/proteins.faa` on disk
(D-004), not into `feature_attribute`, since some are 2 kB and `explain` reads
that table on every call.

### 3. Codon rows sit *inside* the CDS, so mapping them to CDS double-counts — **changed**

`stop_codon` 644-646 lies within `CDS` 644-739; `start_codon` 1838-1840 within
`CDS` 1651-1840. This reversed a vocabulary seed I had written and defended two
decisions earlier. **D-022.** GTF is the opposite case — its stop codon sits
outside the CDS — and that is coordinate reconciliation (D-013 field 5), a
different mechanism that must not be conflated with type mapping.

### 4. funannotate's `-T1` suffix must not be stripped — **changed**

mRNA `PHE_000001-T1` minus `-T1` is `PHE_000001`: the ID of its own parent gene.
True for all 21 transcripts in the fixture. **D-023**: no ID normalisation in
any profile, in either direction.

### 5. AUGUSTUS has no `exon` features — **changed**

Gene, transcript, CDS, intron, start_codon, stop_codon, and nothing else.

**Do not synthesise exons.** With no predicted UTR the exons would coincide with
the CDS segments and duplicate stored information; where AUGUSTUS *does* predict
UTRs, a synthesised exon would be wrong. `stats` reports zero exons for an
AUGUSTUS load, which is the honest answer.

This needs a new profile field — see field 6 below — or "the file had no exons"
and "we dropped the exons" are indistinguishable, which is the ambiguity this
project exists to remove.

### 6. CDS segments share an ID; exon segments do not — **changed**

CDS rows of one transcript repeat a single ID (`g1.t1.cds`,
`PHE_000001-T1.cds`, `cds-XP_011315562.1`); each exon gets its own
(`PHE_000001-T1.exon1`, `exon-XM_011327722.1-1`). Both patterns hold in all
three dialects, so an implementation assuming either rule universally is wrong
half the time. The D-014 discontinuous-feature exemption is exercised by CDS
rows only; a repeated *exon* ID in these dialects is a genuine error.

### 7. funannotate emits minus-strand segments in descending file order

11 of its 21 discontinuous groups. Evidence for the milestone-1 convention that
`segment_index` follows coordinate order, not file order.

### 8. RefSeq hangs mRNA and tRNA off `pseudogene` parents

Present in `ALLOWED_PARENTS` already; now asserted against the real hierarchy
rather than assumed. `gene_biotype` values observed: `protein_coding`, `tRNA`,
`rRNA`, `pseudogene`, `tRNA_pseudogene`.

### 9. RefSeq URL-encodes attribute values; the others do not — **changed**

`strain=PH-1%3B NRRL 31084` — `%3B` is a semicolon that would otherwise split
the attribute. Percent-decoding is mandated by the GFF3 spec and must happen
*after* splitting on `;` and `,`, in every dialect. Before, and you get silent
field-splitting; never, and you get silent data corruption.

### 10. `Dbxref` is a comma-separated list, and funannotate misspells it — **changed**

`Dbxref=GeneID:23558436,GenBank:XM_011327722.1` is two `feature_attribute` rows,
not one. funannotate writes `DBxref` (capital B) and terminates every line with
a trailing `;`. Attribute keys are matched case-insensitively and stored as
written; an empty final field is discarded rather than becoming a key with an
empty value.

### 11. RefSeq's `region` row is region metadata, not a gene feature — **changed**

The first data line carries `Dbxref=taxon:229533` and
`strain=PH-1%3B NRRL 31084`, and `##sequence-region NC_026477.1 1 8033942` gives
the length. The profile uses these to **cross-check** the FASTA-derived
`sequence_region` rows and the user-supplied `ncbi_taxon_id`; a mismatch is a
warning naming both values, never a silent overwrite. D-016 stands — the config
is authoritative, because the file is the thing being validated.

### 12. RefSeq marks partial features, and the flags matter

`partial=true;start_range=.,3250;end_range=5086,.` means the feature runs off the
end of the assembled sequence. Stored as attributes. A partial CDS whose length
is not a multiple of three is expected and must not be flagged as an error,
which a naive validity check would do — and which would also break the AUGUSTUS
oracle in finding 2 if applied indiscriminately across dialects.

### 13. Comment lines are interleaved with data — **changed**

AUGUSTUS writes `# start gene g1` and whole evidence blocks *between* feature
rows, not only in a header. The skip rule is "first character is `#`", applied
identically by the tokenizer and by D-008's naive counting pass. No comment line
in any fixture contains a tab, so a field-count heuristic is unnecessary and the
simple rule is safe.

## Consequence for the profile fields

D-013's five declared fields become five again, but not the same five:

1. **Detection signature** — column 2, plus an optional header regex.
2. ~~ID normalisation~~ — **deleted** (D-023). IDs are stored verbatim.
3. **Parent inference** — explicit `Parent` for all three of these dialects.
   Synthesis is needed only for GTF (milestone 5).
4. **Type mapping** — into the controlled vocabulary, non-lossy (D-017).
5. **Coordinate reconciliation** — a no-op for all three of these. GTF's stop
   codon absorption is the only known case.
6. **Provided-types declaration** *(new, finding 5)* — which feature types this
   dialect is expected to emit. AUGUSTUS declares no `exon`, so zero exons is a
   recorded expectation; a RefSeq file with zero exons is a warning.
7. **Comment extraction** *(new, finding 2)* — optional, dialect-specific. Only
   AUGUSTUS declares any: the protein blocks, plus the load-time length check.
   Hint-support statistics are deliberately not stored in v2.0 — they are
   prediction diagnostics, not annotation, and no query needs them.
