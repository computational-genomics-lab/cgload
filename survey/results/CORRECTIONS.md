# Corrections to the survey's results

## Medaka (GCF_002234675.1), 23 October 2026

**What was wrong.** The survey read an incomplete copy of the medaka RefSeq
GFF3. Its download, `survey/downloads/GCF_002234675.1/genomic.gff` on the lab
server, has 292,477 lines (MD5 `e4475fda2191536c132c447f391cca87`), and its
last line stops in the middle of an attribute value:

    NC_019863.2  Gnomon  CDS  29575426  29575527  .  +  0  ID=cds-XP_004069222.1;...;Name=XP_00406922

So the survey's medaka row described 292,462 rows on 5 of the assembly's 25
sequences. `run_panel.sh` reuses a `genomic.gff` that is already present and
does not test the archive it came from, so a download cut short was surveyed
as if it were whole.

**How it was found.** The benchmark checks that each annotation file has the
MD5 the survey recorded (`benchmark/PROTOCOL.md`, section 2). Medaka's did not
match. NCBI's file, fetched for the benchmark on 21 and 22 October 2026, had
the same MD5 both times (`db602118b6d3d858ed0f09750d3afb86`): 1,399,107 lines,
1,399,050 feature rows on 25 sequences. The other 29 files match the survey's
exactly.

**What was done.** The survey was run again on the complete file, with the same
`cgsurvey.py` (MD5 prefix `e37f86998b8f`, as in the first run) and the same
arguments as `run_panel.sh` gives it. Its row replaced the old one in
`survey.tsv`; `downloaded` is when the file was fetched (the start of the
benchmark preflight's download step). The run is recorded in `provenance.tsv`.

**What changed.** The counts that describe the file: 1,399,050 rows (was
292,462), 1,368,195 parent-child edges (285,969), 25 sequences (5), 1
organellar sequence (0), 587,965 multi-row IDs (122,583), 669 exons directly
under a pseudogene (108). The findings did not change: medaka has no defect of
any kind, and `hits.tsv` has no medaka rows, before or after.
