# Deliberately malformed fixtures

Files here are **not valid input**. They reproduce defects observed in published
provider data, and every one of them makes cgload refuse a load by default.

They live in a subdirectory rather than alongside the real fixtures because
several tests iterate "every fixture" and assert it loads cleanly, parses, round
trips, or reports no envelope violations. A broken file in that set would either
break those tests or, worse, force them to be weakened with per-name exclusions
until they stopped asserting anything.

| File | Defect | Decision |
|---|---|---|
| `trans_spliced_envelope.gff3` | gene envelope excludes its own exons; two levels wrong | D-046, D-055 |
| `inverted_row_envelope.gff3` | row written start > end, with neither number an extent | D-056 |

Both are reduced from *Zea mays* B73 GCF_902167145.1, organellar genes, where
NCBI's GFF3 converter writes a trans-spliced gene's extent from one scattered
segment instead of from all of them. The correct spans are known because NCBI
publishes the same assembly in GenBank format, where they are complete.
