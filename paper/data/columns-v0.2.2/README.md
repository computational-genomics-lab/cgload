# GFF3 columns 2, 6 and 8, row by row (cgload 0.2.2)

`benchmark/run_columns.py` on the 30 genomes of the benchmark, on the lab
server on 20 November 2026, with cgload 0.2.2 (commit `dbb51ce`, tag
`v0.2.2`) on SQLite. For each genome it loaded the file, exported it, and
compared the program (column 2), score (column 6) and phase (column 8) of
every row with the input, using `benchmark/columns.py`. This check is not in
the registered protocol; `benchmark/PROTOCOL.md`, section 12, says why it was
added.

- `columns.tsv`: one row per genome, with the counts.
- `columns/ACCESSION.json`: the counts and up to ten examples of each kind.
- `summary.txt`: what `benchmark/summarise_columns.py` printed.

In all 28,211,625 rows no value of the three columns differed. The 12 rows
whose span differed are the 12 features that `--repair-envelopes` repaired.
`input_score_beyond_six_digits` counts the scores that cgload 0.2.1's export
rounded (D-109): 52, of which 49 in human and 3 in mouse.
`input_source_differs_from_first_row` counts the rows that its export gave
their feature's first-row column 2 (D-110): none.
