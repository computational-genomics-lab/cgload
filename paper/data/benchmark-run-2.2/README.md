# Benchmark results reported in the paper

The results folder of the benchmark run on the lab server, as `benchmark/run.sh`
left it after `benchmark-run-2.2`. `benchmark/PROTOCOL.md` defines every table
and column; its section 12 says which tag each result comes from.

Left out here, and kept on the lab server: the PostgreSQL server log, the
reference databases (`truth/*.sqlite`, which `benchmark/truth.py` rebuilds from
the input files), and the results that later repeats replaced (`superseded/`).

Units. `*_peak_mb` columns are MiB (2^20 bytes): `ru_maxrss` in KiB divided by
1,024. The sizes that `summarise.py` writes to `tables/` (`*_db_mb`,
`data_dir_mb`, `max_leftover_mb`) are MB (10^6 bytes). The paper gives every
size in binary units (KiB, MiB, GiB), so 130.4 MB in `tables/kills.tsv` is the
124 MiB of the paper's Table 6.

`paper/figures.py` makes the paper's tables and figures from this folder.
