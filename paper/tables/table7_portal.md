| Request | Requests | Median, ms | 95th percentile, ms | Slowest, ms | Median, ms (PostgreSQL, human) | 95th percentile, ms (PostgreSQL, human) |
| :---------- | ----: | ----: | ----: | ----: | ------: | ------: |
| Search by the gene's name | 6,000 | 23.0 | 39.5 | 130 | 50.0 | 89.4 |
| Gene page | 6,000 | 20.1 | 45.0 | 837 | 48.3 | 141 |
| Region page, gene and 5 kb either side | 6,000 | 19.3 | 42.2 | 95.8 | 33.2 | 66.7 |
| Protein, as FASTA | 4,830 | 18.5 | 42.4 | 558 | 52.7 | 153 |
| 64 KiB of the annotation file | 6,000 | 3.4 | 21.1 | 54.0 | 3.6 | 32.9 |
| 64 KiB of the assembly file | 6,000 | 3.3 | 11.7 | 52.4 | 3.5 | 3.7 |
