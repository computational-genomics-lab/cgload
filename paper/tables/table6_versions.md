| Software | Version | Used for |
| :-- | :-- | :-- |
| cgload | 0.2.1 | the loader presented; benchmarked at commits 8c999e2 and 9b97178, the same code apart from its version string and comments |
| cgload | 0.2.2 | the check of GFF3 columns 2, 6 and 8, run on 20 November 2026 (tag v0.2.2, commit dbb51ce) |
| cgload | 0.3.0 | the portal (Q9), at commit 6e4c5b5, run on 25 November 2026 |
| Flask, Werkzeug, Jinja2 | 3.1.3, 3.1.9, 3.1.6 | the portal |
| waitress | 3.0.2 | the portal's web server |
| igv.js | 3.8.10 | the portal's genome browser, served with cgload |
| Python | 3.12.15 | cgload, gffutils and the benchmark scripts |
| SQLite | 3.54.1 | cgload and gffutils databases |
| PostgreSQL | 16.15 | cgload and Chado databases |
| SQLAlchemy | 2.1.5 | cgload |
| pyfaidx | 0.9.1 | cgload |
| click | 8.5.0 | cgload |
| psycopg | 3.3.6 | cgload on PostgreSQL |
| HTSlib (tabix, bgzip) | 1.24 | query timing |
| gffutils | 0.14 | loader compared |
| R | 4.5.3 | txdbmaker |
| txdbmaker | 1.6.2 | loader compared |
| GenomicFeatures | 1.62.0 | txdbmaker |
| rtracklayer | 1.70.1 | txdbmaker |
| python-chado | 2.3.9 | loader compared |
| Chado schema | 1.31 (GMOD/Chado commit ef04588) | python-chado's database |
| Sequence Ontology | SO-Ontologies commit 4340c14 | python-chado's database (so.obo) |
| Python, for python-chado | 3.11.17 | python-chado |
| SQLAlchemy, for python-chado | 1.4.54 | python-chado |
| Biopython | 1.81 | python-chado |
| bcbio-gff | 0.7.1 | python-chado |
