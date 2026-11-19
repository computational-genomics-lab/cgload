| Loader | Genome | Kills | State afterwards | Notes |
| :-- | :-- | --: | :-- | :-- |
| cgload, new load | Arabidopsis thaliana | 10 | 10 empty | the repeated load succeeded 10 times; 9 left staging files, up to 124 MiB |
| cgload, new load | Fusarium graminearum | 10 | 10 empty | the repeated load succeeded 10 times; 7 left staging files, up to 35 MiB |
| cgload, `--force` reload | Arabidopsis thaliana | 10 | 10 intact | the earlier load passed verify 10 times; 6 left staging files, up to 116 MiB |
| cgload, `--force` reload | Fusarium graminearum | 10 | 10 intact | the earlier load passed verify 10 times; 6 left staging files, up to 35 MiB |
| gffutils, create_unique | Arabidopsis thaliana | 10 | 10 empty | gffutils opens 0 of the 10 files left |
| gffutils, create_unique | Fusarium graminearum | 10 | 3 complete, 7 empty | gffutils opens 0 of the 10 files left |
| Chado (python-chado) | Fusarium graminearum | 5 | 5 empty | – |
