# The package name, and how to change it

## History

Two earlier names were rejected after collision checks, both discovered late:

| name | collision | found by |
|---|---|---|
| `tessera` | a tissue-segmentation algorithm, bioRxiv 2025, active repo | bioRxiv search |
| `sextant` | PyPI 0.5.0, "CLI tool to work with helm charts"; the console script collides | `pip index versions` |

The current name is **`cgload`** — "comparative genomics load". Verify it
yourself before the first release; see the six checks below.

> **A warning about this document.** An earlier version of it described the
> `sextant` collision. When `rename.py` was run, it rewrote every occurrence of
> the old name — including inside the prose *about* the collision — producing a
> document that claimed the new name was taken by a Helm CLI. A mechanical
> rename cannot distinguish a name in use from a name being discussed. **After
> any rename, read this file rather than trusting the script's verification
> output.**

## Renaming

If the name must change again, use `rename.py` rather than a sed one-liner:

```console
python3 rename.py <newname>
```

## Why a script

Occurrences fall into three classes, in increasing order of how quietly a miss
breaks things.

**Imports and paths** — `src/<name>/`, `from <name>.x import y`, the `[project]
name`, the console script entry point. A miss is a loud `ImportError`.

**Prose** — README, docs, docstrings. A miss is embarrassing, not dangerous.

**On-disk format strings** — and these are the reason for the script:

| constant | value | consequence of a partial rename |
|---|---|---|
| `EXPORT_PRAGMA` | `##<name>-export` | exported files stop identifying themselves; a re-import is misdetected as the dialect the data came from |
| `LINKAGE_ATTRIBUTE` | `<name>_linkage_method` | linkage provenance is written under one name and read under another, so every computed edge is silently recorded as having come from the file |
| vocabulary pragma | `##<name>-vocabulary` | the hash is written but not recognised |

`LINKAGE_ATTRIBUTE` is declared as a literal in three modules that deliberately
do not import one another, so that the exporter and the reader cannot drift into
sharing an assumption. A partial rename leaves them disagreeing — **and every
test still passes**, because each module is internally consistent. That exact
failure has occurred once before.

The script verifies afterwards that no occurrence of the old name survives and
that the constants agree across modules. It reads them rather than importing, so
it works before the package is installed.

## Choosing a name

Six checks, all of them, before committing:

1. **PyPI** — `pip index versions <name>`, and open `pypi.org/project/<name>`.
2. **GitHub** — search the name plus "genomics" and plus "bioinformatics".
3. **bioRxiv** — this is the one that caught the first collision; a preprint
   with no package yet will not appear on PyPI.
4. **Bioconda and Bioconductor** — a tool can exist there and nowhere else.
5. **Google Scholar** — an older tool may have a paper and a dead website.
6. **The console script** — `which <name>` on a machine with common tools
   installed.

Prefer a coined or obscure word. Every evocative single English word in the
navigation, surveying and verification space is taken; two picks drawn from that
well collided in succession.

## Afterwards

```console
pip install -e ".[dev]"
pytest
cd examples/minimal && <name> init && <name> load --config organism.toml
<name> export --assembly-id 1 --out /tmp/rt.gff3 && head -3 /tmp/rt.gff3
```

The last line is the one that catches a class-3 miss: the header must read
`##<name>-export`. A passing test suite does not prove this.
