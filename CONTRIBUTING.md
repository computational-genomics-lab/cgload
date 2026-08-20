# Contributing

Bug reports, dialect reports and pull requests are all welcome.

## The most useful thing you can send

**A file that cgload refuses, or loads wrongly.** Nearly every rule in this
codebase was written from the files available at the time, and nearly every one
has since been corrected by a file that behaved differently. If cgload rejects
something valid, that is a bug and the file is the bug report.

If the file is large, the useful part is small: run the pre-flight scanner and
send its output plus twenty lines around the offending row.

```console
python3 preflight_gff3.py your_annotation.gff3
```

Please include the accession if it is public, so the case can be reproduced.

## Reporting a bug

Include the command you ran, what happened, what you expected, and the output of
`cgload --version`. If a load failed, the error message names a file and a line
— include both.

## Development

```console
git clone https://github.com/computational-genomics-lab/cgload
cd cgload
pip install -e ".[dev]"
pytest
ruff check .
```

Tests must pass and lint must be clean before a pull request is reviewed.

## What a change needs

**A test that fails without it.** For a parser change, that means a fixture: a
small slice of a real file, graph-closed so no feature is missing its parent and
no parent is missing its children. Synthetic fixtures are acceptable when the
case cannot be found in real data, but say so in the file's header.

**A note in `tests/fixtures/PROVENANCE.md`** giving the accession, the region,
and what the fixture is for.

**An entry in `docs/decisions.md`** if the change alters behaviour rather than
fixing an outright mistake. The format is: what was found, where, what changed,
and the measurement that justified it. This file is the reason the codebase is
comprehensible; keeping it current is not optional.

## Adding a dialect

A new annotation producer should be a declared profile in
`src/cgload/parsers/profiles.py`, not a branch in the parser. If a dialect
cannot be expressed as a declaration, that is worth discussing in an issue
first — it usually means the profile model needs a new field rather than the
parser needing a special case.

Every profile needs a fixture from a real file produced by that tool.

## Code style

`ruff` enforces the mechanical parts. Beyond that: comments should record *why*,
particularly the evidence behind a rule. A comment saying what the code does is
usually redundant; one saying which file broke and how is not.

## Conduct

See [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
