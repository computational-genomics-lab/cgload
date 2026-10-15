# Releasing cgload

A release is a GitHub release. Publishing it runs
`.github/workflows/release.yml`, which builds the package, checks it, and
uploads it to PyPI. No password or API token is stored anywhere: PyPI trusts
that one workflow in this repository ("trusted publishing").

## Once: let PyPI trust the workflow

Someone who will own the `cgload` project on PyPI does this.

1. Make an account at <https://pypi.org> and turn on two-factor
   authentication (PyPI requires it).
2. Go to *Your account → Publishing → Add a new pending publisher*, choose
   GitHub, and enter:
   - PyPI project name: `cgload`
   - Owner: `computational-genomics-lab`
   - Repository name: `cgload`
   - Workflow name: `release.yml`
   - Environment name: `pypi`
3. Save. The first release creates the project, and that account owns it.

To rehearse on TestPyPI first, do the same at <https://test.pypi.org> with the
workflow `testpypi.yml` and the environment `testpypi`, then run that workflow
by hand from the Actions tab.

## Each release

1. Set `__version__` in `src/cgload/__init__.py`, move the *Unreleased*
   notes in `CHANGELOG.md` under the new version and date, commit to `main`
   and push. For the first release, also change the README's install section
   to `pip install cgload`. Wait until CI passes on that commit.
2. On GitHub, *Releases → Draft a new release*. Create the tag `vX.Y.Z` on
   `main`, matching `__version__` exactly; the workflow stops if they differ.
   Paste the changelog section as the notes, and publish.
3. Check the *Release* workflow in the Actions tab. When it finishes,
   `pip install cgload==X.Y.Z` works.

A version on PyPI cannot be replaced, only followed by a newer one. If
something is wrong, fix it and release `X.Y.Z+1`.

## Without the workflow

From a clean checkout of the tag, with a PyPI API token kept out of the
repository and out of shell history:

```console
python -m pip install build twine
python -m build
twine check dist/*
twine upload dist/*
```

## Afterwards

- Archive the release on Zenodo for a DOI: upload the source archive from the
  GitHub release page to a new version of the project's Zenodo record. The
  GitHub–Zenodo integration needs an organisation admin, so this step is
  manual.
- A bioconda recipe can then point at the PyPI source archive. It would let
  users whose system Python is older than 3.11 install cgload with one
  command.
