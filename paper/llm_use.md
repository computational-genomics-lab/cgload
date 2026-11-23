# Additional file 5: use of large language models

Springer Nature's policy, which BMC Bioinformatics follows, does not accept a
large language model as an author and asks that its use be documented in the
Methods section of the paper. This note is the detailed record behind the
section "Use of large language models" in the Methods of the paper. Its counts
come from the repository's git history and can be checked there.

## Which model, and how it was used

Claude (Anthropic), from 20 August to 23 November 2026, mainly the model
Claude Opus 5.5. The model worked on a copy of the repository and could not
reach the lab server: every benchmark run was started by the authors on that
server, and the authors sent the results back to the model. Every change the
model made was committed under the name of the first author, and reached the
public repository only when the first author pushed it.

From 14 October 2026 each commit that the model wrote carries a trailer naming
it (`Co-Authored-By: Claude …`). Up to tag `v0.2.2`, the version the paper
describes, 77 of the 90 commits carry one, from 14 October to 19 November 2026;
74 of them name Claude Opus 5.5, two Claude Fable 5.2 and one Claude Opus 5.
Of the 13 commits without a trailer, the 11 from 20 August to 8 October 2026,
among them the census scripts in `survey/`, were also made with Claude's help,
before the trailer was used; the model versions used for them were not
recorded. The other two are small edits by the authors on 10 and 12 November
2026 (the full text of the MIT licence, and the removal of internal review
notes). The commits that the model wrote after `v0.2.2` carry the trailer too.

## What the model did

- **cgload's code and tests.** It wrote and changed code in `src/` and `tests/`,
  including fixes for the problems the benchmark exposed (decisions D-095 to
  D-110 in `docs/decisions.md`), and wrote the documentation in `docs/`, the
  README and the changelog.
- **The census and the benchmark.** It helped write the census scripts
  (`survey/`). It drafted the benchmark protocol (`benchmark/PROTOCOL.md`) with
  the authors, wrote the benchmark scripts and their tests (`benchmark/`), wrote
  the shell scripts the authors ran on the server, and analysed each run's
  results. It found and corrected errors in its own scripts; each is recorded
  with its date in the protocol, section 12.
- **Figures and tables.** It wrote `paper/figures.py`, which draws every figure
  and table of the paper from the committed results. No figure or number was
  made by asking the model for it directly, and the model created or altered no
  image. It also wrote `paper/build.py`.
- **The manuscript.** It drafted the text of the manuscript, from the results
  and under the authors' direction. The authors read the draft against the code
  and the results, revised it, and wrote the declarations.
- **Checking.** It checked every number in the manuscript against the results
  with a script, and searched the literature for related tools. On 11 and 23
  November 2026 the authors asked it to review the draft, and then the
  submission files, for errors, against the journal's guidelines and against
  the repository and its history. It corrected the inconsistencies it found in
  the text, the cover letter, the additional files and `paper/figures.py`; the
  authors reviewed each correction.

## ChatGPT, as reviewer

On 19 November 2026 the authors gave the submission files to ChatGPT (OpenAI)
and asked it for errors and limitations. Every point was checked against the
code and the results before anything changed. The review led to corrections in
the text, to the row-by-row comparison of GFF3 columns 2, 6 and 8
(`benchmark/columns.py`), and through it to two fixes in cgload 0.2.2
(decisions D-109 and D-110).

## Where the model's output is kept

In the public repository, https://github.com/computational-genomics-lab/cgload:
the code it wrote is in the repository's history; its draft of the manuscript
is `paper/manuscript.md`, with every earlier version in the history; and the
figures and tables are the output of `paper/figures.py`.

## How the authors checked the work

The authors ran every benchmark step on their own server, from the scripts in
the repository, and every result in the paper comes from those runs
(`paper/data/benchmark-run-2.2/` and `paper/data/columns-v0.2.2/`). The test
suite (575 tests, and 25 more for the benchmark) runs in continuous integration
on every push, on Linux, macOS and Windows. The authors reviewed the code, the
results and the text, and are responsible for them.
