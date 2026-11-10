# Supplementary note: use of a large language model

The journal follows the ISCB policy for acceptable use of large language models
(LLMs), which asks that any use be declared, with enough detail, in the
supplementary materials. This note is that declaration. Its counts come from
the repository's git history and can be checked there.

<!-- AUTHORS: the facts below are from the git history. The two places marked
"to be completed" need your own statement before submission. -->

## Which model, and how it was used

Claude (Anthropic), mostly the model Claude Opus 5.5, used through Claude
Code in two working sessions linked from the commit messages. The model worked
in a cloud workspace on a copy of the repository. It could not reach the lab
server: every benchmark run was started by the authors on that server, and the
results were sent back to the model by the authors. Every change the model made
was committed under the name of the first author, with a trailer naming the
model (`Co-Authored-By: Claude …`), and reached the public repository only when
the first author pushed it.

Up to tag v0.2.0, 70 of the 81 commits in the repository carry such a
trailer, from 14 October to 10 November 2026. The other 11, from 20 August to 8 October 2026, carry none.
<!-- AUTHORS (to be completed): say whether an LLM helped with those 11. -->

## What the model did

These uses fall under the policy's acceptable uses: help with writing code,
documenting code, and finding inconsistencies.

- **cgload's code and tests.** It wrote and changed code in `src/` and `tests/`,
  including fixes for the problems the benchmark exposed (decisions D-095 to
  D-108 in `docs/decisions.md`), and wrote the documentation in `docs/`, the
  README and the changelog.
- **The benchmark.** It drafted the protocol (`benchmark/PROTOCOL.md`) with the
  authors, wrote the benchmark scripts and their tests (`benchmark/`), wrote the
  shell scripts the authors ran on the server, and analysed each run's results.
  It found and corrected errors in its own scripts; each is recorded with its
  date in the protocol, section 12.
- **Figures and tables.** It wrote `paper/figures.py`, which draws every figure
  and table of the paper from the committed results; no figure or number was
  made by prompting the model directly. It also wrote `paper/build.py` and
  `paper/deposit_zenodo.py`.
- **Checking.** It checked every number in the manuscript against the results
  with a script, and searched the literature for related tools.
- **A working draft of the manuscript.** It drafted a version of the text,
  which held the facts, numbers and structure for the authors to work from.
  <!-- AUTHORS (to be completed): the policy does not allow an LLM to draft the
  paper. State here that the submitted text was written by the authors, once
  that is true, and remove this comment. -->

## How the authors checked the work

The authors ran every benchmark step on their own server, from the scripts in
the repository, and every result in the paper comes from those runs
(`paper/data/benchmark-run-2.2/`). The test suite (573 tests, and 20 more for
the benchmark) runs in continuous integration on every push. The authors reviewed the code and the results
and are responsible for them.
