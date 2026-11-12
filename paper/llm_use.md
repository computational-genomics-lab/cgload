# Additional file 5: use of a large language model

GigaScience allows generative AI tools when their use is disclosed, cited, and
summarised at the end of the paper, with the tool's output kept in an open
repository. This note is the detailed record behind the paper's disclosure
section. Its counts come from the repository's git history and can be checked
there.

<!-- AUTHORS: the facts below are from the git history. The two places marked
"to be completed" need your own statement before submission. -->

## Which model, and how it was used

Claude (Anthropic), mostly the model Claude Opus 5.5, used through Claude
Code in two working sessions that the commit messages link to, from 14 October
to 12 November 2026. The model worked in a cloud workspace on a copy of the
repository. It could not reach the lab server: every benchmark run was started
by the authors on that server, and the authors sent the results back to the
model. Every change the model made was committed under the name of the first
author, with a trailer naming the model (`Co-Authored-By: Claude …`), and
reached the public repository only when the first author pushed it.

Up to tag `v0.2.0`, 70 of the 81 commits in the repository carry such a
trailer, from 14 October to 10 November 2026. The other 11, from 20 August to
8 October 2026, carry none.
<!-- AUTHORS (to be completed): say whether an LLM helped with those 11. -->

## What the model did

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
  and table of the paper from the committed results. No figure or number was
  made by asking the model for it directly, and the model created or altered no
  image. It also wrote `paper/build.py`.
- **Checking.** It checked every number in the manuscript against the results
  with a script, and searched the literature for related tools.
- **The manuscript.** It drafted the text of the manuscript, from the results
  and under the authors' direction.
  <!-- AUTHORS (to be completed): say what you changed in the draft and how you
  checked it, so that this entry describes the text you submit. -->

## Where the model's output is kept

In the public repository, https://github.com/computational-genomics-lab/cgload:
the code it wrote is in the commits that name it; its draft of the manuscript
is `paper/manuscript.md`, with every earlier version in the history; and the
figures and tables are the output of `paper/figures.py`. 

## How the authors checked the work

The authors ran every benchmark step on their own server, from the scripts in
the repository, and every result in the paper comes from those runs
(`paper/data/benchmark-run-2.2/`). The test suite (573 tests, and 20 more for
the benchmark) runs in continuous integration on every push. The authors
reviewed the code, the results and the text, and are responsible for them.
