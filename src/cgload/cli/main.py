"""The ``cgload`` command group.

Design note (P1, fail loudly): every subcommand in the published surface is
registered from the first commit. Commands whose milestone has not landed
exit non-zero with a message naming the milestone. The alternative -- adding
commands to ``--help`` as they are written -- means a user cannot tell the
difference between "this tool does not do that" and "this version does not do
that yet", and a half-built command that returns 0 is exactly the silent
wrongness this project exists to avoid.
"""

from __future__ import annotations

import click

from cgload import __version__
from cgload.cli.commands.explain import explain_command
from cgload.cli.commands.export import export_command
from cgload.cli.commands.init import init_command
from cgload.cli.commands.load import load_command
from cgload.cli.commands.stats import stats_command
from cgload.cli.commands.verify import verify_command

# Milestone numbers refer to docs/roadmap.md. Update the entry, and delete it
# from this table, in the same commit that implements the command.
PENDING: dict[str, tuple[int, str]] = {
    "query": (9, "region, annotation and KO queries"),
    "remove": (9, "transactional removal of an assembly"),
}


class PendingCommand(click.Command):
    """A command that is part of the published surface but not yet built."""

    def __init__(self, name: str, milestone: int, summary: str) -> None:
        self.milestone = milestone
        super().__init__(
            name=name,
            callback=self._fail,
            help=f"[milestone {milestone}] {summary}",
            short_help=f"(not yet built \u2014 milestone {milestone})",
            context_settings={"ignore_unknown_options": True},
            params=[click.Argument(["args"], nargs=-1, type=click.UNPROCESSED)],
        )

    def _fail(self, args: tuple[str, ...]) -> None:  # noqa: ARG002
        click.echo(
            f"cgload {self.name}: not implemented in {__version__}.\n"
            f"Scheduled for milestone {self.milestone}; see docs/roadmap.md.",
            err=True,
        )
        raise click.exceptions.Exit(2)


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, "-V", "--version", prog_name="cgload")
def cli() -> None:
    """Load genome annotations into a queryable relational database.

    cgload reconciles the annotation dialects that different pipelines
    produce, stores features in a documented relational schema, keeps
    sequence in indexed FASTA on disk, and verifies after every load that
    nothing was silently dropped.
    """


cli.add_command(init_command)
cli.add_command(load_command)
cli.add_command(verify_command)
cli.add_command(explain_command)
cli.add_command(export_command)
cli.add_command(stats_command)

for _name, (_milestone, _summary) in PENDING.items():
    cli.add_command(PendingCommand(_name, _milestone, _summary))
