"""Milestone 0 tests.

There is no functionality yet, so these assert the two things that can be
wrong at this stage: that the package installs and its entry point runs, and
that an unbuilt command fails loudly instead of appearing to succeed.
"""

import subprocess
import sys

import click
import pytest
from cgload import __version__
from cgload.cli.main import PENDING, PendingCommand, cli
from click.testing import CliRunner


@pytest.fixture
def run():
    return CliRunner()


def test_help_exits_zero(run):
    result = run.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "cgload reconciles" in result.output


def test_version_flag_reports_package_version(run):
    result = run.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help_lists_every_v1_command(run):
    """The published surface must be visible from the first release.

    If a command is added to the CLI without a roadmap entry, or a roadmap
    entry is deleted without the command being implemented, this fails.
    """
    expected = {
        "init", "load", "verify", "explain",
        "stats", "query", "export", "remove",
        "index", "serve",
    }
    result = run.invoke(cli, ["--help"])
    listed = set(cli.commands)
    assert listed == expected
    for name in expected:
        assert name in result.output


@pytest.fixture
def pending():
    """A group holding one stub. Every published command is built now, so the
    stub mechanism (P1) is exercised on a command made for the test."""
    group = click.Group()
    group.add_command(PendingCommand("someday", 12, "not built yet"))
    return group


def test_every_published_command_is_built():
    assert PENDING == {}


def test_unimplemented_commands_exit_nonzero(run, pending):
    """P1: a command that does not work must not return success."""
    result = run.invoke(pending, ["someday"])
    assert result.exit_code == 2


def test_unimplemented_commands_name_their_milestone(run, pending):
    result = run.invoke(pending, ["someday"])
    assert "not implemented" in result.output
    assert "milestone 12" in result.output


def test_unimplemented_commands_ignore_arguments(run, pending):
    """Stubs must not mask a real failure as an argument-parsing error."""
    result = run.invoke(pending, ["someday", "--config", "nonexistent.toml", "--force"])
    assert result.exit_code == 2
    assert "not implemented" in result.output


def test_unknown_command_is_rejected(run):
    result = run.invoke(cli, ["definitely-not-a-command"])
    assert result.exit_code != 0


def test_module_entry_point_matches_script(run):
    """`python -m cgload` must behave like the installed `cgload` script."""
    completed = subprocess.run(
        [sys.executable, "-m", "cgload", "--version"],
        capture_output=True, text=True, check=True,
    )
    assert __version__ in completed.stdout


def test_version_is_pep440_parsable():
    from importlib.metadata import version as installed_version

    assert installed_version("cgload") == __version__


#: Flags of other programs that cgload's messages legitimately mention.
OTHER_PROGRAMS_FLAGS = {
    "--csi",  # tabix --csi, for a sequence too long for a .tbi index (D-082)
}


def test_every_flag_a_message_names_exists():
    """A message must not send the user to a flag that does not exist.

    Two did: a detection error said "Pass --dialect", a flag no command has
    (the choice lives in the config file), and an earlier one pointed at
    `cgload remap`. This reads every string literal in the package, f-strings
    included, and checks each `--flag` in it against the options the commands
    actually declare.
    """
    import ast
    import re
    from pathlib import Path

    import cgload

    declared = {"--help", "--version"}
    for command in cli.commands.values():
        for param in command.params:
            declared.update(opt for opt in param.opts + param.secondary_opts
                            if opt.startswith("--"))

    named: dict[str, str] = {}
    for path in Path(cgload.__file__).parent.rglob("*.py"):
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                for flag in re.findall(r"(?<![\w-])--[a-z][a-z-]*[a-z]", node.value):
                    named.setdefault(flag, path.name)

    unknown = {flag: where for flag, where in named.items()
               if flag not in declared | OTHER_PROGRAMS_FLAGS}
    assert not unknown, f"messages name flags no command declares: {unknown}"


def test_serve_says_where_the_portal_opens():
    """Every interface by default (D-115), with the addresses to use; a refused
    connection from another computer would otherwise say nothing."""
    from cgload.cli.commands.serve import where_it_listens

    lines = where_it_listens("0.0.0.0", 8000, "", ["10.0.0.239"])
    assert "http://127.0.0.1:8000/" in lines[0]
    assert "or from another computer:   http://10.0.0.239:8000/" in lines[1]
    assert any("--host 127.0.0.1" in line for line in lines)
    assert any("HTTPS" in line for line in lines) and lines[-1].endswith("Ctrl-C.")
    unknown = where_it_listens("0.0.0.0", 8080, "/genomes", [])
    assert "http://<this computer's address>:8080/genomes" in unknown[1]
    private = where_it_listens("127.0.0.1", 8000)
    assert "http://127.0.0.1:8000/" in private[0] and "leave out --host" in private[1]
    assert "http://10.0.0.239:8000/" in where_it_listens("10.0.0.239", 8000)[0]
    assert "http://[fe80::1]:8000/" in where_it_listens("fe80::1", 8000)[0]
    assert "http://[::1]:8000/" in where_it_listens("::1", 8000)[0]


def test_the_addresses_offered_are_not_loopback():
    import ipaddress

    from cgload.cli.commands.serve import network_addresses

    for address in network_addresses():
        assert not ipaddress.ip_address(address).is_loopback


def _listening() -> tuple:
    import socket

    taken = socket.socket()
    taken.bind(("127.0.0.1", 0))
    taken.listen(1)
    return taken, taken.getsockname()[1]


def test_a_taken_default_port_moves_to_the_next_free_one():
    pytest.importorskip("waitress")
    from cgload.cli.commands.serve import open_server

    taken, port = _listening()
    try:
        server, used = open_server(lambda environ, start: [], "127.0.0.1", port,
                                   port_given=False, threads=1, url_prefix="")
        server.close()
        assert used > port
    finally:
        taken.close()


def test_a_taken_port_that_was_asked_for_is_refused_in_words():
    pytest.importorskip("waitress")
    from cgload.cli.commands.serve import PortalStartError, open_server

    taken, port = _listening()
    try:
        with pytest.raises(PortalStartError) as refused:
            open_server(lambda environ, start: [], "127.0.0.1", port, port_given=True,
                        threads=1, url_prefix="")
    finally:
        taken.close()
    assert f"port {port} is already in use" in str(refused.value)
    assert f"--port {port + 1}" in str(refused.value)


def test_other_start_up_errors_are_one_sentence():
    import errno
    import socket

    from cgload.cli.commands.serve import why_not

    assert "administrator rights" in why_not(OSError(errno.EACCES, "denied"), "0.0.0.0", 80)
    assert "no address 10.9.9.9" in why_not(OSError(errno.EADDRNOTAVAIL, "x"), "10.9.9.9", 8000)
    assert "can look up" in why_not(socket.gaierror(-2, "unknown"), "nowhere.invalid", 8000)


def test_serve_refuses_a_taken_port_without_a_traceback(tmp_path):
    """The whole command, up to the port: a clean message and exit status 1."""
    pytest.importorskip("flask")
    from cgload.cli.main import cli
    from click.testing import CliRunner

    db = tmp_path / "cgload.db"
    runner = CliRunner()
    assert runner.invoke(cli, ["init", "--db-url", f"sqlite:///{db}"]).exit_code == 0
    taken, port = _listening()
    try:
        result = runner.invoke(cli, ["serve", "--db-url", f"sqlite:///{db}", "--host",
                                     "127.0.0.1", "--port", str(port)])
    finally:
        taken.close()
    assert result.exit_code == 1
    assert f"port {port} is already in use" in result.output
    assert "Traceback" not in result.output
