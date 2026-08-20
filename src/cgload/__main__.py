"""Enable ``python -m cgload``."""

from __future__ import annotations

from cgload.cli.main import cli


def main() -> None:
    cli(prog_name="cgload")


if __name__ == "__main__":
    main()
