#!/usr/bin/env python3
"""Run one command and record its wall time, CPU time and peak memory.

    python3 benchmark/measure.py --json OUT.json [--timeout S] [--stdout F] [--stderr F] -- CMD...

Wall time is measured around the process; CPU time and peak resident memory
come from the kernel's accounting of the child (getrusage), so they cover the
command and anything it waited for. A command still running after --timeout
seconds is killed with its whole process group and recorded as timed out (exit
code 124). Each measurement runs in a fresh copy of this program, so one
measurement's peak cannot leak into the next.
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import signal
import subprocess
import sys
import time
from datetime import UTC, datetime
from pathlib import Path


def run(command: list[str], *, timeout: float | None, stdout: Path | None,
        stderr: Path | None) -> dict:
    out = open(stdout, "w") if stdout else None  # noqa: SIM115 - closed below
    err = open(stderr, "w") if stderr else None  # noqa: SIM115 - closed below
    started_at = datetime.now(UTC).strftime("%Y-%m-%dT%H:%M:%SZ")
    started = time.monotonic()
    timed_out = False
    try:
        try:
            process = subprocess.Popen(command, stdout=out, stderr=err, start_new_session=True)
        except OSError as error:
            # The program is missing or cannot run: record it like any other failure.
            if err:
                err.write(f"measure.py: could not start {command[0]}: {error}\n")
            code = 127
        else:
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                timed_out = True
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                code = 124
    finally:
        for handle in (out, err):
            if handle:
                handle.close()
    wall = time.monotonic() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return {
        "command": command,
        "started": started_at,
        "exit_code": code,
        "timed_out": timed_out,
        "wall_s": round(wall, 3),
        "cpu_s": round(usage.ru_utime + usage.ru_stime, 3),
        # ru_maxrss is in kilobytes on Linux.
        "peak_rss_mb": round(usage.ru_maxrss / 1024, 1),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--timeout", type=float, default=None)
    parser.add_argument("--stdout", type=Path)
    parser.add_argument("--stderr", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("no command given")
    result = run(command, timeout=args.timeout, stdout=args.stdout, stderr=args.stderr)
    args.json.write_text(json.dumps(result) + "\n")
    return result["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
