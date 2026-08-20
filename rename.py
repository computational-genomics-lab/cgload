#!/usr/bin/env python3
"""Rename the package, including the strings that fail silently if missed.

    python3 rename.py loxodrome

Three classes of occurrence, in increasing order of how quietly they break:

1. **Imports and paths.** `src/<name>/`, `from <name>.x import y`,
   `[project] name`, the console script. A miss here is a loud ImportError.
2. **Prose.** README, docs, docstrings. A miss is embarrassing, not dangerous.
3. **On-disk format strings.** `##<name>-export`, `##<name>-vocabulary`,
   `<name>_linkage_method`, the dialect profile name. **A miss here leaves the
   test suite passing while exported files stop identifying themselves**, or
   while linkage provenance is silently discarded on re-import.

Class 3 is why this is a script rather than a sed one-liner. `LINKAGE_ATTRIBUTE`
is declared as a literal in three modules that deliberately do not import each
other, so a partial rename leaves them disagreeing -- which is exactly the bug
that erased linkage provenance once before.

The script verifies afterwards that no occurrence of the old name survives and
that the three constants agree.
"""

from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

OLD = "cgload"

#: NOTE: a rename rewrites this file too, so OLD above is always the *current*
#: name once the script has run. That is correct -- the next rename starts from
#: here -- but it means the constant is not a record of history. See
#: docs/RENAMING.md for which names have been used and why they were rejected.

TEXT_SUFFIXES = {
    ".py", ".md", ".toml", ".yml", ".yaml", ".cfg", ".txt", ".ini",
    ".gff3", ".gbff", ".json", ".in",
}

SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", ".ruff_cache", "build", "dist"}

#: Extensionless files that still carry the name. `.gitignore` was missed by an
#: earlier version of this script -- it has no suffix, so a suffix allowlist
#: skips it silently.
EXTENSIONLESS = {"LICENSE", ".gitignore", ".gitattributes", "Makefile", "Dockerfile"}


def substitutions(new: str) -> list[tuple[str, str]]:
    """Ordered longest-first so no replacement eats another's prefix."""
    return [
        (f"##{OLD}-export", f"##{new}-export"),
        (f"##{OLD}-vocabulary", f"##{new}-vocabulary"),
        (f"{OLD}_linkage_method", f"{new}_linkage_method"),
        (OLD.upper(), new.upper()),
        (OLD.capitalize(), new.capitalize()),
        (OLD, new),
    ]


def rewrite(path: Path, pairs: list[tuple[str, str]]) -> bool:
    try:
        text = path.read_text(encoding="utf-8")
    except (UnicodeDecodeError, OSError):
        return False
    original = text
    for old, new in pairs:
        text = text.replace(old, new)
    if text != original:
        path.write_text(text, encoding="utf-8")
        return True
    return False


def main(new: str) -> int:
    if not re.fullmatch(r"[a-z][a-z0-9_]{2,30}", new):
        print(f"'{new}' is not a usable package name: lowercase, starts with a "
              "letter, letters/digits/underscore only.")
        return 2

    root = Path(__file__).resolve().parent
    pairs = substitutions(new)

    # 1. the package directory, before rewriting anything that names it
    package = root / "src" / OLD
    if package.is_dir():
        shutil.move(str(package), str(root / "src" / new))
        print(f"  moved src/{OLD}/ -> src/{new}/")

    # 2. every text file
    changed = 0
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        if any(part in SKIP_DIRS for part in path.parts):
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name not in EXTENSIONLESS:
            continue
        if path.name == Path(__file__).name:
            continue
        if rewrite(path, pairs):
            changed += 1
    print(f"  rewrote {changed} file(s)")

    # 3. verification -- the part that matters
    print("\nverifying:")
    survivors: list[str] = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or any(p in SKIP_DIRS for p in path.parts):
            continue
        if path.name == Path(__file__).name:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if re.search(OLD, text, re.IGNORECASE):
            for number, line in enumerate(text.splitlines(), 1):
                if re.search(OLD, line, re.IGNORECASE):
                    survivors.append(f"    {path.relative_to(root)}:{number}: {line.strip()[:90]}")

    if survivors:
        print(f"  {len(survivors)} surviving reference(s) to '{OLD}':")
        print("\n".join(survivors[:20]))
        if len(survivors) > 20:
            print(f"    ... and {len(survivors) - 20} more")
    else:
        print(f"  no surviving reference to '{OLD}'")

    # 4. the three constants must agree across modules that do not import
    #    each other. Checked by reading, not importing, so this works before
    #    the package is installed.
    print("\n  cross-module constants:")
    found: dict[str, set[str]] = {}
    for name in ("LINKAGE_ATTRIBUTE", "EXPORT_PRAGMA"):
        values: set[str] = set()
        for path in (root / "src" / new).rglob("*.py"):
            for match in re.finditer(rf'^{name}\s*=\s*"([^"]+)"', path.read_text(), re.M):
                values.add(match.group(1))
        found[name] = values
        state = "agree" if len(values) <= 1 else "DISAGREE"
        print(f"    {name:20s} {state}: {sorted(values)}")

    disagreement = any(len(v) > 1 for v in found.values())
    if disagreement:
        print("\n  A disagreement here is the silent failure this script exists "
              "to prevent: exports stop identifying themselves, or linkage "
              "provenance is discarded on re-import, while every test passes.")

    print("\nnext:")
    print("  pip install -e '.[dev]' && pytest")
    print(f"  cd examples/minimal && {new} init && {new} load --config organism.toml")
    print(f"  {new} export --assembly-id 1 --out /tmp/rt.gff3 && head -3 /tmp/rt.gff3")
    print(f"    ^ the export header must read ##{new}-export")

    return 1 if (survivors or disagreement) else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
