"""Reading the assembly, and building its index.

pyfaidx writes a samtools-compatible ``.fai`` (D-004), which is what gives
O(1) access to any region and what a genome browser reads directly.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from pyfaidx import Fasta, FastaIndexingError


class AssemblyReadError(RuntimeError):
    """The assembly is unreadable, empty, or has duplicate sequence names."""


@dataclass(frozen=True)
class Region:
    """One scaffold, contig or chromosome."""

    source_id: str
    length: int


def index_and_read(path: Path) -> tuple[Path, list[Region]]:
    """Write ``<path>.fai`` beside the assembly and return its regions.

    Errors rather than warns on: an empty file, a duplicate sequence name, or
    a zero-length record. Each of those loads "successfully" in a permissive
    reader and produces a database that is quietly missing or double-counting
    a scaffold.
    """
    path = Path(path)
    try:
        with Fasta(str(path), duplicate_action="stop", rebuild=True) as handle:
            regions = [Region(source_id=name, length=len(record))
                       for name, record in handle.items()]
    except FastaIndexingError as exc:
        raise AssemblyReadError(f"{path.name}: {exc}") from exc
    except ValueError as exc:  # pyfaidx raises this for duplicate names
        raise AssemblyReadError(f"{path.name}: {exc}") from exc

    if not regions:
        raise AssemblyReadError(f"{path.name}: contains no sequences")

    empty = [r.source_id for r in regions if r.length == 0]
    if empty:
        raise AssemblyReadError(
            f"{path.name}: {len(empty)} zero-length sequence(s), first is {empty[0]!r}"
        )

    names = [r.source_id for r in regions]
    if len(set(names)) != len(names):
        seen: set[str] = set()
        dupe = next(n for n in names if n in seen or seen.add(n))
        raise AssemblyReadError(f"{path.name}: duplicate sequence name {dupe!r}")

    index = path.with_suffix(path.suffix + ".fai")
    if not index.exists():
        raise AssemblyReadError(f"{path.name}: index was not written")
    return index, regions


def iter_regions(path: Path) -> Iterator[Region]:
    _, regions = index_and_read(path)
    yield from regions
