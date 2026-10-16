"""The benchmark's common record of "what a genome annotation holds".

Every loader under test is reduced to the same small SQLite table, and so is the
input file itself, so that one comparison serves them all. This module defines
that table and the few rules needed to fill it the same way from any source.

It uses the Python standard library only and imports nothing from cgload or from
any other loader under test. The reference count built from it (``truth.py``)
must not share code with the programs it judges; that is the property the whole
comparison rests on.

A *feature* is the set of rows that share one ``ID`` on one file (GFF3's rule for
a feature written on several rows), or a single row that has no ``ID``. Its
fields:

* ``key``: the ID, or ``anon:<digest>`` for a row without one (see ``anon_key``).
* ``seqid``, ``type``: columns 1 and 3, as written.
* ``strand``: every strand its rows give, sorted and comma-joined.
* ``start``, ``end``: the smallest and largest coordinate over all its rows.
* ``segs``: every row's ``start-end``, sorted, comma-joined.
* ``parents``: the IDs named by ``Parent``, sorted and comma-joined.
* ``attrs``: a digest of every other attribute: keys sorted, values sorted,
  percent-escapes decoded, rows of one feature merged.
* ``records``: how many separate records a loader keeps for the feature. One is
  normal; more means the loader split one feature into several.
"""

from __future__ import annotations

import gzip
import hashlib
import sqlite3
from collections.abc import Iterable, Iterator
from pathlib import Path
from urllib.parse import unquote

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS sequence (
    seqid TEXT PRIMARY KEY, length INTEGER, circular INTEGER, genome TEXT
);
CREATE TABLE IF NOT EXISTS feature (
    key TEXT PRIMARY KEY,
    seqid TEXT, type TEXT, strand TEXT,
    start INTEGER, "end" INTEGER,
    nseg INTEGER, segs TEXT,
    parents TEXT,
    attrs TEXT,
    records INTEGER,
    flags TEXT
);
CREATE TABLE IF NOT EXISTS part (
    key TEXT, part TEXT, segs TEXT, PRIMARY KEY (key, part)
);
"""

#: The `genome=` values of a region row that mean an organelle. Our own list,
#: written from the INSDC /organelle vocabulary, not taken from any loader.
ORGANELLE_GENOMES = frozenset({
    "mitochondrion", "chloroplast", "plastid", "apicoplast", "chromatophore",
    "chromoplast", "cyanelle", "hydrogenosome", "kinetoplast", "leucoplast",
    "mitosome", "nucleomorph", "proplastid",
})

#: Hard cases, in the order the report lists them. Each name is a flag that
#: ``truth.py`` sets on a feature; the protocol defines each one.
HARD_CASES = (
    "multi_row",            # written on two or more rows sharing one ID
    "overlapping_rows",     # two of its rows overlap (programmed frameshifts)
    "trans_spliced",        # a row carries exception=trans-splicing
    "strand_unknown",       # strand '?'
    "crosses_origin",       # on a circular sequence, ends past its length
    "past_end",             # wholly or partly past the end of its sequence otherwise
    "inverted",             # a row with start > end
    "escapes_parent",       # not inside the span of a parent it names
    "exon_under_pseudogene",  # an exon whose parent is a pseudogene
    "multi_parent",         # names two or more parents
    "dangling_parent",      # names a parent the file does not contain
    "no_id",                # a row with no ID
    "multi_valued",         # an attribute with two or more values
    "rows_disagree",        # its rows carry different attributes
    "organelle",            # on a mitochondrial or plastid sequence
)

_SEP_VALUE = "\x1f"
_SEP_KEY = "\x1e"


def open_text(path: Path):
    """Open plain or gzipped text, by magic bytes."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def parse_attributes(column: str) -> dict[str, list[str]]:
    """Column 9 as key -> values, escapes decoded after splitting on commas."""
    found: dict[str, list[str]] = {}
    if column in ("", "."):
        return found
    for chunk in column.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        key, sep, value = chunk.partition("=")
        if not sep:
            continue
        key = unquote(key.strip())
        found.setdefault(key, []).extend(unquote(v) for v in value.split(","))
    return found


def merge_attributes(target: dict[str, set[str]], attributes: dict[str, list[str]]) -> None:
    for key, values in attributes.items():
        if key in ("ID", "Parent"):
            continue
        target.setdefault(key, set()).update(values)


def attrs_digest(attributes: dict[str, Iterable[str]]) -> str:
    """One digest for a feature's attributes: keys sorted, values sorted."""
    canonical = _SEP_KEY.join(
        key + "=" + _SEP_VALUE.join(sorted(set(values)))
        for key, values in sorted(attributes.items())
        if key not in ("ID", "Parent")
    )
    return hashlib.sha1(canonical.encode("utf-8"), usedforsecurity=False).hexdigest()


def anon_key(seqid: str, ftype: str, start: int, end: int, strand: str, parents: str) -> str:
    """A stable key for a row with no ID, from what the row says."""
    text = "\t".join((seqid, ftype, str(start), str(end), strand, parents))
    return "anon:" + hashlib.sha1(text.encode("utf-8"), usedforsecurity=False).hexdigest()[:24]


def segs_text(segments: Iterable[tuple[int, int]]) -> str:
    return ",".join(f"{start}-{end}" for start, end in sorted(segments))


def parse_segs(text: str) -> list[tuple[int, int]]:
    if not text:
        return []
    out = []
    for piece in text.split(","):
        start, _, end = piece.partition("-")
        out.append((int(start), int(end)))
    return out


def iter_rows(path: Path) -> Iterator[tuple[int, list[str]]]:
    """(line number, columns) for every GFF3 data row before any ##FASTA."""
    with open_text(path) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith("##FASTA"):
                return
            if line.startswith("#") or not line.strip():
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) == 9:
                yield number, columns


def create(path: Path) -> sqlite3.Connection:
    """A new, empty record database at ``path`` (an old one is replaced)."""
    if path.exists():
        path.unlink()
    connection = sqlite3.connect(path)
    connection.executescript(SCHEMA)
    connection.execute("PRAGMA journal_mode = OFF")
    connection.execute("PRAGMA synchronous = OFF")
    return connection


def write_features(connection: sqlite3.Connection, rows: Iterable[tuple]) -> int:
    """Insert (key, seqid, type, strand, start, end, nseg, segs, parents, attrs,
    records, flags) tuples; return how many."""
    count = 0
    batch = []
    for row in rows:
        batch.append(row)
        if len(batch) >= 50_000:
            connection.executemany("INSERT INTO feature VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", batch)
            count += len(batch)
            batch.clear()
    if batch:
        connection.executemany("INSERT INTO feature VALUES (?,?,?,?,?,?,?,?,?,?,?,?)", batch)
        count += len(batch)
    return count


def set_meta(connection: sqlite3.Connection, **values) -> None:
    connection.executemany(
        "INSERT OR REPLACE INTO meta VALUES (?, ?)",
        [(key, str(value)) for key, value in values.items()],
    )


def locked_append(path: Path, header: str, lines: Iterable[str]) -> None:
    """Append lines to a shared table, writing the header first if the table is
    new. Runners work in parallel, so the file is locked while it is written."""
    import fcntl

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            if handle.tell() == 0:
                handle.write(header.rstrip("\n") + "\n")
            for line in lines:
                handle.write(line.rstrip("\n") + "\n")
            handle.flush()
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def unique_anon(key: str, seen: dict[str, int]) -> str:
    """Number identical rows without an ID, so each keeps a key of its own."""
    count = seen.get(key, 0) + 1
    seen[key] = count
    return key if count == 1 else f"{key}#{count}"
