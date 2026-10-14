"""BGZF compression and a tabix (.tbi) index, written without external programs.

The README promises an annotation in the data directory that JBrowse 2, IGV and
tabix can read straight from disk: bgzip-compressed, tabix-indexed GFF3. Nothing
wrote it (AUDIT-10). Calling `bgzip` and `tabix` would add two programs to
install and a subprocess per load, which is where D-005 says genome tools pick
up their attack surface; pysam would add a compiled dependency with no Windows
wheels. Both formats are short published specifications (SAMv1 section 4.1 for
BGZF, the tabix section of the SAM/BAM index specification for .tbi), so they
are written here with ``zlib`` and ``struct`` (D-082).

Only what the data directory needs is implemented: writing, the GFF preset,
coordinates below 2**29 (the TBI format's limit; a longer sequence needs a CSI
index, which is not written, and the caller says so).
"""

from __future__ import annotations

import struct
import zlib
from pathlib import Path
from typing import BinaryIO

#: Uncompressed bytes per BGZF block. htslib uses 0xff00 so that even
#: incompressible input fits in the 64 KiB a block may occupy.
BLOCK_INPUT = 0xFF00

#: The empty block that ends every BGZF file (SAMv1 4.1.2).
EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

#: TBI's binning scheme covers positions below 2**29.
TBI_MAX_POSITION = 1 << 29

_MIN_SHIFT = 14  # 16 kb linear-index windows
_PSEUDO_BIN = 37450  # per-reference summary bin written by htslib


class IndexLimitError(ValueError):
    """A coordinate is beyond what a .tbi index can address."""


def _block(data: bytes, level: int) -> bytes:
    compressor = zlib.compressobj(level, zlib.DEFLATED, -15)
    payload = compressor.compress(data) + compressor.flush()
    size = 18 + len(payload) + 8
    header = (
        b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff"
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<HH", 2, size - 1)
    )
    trailer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + payload + trailer


class BgzfWriter:
    """Write BGZF, and report the virtual offset of every byte written.

    A virtual offset is the compressed offset of a block's start shifted left
    16 bits, plus the offset within that block's uncompressed data. A tabix
    index is a map from coordinates to virtual offsets.
    """

    def __init__(self, handle: BinaryIO, *, level: int = 6) -> None:
        self._handle = handle
        self._level = level
        self._buffer = bytearray()
        self._block_start = 0

    def tell(self) -> int:
        return (self._block_start << 16) | len(self._buffer)

    def write(self, data: bytes) -> None:
        view = memoryview(data)
        while view:
            room = BLOCK_INPUT - len(self._buffer)
            self._buffer += view[:room]
            view = view[room:]
            if len(self._buffer) == BLOCK_INPUT:
                self._flush()

    def _flush(self) -> None:
        if not self._buffer:
            return
        block = _block(bytes(self._buffer), self._level)
        self._handle.write(block)
        self._block_start += len(block)
        self._buffer.clear()

    def close(self) -> None:
        self._flush()
        self._handle.write(EOF_BLOCK)


def reg2bin(begin: int, end: int) -> int:
    """The smallest bin holding [begin, end), 0-based half-open (SAMv1 5.3)."""
    end -= 1
    if begin >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (begin >> 14)
    if begin >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (begin >> 17)
    if begin >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (begin >> 20)
    if begin >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (begin >> 23)
    if begin >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (begin >> 26)
    return 0


class _Reference:
    def __init__(self, name: str) -> None:
        self.name = name
        self.bins: dict[int, list[list[int]]] = {}
        self.linear: list[int] = []
        self.first = -1
        self.last = 0
        self.records = 0
        self.previous_begin = -1


class TabixIndex:
    """Build a .tbi index for GFF3 records added in sorted order.

    The GFF preset of tabix: sequence in column 1, start in 4, end in 5, 1-based
    closed coordinates, lines starting with ``#`` skipped.
    """

    def __init__(self) -> None:
        self._references: list[_Reference] = []
        self._by_name: dict[str, _Reference] = {}

    def add(self, seqid: str, start: int, end: int, begin_offset: int, end_offset: int) -> None:
        """One record: 1-based closed [start, end], and its virtual offsets."""
        if end > TBI_MAX_POSITION:
            raise IndexLimitError(
                f"{seqid}:{end:,} is beyond the {TBI_MAX_POSITION:,} positions a .tbi "
                f"index can address"
            )
        reference = self._by_name.get(seqid)
        if reference is None:
            reference = _Reference(seqid)
            self._references.append(reference)
            self._by_name[seqid] = reference
        elif reference is not self._references[-1]:
            raise ValueError(f"records for {seqid} are not contiguous")

        begin = start - 1
        if begin < reference.previous_begin:
            raise ValueError(f"records for {seqid} are not sorted by start")
        reference.previous_begin = begin

        chunks = reference.bins.setdefault(reg2bin(begin, end), [])
        if chunks and chunks[-1][1] == begin_offset:
            chunks[-1][1] = end_offset
        else:
            chunks.append([begin_offset, end_offset])

        first_window = begin >> _MIN_SHIFT
        last_window = (end - 1) >> _MIN_SHIFT
        if len(reference.linear) <= last_window:
            reference.linear.extend([0] * (last_window + 1 - len(reference.linear)))
        for window in range(first_window, last_window + 1):
            if reference.linear[window] == 0:
                reference.linear[window] = begin_offset

        if reference.first < 0:
            reference.first = begin_offset
        reference.last = end_offset
        reference.records += 1

    def write(self, path: Path) -> None:
        names = b"".join(ref.name.encode() + b"\x00" for ref in self._references)
        body = bytearray(b"TBI\x01")
        body += struct.pack("<8i", len(self._references), 0, 1, 4, 5, ord("#"), 0, len(names))
        body += names
        for ref in self._references:
            linear = list(ref.linear)
            for window in range(1, len(linear)):
                if linear[window] == 0:
                    linear[window] = linear[window - 1]
            body += struct.pack("<i", len(ref.bins) + 1)
            for number in sorted(ref.bins):
                chunks = ref.bins[number]
                body += struct.pack("<Ii", number, len(chunks))
                for begin_offset, end_offset in chunks:
                    body += struct.pack("<QQ", begin_offset, end_offset)
            body += struct.pack("<Ii", _PSEUDO_BIN, 2)
            body += struct.pack("<QQQQ", max(ref.first, 0), ref.last, ref.records, 0)
            body += struct.pack("<i", len(linear))
            body += struct.pack(f"<{len(linear)}Q", *linear)
        body += struct.pack("<Q", 0)
        with Path(path).open("wb") as handle:
            writer = BgzfWriter(handle)
            writer.write(bytes(body))
            writer.close()
