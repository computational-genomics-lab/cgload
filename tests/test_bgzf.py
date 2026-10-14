"""The data directory's annotation: BGZF and a tabix index, checked independently.

`cgload.bgzf` writes both formats itself (D-082). These tests read them back
with code written here from the specifications, which shares nothing with the
writer, and, where pysam is installed, with htslib's own reader.
"""

from __future__ import annotations

import gzip
import random
import struct
import zlib
from pathlib import Path

import pytest
from cgload.bgzf import BgzfWriter, IndexLimitError, TabixIndex
from cgload.export import write_indexed_gff3

# ---------------------------------------------------------------------------
# an independent reader
# ---------------------------------------------------------------------------


def _blocks(raw: bytes) -> dict[int, bytes]:
    """Compressed offset of each BGZF block -> its uncompressed data."""
    blocks, offset = {}, 0
    while offset < len(raw):
        assert raw[offset : offset + 4] == b"\x1f\x8b\x08\x04", "not a BGZF block"
        size = struct.unpack("<H", raw[offset + 16 : offset + 18])[0] + 1
        blocks[offset] = zlib.decompress(raw[offset + 18 : offset + size - 8], -15)
        offset += size
    return blocks


def _read_tbi(path: Path):
    data = gzip.decompress(path.read_bytes())
    assert data[:4] == b"TBI\x01"
    n_ref, fmt, col_seq, col_beg, col_end, meta, skip, l_nm = struct.unpack("<8i", data[4:36])
    names = data[36 : 36 + l_nm].split(b"\x00")[:-1]
    offset = 36 + l_nm
    references = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack("<i", data[offset : offset + 4])
        offset += 4
        bins = {}
        for _ in range(n_bin):
            number, n_chunk = struct.unpack("<Ii", data[offset : offset + 8])
            offset += 8
            chunks = [
                struct.unpack("<QQ", data[offset + 16 * i : offset + 16 * i + 16])
                for i in range(n_chunk)
            ]
            offset += 16 * n_chunk
            bins[number] = chunks
        (n_intv,) = struct.unpack("<i", data[offset : offset + 4])
        offset += 4
        linear = list(struct.unpack(f"<{n_intv}Q", data[offset : offset + 8 * n_intv]))
        offset += 8 * n_intv
        references.append((bins, linear))
    header = (fmt, col_seq, col_beg, col_end, chr(meta), skip)
    return header, [name.decode() for name in names], references


def _reg2bins(begin: int, end: int) -> list[int]:
    end -= 1
    found = [0]
    for shift, first in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        found.extend(range(first + (begin >> shift), first + (end >> shift) + 1))
    return found


def _query(gz: Path, seqid: str, begin: int, end: int) -> list[str]:
    """Rows overlapping 0-based half-open [begin, end), found through the index."""
    _header, names, references = _read_tbi(gz.with_name(gz.name + ".tbi"))
    bins, linear = references[names.index(seqid)]
    floor = linear[min(begin >> 14, len(linear) - 1)] if linear else 0

    blocks = _blocks(gz.read_bytes())
    position, total, parts = {}, 0, []
    for offset in sorted(blocks):
        position[offset] = total
        total += len(blocks[offset])
        parts.append(blocks[offset])
    text = b"".join(parts)

    def absolute(virtual: int) -> int:
        return position[virtual >> 16] + (virtual & 0xFFFF)

    found = set()
    for number in _reg2bins(begin, end):
        for start_offset, end_offset in bins.get(number, []):
            if end_offset <= floor:
                continue
            for line in text[absolute(start_offset) : absolute(end_offset)].decode().splitlines():
                columns = line.split("\t")
                if columns[0] == seqid and int(columns[3]) - 1 < end and int(columns[4]) > begin:
                    found.add(line)
    return sorted(found)


# ---------------------------------------------------------------------------
# features shaped like the parser's output
# ---------------------------------------------------------------------------


class _Segment:
    def __init__(self, start: int, end: int) -> None:
        self.start, self.end, self.phase, self.score = start, end, None, None


class _Feature:
    def __init__(self, seqid: str, identifier: str, spans: list[tuple[int, int]]) -> None:
        self.seqid, self.source_id = seqid, identifier
        self.segments = [_Segment(a, b) for a, b in spans]
        self.source_type = self.feature_type = "gene"
        self.strand, self.source_program = "+", "test"
        self.parent_source_id = self.linkage_method = None
        self.attributes = {"Note": [identifier]}


def _random_features(seed: int = 3):
    rnd = random.Random(seed)
    regions = [("chr1", 5_000_000), ("chr2", 40_000), ("chr3", 90_000_000)]
    features = []
    for name, length in regions:
        for number in range(1500):
            start = rnd.randint(1, length - 1)
            end = min(length, start + rnd.choice([10, 900, 20_000, 400_000]))
            spans = [(start, end)]
            if end - start > 100 and rnd.random() < 0.3:
                middle = (start + end) // 2
                spans = [(start, middle - 10), (middle + 10, end)]
            features.append(_Feature(name, f"{name}-{number}", spans))
    return regions, features


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------


def test_bgzf_is_gzip_and_ends_with_the_eof_block(tmp_path):
    payload = b"".join(f"line {n}\n".encode() for n in range(40_000))
    path = tmp_path / "x.gz"
    with path.open("wb") as handle:
        writer = BgzfWriter(handle)
        writer.write(payload)
        writer.close()
    raw = path.read_bytes()
    assert gzip.decompress(raw) == payload
    assert raw.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
    assert len(_blocks(raw)) > 2


def test_virtual_offsets_point_at_the_bytes_written(tmp_path):
    path = tmp_path / "x.gz"
    marks = []
    with path.open("wb") as handle:
        writer = BgzfWriter(handle)
        for number in range(30_000):
            marks.append((writer.tell(), f"record {number}\n"))
            writer.write(marks[-1][1].encode())
        writer.close()
    blocks = _blocks(path.read_bytes())
    for offset, text in random.Random(1).sample(marks, 500):
        block, within = offset >> 16, offset & 0xFFFF
        joined = blocks[block] + b"".join(b for o, b in sorted(blocks.items()) if o > block)
        assert joined[within : within + len(text)].decode() == text


def test_index_queries_return_exactly_the_overlapping_rows(tmp_path):
    regions, features = _random_features()
    gz = tmp_path / "a.gff3.gz"
    assert write_indexed_gff3(features, gz, regions=regions, version="test")
    header, names, _refs = _read_tbi(gz.with_name(gz.name + ".tbi"))
    assert header == (0, 1, 4, 5, "#", 0)
    assert names == ["chr1", "chr2", "chr3"]

    rows = [line for line in gzip.decompress(gz.read_bytes()).decode().splitlines()
            if not line.startswith("#")]
    rnd = random.Random(9)
    for _ in range(150):
        name, length = rnd.choice(regions)
        begin = rnd.randint(0, length - 1)
        end = min(length, begin + rnd.choice([1, 500, 50_000, 2_000_000]))
        expected = sorted(
            line for line in rows
            if line.split("\t")[0] == name
            and int(line.split("\t")[3]) - 1 < end and int(line.split("\t")[4]) > begin
        )
        assert _query(gz, name, begin, end) == expected


def test_htslib_reads_the_index_when_pysam_is_installed(tmp_path):
    pysam = pytest.importorskip("pysam")
    regions, features = _random_features(seed=5)
    gz = tmp_path / "a.gff3.gz"
    assert write_indexed_gff3(features, gz, regions=regions, version="test")
    rows = [line for line in gzip.decompress(gz.read_bytes()).decode().splitlines()
            if not line.startswith("#")]
    tabix = pysam.TabixFile(str(gz))
    rnd = random.Random(11)
    for _ in range(300):
        name, length = rnd.choice(regions)
        begin = rnd.randint(0, length - 1)
        end = min(length, begin + rnd.choice([1, 500, 50_000, 2_000_000]))
        expected = sorted(
            line for line in rows
            if line.split("\t")[0] == name
            and int(line.split("\t")[3]) - 1 < end and int(line.split("\t")[4]) > begin
        )
        assert sorted(tabix.fetch(name, begin, end)) == expected


def test_rows_must_arrive_sorted_and_contiguous():
    index = TabixIndex()
    index.add("chr1", 100, 200, 0, 10)
    with pytest.raises(ValueError, match="not sorted"):
        index.add("chr1", 50, 60, 10, 20)
    index.add("chr2", 5, 6, 20, 30)
    with pytest.raises(ValueError, match="not contiguous"):
        index.add("chr1", 300, 400, 30, 40)


def test_a_feature_on_an_unnamed_sequence_is_refused_not_dropped(tmp_path):
    features = [_Feature("chr1", "g1", [(10, 20)]), _Feature("chrX", "g2", [(10, 20)])]
    gz = tmp_path / "a.gff3.gz"
    with pytest.raises(ValueError, match="g2 is on chrX"):
        write_indexed_gff3(features, gz, regions=[("chr1", 1000)], version="t")


def test_a_coordinate_beyond_the_tbi_limit_writes_the_file_without_an_index(tmp_path):
    features = [_Feature("big", "g1", [(600_000_000, 600_000_100)])]
    gz = tmp_path / "big.gff3.gz"
    assert not write_indexed_gff3(features, gz, regions=[("big", 700_000_000)], version="t")
    assert gz.exists()
    assert not gz.with_name(gz.name + ".tbi").exists()
    with pytest.raises(IndexLimitError):
        TabixIndex().add("big", 600_000_000, 600_000_100, 0, 1)
