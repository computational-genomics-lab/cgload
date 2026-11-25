"""Sequence for the portal: the stored assembly, spliced segments, proteins.

The database holds no sequence (D-004). Every base the portal shows is read
from the assembly that ``load`` copied into the data directory and indexed, so
a sequence on a page is the sequence the annotation was loaded against.

Translation follows NCBI's genetic codes (D-113). A protein shown by the portal
is translated by cgload from the CDS rows; the GFF3 files cgload reads carry no
protein sequence to compare with, so the page says where a file declares an
exception (``transl_except``, ``exception``) that a plain translation does not
apply.
"""

from __future__ import annotations

import threading
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

from pyfaidx import Fasta

_COMPLEMENT = str.maketrans("ACGTURYKMBDHVNacgturykmbdhvn", "TGCAAYRMKVHDBNtgcaayrmkvhdbn")


def reverse_complement(sequence: str) -> str:
    """The reverse complement, IUPAC ambiguity codes included; case is kept."""
    return sequence.translate(_COMPLEMENT)[::-1]


#: NCBI's genetic codes, by table number: the amino acid of each of the 64
#: codons in TCAG order (TTT, TTC, TTA, TTG, TCT, ...), the codons that may
#: start a protein, and the table's name. From NCBI's gc.prt; ``*`` is a stop.
GENETIC_CODES: dict[int, tuple[str, str, str]] = {
    1: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M---------------M---------------M----------------------------",
        "Standard"),
    2: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "--------------------------------MMMM---------------M------------",
        "Vertebrate Mitochondrial"),
    3: ("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "----------------------------------MM---------------M------------",
        "Yeast Mitochondrial"),
    4: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM---------------M------------MMMM---------------M------------",
        "Mold, Protozoan and Coelenterate Mitochondrial"),
    5: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "---M----------------------------MMMM---------------M------------",
        "Invertebrate Mitochondrial"),
    6: ("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "-----------------------------------M----------------------------",
        "Ciliate, Dasycladacean and Hexamita Nuclear"),
    9: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "-----------------------------------M---------------M------------",
        "Echinoderm and Flatworm Mitochondrial"),
    10: ("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Euplotid Nuclear"),
    11: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M---------------M------------MMMM---------------M------------",
         "Bacterial, Archaeal and Plant Plastid"),
    12: ("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-------------------M---------------M----------------------------",
         "Alternative Yeast Nuclear"),
    13: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
         "---M------------------------------MM---------------M------------",
         "Ascidian Mitochondrial"),
    14: ("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Alternative Flatworm Mitochondrial"),
    15: ("FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Blepharisma Macronuclear"),
    16: ("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Chlorophycean Mitochondrial"),
    21: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
         "-----------------------------------M---------------M------------",
         "Trematode Mitochondrial"),
    22: ("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Scenedesmus obliquus Mitochondrial"),
    23: ("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "--------------------------------M--M---------------M------------",
         "Thraustochytrium Mitochondrial"),
    24: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
         "---M---------------M---------------M---------------M------------",
         "Rhabdopleuridae Mitochondrial"),
    25: ("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M-------------------------------M---------------M------------",
         "Candidate Division SR1 and Gracilibacteria"),
    26: ("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-------------------M---------------M----------------------------",
         "Pachysolen tannophilus Nuclear"),
    27: ("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Karyorelict Nuclear"),
    28: ("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Condylostoma Nuclear"),
    29: ("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Mesodinium Nuclear"),
    30: ("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Peritrich Nuclear"),
    31: ("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "-----------------------------------M----------------------------",
         "Blastocrithidia Nuclear"),
    32: ("FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M---------------M------------MMMM---------------M------------",
         "Balanophoraceae Plastid"),
    33: ("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
         "---M---------------M---------------M---------------M------------",
         "Cephalodiscidae Mitochondrial"),
}

_BASES = "TCAG"
_CODON_INDEX = {
    first + second + third: index
    for index, (first, second, third) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES
    )
}


def translate(cds: str, table: int = 1, *, from_start: bool = True) -> str:
    """Translate a coding sequence that begins with a whole codon.

    A codon with any base other than A, C, G or T becomes ``X``; a final
    incomplete codon is dropped. With ``from_start``, the first codon is read
    as ``M`` when the table lists it as a start codon, as NCBI does for a CDS
    that is complete at its 5' end.
    """
    if table not in GENETIC_CODES:
        raise ValueError(f"no NCBI genetic code {table}")
    amino, starts, _ = GENETIC_CODES[table]
    sequence = cds.upper().replace("U", "T")
    protein = []
    for position in range(0, len(sequence) - 2, 3):
        index = _CODON_INDEX.get(sequence[position : position + 3])
        if index is None:
            protein.append("X")
        elif position == 0 and from_start and starts[index] == "M":
            protein.append("M")
        else:
            protein.append(amino[index])
    return "".join(protein)


@dataclass(frozen=True)
class Piece:
    """One row of a feature: 1-based, inclusive, on its own strand."""

    start: int
    end: int
    strand: str
    phase: int | None = None


def transcript_order(pieces: Sequence[Piece],
                     file_order: Sequence[int] | None = None) -> list[Piece]:
    """The rows of one feature in the order they are transcribed.

    On one strand that is ascending coordinates on ``+`` and descending on
    ``-``, whatever order the file wrote them in. A trans-spliced feature has
    rows on both strands (D-068, D-096), and only the file says their order, so
    ``file_order`` (each row's input line) decides it there.
    """
    strands = {piece.strand for piece in pieces}
    if len(strands) > 1 and file_order is not None:
        return [piece for _, piece in sorted(zip(file_order, pieces, strict=True),
                                             key=lambda pair: pair[0])]
    ordered = sorted(pieces, key=lambda piece: (piece.start, piece.end))
    if strands == {"-"}:
        ordered.reverse()
    return ordered


class AssemblySequence:
    """The stored assembly of one genome, read through its ``.fai`` index.

    pyfaidx reads with seek and read on one file handle, so reads are
    serialised with a lock: the portal serves several requests at once.
    """

    def __init__(self, path: Path, *, circular: Iterable[str] = ()) -> None:
        self.path = Path(path)
        self.circular = frozenset(circular)
        self._lock = threading.Lock()
        self._fasta: Fasta | None = None

    def _handle(self) -> Fasta:
        if self._fasta is None:
            # rebuild=False: the portal reads the index load wrote, never writes one.
            self._fasta = Fasta(str(self.path), rebuild=False, as_raw=True,
                                sequence_always_upper=False)
        return self._fasta

    def length(self, seqid: str) -> int:
        with self._lock:
            return len(self._handle()[seqid])

    def fetch(self, seqid: str, start: int, end: int) -> str:
        """Bases ``start`` to ``end``, 1-based and inclusive.

        On a circular sequence a row may end past the sequence's length, as
        GFF3 writes a row that crosses the origin (D-081); the bases past the
        end are read from its beginning.
        """
        if start < 1 or end < start:
            raise ValueError(f"cannot read {seqid}:{start}-{end}")
        with self._lock:
            record = self._handle()[seqid]
            length = len(record)
            if end <= length:
                return str(record[start - 1 : end])
            if seqid not in self.circular or end - start + 1 > length:
                raise ValueError(f"{seqid}:{start}-{end} runs past the end of {seqid}")
            if start > length:
                return str(record[start - length - 1 : end - length])
            return str(record[start - 1 : length]) + str(record[0 : end - length])

    def piece(self, seqid: str, piece: Piece) -> str:
        bases = self.fetch(seqid, piece.start, piece.end)
        return reverse_complement(bases) if piece.strand == "-" else bases

    def spliced(self, seqid: str, pieces: Sequence[Piece]) -> str:
        """The pieces joined in the order given, each read on its own strand."""
        return "".join(self.piece(seqid, piece) for piece in pieces)

    def close(self) -> None:
        with self._lock:
            if self._fasta is not None:
                self._fasta.close()
                self._fasta = None


def coding_sequence(sequence: AssemblySequence, seqid: str, ordered: Sequence[Piece]) -> str:
    """A CDS's bases from its first whole codon, ``ordered`` in transcript order.

    The first row's phase says how many bases precede that codon (GFF3): zero
    for a CDS complete at its 5' end, one or two where the file says the CDS
    begins mid-codon.
    """
    bases = sequence.spliced(seqid, ordered)
    skip = (ordered[0].phase or 0) if ordered else 0
    return bases[skip:]


def fasta_text(header: str, bases: str, width: int = 60) -> str:
    lines = [f">{header}"]
    lines.extend(bases[index : index + width] for index in range(0, len(bases), width))
    return "\n".join(lines) + "\n"
