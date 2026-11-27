"""BLAST against the genomes of a cgload database (D-116).

Two halves.

``build`` (``cgload blastdb``) writes two BLAST+ databases into each assembly's
data directory, under ``blast/``: ``genome``, the assembly's sequences, and
``proteins``, every CDS feature translated as the portal translates it. Each is
checked by counting, as a load is: the genome database must hold every
sequence of the assembly and the protein database one protein for every CDS
feature, and BLAST+'s own count of what it stored must agree. The counts and
the outcome go to ``blast/manifest.json``. A database is replaced only once
its successor has passed; one that fails is not kept, and neither is its
predecessor, which may name rows of an earlier load. One build of an
assembly runs at a time: a second is refused at once and changes nothing,
since its clean-up would remove what the first is building (D-117).

``Blast`` (``cgload serve``) runs a search, and is built so that a visitor's
text never reaches anything but BLAST's input:

* BLAST+ is started with a list of arguments, never through a shell, and the
  programs are found when the portal starts, never from a request.
* The query reaches BLAST on its standard input, as one FASTA record whose
  header is always ``>query``; it is never an argument or a file name, and no
  text from a visitor or from the annotation files is an argument either.
* The databases are those the manifests list, chosen by assembly id.
* The program, task, E-value and number of hits each come from a fixed list;
  the query is one sequence, of the letters of its alphabet, at most
  10,000 long.
* Only a few searches run at once, as many again may wait for a turn, and
  any more are turned away at once, so searches cannot take every thread the
  pages need. Each is stopped after a time limit, and what BLAST writes to
  its error stream is logged, not shown.
"""

from __future__ import annotations

import contextlib
import json
import logging
import os
import re
import shutil
import subprocess  # BLAST+ only, with argument lists (module docstring)  # nosec B404
import sys
import tempfile
import threading
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy import and_, func, select

from cgload.db.export_source import StoredFeature, features_where
from cgload.db.schema import feature, sequence_region

log = logging.getLogger("cgload.portal.blast")

#: program: (the alphabet of its query, the database it searches)
PROGRAMS = {
    "blastn": ("nucleotide", "genome"),
    "tblastn": ("protein", "genome"),
    "blastp": ("protein", "proteins"),
    "blastx": ("nucleotide", "proteins"),
}
#: blastn's tasks; the first is the default (fast, for near-identical sequence).
TASKS = ("megablast", "blastn")
EVALUES = ("10", "1", "0.1", "0.01", "0.001", "1e-5", "1e-10", "1e-20", "1e-50")
DEFAULT_EVALUE = "0.001"
HIT_CHOICES = (10, 25, 50, 100)
DEFAULT_HITS = 25
MAX_HSPS = 3
MAX_QUERY = 10_000
MIN_QUERY = 10
MAX_OUTPUT = 32 * 1024 * 1024
NUCLEOTIDE_LETTERS = frozenset("ACGTURYKMSWBDHVN")
PROTEIN_LETTERS = frozenset("ABCDEFGHIKLMNOPQRSTUVWXYZ*")
BINARIES = ("blastn", "tblastn", "blastp", "blastx", "makeblastdb")
OUTFMT = ("6 qseqid stitle pident length mismatch gapopen qstart qend sstart send "
          "evalue bitscore qlen slen sstrand qseq sseq")
MANIFEST = "manifest.json"
GENOME_TITLE = re.compile(r"^g(\d+) (\S+)$")
PROTEIN_TITLE = re.compile(r"^p(\d+)\.(\d+)(?: (.*))?$")


class QueryError(ValueError):
    """A search the portal will not run, and why, in words for its user."""


class BlastBusy(RuntimeError):
    """Every BLAST slot was taken for as long as a search may wait."""


class BuildBusy(RuntimeError):
    """Another ``cgload blastdb`` holds this assembly's build lock."""


class BlastFailed(RuntimeError):
    """BLAST did not finish a search: too slow, or an error of its own."""


# --- the request --------------------------------------------------------------------


@dataclass(frozen=True)
class SearchRequest:
    program: str
    query: str
    assembly_ids: tuple[int, ...]
    evalue: str = DEFAULT_EVALUE
    max_hits: int = DEFAULT_HITS
    task: str | None = None
    title: str = "query"


def clean_query(text: str, alphabet: str) -> tuple[str, str]:
    """The sequence of a pasted query, and its FASTA header's first word if it
    had one. Accepts one FASTA record or bare sequence, with spaces, line breaks
    and position numbers, as GenBank prints sequence."""
    lines = text.replace("\r", "\n").split("\n")
    headers = [line for line in lines if line.startswith(">")]
    if len(headers) > 1:
        raise QueryError(f"One sequence at a time, please: this has {len(headers)} "
                         "FASTA records.")
    title = headers[0][1:].split(maxsplit=1)[0][:60] if headers and headers[0][1:].split() \
        else "query"
    sequence = re.sub(r"[\s\d]", "", "".join(line for line in lines
                                             if not line.startswith(">"))).upper()
    if alphabet == "nucleotide":
        sequence = sequence.replace("-", "").replace("U", "T")
    else:
        sequence = sequence.replace("-", "").rstrip("*")
    if not sequence:
        raise QueryError("There is no sequence to search.")
    if len(sequence) > MAX_QUERY:
        raise QueryError(f"The sequence is {len(sequence):,} letters long; "
                         f"a search takes at most {MAX_QUERY:,}.")
    allowed = NUCLEOTIDE_LETTERS if alphabet == "nucleotide" else PROTEIN_LETTERS
    wrong = sorted(set(sequence) - allowed)
    if wrong:
        hint = (" It looks like a protein: choose blastp or tblastn."
                if alphabet == "nucleotide" and set(sequence) <= PROTEIN_LETTERS else "")
        letters = ("A, C, G, T and the IUPAC codes" if alphabet == "nucleotide"
                   else "the letters of the amino acids")
        raise QueryError(f"A {alphabet} sequence is {letters}; this one also has "
                         f"{' '.join(wrong[:12])}{' …' if len(wrong) > 12 else ''}.{hint}")
    if len(sequence) < MIN_QUERY:
        raise QueryError(f"The sequence is {len(sequence)} letters long; "
                         f"a search needs at least {MIN_QUERY}.")
    title = re.sub(r"[^A-Za-z0-9._:|-]", "_", title) or "query"
    return sequence, title


def make_request(form: dict, available) -> SearchRequest:
    """A search from what a form or a JSON body sent, every value checked.
    ``available(program)`` gives the assemblies that program can search."""
    def one(name: str, default: str = "") -> str:
        value = form.get(name, default)
        if isinstance(value, list):
            value = value[0] if value else default
        return str(value) if value is not None else default

    program = one("program", "blastp")
    if program not in PROGRAMS:
        raise QueryError(f"The program is one of {', '.join(PROGRAMS)}.")
    task = None
    if program == "blastn":  # a form sends a task whatever the program; only blastn has one
        task = one("task") or TASKS[0]
        if task not in TASKS:
            raise QueryError(f"blastn's task is one of {', '.join(TASKS)}.")
    evalue = one("evalue", DEFAULT_EVALUE)
    if evalue not in EVALUES:
        raise QueryError(f"The E-value is one of {', '.join(EVALUES)}.")
    hits = one("max_hits", str(DEFAULT_HITS))
    if not hits.isdigit() or int(hits) not in HIT_CHOICES:
        raise QueryError(f"The number of hits is one of {', '.join(map(str, HIT_CHOICES))}.")
    available = sorted(set(available(program)))
    raw = form.get("genome", "all")
    chosen = raw if isinstance(raw, list) else [raw]
    if chosen in (["all"], [], [""]):
        ids = available
    else:
        ids = []
        for value in chosen:
            text = str(value)
            if not text.isdigit() or int(text) not in available:
                raise QueryError(f"{program} cannot search that genome here; choose another, "
                                 "or all of them.")
            ids.append(int(text))
    if not ids:
        raise QueryError(f"No genome here has a database {program} can search.")
    raw_query = one("query")
    if len(raw_query) > 3 * MAX_QUERY:
        raise QueryError(f"The query is too long: a search takes at most {MAX_QUERY:,} letters.")
    sequence, title = clean_query(raw_query, PROGRAMS[program][0])
    return SearchRequest(program=program, query=sequence, assembly_ids=tuple(sorted(set(ids))),
                         evalue=evalue, max_hits=int(hits), task=task, title=title)


# --- the result -----------------------------------------------------------------------


@dataclass
class Hsp:
    pident: float
    length: int
    mismatches: int
    gap_opens: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    strand: str | None
    qseq: str
    sseq: str

    def match_line(self) -> str:
        """``|`` under each identical pair, as BLAST marks identities."""
        return "".join("|" if a == b and a != "-" else " "
                       for a, b in zip(self.qseq, self.sseq, strict=False))


    def blocks(self, program: str, width: int = 60
               ) -> list[tuple[int, str, str, str, int, int, int]]:
        """The alignment in rows of ``width``: query start, query, match line,
        subject, subject start, query end, subject end; positions step by three
        on the side BLAST translated."""
        q_step = 3 if program == "blastx" else 1
        s_step = 3 if program == "tblastn" else 1
        q_dir = 1 if self.qend >= self.qstart else -1
        s_dir = 1 if self.send >= self.sstart else -1
        match = self.match_line()
        q_at, s_at = self.qstart, self.sstart
        rows = []
        for begin in range(0, len(self.qseq), width):
            q, m, t = (self.qseq[begin:begin + width], match[begin:begin + width],
                       self.sseq[begin:begin + width])
            q_letters = sum(1 for letter in q if letter != "-")
            s_letters = sum(1 for letter in t if letter != "-")
            q_end = q_at + q_dir * (q_letters * q_step - 1) if q_letters else q_at
            s_end = s_at + s_dir * (s_letters * s_step - 1) if s_letters else s_at
            rows.append((q_at, q, m, t, s_at, q_end, s_end))
            q_at = q_end + q_dir if q_letters else q_at
            s_at = s_end + s_dir if s_letters else s_at
        return rows


@dataclass
class Hit:
    kind: str  # "genome" or "proteins"
    assembly_id: int
    subject: str  # a seqid, or a feature's row
    description: str
    length: int
    hsps: list[Hsp] = field(default_factory=list)

    @property
    def best(self) -> Hsp:
        return self.hsps[0]

    @property
    def name(self) -> str:
        """A sequence's name, or a protein's: its feature's name, as on its page."""
        if self.kind == "genome":
            return self.subject
        return self.description.split(" ", 1)[0] or f"row {self.subject}"

    @property
    def product(self) -> str:
        """A protein's product, as its CDS states it; empty for a genome hit."""
        if self.kind == "genome":
            return ""
        parts = self.description.split(" ", 1)
        return parts[1] if len(parts) > 1 else ""


@dataclass
class SearchResult:
    request: SearchRequest
    hits: list[Hit]
    seconds: float
    databases: int


def span_of(hit: Hit, within: int) -> tuple[int, int]:
    """The subject span of a hit's best HSP, widened to its other HSPs on the
    same strand that lie within ``within`` bases of it."""
    low, high = sorted((hit.best.sstart, hit.best.send))
    first_low, first_high = low, high
    for hsp in hit.hsps[1:]:
        a, b = sorted((hsp.sstart, hsp.send))
        if hsp.strand == hit.best.strand and a >= first_low - within and b <= first_high + within:
            low, high = min(low, a), max(high, b)
    return low, high


NCBI_BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
NCBI_ACCESSION = re.compile(r"[A-Z]{1,4}_?\d{5,9}(\.\d+)?")
NCBI_MAX_LETTERS = 4000


def ncbi_url(program: str, sequence: str, accession: str | None = None) -> str | None:
    """A link that opens NCBI's BLAST with this query filled in, to search
    NCBI's own databases: by accession when there is one, else by sequence
    when it is short enough to travel in a link. Nothing is sent from here."""
    from urllib.parse import urlencode

    if program not in PROGRAMS:
        return None
    if accession and NCBI_ACCESSION.fullmatch(accession):
        query = accession
    elif 0 < len(sequence) <= NCBI_MAX_LETTERS:
        query = sequence
    else:
        return None
    return NCBI_BLAST + "?" + urlencode({"PAGE_TYPE": "BlastSearch", "PROGRAM": program,
                                         "QUERY": query, "LINK_LOC": "blasthome"})


def parse(text: str, kind: str) -> list[Hit]:
    """BLAST's tabular output (``OUTFMT``) as hits, in BLAST's order."""
    hits: dict[tuple[int, str], Hit] = {}
    for line in text.splitlines():
        cells = line.split("\t")
        if len(cells) != 17:
            continue
        title = cells[1]
        if kind == "genome":
            match = GENOME_TITLE.match(title)
            if not match:
                continue
            assembly_id, subject, description = int(match.group(1)), match.group(2), ""
        else:
            match = PROTEIN_TITLE.match(title)
            if not match:
                continue
            assembly_id, subject = int(match.group(1)), match.group(2)
            description = match.group(3) or ""
        try:
            hsp = Hsp(pident=float(cells[2]), length=int(cells[3]), mismatches=int(cells[4]),
                      gap_opens=int(cells[5]), qstart=int(cells[6]), qend=int(cells[7]),
                      sstart=int(cells[8]), send=int(cells[9]), evalue=float(cells[10]),
                      bitscore=float(cells[11]),
                      strand=cells[14] if cells[14] in ("plus", "minus") else None,
                      qseq=cells[15], sseq=cells[16])
        except ValueError:
            continue
        key = (assembly_id, subject)
        if key not in hits:
            hits[key] = Hit(kind=kind, assembly_id=assembly_id, subject=subject,
                            description=description, length=int(cells[13]))
        hits[key].hsps.append(hsp)
    return list(hits.values())


# --- running it -----------------------------------------------------------------------


def find_binaries(directory: Path | None = None) -> dict[str, str]:
    """The BLAST+ programs this computer has: in ``directory`` if given, else
    on PATH or beside the Python running cgload, as in a conda environment
    holding both (a service's PATH may not name it). Resolved once, when the
    portal starts; never from a request."""
    places = ([str(directory)] if directory is not None
              else [os.environ.get("PATH", ""), str(Path(sys.executable).parent)])
    found = {}
    for name in BINARIES:
        for place in places:
            path = shutil.which(name, path=place)
            if path:
                found[name] = str(Path(path).resolve())
                break
    return found


def version_of(binary: str) -> str:
    try:
        done = subprocess.run([binary, "-version"], capture_output=True, text=True,  # noqa: S603  # nosec B603
                              timeout=30, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    first = done.stdout.strip().splitlines()
    return first[0].split(":", 1)[-1].strip() if first else "unknown"


@dataclass
class Databases:
    """One assembly's BLAST databases, as its manifest describes them."""

    assembly_id: int
    directory: Path
    genome: bool
    proteins: bool
    manifest: dict


def databases_of(genomes: Iterable) -> dict[int, Databases]:
    """The assemblies whose BLAST databases passed their check and still match
    the load: a manifest from before a reload names rows that no longer exist."""
    found = {}
    for genome in genomes:
        if genome.directory is None:
            continue
        directory = genome.directory / "blast"
        try:
            manifest = json.loads((directory / MANIFEST).read_text())
        except (OSError, ValueError):
            continue
        if manifest.get("assembly_id") != genome.assembly_id:
            continue
        loaded = genome.loaded_at.isoformat() if genome.loaded_at else None
        if manifest.get("loaded_at") != loaded:
            log.warning("BLAST databases of %s are older than its load: run cgload blastdb",
                        genome.label)
            continue
        have = {kind: manifest.get(kind, {}).get("check") == "PASS"
                and any(directory.glob(f"{kind}.*"))
                for kind in ("genome", "proteins")}
        if any(have.values()):
            found[genome.assembly_id] = Databases(genome.assembly_id, directory,
                                                  have["genome"], have["proteins"], manifest)
    return found


class Blast:
    """Searches against the portal's genomes, a few at a time."""

    def __init__(self, portal, *, binaries: dict[str, str] | None = None, jobs: int = 2,
                 threads: int = 4, timeout: float = 120.0, wait: float = 30.0) -> None:
        self.portal = portal
        self.binaries = find_binaries() if binaries is None else binaries
        if jobs < 1:
            raise ValueError("BLAST needs at least one job; leave it out to turn it off")
        self.jobs, self.threads, self.timeout, self.wait = jobs, threads, timeout, wait
        self._slots = threading.BoundedSemaphore(jobs)  # searches running
        self._admitted = threading.BoundedSemaphore(2 * jobs)  # running or waiting
        self.refresh()

    def refresh(self) -> None:
        self.databases = databases_of(self.portal.genomes())

    @property
    def programs(self) -> list[str]:
        return [name for name in PROGRAMS if name in self.binaries]

    @property
    def ready(self) -> bool:
        return bool(self.programs) and bool(self.databases)

    def searchable(self, program: str) -> list[int]:
        kind = PROGRAMS[program][1]
        return [aid for aid, dbs in sorted(self.databases.items()) if getattr(dbs, kind)]

    def search(self, request: SearchRequest) -> SearchResult:
        if request.program not in self.binaries:
            raise QueryError(f"{request.program} is not installed on this server.")
        kind = PROGRAMS[request.program][1]
        targets = [self.databases[aid].directory / kind for aid in request.assembly_ids
                   if aid in self.databases and getattr(self.databases[aid], kind)]
        if not targets:
            raise QueryError("None of the genomes chosen has that kind of BLAST database.")
        if not self._admitted.acquire(blocking=False):
            raise BlastBusy("as many BLAST searches as the portal allows are running or waiting")
        try:
            if not self._slots.acquire(timeout=self.wait):
                raise BlastBusy("the portal is running as many BLAST searches as it allows")
            started = time.monotonic()
            try:
                text = "".join(self._run(request, group) for group in _groups(targets))
            finally:
                self._slots.release()
        finally:
            self._admitted.release()
        hits = parse(text, kind)
        hits.sort(key=lambda hit: (hit.best.evalue, -hit.best.bitscore))
        return SearchResult(request, hits[:request.max_hits], time.monotonic() - started,
                            len(targets))

    def command(self, request: SearchRequest, names: list[str]) -> list[str]:
        """The argument list for one run of BLAST+; the query is not in it."""
        args = [self.binaries[request.program], "-db", " ".join(names), "-outfmt", OUTFMT,
                "-evalue", request.evalue, "-max_target_seqs", str(request.max_hits),
                "-max_hsps", str(MAX_HSPS), "-num_threads", str(self.threads)]
        if request.task:
            args += ["-task", request.task]
        return args

    def _run(self, request: SearchRequest, group: tuple[Path, list[str]]) -> str:
        cwd, names = group
        fasta = f">query\n{request.query}\n"
        environment = {name: value for name, value in os.environ.items()
                       if name not in ("BLASTDB", "NCBI")}
        try:
            done = subprocess.run(  # noqa: S603 - an argument list of fixed choices; no shell  # nosec B603
                self.command(request, names), input=fasta, capture_output=True, text=True,
                timeout=self.timeout, cwd=cwd, check=False, env=environment)
        except subprocess.TimeoutExpired:
            raise BlastFailed(f"The search took longer than {self.timeout:.0f} s and was "
                              "stopped. Search one genome rather than all of them, or a "
                              "shorter sequence.") from None
        except OSError as exc:
            log.error("BLAST could not start: %s", exc)
            raise BlastFailed("BLAST could not start on this server.") from None
        if done.returncode != 0:
            log.warning("BLAST exited with %s: %s", done.returncode, done.stderr.strip()[:2000])
            raise BlastFailed("BLAST could not search this sequence.")
        return done.stdout[:MAX_OUTPUT]


def _groups(targets: list[Path]) -> list[tuple[Path, list[str]]]:
    """How to name the databases (``<assembly>/blast/<kind>``) to BLAST+,
    which takes a space-separated list and refuses paths with spaces: relative
    to the folder that holds their assemblies, in one run, when none of those
    names has a space (cgload's own folder names never do); otherwise one run
    per database, from its assembly's folder. Every name keeps a folder in it,
    ``blast/`` at least, so BLAST+ reads it from that folder and never looks
    for a database of that name in the places BLASTDB names."""
    homes = [path.parent.parent for path in targets]
    common = Path(os.path.commonpath([str(home) for home in homes]))
    names = [path.relative_to(common).as_posix() for path in targets]
    if not any(re.search(r"\s", name) for name in names):
        return [(common, names)]
    return [(path.parent.parent, [f"blast/{path.name}"]) for path in targets]


# --- building the databases -----------------------------------------------------------


@dataclass
class BuildResult:
    assembly_id: int
    label: str
    check: str  # PASS or FAIL
    sequences: int = 0
    proteins: int = 0
    genome_letters: int = 0
    protein_letters: int = 0
    cds_features: int = 0
    with_internal_stop: int = 0
    seconds: float = 0.0
    problems: list[str] = field(default_factory=list)


def proteins_of(portal, genome) -> Iterator[tuple[StoredFeature, str]]:
    """Every CDS feature of an assembly with its protein, as the portal
    translates it, one sequence region at a time."""
    from cgload.portal.data import _partial_start, _pieces, file_order, genetic_code
    from cgload.portal.sequence import coding_sequence, transcript_order, translate

    handle = portal.sequence(genome)
    with portal.engine.connect() as connection:
        regions = connection.execute(
            select(sequence_region.c.sequence_region_id)
            .where(sequence_region.c.assembly_id == genome.assembly_id)
            .order_by(sequence_region.c.sequence_region_id)).scalars().all()
        for region in regions:
            items = features_where(connection, and_(
                feature.c.assembly_id == genome.assembly_id,
                feature.c.sequence_region_id == region,
                feature.c.feature_type == "CDS"))
            for item in items:
                ordered = transcript_order(_pieces(item), file_order(connection, item))
                cds = coding_sequence(handle, item.seqid, ordered)
                protein = translate(cds, genetic_code(item),
                                    from_start=(ordered[0].phase or 0) == 0
                                    and not _partial_start(item))
                yield item, protein


def _title(text: str) -> str:
    """One line of printable text for a FASTA header: tabs and line breaks
    would break BLAST's tabular output, which this text reaches."""
    return re.sub(r"\s+", " ", re.sub(r"[\x00-\x1f\x7f]", " ", text)).strip()[:200]


def _genome_records(genome, result: BuildResult) -> Iterator[str]:
    """The assembly as FASTA, each header ``g<assembly> <seqid>``."""
    with genome.fasta.open() as handle:
        for line in handle:
            if line.startswith(">"):
                result.sequences += 1
                yield f">g{genome.assembly_id} {line[1:].split()[0]}\n"
            else:
                result.genome_letters += len(line.strip())
                yield line


def _protein_records(portal, genome, result: BuildResult) -> Iterator[str]:
    """Every protein as FASTA, each header ``p<assembly>.<row> <name> <product>``."""
    from cgload.portal.data import display_name

    for item, protein in proteins_of(portal, genome):
        result.proteins += 1
        if "*" in protein.rstrip("*"):
            result.with_internal_stop += 1
        # The name is one word, so that a hit's name and product can be told
        # apart again (Hit.name); a name with spaces keeps them as underscores.
        name = re.sub(r"\s+", "_", display_name(item).strip()) or f"row{item.row_id}"
        product = " ".join(item.attributes.get("product", []))
        yield f">p{genome.assembly_id}.{item.row_id} {_title(f'{name} {product}')}\n"
        body = protein.rstrip("*") or "X"
        result.protein_letters += len(body)
        for begin in range(0, len(body), 60):
            yield body[begin:begin + 60] + "\n"


LOCK = "blast.lock"


@contextlib.contextmanager
def building(directory: Path) -> Iterator[None]:
    """Hold the build lock of the assembly in ``directory`` for as long as the
    block runs, or raise BuildBusy at once if another process holds it. The
    lock is the operating system's on ``blast.lock``, so it goes with the
    process that held it, however that process ends (D-117)."""
    handle = (directory / LOCK).open("a+")
    try:
        try:
            if os.name == "nt":
                import msvcrt

                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            raise BuildBusy("another cgload blastdb is building this assembly's BLAST "
                            "databases now; this one changed nothing. Run it again when "
                            "that one has finished.") from None
        try:
            yield
        finally:
            if os.name == "nt":
                import msvcrt

                with contextlib.suppress(OSError):
                    msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
    finally:
        handle.close()


def build(portal, genomes, makeblastdb: str, *, log_line=print) -> list[BuildResult]:
    """Build and check the BLAST databases of ``genomes``; see the module docstring."""
    results = []
    tool = version_of(makeblastdb)
    for genome in genomes:
        result = BuildResult(genome.assembly_id, genome.label, "FAIL")
        results.append(result)
        if genome.directory is None or genome.fasta is None:
            result.problems.append("its assembly file is not where the database says; "
                                   "give --data-dir")
        else:
            try:
                with building(genome.directory):
                    _build_one(portal, genome, makeblastdb, tool, result)
            except BuildBusy as exc:
                result.problems.append(str(exc))
        log_line(f"assembly {genome.assembly_id} ({genome.label}): "
                 f"{_count(result.sequences, 'sequence')}, "
                 f"{_count(result.proteins, 'protein')}: {result.check}"
                 + "".join(f"\n  {problem}" for problem in result.problems))
    return results


def _build_one(portal, genome, makeblastdb: str, tool: str, result: BuildResult) -> None:
    """One assembly's build, under its lock: blast.new, checked, then swapped in."""
    from cgload import __version__

    started = time.monotonic()
    final = genome.directory / "blast"
    work = genome.directory / "blast.new"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir()
    try:
        with portal.engine.connect() as connection:
            expected_sequences = connection.execute(
                select(func.count()).select_from(sequence_region)
                .where(sequence_region.c.assembly_id == genome.assembly_id)).scalar_one()
            expected_cds = connection.execute(
                select(func.count()).select_from(feature).where(
                    feature.c.assembly_id == genome.assembly_id,
                    feature.c.feature_type == "CDS",
                    feature.c.segment_index == 0)).scalar_one()

        added_genome = _makeblastdb(makeblastdb, "nucl", work, "genome",
                                    f"cgload assembly {genome.assembly_id} genome",
                                    _genome_records(genome, result))
        added_proteins = (_makeblastdb(makeblastdb, "prot", work, "proteins",
                                       f"cgload assembly {genome.assembly_id} proteins",
                                       _protein_records(portal, genome, result))
                          if expected_cds else 0)
        result.cds_features = expected_cds
        if not (result.sequences == expected_sequences == added_genome):
            result.problems.append(
                f"genome: {expected_sequences:,} sequences in the database, "
                f"{result.sequences:,} read from the FASTA, {added_genome:,} in BLAST")
        if not (result.proteins == expected_cds == added_proteins):
            result.problems.append(
                f"proteins: {expected_cds:,} CDS features, {result.proteins:,} "
                f"translated, {added_proteins:,} in BLAST")
        result.check = "FAIL" if result.problems else "PASS"
        result.seconds = time.monotonic() - started
        manifest = {
            "assembly_id": genome.assembly_id, "label": genome.label,
            "loaded_at": genome.loaded_at.isoformat() if genome.loaded_at else None,
            "cgload": __version__, "makeblastdb": tool,
            "built_at": datetime.now(UTC).isoformat(timespec="seconds"),
            "genome": {"sequences": result.sequences, "in_database": expected_sequences,
                       "in_blast": added_genome, "letters": result.genome_letters,
                       "check": result.check
                       if not any(p.startswith("genome") for p in result.problems)
                       else "FAIL"},
            "proteins": {"cds_features": expected_cds, "translated": result.proteins,
                         "in_blast": added_proteins, "letters": result.protein_letters,
                         "with_internal_stop": result.with_internal_stop,
                         "check": result.check
                         if not any(p.startswith("proteins") for p in result.problems)
                         else "FAIL"},
        }
        (work / MANIFEST).write_text(json.dumps(manifest, indent=1) + "\n")
    except (OSError, subprocess.SubprocessError, ValueError, KeyError) as exc:
        result.problems.append(str(exc))
        result.check = "FAIL"
    if result.check == "PASS":
        old = genome.directory / "blast.old"
        shutil.rmtree(old, ignore_errors=True)
        if final.exists():
            final.rename(old)
        work.rename(final)
        shutil.rmtree(old, ignore_errors=True)
    else:
        # A failed build leaves nothing, not even the old databases: they may
        # name rows of an earlier load.
        shutil.rmtree(work, ignore_errors=True)
        shutil.rmtree(final, ignore_errors=True)


def _count(number: int, word: str) -> str:
    return f"{number:,} {word}{'' if number == 1 else 's'}"


def _makeblastdb(binary: str, dbtype: str, directory: Path, name: str, title: str,
                 records: Iterator[str]) -> int:
    """Stream ``records`` into makeblastdb; the number of sequences it added.

    Its output goes to files, not pipes: makeblastdb can write more warnings
    than a pipe holds while it is still being fed, and both would wait."""
    with (tempfile.TemporaryFile("w+") as out, tempfile.TemporaryFile("w+") as err):
        process = subprocess.Popen(  # noqa: S603 - a fixed argument list; names written here  # nosec B603
            [binary, "-in", "-", "-dbtype", dbtype, "-out", name, "-title", title,
             "-blastdb_version", "5"],
            stdin=subprocess.PIPE, stdout=out, stderr=err, text=True, cwd=directory)
        try:
            for record in records:
                process.stdin.write(record)
        except BrokenPipeError:
            pass  # makeblastdb stopped reading; its exit status says why
        finally:
            with contextlib.suppress(BrokenPipeError):
                process.stdin.close()
        process.wait()
        out.seek(0)
        err.seek(0)
        stdout, stderr = out.read(), err.read()
    if process.returncode != 0:
        said = " ".join(line.strip() for line in stdout.splitlines()[-3:] if line.strip())
        raise OSError(f"makeblastdb failed for {name}: {stderr.strip()[:500]}"
                      + (f" (its output ended: {said[:300]})" if said else ""))
    added = re.findall(r"added (\d+) sequences", stdout)
    return sum(int(count) for count in added)
