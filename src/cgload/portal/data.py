"""What the portal reads, and from where (D-111).

Every question the portal asks of the database goes through ``cgload.query`` or
the search index (``cgload.portal.search``), so a page shows what ``cgload
query`` would print for the same question, and the waits measured for those
queries in the benchmark are the waits behind a page. Sequence comes from the
assembly in the data directory (``cgload.portal.sequence``), and the load's own
receipt, ``manifest.json``, from beside it.

Nothing here writes. ``read_only_engine`` opens SQLite read-only and asks
PostgreSQL for read-only transactions, so a fault in the portal cannot change
the database it serves.
"""

from __future__ import annotations

import json
import threading
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from sqlalchemy import Connection, Engine, and_, create_engine, event, func, select
from sqlalchemy.engine import make_url

from cgload.assembly import ASSEMBLY_FILENAME
from cgload.db.counts import feature_count_expression
from cgload.db.engine import validated_url
from cgload.db.export_source import StoredFeature, features_where
from cgload.db.schema import (
    annotation_run,
    assembly,
    feature,
    organism,
    protein_annotation,
    sequence_region,
    source_file,
    strain,
)
from cgload.portal import search as search_index
from cgload.portal.sequence import (
    GENETIC_CODES,
    AssemblySequence,
    Piece,
    coding_sequence,
    transcript_order,
    translate,
)
from cgload.query import QueryError, _rows_in_region, by_id, resolve_region

ANNOTATION_FILENAME = "annotation.gff3.gz"
#: A page lists at most this many features of a region; the browser shows all.
REGION_TABLE_LIMIT = 500
#: The longest region a page lists features for, and the longest sequence a
#: download writes; the whole assembly is a file download of its own.
REGION_LIMIT = 5_000_000
#: Features of one family drawn or listed in full; a bigger family is listed
#: up to this many, with a download of the whole.
FAMILY_LIMIT = 5000
#: Types that are parts of a transcript, drawn inside it.
EXON_TYPES = ("exon", "five_prime_UTR", "three_prime_UTR")
TYPE_ORDER = {"gene": 0, "pseudogene": 0, "mRNA": 1, "ncRNA": 1, "tRNA": 1, "rRNA": 1,
              "CDS": 2, "exon": 3}


def read_only_engine(db_url: str, *, from_environment: bool = False) -> Engine:
    """An engine that cannot write: SQLite opened in read-only mode with
    ``query_only`` on, and PostgreSQL in read-only transactions."""
    url = validated_url(db_url, from_environment=from_environment)
    backend = url.get_backend_name()
    if backend == "sqlite":
        path = Path(url.database or "")
        if not url.database or url.database == ":memory:" or not path.exists():
            raise FileNotFoundError(f"no SQLite database at {path.resolve()}")
        # A file: URI with mode=ro, so SQLite itself refuses to write.
        engine = create_engine(
            make_url("sqlite:///" + path.resolve().as_uri()).set(
                query={"mode": "ro", "uri": "true"}),
            future=True,
            connect_args={"check_same_thread": False},
        )

        @event.listens_for(engine, "connect")
        def _read_only(dbapi_connection, _record):  # pragma: no cover - trivial
            cursor = dbapi_connection.cursor()
            cursor.execute("PRAGMA query_only=ON")
            cursor.close()

        return engine
    if backend == "postgresql":
        return create_engine(
            url, future=True, pool_pre_ping=True,
            connect_args={"options": "-c default_transaction_read_only=on"},
        )
    raise ValueError(f"the portal serves SQLite and PostgreSQL databases, not {backend}")


@dataclass
class Genome:
    """One loaded assembly, and where its files are."""

    assembly_id: int
    genus: str
    species: str
    strain: str
    name: str
    version: str | None
    accession: str | None
    cgload_version: str
    loaded_at: datetime | None
    directory: Path | None
    fasta: Path | None
    circular: frozenset[str]
    sequences: int
    indexed: dict | None = None
    features: int | None = None
    rows: int | None = None

    @property
    def organism(self) -> str:
        return f"{self.genus} {self.species}"

    @property
    def label(self) -> str:
        parts = [self.organism]
        if self.strain and self.strain != "reference":
            parts.append(self.strain)
        parts.append(self.name)
        return " · ".join(parts)

    def path(self, name: str) -> Path | None:
        """A file of this assembly's directory that the portal serves, if present."""
        if self.directory is None:
            return None
        allowed = {
            "assembly.fa": self.fasta,
            "assembly.fa.fai": None if self.fasta is None else Path(str(self.fasta) + ".fai"),
            ANNOTATION_FILENAME: self.directory / ANNOTATION_FILENAME,
            ANNOTATION_FILENAME + ".tbi": self.directory / (ANNOTATION_FILENAME + ".tbi"),
        }
        path = allowed.get(name)
        return path if path is not None and path.is_file() else None

    def receipt(self) -> dict | None:
        """The load's manifest: files, checksums, checks and repairs."""
        if self.directory is None:
            return None
        manifest = self.directory / "manifest.json"
        try:
            return json.loads(manifest.read_text())
        except (OSError, ValueError):
            return None


@dataclass
class FamilyView:
    """A feature with the whole family it belongs to."""

    genome: Genome
    selected: StoredFeature
    top: StoredFeature
    members: list[StoredFeature]
    children: dict[int, list[StoredFeature]]
    truncated: bool
    annotations: dict[int, list[dict]] = field(default_factory=dict)

    def depth_first(self) -> list[tuple[int, StoredFeature]]:
        """(depth, feature) from the top down, children in position order."""
        ordered: list[tuple[int, StoredFeature]] = []
        seen: set[int] = set()
        stack = [(0, self.top)]
        while stack:
            depth, item = stack.pop()
            if item.row_id in seen:
                continue
            seen.add(item.row_id)
            ordered.append((depth, item))
            stack.extend((depth + 1, child)
                         for child in reversed(self.children.get(item.row_id, [])))
        ordered.extend((0, item) for item in self.members if item.row_id not in seen)
        return ordered


class PortalData:
    """The portal's read-only view of one database and its data directory."""

    def __init__(self, engine: Engine, *, data_dir: Path | None = None) -> None:
        self.engine = engine
        self.data_dir = Path(data_dir).resolve() if data_dir is not None else None
        self._sequences: dict[int, AssemblySequence] = {}
        self._lock = threading.Lock()
        self._counts: dict[int, tuple[int, int]] = {}
        self._genomes: dict[int, Genome] = {}
        self.refresh()

    # -- genomes -----------------------------------------------------------

    def refresh(self) -> None:
        with self.engine.connect() as connection:
            indexed = search_index.indexed_assemblies(connection)
            rows = connection.execute(
                select(assembly, organism.c.genus, organism.c.species,
                       strain.c.name.label("strain_name"))
                .select_from(assembly.join(strain).join(organism))
                .order_by(organism.c.genus, organism.c.species, assembly.c.assembly_id)
            ).all()
            regions = connection.execute(
                select(sequence_region.c.assembly_id, sequence_region.c.source_id,
                       sequence_region.c.circular)
            ).all()
        sequences: dict[int, int] = defaultdict(int)
        circular: dict[int, set[str]] = defaultdict(set)
        for region in regions:
            sequences[region.assembly_id] += 1
            if region.circular:
                circular[region.assembly_id].add(region.source_id)
        genomes = {}
        for row in rows:
            directory, fasta = self._files(row.data_root, row.fasta_path)
            genomes[row.assembly_id] = Genome(
                assembly_id=row.assembly_id, genus=row.genus, species=row.species,
                strain=row.strain_name, name=row.name, version=row.version,
                accession=row.accession, cgload_version=row.cgload_version,
                loaded_at=row.loaded_at, directory=directory, fasta=fasta,
                circular=frozenset(circular[row.assembly_id]),
                sequences=sequences[row.assembly_id], indexed=indexed.get(row.assembly_id),
            )
            if row.assembly_id in self._counts:
                genomes[row.assembly_id].features, genomes[row.assembly_id].rows = (
                    self._counts[row.assembly_id])
        self._genomes = genomes

    def _files(self, data_root: str | None, fasta_path: str | None
               ) -> tuple[Path | None, Path | None]:
        if not fasta_path:
            return None, None
        root = self.data_dir or (Path(data_root) if data_root else None)
        if root is None:
            return None, None
        fasta = (root / fasta_path).resolve()
        if not fasta.is_file():
            return None, None
        return fasta.parent, fasta

    def count_features(self) -> None:
        """Count features and rows of every assembly once; slow on large
        databases, so the server runs it in the background at start-up."""
        with self.engine.connect() as connection:
            for row in connection.execute(
                select(feature.c.assembly_id, feature_count_expression().label("features"),
                       func.count().label("rows"))
                .group_by(feature.c.assembly_id)
            ):
                self._counts[row.assembly_id] = (row.features, row.rows)
                genome = self._genomes.get(row.assembly_id)
                if genome is not None:
                    genome.features, genome.rows = row.features, row.rows

    def genomes(self) -> list[Genome]:
        return list(self._genomes.values())

    def genome(self, assembly_id: int) -> Genome:
        try:
            return self._genomes[assembly_id]
        except KeyError:
            raise LookupError(f"no assembly with id {assembly_id}") from None

    def sequence(self, genome: Genome) -> AssemblySequence:
        if genome.fasta is None:
            raise FileNotFoundError(
                f"the assembly file of {genome.label} is not where the database says")
        with self._lock:
            handle = self._sequences.get(genome.assembly_id)
            if handle is None:
                handle = AssemblySequence(genome.fasta, circular=genome.circular)
                self._sequences[genome.assembly_id] = handle
            return handle

    def close(self) -> None:
        for handle in self._sequences.values():
            handle.close()
        self.engine.dispose()

    def counts_by_type(self, assembly_id: int) -> list[tuple[str, int]]:
        with self.engine.connect() as connection:
            return [(row[0], row[1]) for row in connection.execute(
                select(feature.c.feature_type, feature_count_expression())
                .where(feature.c.assembly_id == assembly_id)
                .group_by(feature.c.feature_type)
                .order_by(feature_count_expression().desc()))]

    def longest_sequences(self, assembly_id: int, limit: int = 25) -> list[tuple[str, int]]:
        with self.engine.connect() as connection:
            return [(row.source_id, row.length) for row in connection.execute(
                select(sequence_region.c.source_id, sequence_region.c.length)
                .where(sequence_region.c.assembly_id == assembly_id)
                .order_by(sequence_region.c.length.desc(), sequence_region.c.source_id)
                .limit(limit))]

    def analyses(self, assembly_id: int) -> list[dict]:
        with self.engine.connect() as connection:
            return [dict(row._mapping) for row in connection.execute(
                select(annotation_run.c.tool, annotation_run.c.tool_version,
                       annotation_run.c.query_count, annotation_run.c.matched_count,
                       annotation_run.c.unmatched_count)
                .where(annotation_run.c.assembly_id == assembly_id)
                .order_by(annotation_run.c.annotation_run_id))]

    def source_files(self, assembly_id: int) -> list[dict]:
        """The files a load read, by name. The database keeps the full path on
        the machine that loaded them; a public page shows only the name."""
        with self.engine.connect() as connection:
            rows = connection.execute(
                select(source_file.c.path, source_file.c.sha256, source_file.c.file_format,
                       source_file.c.detected_dialect)
                .where(source_file.c.assembly_id == assembly_id)
                .order_by(source_file.c.source_file_id)).all()
        return [{"name": file_name(row.path), "sha256": row.sha256,
                 "file_format": row.file_format, "detected_dialect": row.detected_dialect}
                for row in rows]

    # -- search ------------------------------------------------------------

    def search(self, genome: Genome, text: str, *, page: int = 1, per_page: int = 50,
               limit: int = search_index.MAX_MATCHES):
        """One page of matches, best first, each with the top of its family.

        Matches are ordered by how they matched, then genes before transcripts
        before CDS, then by position; only the page shown is read in full.
        """
        with self.engine.connect() as connection:
            result = search_index.search(connection, genome.assembly_id, text,
                                         indexed=genome.indexed is not None, limit=limit)
            ids = [match.feature_id for match in result.matches]
            placed: dict[int, tuple] = {}
            for begin in range(0, len(ids), 400):
                for row in connection.execute(
                    select(feature.c.feature_id, feature.c.feature_type,
                           sequence_region.c.source_id, feature.c.start)
                    .join(sequence_region,
                          sequence_region.c.sequence_region_id == feature.c.sequence_region_id)
                    .where(feature.c.feature_id.in_(ids[begin : begin + 400]))
                ):
                    placed[row.feature_id] = (TYPE_ORDER.get(row.feature_type, 9),
                                              row.source_id, row.start)
            ordered = sorted((match for match in result.matches if match.feature_id in placed),
                             key=lambda match: (match.rank, placed[match.feature_id]))
            first = (max(page, 1) - 1) * per_page
            shown = ordered[first : first + per_page]
            features = self._first_rows(connection, genome.assembly_id,
                                        [match.feature_id for match in shown])
            tops = _tops(connection, genome.assembly_id, list(features.values()))
        rows = [(match, features[match.feature_id], tops.get(match.feature_id))
                for match in shown if match.feature_id in features]
        return result, rows, len(ordered)

    @staticmethod
    def _first_rows(connection: Connection, assembly_id: int, row_ids: list[int]
                    ) -> dict[int, StoredFeature]:
        """Whole features, every row of each, by the row id of their first row."""
        return _whole_features(connection, assembly_id, row_ids)

    # -- one feature and its family ------------------------------------------

    def row_of_id(self, genome: Genome, source_id: str) -> int:
        with self.engine.connect() as connection:
            row = connection.execute(
                select(feature.c.feature_id).where(
                    feature.c.assembly_id == genome.assembly_id,
                    feature.c.source_id == source_id, feature.c.segment_index == 0)
                .order_by(feature.c.feature_id).limit(1)).first()
        if row is None:
            raise LookupError(f"no feature with ID {source_id!r} in {genome.label}")
        return row.feature_id

    def family(self, genome: Genome, row_id: int) -> FamilyView:
        with self.engine.connect() as connection:
            row = connection.execute(
                select(feature.c.feature_id, feature.c.source_id, feature.c.parent_id,
                       feature.c.assembly_id)
                .where(feature.c.feature_id == row_id)).first()
            if row is None or row.assembly_id != genome.assembly_id:
                raise LookupError(f"no feature {row_id} in {genome.label}")
            if row.source_id is not None:
                # Any row of a feature written on several rows leads to its first row.
                row_id = connection.execute(
                    select(func.min(feature.c.feature_id)).where(
                        feature.c.assembly_id == genome.assembly_id,
                        feature.c.source_id == row.source_id,
                        feature.c.segment_index == 0)).scalar_one() or row_id
            top_row = _top_row(connection, row_id)
            top = connection.execute(
                select(feature.c.source_id).where(feature.c.feature_id == top_row)).first()
            if top is not None and top.source_id is not None:
                members = by_id(connection, genome.assembly_id, top.source_id,
                                descendants=True).features
            else:
                members = features_where(connection, feature.c.feature_id == top_row)
            truncated = len(members) > FAMILY_LIMIT
            by_row = {item.row_id: item for item in members}
            selected = by_row.get(row_id)
            if selected is None:  # an anonymous row below an anonymous top
                selected = features_where(connection, feature.c.feature_id == row_id)[0]
                members.append(selected)
                by_row[row_id] = selected
            annotations = _annotations(connection, list(by_row))
        top_item = by_row.get(top_row, selected)
        children: dict[int, list[StoredFeature]] = defaultdict(list)
        for item in members:
            if item.parent_row_id is not None and item.parent_row_id in by_row:
                children[item.parent_row_id].append(item)
        for kids in children.values():
            kids.sort(key=lambda item: (TYPE_ORDER.get(item.feature_type, 9), item.start,
                                        item.end, item.source_id or ""))
        return FamilyView(genome, selected, top_item, members, dict(children), truncated,
                          annotations)

    # -- regions -------------------------------------------------------------

    def region(self, genome: Genome, text: str):
        """The region asked for, and the top-level features with a row in it.

        A parent's rows enclose its children's (``load`` checks it), so the
        top-level features overlapping a region are every family that has a
        row there; their members are a click away.
        """
        with self.engine.connect() as connection:
            region = resolve_region(connection, genome.assembly_id, text)
            if region.end - region.start + 1 > REGION_LIMIT:
                return region, [], True
            rows = _rows_in_region(region)
            tops = select(feature.c.feature_id).where(
                feature.c.feature_id.in_(rows), feature.c.parent_id.is_(None))
            ids = sorted(set(connection.execute(tops.limit(REGION_TABLE_LIMIT * 4))
                             .scalars()))
            features = _whole_features(connection, genome.assembly_id, ids)
        listed = sorted({item.row_id: item for item in features.values()}.values(),
                        key=lambda item: (item.start, item.end, item.source_id or ""))
        return region, listed[:REGION_TABLE_LIMIT], len(listed) > REGION_TABLE_LIMIT

    # -- sequences -----------------------------------------------------------

    def sequences_of(self, family: FamilyView, item: StoredFeature) -> dict[str, str]:
        """The sequences a page offers for one feature of a family.

        ``genomic``: the feature's whole span, on its strand. ``spliced``: its
        exons (or, without exons, its own rows) joined in transcript order.
        ``cds``: for a CDS, its rows from the first whole codon; ``protein``: that
        translated with the file's ``transl_table``, or the standard code.
        """
        handle = self.sequence(family.genome)
        result: dict[str, str] = {}
        whole = Piece(item.start, item.end, "-" if item.strand == "-" else "+")
        result["genomic"] = handle.piece(item.seqid, whole)
        if item.feature_type == "CDS":
            with self.engine.connect() as connection:
                lines = file_order(connection, item)
            ordered = transcript_order(_pieces(item), lines)
            result["cds"] = coding_sequence(handle, item.seqid, ordered)
            result["protein"] = translate(
                result["cds"], genetic_code(item),
                from_start=(ordered[0].phase or 0) == 0 and not _partial_start(item))
            return result
        exons = [child for child in family.children.get(item.row_id, [])
                 if child.feature_type == "exon"]
        pieces = [piece for exon in exons for piece in _pieces(exon)] or _pieces(item)
        if len(pieces) > 1 or exons:
            result["spliced"] = handle.spliced(item.seqid, transcript_order(pieces))
        return result


def display_name(item: StoredFeature | None) -> str:
    """What a feature is called on a page: its Name; for a gene, its gene
    symbol or locus tag; otherwise its own ID, which a transcript or a CDS has
    and which its gene symbol would not tell apart from its siblings."""
    if item is None:
        return ""
    names = item.attributes.get("Name")
    if names:
        return names[0]
    keys = ("gene", "locus_tag")
    if item.feature_type in ("gene", "pseudogene"):
        for key in keys:
            if item.attributes.get(key):
                return item.attributes[key][0]
    if item.source_id:
        return item.source_id
    for key in keys:
        if item.attributes.get(key):
            return item.attributes[key][0]
    return f"{item.source_type} {item.seqid}:{item.start:,}"


def genetic_code(item: StoredFeature) -> int:
    """The NCBI table a CDS names in ``transl_table``, else the standard code."""
    for value in item.attributes.get("transl_table", []):
        try:
            table = int(value)
        except ValueError:
            continue
        if table in GENETIC_CODES:
            return table
    return 1


def file_name(path: str) -> str:
    """The last part of a path written on any system."""
    return path.replace("\\", "/").rstrip("/").rsplit("/", 1)[-1]


def exceptions(item: StoredFeature) -> list[str]:
    """What the file says a plain translation of this CDS does not account for."""
    notes = []
    for key in ("transl_except", "exception", "ribosomal_slippage"):
        for value in item.attributes.get(key, []):
            notes.append(f"{key}={value}")
    return notes


def _partial_start(item: StoredFeature) -> bool:
    start_range = item.attributes.get("start_range", [])
    end_range = item.attributes.get("end_range", [])
    # NCBI marks a CDS incomplete at its 5' end with start_range=.,N on the plus
    # strand and end_range=N,. on the minus strand.
    if item.strand == "-":
        return any(value.endswith(",.") for value in end_range)
    return any(value.startswith(".,") for value in start_range)


def _pieces(item: StoredFeature) -> list[Piece]:
    return [Piece(segment.start, segment.end, segment.strand or item.strand, segment.phase)
            for segment in item.segments]


def _top_row(connection: Connection, row_id: int) -> int:
    """The first row at the top of a feature's parent chain."""
    seen = {row_id}
    current = row_id
    while True:
        parent = connection.execute(
            select(feature.c.parent_id).where(feature.c.feature_id == current)
        ).scalar_one_or_none()
        if parent is None or parent in seen:
            return current
        seen.add(parent)
        current = parent


def _whole_features(connection: Connection, assembly_id: int, row_ids: list[int]
                    ) -> dict[int, StoredFeature]:
    """Every row of each feature whose first row is in ``row_ids``, by that row."""
    wanted = set(row_ids)
    found: dict[int, StoredFeature] = {}
    ordered = sorted(wanted)
    for begin in range(0, len(ordered), 400):
        batch = ordered[begin : begin + 400]
        names = sorted({name for (name,) in connection.execute(
            select(feature.c.source_id).where(feature.c.feature_id.in_(batch),
                                               feature.c.source_id.is_not(None)))})
        for start in range(0, len(names), 400):
            for item in features_where(connection, and_(
                    feature.c.assembly_id == assembly_id,
                    feature.c.source_id.in_(names[start : start + 400]))):
                if item.row_id in wanted:
                    found[item.row_id] = item
        for item in features_where(connection, and_(feature.c.feature_id.in_(batch),
                                                    feature.c.source_id.is_(None))):
            found[item.row_id] = item
    return found


def _tops(connection: Connection, assembly_id: int, items: list[StoredFeature]
          ) -> dict[int, StoredFeature]:
    """Each feature's top-level feature, walked a level at a time for all at once."""
    parent_of = {item.row_id: item.parent_row_id for item in items}
    missing = {parent for parent in parent_of.values()
               if parent is not None and parent not in parent_of}
    while missing:
        wanted = sorted(missing)
        for begin in range(0, len(wanted), 400):
            for row in connection.execute(
                select(feature.c.feature_id, feature.c.parent_id)
                .where(feature.c.feature_id.in_(wanted[begin : begin + 400]))
            ):
                parent_of[row.feature_id] = row.parent_id
        for row_id in wanted:
            parent_of.setdefault(row_id, None)
        missing = {parent for parent in parent_of.values()
                   if parent is not None and parent not in parent_of}
    top_of: dict[int, int] = {}
    for item in items:
        node, seen = item.row_id, set()
        while parent_of.get(node) is not None and node not in seen:
            seen.add(node)
            node = parent_of[node]
        top_of[item.row_id] = node
    tops = _whole_features(connection, assembly_id, sorted(set(top_of.values())))
    return {row: tops[top] for row, top in top_of.items() if top in tops}


def file_order(connection: Connection, item: StoredFeature) -> list[int] | None:
    """Each row's input line, in the order of ``item.segments``, for a feature
    whose rows lie on both strands; None for any other feature."""
    if len({segment.strand for segment in item.segments}) < 2:
        return None
    lines = {
        (row.start, row.end, row.strand): row.source_line
        for row in connection.execute(
            select(feature.c.start, feature.c.end, feature.c.strand, feature.c.source_line)
            .where(feature.c.feature_id.in_([segment.line_number
                                             for segment in item.segments])))
    }
    order = [lines.get((segment.start, segment.end, segment.strand)) for segment in item.segments]
    return None if any(line is None for line in order) else order


def _annotations(connection: Connection, row_ids: list[int]) -> dict[int, list[dict]]:
    found: dict[int, list[dict]] = defaultdict(list)
    for begin in range(0, len(row_ids), 400):
        batch = row_ids[begin : begin + 400]
        for row in connection.execute(
            select(protein_annotation.c.feature_id, protein_annotation.c.analysis,
                   protein_annotation.c.accession, protein_annotation.c.description,
                   protein_annotation.c.protein_start, protein_annotation.c.protein_end,
                   protein_annotation.c.evalue, annotation_run.c.tool)
            .join(annotation_run,
                  annotation_run.c.annotation_run_id == protein_annotation.c.annotation_run_id)
            .where(protein_annotation.c.feature_id.in_(batch))
            .order_by(protein_annotation.c.feature_id, annotation_run.c.tool,
                      protein_annotation.c.analysis, protein_annotation.c.accession)
        ):
            found[row.feature_id].append(dict(row._mapping))
    return dict(found)


__all__ = [
    "ANNOTATION_FILENAME", "ASSEMBLY_FILENAME", "FamilyView", "Genome", "PortalData",
    "QueryError", "exceptions", "genetic_code", "read_only_engine",
]
