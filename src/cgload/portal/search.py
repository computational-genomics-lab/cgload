"""The portal's search index: every name a feature is known by (D-112).

The schema indexes a feature's ID, but not its attribute values (D-101), and
people search for a gene by its name, a synonym, a locus tag, a protein ID or a
word of its product. ``cgload index`` therefore writes one row per name into a
table of its own, ``search_term``, and the portal looks names up there with an
index, on SQLite and PostgreSQL alike.

What is indexed. Every feature except the parts of another feature (exons,
UTRs, introns, start and stop codons), whose names repeat their parent's and
which a user reaches through it. For each such feature: its ID; the whole value
of each attribute in ``NAME_KEYS``, ``XREF_KEYS`` and ``TEXT_KEYS``; for a
cross-reference, also the part after its database prefix (``GeneID:672`` is
found by ``672``); and for free text, each word.

The index is checked after it is built, the way ``load`` checks a load: by
counts that do not go through the code that built it. For each key, the
features carrying that key are counted in the attribute table, and the features
with a term of that key in the index; the two must be equal (``check_index``).
A search that cannot find a feature whose name is in the file is the failure
the check is there to catch.

The tables carry no foreign keys, so a database without the portal never needs
them, and removing an assembly deletes its rows here too
(``delete_index_rows``).
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from sqlalchemy import (
    Column,
    Connection,
    DateTime,
    Engine,
    Index,
    Integer,
    MetaData,
    String,
    Table,
    and_,
    delete,
    func,
    insert,
    inspect,
    intersect,
    select,
)

from cgload import __version__
from cgload.db.schema import (
    NAMING_CONVENTION,
    annotation_run,
    assembly,
    feature,
    feature_attribute,
    protein_annotation,
)
from cgload.query import QueryError, term_variants

portal_metadata = MetaData(naming_convention=NAMING_CONVENTION)

#: PostgreSQL compares text by the database's locale unless told otherwise, and
#: a prefix search written as a range (``term >= 'abc' AND term < 'abc￿'``)
#: is only right under a byte-wise order. SQLite's default order is byte-wise.
_TERM_TYPE = String(200).with_variant(String(200, collation="C"), "postgresql")
TERM_LENGTH = 200

search_term = Table(
    "search_term",
    portal_metadata,
    Column("search_term_id", Integer, primary_key=True),
    Column("assembly_id", Integer, nullable=False),
    # Lower case, at most TERM_LENGTH characters.
    Column("term", _TERM_TYPE, nullable=False),
    # The feature's first row, as `feature_attribute` and `protein_annotation` use.
    Column("feature_id", Integer, nullable=False),
    # The attribute key the term came from ("Name", "Dbxref", "ID", ...), with
    # ":word" or ":suffix" for a term taken from part of the value.
    Column("field", String(64), nullable=False),
)
Index("ix_search_term_assembly_term", search_term.c.assembly_id, search_term.c.term)
Index("ix_search_term_feature", search_term.c.feature_id)

search_index = Table(
    "search_index",
    portal_metadata,
    Column("assembly_id", Integer, primary_key=True, autoincrement=False),
    Column("cgload_version", String(64), nullable=False),
    Column("built_at", DateTime(timezone=True), nullable=False, server_default=func.now()),
    Column("features", Integer, nullable=False),
    Column("terms", Integer, nullable=False),
)

#: Attribute keys whose whole value is a name somebody may type.
NAME_KEYS: tuple[str, ...] = (
    "Name", "gene", "gene_name", "gene_id", "gene_synonym", "locus_tag", "old_locus_tag",
    "protein_id", "transcript_id", "Alias", "standard_name",
)
#: Cross-references: the whole value, and the part after the database prefix.
XREF_KEYS: tuple[str, ...] = ("Dbxref", "db_xref", "Ontology_term")
#: Free text: the whole value, and each word.
TEXT_KEYS: tuple[str, ...] = ("product", "description")
INDEXED_KEYS: tuple[str, ...] = NAME_KEYS + XREF_KEYS + TEXT_KEYS

#: Parts of a feature, found through the feature they belong to.
PART_TYPES: tuple[str, ...] = (
    "exon", "five_prime_UTR", "three_prime_UTR", "intron", "start_codon", "stop_codon",
)

#: Words too common to find anything by.
STOPWORDS = frozenset({
    "a", "an", "and", "as", "at", "by", "for", "from", "in", "is", "of", "on", "or", "the",
    "to", "with",
})
_WORD = re.compile(r"[0-9a-z]+")

#: How many rows go to the database in one INSERT.
_BATCH = 5000


def normal(text: str) -> str:
    """The form a term is stored and looked up in."""
    return text.strip().lower()[:TERM_LENGTH]


def words(text: str) -> list[str]:
    """The words of free text, as the index stores them."""
    return [word for word in _WORD.findall(text.lower())
            if len(word) >= 2 and word not in STOPWORDS]


def terms_of(key: str, value: str) -> Iterator[tuple[str, str]]:
    """(term, field) pairs for one attribute value."""
    whole = normal(value)
    if not whole:
        return
    yield whole, key
    if key in XREF_KEYS and ":" in whole:
        suffix = normal(whole.split(":", 1)[1])
        if suffix:
            yield suffix, f"{key}:suffix"
    elif key in TEXT_KEYS:
        for word in words(whole):
            if word != whole:
                yield word, f"{key}:word"


def _indexable(assembly_id: int):
    return and_(
        feature.c.assembly_id == assembly_id,
        feature.c.segment_index == 0,
        feature.c.feature_type.not_in(PART_TYPES),
    )


def has_index_tables(connection: Connection) -> bool:
    return inspect(connection).has_table("search_term") and inspect(connection).has_table(
        "search_index")


def indexed_assemblies(connection: Connection) -> dict[int, dict]:
    """Assemblies with a search index, and when and by which cgload it was built."""
    if not has_index_tables(connection):
        return {}
    return {
        row.assembly_id: dict(row._mapping)
        for row in connection.execute(select(search_index))
    }


def delete_index_rows(connection: Connection, assembly_id: int) -> None:
    """Delete one assembly's index, inside the caller's transaction."""
    if not has_index_tables(connection):
        return
    connection.execute(delete(search_term).where(search_term.c.assembly_id == assembly_id))
    connection.execute(delete(search_index).where(search_index.c.assembly_id == assembly_id))


@dataclass(frozen=True)
class KeyCount:
    key: str
    in_attributes: int
    in_index: int

    @property
    def ok(self) -> bool:
        return self.in_attributes == self.in_index


@dataclass(frozen=True)
class IndexReport:
    assembly_id: int
    features: int
    terms: int
    counts: tuple[KeyCount, ...]

    @property
    def ok(self) -> bool:
        return all(count.ok for count in self.counts)

    def failures(self) -> list[KeyCount]:
        return [count for count in self.counts if not count.ok]


#: Feature rows read per step while building: a step's rows are read whole
#: before its terms are written, so no cursor stays open across an INSERT.
_WINDOW = 20000


def _term_rows(connection: Connection, assembly_id: int, low: int, high: int) -> list[dict]:
    """The (term, feature, field) rows of the features with row ids in
    ``[low, high]``, without repeats within a feature."""
    window = feature.c.feature_id.between(low, high)
    rows = []
    for row in connection.execute(
        select(feature.c.feature_id, feature.c.source_id)
        .where(_indexable(assembly_id), window, feature.c.source_id.is_not(None))
    ).all():
        term = normal(row.source_id)
        if term:
            rows.append({"assembly_id": assembly_id, "term": term,
                         "feature_id": row.feature_id, "field": "ID"})
    current, seen = None, set()
    for row in connection.execute(
        select(feature_attribute.c.feature_id, feature_attribute.c.key,
               feature_attribute.c.value)
        .join(feature, feature.c.feature_id == feature_attribute.c.feature_id)
        .where(_indexable(assembly_id), feature_attribute.c.feature_id.between(low, high),
               feature_attribute.c.key.in_(INDEXED_KEYS))
        .order_by(feature_attribute.c.feature_id, feature_attribute.c.feature_attribute_id)
    ).all():
        if row.feature_id != current:
            current, seen = row.feature_id, set()
        for term, field in terms_of(row.key, row.value):
            if (term, field) not in seen:
                seen.add((term, field))
                rows.append({"assembly_id": assembly_id, "term": term,
                             "feature_id": row.feature_id, "field": field})
    return rows


def build_index(engine: Engine, assembly_ids: Iterable[int] | None = None) -> list[IndexReport]:
    """Build and check the index of each assembly, one transaction each.

    An assembly whose index fails its check is rolled back to no index at all,
    as a failed load is, and the report says which key fell short.
    """
    portal_metadata.create_all(engine)
    with engine.connect() as connection:
        known = list(connection.execute(
            select(assembly.c.assembly_id).order_by(assembly.c.assembly_id)).scalars())
    wanted = known if assembly_ids is None else list(assembly_ids)
    missing = sorted(set(wanted) - set(known))
    if missing:
        raise LookupError(f"no assembly with id {missing[0]}")
    reports = []
    for assembly_id in wanted:
        with engine.connect() as connection:
            transaction = connection.begin()
            try:
                delete_index_rows(connection, assembly_id)
                low, high = connection.execute(
                    select(func.min(feature.c.feature_id), func.max(feature.c.feature_id))
                    .where(feature.c.assembly_id == assembly_id)).one()
                terms = 0
                if low is not None:
                    for start in range(low, high + 1, _WINDOW):
                        rows = _term_rows(connection, assembly_id, start, start + _WINDOW - 1)
                        for begin in range(0, len(rows), _BATCH):
                            connection.execute(insert(search_term), rows[begin:begin + _BATCH])
                        terms += len(rows)
                report = check_index(connection, assembly_id, terms=terms)
                if report.ok:
                    connection.execute(insert(search_index).values(
                        assembly_id=assembly_id, cgload_version=__version__,
                        features=report.features, terms=terms))
                    transaction.commit()
                else:
                    transaction.rollback()
            except BaseException:
                transaction.rollback()
                raise
        reports.append(report)
    return reports


def check_index(connection: Connection, assembly_id: int, *, terms: int | None = None
                ) -> IndexReport:
    """Count, for each key, the features that carry it and the features the
    index can find by it. Pure SQL counts on both sides, so a fault in how
    terms are made cannot hide itself."""
    counts = []
    with_id = select(func.count(func.distinct(feature.c.feature_id))).where(
        _indexable(assembly_id), feature.c.source_id.is_not(None),
        func.trim(feature.c.source_id) != "")
    counts.append(KeyCount("ID", _count(connection, with_id), _indexed(connection, assembly_id,
                                                                        "ID")))
    for key in INDEXED_KEYS:
        carrying = (
            select(func.count(func.distinct(feature_attribute.c.feature_id)))
            .select_from(feature_attribute.join(
                feature, feature.c.feature_id == feature_attribute.c.feature_id))
            .where(_indexable(assembly_id), feature_attribute.c.key == key,
                   func.trim(feature_attribute.c.value) != "")
        )
        counts.append(KeyCount(key, _count(connection, carrying),
                               _indexed(connection, assembly_id, key)))
    features = _count(connection, select(func.count(func.distinct(search_term.c.feature_id)))
                      .where(search_term.c.assembly_id == assembly_id))
    if terms is None:
        terms = _count(connection, select(func.count()).select_from(search_term)
                       .where(search_term.c.assembly_id == assembly_id))
    return IndexReport(assembly_id, features, terms, tuple(counts))


def _count(connection: Connection, statement) -> int:
    return int(connection.execute(statement).scalar_one() or 0)


def _indexed(connection: Connection, assembly_id: int, field: str) -> int:
    return _count(connection, select(func.count(func.distinct(search_term.c.feature_id)))
                  .where(search_term.c.assembly_id == assembly_id, search_term.c.field == field))


# --------------------------------------------------------------------------
# Searching
# --------------------------------------------------------------------------

#: How a match was made, best first.
EXACT_ID, EXACT_NAME, PREFIX, ALL_WORDS, FUNCTIONAL = range(5)
RANK_LABELS = {
    EXACT_ID: "ID",
    EXACT_NAME: "name",
    PREFIX: "name begins with",
    ALL_WORDS: "every word",
    FUNCTIONAL: "functional annotation",
}
#: A search stops collecting after this many features.
MAX_MATCHES = 1000
MAX_QUERY_LENGTH = 200
#: Above every character, so ``term < prefix + _TOP`` bounds a prefix.
_TOP = "\U0010ffff"


@dataclass(frozen=True)
class Match:
    feature_id: int
    rank: int
    field: str


@dataclass(frozen=True)
class SearchResult:
    query: str
    matches: tuple[Match, ...]
    #: True when MAX_MATCHES stopped the search before it ran out of matches.
    truncated: bool
    #: False when the assembly has no index: only IDs and functional terms are searched.
    indexed: bool


def _where_words_are(connection: Connection, assembly_id: int, tokens: list[str],
                     ids: list[int]) -> list[tuple[int, str]]:
    """Each feature with the attribute keys its words were found in ("product",
    or "description, product" when they are spread over two)."""
    keys: dict[int, set[str]] = {feature_id: set() for feature_id in ids}
    for start in range(0, len(ids), 500):
        rows = connection.execute(
            select(search_term.c.feature_id, search_term.c.field).where(
                search_term.c.assembly_id == assembly_id, search_term.c.term.in_(tokens),
                search_term.c.feature_id.in_(ids[start:start + 500])))
        for feature_id, field in rows:
            keys[feature_id].add(field.split(":", 1)[0])
    return [(feature_id, ", ".join(sorted(keys[feature_id])) or "words") for feature_id in ids]


def search(connection: Connection, assembly_id: int, text: str, *, indexed: bool = True,
           limit: int = MAX_MATCHES) -> SearchResult:
    """Features of one assembly that ``text`` names, best matches first.

    In order: the feature whose ID is exactly ``text``; features with a name or
    cross-reference equal to it, in any case; names that begin with it; with
    several words, features whose text has every one; and features that
    eggNOG-mapper or InterProScan annotated with it as an accession.
    """
    query = text.strip()
    if not query:
        raise QueryError("the search is empty")
    if len(query) > MAX_QUERY_LENGTH:
        raise QueryError(f"a search is at most {MAX_QUERY_LENGTH} characters")
    found: dict[int, Match] = {}

    def add(rows, rank: int, field: str | None = None) -> bool:
        for row in rows:
            feature_id, row_field = (row[0], row[1]) if field is None else (row[0], field)
            if feature_id not in found:
                found[feature_id] = Match(feature_id, rank, row_field)
                if len(found) >= limit:
                    return True
        return False

    exact = select(feature.c.feature_id).where(
        feature.c.assembly_id == assembly_id, feature.c.source_id == query,
        feature.c.segment_index == 0).limit(limit)
    if add(connection.execute(exact), EXACT_ID, "ID"):
        return SearchResult(query, tuple(found.values()), True, indexed)

    lowered = normal(query)
    if indexed:
        same = (select(search_term.c.feature_id, search_term.c.field)
                .where(search_term.c.assembly_id == assembly_id, search_term.c.term == lowered)
                .order_by(search_term.c.feature_id).limit(limit))
        if add(connection.execute(same), EXACT_NAME):
            return SearchResult(query, tuple(found.values()), True, indexed)
        if len(lowered) >= 2:
            begins = (select(search_term.c.feature_id, search_term.c.field)
                      .where(search_term.c.assembly_id == assembly_id,
                             search_term.c.term > lowered,
                             search_term.c.term < lowered + _TOP)
                      .order_by(search_term.c.term, search_term.c.feature_id).limit(limit))
            if add(connection.execute(begins), PREFIX):
                return SearchResult(query, tuple(found.values()), True, indexed)
        tokens = sorted(set(words(lowered)))
        if len(tokens) >= 2:
            each = [select(search_term.c.feature_id).where(
                search_term.c.assembly_id == assembly_id, search_term.c.term == token)
                for token in tokens]
            every = intersect(*each).subquery()
            ids = list(connection.execute(select(every.c.feature_id)
                                          .order_by(every.c.feature_id).limit(limit)).scalars())
            if add(_where_words_are(connection, assembly_id, tokens, ids), ALL_WORDS):
                return SearchResult(query, tuple(found.values()), True, indexed)

    try:
        variants = term_variants(query)
    except QueryError:
        variants = []
    if variants:
        annotated = (
            select(protein_annotation.c.feature_id, protein_annotation.c.analysis)
            .join(annotation_run,
                  annotation_run.c.annotation_run_id == protein_annotation.c.annotation_run_id)
            .where(annotation_run.c.assembly_id == assembly_id,
                   protein_annotation.c.accession.in_(variants))
            .order_by(protein_annotation.c.feature_id).limit(limit)
        )
        if add(connection.execute(annotated), FUNCTIONAL):
            return SearchResult(query, tuple(found.values()), True, indexed)
    return SearchResult(query, tuple(found.values()), False, indexed)
