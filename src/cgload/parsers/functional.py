"""eggNOG-mapper and InterProScan output (milestone 7).

Deliberately not dialect profiles (D-013): these are per-protein tables, not
feature graphs. They have no parents, no segments and no coordinates on the
genome, so the profile abstraction has nothing to offer them and forcing them
into it would bend it for every real dialect.

What they do have in common with everything else here: an identifier that has to
be matched against something already loaded, and a match rate that must be
reported rather than assumed.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path


class FunctionalParseError(ValueError):
    """A functional annotation file that cannot be read. Names file and line."""

    def __init__(self, path: Path, line_number: int, message: str) -> None:
        super().__init__(f"{path.name}:{line_number}: {message}")
        self.path = path
        self.line_number = line_number


@dataclass(frozen=True)
class Hit:
    """One assignment: this query protein has this accession."""

    query_id: str
    analysis: str
    accession: str
    description: str | None = None
    protein_start: int | None = None
    protein_end: int | None = None
    score: float | None = None
    evalue: float | None = None


@dataclass
class FunctionalFile:
    tool: str
    tool_version: str | None
    hits: list[Hit]
    #: Every distinct query in the file, including those that produced no hit.
    #: A query with no annotation is still a query the tool was asked about, and
    #: the match rate must be computed against all of them.
    queries: set[str]


# --------------------------------------------------------------------------
# eggNOG-mapper
# --------------------------------------------------------------------------

#: Columns of `*.emapper.annotations` that carry accessions worth a row, mapped
#: to the analysis name recorded against each. Column order is read from the
#: header rather than assumed, because it has changed between eggNOG versions and
#: a positional reader would silently mis-assign every field after the change.
EGGNOG_TERM_COLUMNS: dict[str, str] = {
    # GO terms dominate every other analysis, and they are not independent
    # findings: eggNOG emits the full inherited set, so a protein annotated as a
    # specific hydrolase also carries "hydrolase activity", "catalytic activity"
    # and every ancestor up the ontology. Measured on real emapper-2.1.4 output:
    # 286,671 of 371,409 stored findings are GO terms, 89.7 per annotated protein,
    # 570 on the deepest one. Anything that reports a count to a user -- `stats`,
    # `query`, the paper -- must say "ontology terms including inherited
    # ancestors" rather than implying 90 independent functions (D-045).
    "GOs": "GO",
    "EC": "EC",
    "KEGG_ko": "KEGG_ko",
    "KEGG_Pathway": "KEGG_Pathway",
    "KEGG_Module": "KEGG_Module",
    "KEGG_Reaction": "KEGG_Reaction",
    "BRITE": "BRITE",
    "CAZy": "CAZy",
    "BiGG_Reaction": "BiGG",
    "PFAMs": "Pfam",
    "eggNOG_OGs": "eggNOG_OG",
    "COG_category": "COG",
}

#: Placeholders eggNOG writes for "nothing here". Treating any of these as an
#: accession would create thousands of rows for a term named "-".
EMPTY_VALUES = frozenset({"", "-", "NA", "na", "None", "."})


def parse_eggnog(path: Path) -> FunctionalFile:
    """Read an `.emapper.annotations` file.

    The header line begins with `#query`; earlier `##` lines carry the version.
    Values within a column are comma-separated, so one protein yields many rows.
    """
    hits: list[Hit] = []
    queries: set[str] = set()
    header: list[str] | None = None
    version: str | None = None

    for line_number, line in _lines(path):
        if line.startswith("##"):
            if "emapper" in line and version is None:
                version = line.lstrip("#").strip()
            continue
        if line.startswith("#"):
            header = line.lstrip("#").rstrip().split("\t")
            header[0] = header[0].strip()
            continue

        if header is None:
            raise FunctionalParseError(
                path,
                line_number,
                "data appears before the '#query' header line, so the columns "
                "cannot be identified. eggNOG column order has changed between "
                "versions and cgload will not guess it.",
            )

        columns = line.rstrip("\n").split("\t")
        if len(columns) != len(header):
            raise FunctionalParseError(
                path,
                line_number,
                f"{len(columns)} columns but the header declares {len(header)}",
            )

        row = dict(zip(header, columns, strict=True))
        query = row.get("query") or row.get("query_name")
        if not query:
            raise FunctionalParseError(path, line_number, "no query identifier")
        queries.add(query)

        evalue = _float(row.get("evalue"))
        score = _float(row.get("score"))
        description = _clean(row.get("Description"))

        for column, analysis in EGGNOG_TERM_COLUMNS.items():
            raw = _clean(row.get(column))
            if raw is None:
                continue
            for accession in _split_terms(raw, analysis):
                hits.append(
                    Hit(
                        query_id=query,
                        analysis=analysis,
                        accession=accession,
                        description=description if analysis == "eggNOG_OG" else None,
                        score=score if analysis == "eggNOG_OG" else None,
                        evalue=evalue if analysis == "eggNOG_OG" else None,
                    )
                )

    if header is None:
        raise FunctionalParseError(path, 0, "no '#query' header line found")

    return FunctionalFile(tool="eggnog", tool_version=version, hits=hits, queries=queries)


#: KEGG publishes every pathway under two identifiers -- `ko00362` (the
#: ortholog reference pathway) and `map00362` (the same pathway's map). eggNOG
#: emits both in `KEGG_Pathway`. Measured on real emapper-2.1.4 output
#: (*Ustilago maydis*, 6,509 proteins): 2,253 rows carry pathway data and in
#: **all 2,253** the `ko` set and the `map` set are identical -- 8,613 of 17,226
#: emitted rows, exactly half, are the same pathway written twice.
#:
#: Storing both doubles every pathway count and forces every downstream query to
#: know the convention. The `ko` form is kept because it is the identifier eggNOG
#: uses in `KEGG_ko` and the one KEGG documents as the ortholog pathway; a `map`
#: entry with no `ko` counterpart is kept rather than dropped, since that case
#: does not occur in the measured file and silently discarding it would be the
#: kind of assumption this project exists to avoid. Recorded as D-031.
def canonical_kegg_pathway(term: str) -> str:
    """Return the stored form of a KEGG pathway identifier.

    KEGG publishes every pathway twice: `ko00362` (ortholog reference pathway)
    and `map00362` (the map). cgload stores the `ko` form (D-031), but
    **`map00362` is the form that appears in KEGG's own URLs and in most
    published papers**, so it is what a user will paste into a search box. Any
    code that accepts a pathway identifier from a human must send it through
    here first, or a pathway that is present in the database returns nothing and
    the user concludes it is missing.

    This lives beside the parser rather than in the query layer because the
    ko/map equivalence is a fact about the source format, and a second copy of
    it in a different module would drift. Recorded as D-032.
    """
    stripped = term.strip()
    if stripped.startswith("map") and stripped[3:].isdigit():
        return "ko" + stripped[3:]
    return stripped


def _dedupe_kegg_pathways(items: list[str]) -> list[str]:
    ko_numbers = {item[2:] for item in items if item.startswith("ko")}
    kept: list[str] = []
    for item in items:
        if item.startswith("map") and item[3:] in ko_numbers:
            continue
        kept.append(item)
    return kept


def _split_terms(raw: str, analysis: str) -> Iterator[str]:
    """Split a multi-valued eggNOG cell.

    `COG_category` is the exception and the reason this is a function rather than
    a `split(",")`: it is a run of single letters with no separator at all, so
    `EGP` is three categories. Splitting it on commas yields one bogus accession
    named `EGP`, and a query for COG `E` would then miss it.
    """
    if analysis == "COG":
        for character in raw.strip():
            if character.isalpha():
                yield character
        return
    parts = [part.strip() for part in raw.split(",")]
    if analysis == "KEGG_Pathway":
        parts = _dedupe_kegg_pathways(parts)
    for cleaned in parts:
        if cleaned and cleaned not in EMPTY_VALUES:
            # eggNOG writes an OG as `296A5@1|root`. The taxonomic scope is kept:
            # `296A5@1|root` and `2RD8W@2759|Eukaryota` are different orthologous
            # groups assigned at different levels, so truncating at `@` would
            # merge distinct assignments. Stored verbatim.
            yield cleaned


# --------------------------------------------------------------------------
# InterProScan
# --------------------------------------------------------------------------

#: The 11 mandatory TSV columns, in order. Columns 12-15 (InterPro accession and
#: description, GO terms, pathways) are present only with the relevant flags, so
#: they are read when there and not required.
INTERPRO_COLUMNS = (
    "protein_accession",
    "md5",
    "length",
    "analysis",
    "signature_accession",
    "signature_description",
    "start",
    "stop",
    "score",
    "status",
    "date",
)


def parse_interproscan(path: Path) -> FunctionalFile:
    """Read an InterProScan TSV.

    Headerless and positional, unlike eggNOG, so the column count is the only
    validation available. One protein yields one row per matched signature, and
    the same signature may match twice at different positions -- which is why
    ``protein_start`` is part of the uniqueness constraint on the table.
    """
    hits: list[Hit] = []
    queries: set[str] = set()

    for line_number, line in _lines(path):
        if line.startswith("#"):
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) < len(INTERPRO_COLUMNS):
            raise FunctionalParseError(
                path,
                line_number,
                f"{len(columns)} columns; InterProScan TSV has at least "
                f"{len(INTERPRO_COLUMNS)}",
            )

        row = dict(zip(INTERPRO_COLUMNS, columns[: len(INTERPRO_COLUMNS)], strict=True))
        extra = columns[len(INTERPRO_COLUMNS) :]
        query = row["protein_accession"].strip()
        if not query:
            raise FunctionalParseError(path, line_number, "no protein accession")
        queries.add(query)

        start = _integer(row["start"], path, line_number, "start")
        stop = _integer(row["stop"], path, line_number, "stop")
        if start is not None and stop is not None and stop < start:
            raise FunctionalParseError(
                path, line_number, f"stop {stop} precedes start {start}"
            )

        hits.append(
            Hit(
                query_id=query,
                analysis=row["analysis"].strip() or "unknown",
                accession=row["signature_accession"].strip(),
                description=_clean(row["signature_description"]),
                protein_start=start,
                protein_end=stop,
                # InterProScan writes the e-value in the score column for
                # HMM-based analyses. Recorded as `evalue` when it looks like one,
                # because reporting 1.2e-40 as a "score" inverts the ordering a
                # user would sort by.
                score=None if _is_evalue(row["score"]) else _float(row["score"]),
                evalue=_float(row["score"]) if _is_evalue(row["score"]) else None,
            )
        )

        # Column 12 is the InterPro entry the signature maps into. A separate
        # assignment, not a property of the signature hit: two signatures from
        # different analyses commonly map to one InterPro entry, and a user asking
        # "which proteins have IPR003599" should get one row per protein.
        if extra:
            interpro = _clean(extra[0])
            if interpro and interpro.startswith("IPR"):
                hits.append(
                    Hit(
                        query_id=query,
                        analysis="InterPro",
                        accession=interpro,
                        description=_clean(extra[1]) if len(extra) > 1 else None,
                        protein_start=start,
                        protein_end=stop,
                    )
                )
        if len(extra) > 2:
            for accession in _split_terms(_clean(extra[2]) or "", "GO"):
                hits.append(
                    Hit(query_id=query, analysis="GO", accession=accession.split("(")[0])
                )

    if not queries:
        raise FunctionalParseError(path, 0, "contains no rows")

    return FunctionalFile(tool="interproscan", tool_version=None, hits=hits, queries=queries)


# --------------------------------------------------------------------------
# Shared
# --------------------------------------------------------------------------


def _lines(path: Path) -> Iterator[tuple[int, str]]:
    with Path(path).open("r", encoding="utf-8", errors="strict") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.strip():
                yield line_number, line


def _clean(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return None if stripped in EMPTY_VALUES else stripped


def _float(value: str | None) -> float | None:
    cleaned = _clean(value)
    if cleaned is None:
        return None
    try:
        return float(cleaned)
    except ValueError:
        return None


def _is_evalue(value: str | None) -> bool:
    number = _float(value)
    return number is not None and 0 < number < 0.01


def _integer(value: str, path: Path, line_number: int, what: str) -> int | None:
    cleaned = _clean(value)
    if cleaned is None:
        return None
    try:
        return int(cleaned)
    except ValueError as exc:
        raise FunctionalParseError(
            path, line_number, f"{what} {cleaned!r} is not an integer"
        ) from exc


def detect_functional_format(path: Path) -> str:
    """Decide by content, never by name (D-035).

    eggNOG output carries a `#query` header; InterProScan is headerless with at
    least 11 tab-separated columns. A file matching neither is refused rather
    than guessed at.
    """
    for _, line in _lines(path):
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            if "query" in line.lower():
                return "eggnog"
            continue
        if len(line.rstrip("\n").split("\t")) >= len(INTERPRO_COLUMNS):
            return "interproscan"
        break
    raise FunctionalParseError(
        path,
        1,
        "not recognised as eggNOG-mapper output (no '#query' header) or as an "
        "InterProScan TSV (fewer than 11 tab-separated columns)",
    )


def parse_functional(path: Path, *, tool: str | None = None) -> FunctionalFile:
    resolved = tool or detect_functional_format(path)
    if resolved == "eggnog":
        return parse_eggnog(path)
    if resolved == "interproscan":
        return parse_interproscan(path)
    raise FunctionalParseError(path, 0, f"unknown tool {resolved!r}")
