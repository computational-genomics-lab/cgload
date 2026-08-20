"""Organism configuration, in TOML (D-010).

Shape, matching D-014's per-file array of tables::

    [organism]
    genus = "Phytophthora"
    species = "megasperma"
    ncbi_taxon_id = 4788        # opaque; cgload resolves nothing (D-016)
    strain = "CJ26"             # optional; defaults to "reference"

    [assembly]
    name = "Phyme_CJ26"
    version = "1.0"                    # optional
    accession = "GCF_000149735.1"      # optional; keys the NCBI cross-check
    fasta = "genomes/CJ26.fa"

    [[files.annotation]]        # read from milestone 3
    path = "genomes/CJ26.gff3"

    [[files.functional]]       # eggNOG or InterProScan output (milestone 7)
    path = "annot/CJ26.emapper.annotations"
    # minimum_match_rate = 1.0 # lower it to load a subset deliberately
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field
from pathlib import Path


class ConfigError(ValueError):
    """The config is missing something, or has something it should not."""


#: Extensions that mean GenBank flat file. Advisory only: the format is settled
#: by content, because the M. aeruginosa download arrived named `.gbff` while
#: actually being gzip (D-035).
GENBANK_SUFFIXES = frozenset({".gb", ".gbk", ".gbff", ".genbank"})


@dataclass(frozen=True)
class AnnotationFile:
    path: Path
    dialect: str | None = None
    id_prefix: str | None = None

    @property
    def looks_like_genbank(self) -> bool:
        name = self.path.name.lower().removesuffix(".gz")
        return any(name.endswith(suffix) for suffix in GENBANK_SUFFIXES)


@dataclass(frozen=True)
class FunctionalFileEntry:
    """One eggNOG or InterProScan output to attach after the genome loads."""

    path: Path
    tool: str | None = None
    #: Fraction of query proteins that must match a feature. Default 1.0: a
    #: half-matching functional file is the classic silent partial success, so
    #: loading a subset has to be asked for (D-040).
    minimum_match_rate: float = 1.0


@dataclass(frozen=True)
class OrganismConfig:
    genus: str
    species: str
    strain: str
    assembly_name: str
    fasta: Path
    ncbi_taxon_id: int | None = None
    assembly_version: str | None = None
    accession: str | None = None
    annotation: tuple[AnnotationFile, ...] = field(default_factory=tuple)
    functional: tuple[FunctionalFileEntry, ...] = field(default_factory=tuple)
    source: Path | None = None


def _require(table: dict, key: str, where: str) -> object:
    if key not in table:
        raise ConfigError(f"[{where}] is missing required key {key!r}")
    return table[key]


def load_config(path: Path) -> OrganismConfig:
    path = Path(path)
    if not path.exists():
        raise ConfigError(f"{path} does not exist")
    with path.open("rb") as fh:
        raw = tomllib.load(fh)

    base = path.parent
    org = raw.get("organism")
    asm = raw.get("assembly")
    if not isinstance(org, dict):
        raise ConfigError("missing [organism] table")
    if not isinstance(asm, dict):
        raise ConfigError("missing [assembly] table")

    taxon = org.get("ncbi_taxon_id")
    if taxon is not None and not isinstance(taxon, int):
        raise ConfigError("ncbi_taxon_id must be an integer identifier (D-016)")

    annotation = tuple(
        AnnotationFile(
            path=(base / _require(entry, "path", "[[files.annotation]]")).resolve(),
            dialect=entry.get("dialect"),
            id_prefix=entry.get("id_prefix"),
        )
        for entry in raw.get("files", {}).get("annotation", [])
    )

    functional = tuple(
        FunctionalFileEntry(
            path=(base / _require(entry, "path", "[[files.functional]]")).resolve(),
            tool=entry.get("tool"),
            minimum_match_rate=float(entry.get("minimum_match_rate", 1.0)),
        )
        for entry in raw.get("files", {}).get("functional", [])
    )

    return OrganismConfig(
        genus=str(_require(org, "genus", "organism")),
        species=str(_require(org, "species", "organism")),
        strain=str(org.get("strain") or "reference"),
        ncbi_taxon_id=taxon,
        assembly_name=str(_require(asm, "name", "assembly")),
        assembly_version=(str(asm["version"]) if "version" in asm else None),
        accession=(str(asm["accession"]) if "accession" in asm else None),
        fasta=(base / _require(asm, "fasta", "assembly")).resolve(),
        annotation=annotation,
        functional=functional,
        source=path.resolve(),
    )
