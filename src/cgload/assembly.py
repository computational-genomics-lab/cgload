"""Registering an assembly: rows and files, together or not at all.

Milestone 2 builds this as an API rather than a command. ``load`` stays a
stub exiting 2 until milestone 4 wires annotation into it, because D-011 says
a command either does its whole job or says plainly that it does not. A
half-built ``load`` that stored sequence and silently ignored annotation is
exactly the silent partial success P1 exists to prevent.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from sqlalchemy import Engine, insert, select

from cgload import __version__
from cgload.config import OrganismConfig
from cgload.db import vocabulary as vocab
from cgload.db.init import assert_vocabulary_current
from cgload.db.schema import assembly, organism, sequence_region, source_file, strain
from cgload.fasta import index_and_read
from cgload.storage import staged_load

ASSEMBLY_FILENAME = "assembly.fa"


class AssemblyExistsError(RuntimeError):
    """This organism/strain/assembly triple is already in the database."""


@dataclass(frozen=True)
class AssemblyResult:
    assembly_id: int
    region_count: int
    total_length: int
    data_directory: Path
    fasta_sha256: str


def register_assembly(
    engine: Engine,
    config: OrganismConfig,
    data_root: Path,
    *,
    copy: bool = True,
    force: bool = False,
) -> AssemblyResult:
    """Copy the assembly in, index it, and record its regions."""
    assert_vocabulary_current(engine)

    if not config.fasta.exists():
        raise FileNotFoundError(f"assembly not found: {config.fasta}")

    with staged_load(
        engine,
        data_root,
        organism=f"{config.genus}_{config.species}",
        strain=config.strain,
        assembly=config.assembly_name,
        force=force,
    ) as load:
        conn = load.connection

        staged_fasta, fasta_hash = load.place(config.fasta, ASSEMBLY_FILENAME, copy=copy)
        index_path, regions = index_and_read(staged_fasta)

        organism_id = _upsert_organism(conn, config)
        strain_id = _upsert_strain(conn, organism_id, config.strain)

        existing = conn.execute(
            select(assembly.c.assembly_id).where(
                assembly.c.strain_id == strain_id,
                assembly.c.name == config.assembly_name,
                assembly.c.version.is_(config.assembly_version)
                if config.assembly_version is None
                else assembly.c.version == config.assembly_version,
            )
        ).scalar_one_or_none()
        if existing is not None and not force:
            raise AssemblyExistsError(
                f"assembly {config.assembly_name!r} for strain {config.strain!r} is "
                f"already loaded (assembly_id={existing}). Use `cgload remove`, or "
                f"load it under a different version."
            )

        assembly_id = conn.execute(
            insert(assembly).values(
                strain_id=strain_id,
                name=config.assembly_name,
                version=config.assembly_version,
                accession=config.accession,
                fasta_path=load.relative(staged_fasta),
                fasta_sha256=fasta_hash,
                vocabulary_version=vocab.VOCABULARY_VERSION,
                vocabulary_hash=vocab.content_hash(),
                cgload_version=__version__,
                loaded_at=datetime.now(UTC),
            )
        ).inserted_primary_key[0]

        conn.execute(
            insert(source_file).values(
                assembly_id=assembly_id,
                path=str(config.fasta),
                sha256=fasta_hash,
                file_format="fasta",
            )
        )

        conn.execute(
            insert(sequence_region),
            [
                {
                    "assembly_id": assembly_id,
                    "source_id": region.source_id,
                    "length": region.length,
                }
                for region in regions
            ],
        )

        load.manifest.update(
            {
                "cgload_version": __version__,
                "vocabulary_version": vocab.VOCABULARY_VERSION,
                "vocabulary_hash": vocab.content_hash(),
                "organism": f"{config.genus} {config.species}",
                "strain": config.strain,
                "assembly": config.assembly_name,
                "ncbi_taxon_id": config.ncbi_taxon_id,
                "accession": config.accession,
                "loaded_at": datetime.now(UTC).isoformat(),
                "files": [
                    {
                        "role": "assembly",
                        "source": str(config.fasta),
                        "stored": load.relative(staged_fasta),
                        "sha256": fasta_hash,
                        "mode": "copy" if copy else "link",
                    },
                    {"role": "index", "stored": load.relative(index_path)},
                ],
                "regions": len(regions),
                "total_length": sum(r.length for r in regions),
            }
        )

        return AssemblyResult(
            assembly_id=assembly_id,
            region_count=len(regions),
            total_length=sum(r.length for r in regions),
            data_directory=load.target,
            fasta_sha256=fasta_hash,
        )


def _upsert_organism(conn, config: OrganismConfig) -> int:
    found = conn.execute(
        select(organism.c.organism_id).where(
            organism.c.genus == config.genus, organism.c.species == config.species
        )
    ).scalar_one_or_none()
    if found is not None:
        return found
    return conn.execute(
        insert(organism).values(
            genus=config.genus, species=config.species, ncbi_taxon_id=config.ncbi_taxon_id
        )
    ).inserted_primary_key[0]


def _upsert_strain(conn, organism_id: int, name: str) -> int:
    found = conn.execute(
        select(strain.c.strain_id).where(
            strain.c.organism_id == organism_id, strain.c.name == name
        )
    ).scalar_one_or_none()
    if found is not None:
        return found
    return conn.execute(
        insert(strain).values(organism_id=organism_id, name=name)
    ).inserted_primary_key[0]
