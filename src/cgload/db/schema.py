"""The relational schema, as SQLAlchemy Core metadata.

Read this file to know what cgload stores. That is a project requirement,
not a nicety, and it is why this is Core rather than the ORM (D-002).

Conventions that hold everywhere below:

* **Coordinates are 1-based and inclusive**, the GFF3 convention, on the
  forward strand regardless of feature strand. GenBank locations arrive from
  Biopython 0-based half-open and are converted in the parser, never here.
* **``source_*`` columns are verbatim input.** Nothing rewrites them. They are
  what ``verify`` counts.
* **One row per segment.** A discontinuous feature -- a multi-exon CDS sharing
  one ID across its segments -- is N rows with the same ``source_id`` and
  ascending ``segment_index``. The feature-level count that ``verify``
  compares against the source is ``COUNT(DISTINCT source_id)``, never
  ``COUNT(*)`` (D-008).
"""

from __future__ import annotations

from sqlalchemy import (
    CheckConstraint,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Index,
    Integer,
    MetaData,
    String,
    Table,
    Text,
    UniqueConstraint,
    func,
)

from cgload.db.vocabulary import LINKAGE_METHODS

# Explicit constraint naming: SQLite will otherwise emit anonymous constraints
# that no ALTER can address later, and error messages that name nothing.
NAMING_CONVENTION = {
    "ix": "ix_%(table_name)s_%(column_0_N_name)s",
    "uq": "uq_%(table_name)s_%(column_0_N_name)s",
    "ck": "ck_%(table_name)s_%(constraint_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s",
    "pk": "pk_%(table_name)s",
}

metadata = MetaData(naming_convention=NAMING_CONVENTION)

_LINKAGE_IN = ", ".join(f"'{m}'" for m in LINKAGE_METHODS)


# --------------------------------------------------------------------------
# Taxonomy. D-016: ncbi_taxon_id is an opaque integer supplied by the user.
# cgload resolves nothing, fetches nothing and vendors nothing.
# --------------------------------------------------------------------------

organism = Table(
    "organism",
    metadata,
    Column("organism_id", Integer, primary_key=True),
    Column("ncbi_taxon_id", Integer, nullable=True),
    Column("genus", String(128), nullable=False),
    Column("species", String(128), nullable=False),
    Column("common_name", String(256), nullable=True),
    UniqueConstraint("genus", "species", name="uq_organism_genus_species"),
)

# A strain always exists. Where the source names none, the loader creates one
# named 'reference', so downstream code never branches on NULL and `stats`
# never has a category that means two different things.
strain = Table(
    "strain",
    metadata,
    Column("strain_id", Integer, primary_key=True),
    Column("organism_id", Integer, ForeignKey("organism.organism_id"), nullable=False),
    Column("name", String(256), nullable=False),
    UniqueConstraint("organism_id", "name", name="uq_strain_organism_name"),
)


# --------------------------------------------------------------------------
# Assembly. Sequence lives on disk (D-004); this holds the reference and the
# provenance needed to detect that the files moved underneath us.
# --------------------------------------------------------------------------

assembly = Table(
    "assembly",
    metadata,
    Column("assembly_id", Integer, primary_key=True),
    Column("strain_id", Integer, ForeignKey("strain.strain_id"), nullable=False),
    Column("name", String(256), nullable=False),
    Column("version", String(64), nullable=True),
    # The NCBI assembly accession, e.g. GCF_000149735.1. Opaque, user-supplied,
    # never resolved (D-016). `verify`'s optional NCBI cross-check keys on this:
    # without it stored counts have no published figure to be compared against.
    Column("accession", String(64), nullable=True),
    # Populated at milestone 2. Relative to the data directory, never absolute:
    # an absolute path makes the database non-portable, which defeats D-002.
    Column("fasta_path", Text, nullable=True),
    Column("fasta_sha256", String(64), nullable=True),
    # Which vocabulary this assembly's feature_type values were derived under.
    # A load against a different hash is refused until `cgload remap` runs.
    Column("vocabulary_version", Integer, nullable=False),
    Column("vocabulary_hash", String(64), nullable=False),
    Column("cgload_version", String(64), nullable=False),
    Column("loaded_at", DateTime(timezone=True), nullable=False, server_default=func.now()),
    UniqueConstraint("strain_id", "name", "version", name="uq_assembly_strain_name_version"),
)

# One row per input file (D-014). Provenance is per file because an organism's
# annotation may arrive as several files of different origin.
source_file = Table(
    "source_file",
    metadata,
    Column("source_file_id", Integer, primary_key=True),
    Column("assembly_id", Integer, ForeignKey("assembly.assembly_id"), nullable=False),
    Column("path", Text, nullable=False),
    Column("sha256", String(64), nullable=False),
    Column("file_format", String(32), nullable=False),  # gff3 | gtf | genbank | fasta
    Column("detected_dialect", String(64), nullable=True),
    Column("detection_confidence", String(16), nullable=True),  # high | low | overridden
    Column("id_prefix", String(64), nullable=True),
    UniqueConstraint("assembly_id", "path", name="uq_source_file_assembly_path"),
)

# Scaffolds, contigs, chromosomes. `length` is authoritative from the .fai.
sequence_region = Table(
    "sequence_region",
    metadata,
    Column("sequence_region_id", Integer, primary_key=True),
    Column("assembly_id", Integer, ForeignKey("assembly.assembly_id"), nullable=False),
    Column("source_id", String(256), nullable=False),  # column 1 of GFF3, verbatim
    Column("length", Integer, nullable=False),
    Column("region_type", String(64), nullable=True),  # chromosome | scaffold | contig
    UniqueConstraint("assembly_id", "source_id", name="uq_sequence_region_assembly_source_id"),
    CheckConstraint("length > 0", name="length_positive"),
)


# --------------------------------------------------------------------------
# The controlled vocabulary, as data (D-017). Seeded at init; its hash is
# recorded on every assembly.
# --------------------------------------------------------------------------

vocabulary = Table(
    "vocabulary",
    metadata,
    Column("term", String(128), primary_key=True),
    Column("canonical_term", String(128), nullable=False),
    Column("tier", String(1), nullable=False),
    CheckConstraint("tier IN ('A', 'B')", name="tier_known"),
)

vocabulary_meta = Table(
    "vocabulary_meta",
    metadata,
    Column("vocabulary_meta_id", Integer, primary_key=True),
    Column("version", Integer, nullable=False),
    Column("content_hash", String(64), nullable=False),
    Column("seeded_at", DateTime(timezone=True), nullable=False, server_default=func.now()),
    CheckConstraint("vocabulary_meta_id = 1", name="singleton"),
)


# --------------------------------------------------------------------------
# Features.
# --------------------------------------------------------------------------

feature = Table(
    "feature",
    metadata,
    Column("feature_id", Integer, primary_key=True),
    Column("assembly_id", Integer, ForeignKey("assembly.assembly_id"), nullable=False),
    Column(
        "sequence_region_id",
        Integer,
        ForeignKey("sequence_region.sequence_region_id"),
        nullable=False,
    ),
    Column("source_file_id", Integer, ForeignKey("source_file.source_file_id"), nullable=False),
    # Column 9 ID after id_prefix, verbatim otherwise. Repeated across the
    # segments of one discontinuous feature.
    #
    # NULL when the input row carried no ID at all -- AUGUSTUS emits intron,
    # start_codon and stop_codon rows with only a Parent. Such a row is its own
    # feature: with no ID there is nothing to group segments by, so it cannot be
    # part of a discontinuous feature. `verify` therefore counts
    # COUNT(DISTINCT source_id) + COUNT(*) over the NULL rows; see D-024, which
    # explains why synthesising an ID here would break D-008's independence.
    Column("source_id", String(512), nullable=True),
    # Column 3 exactly as it appeared. Never rewritten. `verify` counts this.
    Column("source_type", String(128), nullable=False),
    # After synonym mapping. Linkage and `query` use this. Equal to
    # source_type for any term absent from the vocabulary.
    Column("feature_type", String(128), nullable=False),
    Column("start", Integer, nullable=False),  # 1-based inclusive
    Column("end", Integer, nullable=False),  # 1-based inclusive
    Column("strand", String(1), nullable=False),
    Column("phase", Integer, nullable=True),
    Column("score", Float, nullable=True),
    # Ascending coordinate order, strand stored separately. Not file order:
    # some pipelines emit minus-strand segments 3'->5'.
    Column("segment_index", Integer, nullable=False, server_default="0"),
    Column("parent_id", Integer, ForeignKey("feature.feature_id"), nullable=True),
    Column("linkage_method", String(32), nullable=True),
    Column("source_program", String(128), nullable=True),  # column 2, verbatim
    UniqueConstraint(
        "sequence_region_id",
        "source_id",
        "start",
        "end",
        name="uq_feature_region_source_id_span",
    ),
    CheckConstraint("start >= 1", name="start_positive"),
    CheckConstraint("end >= start", name="span_ordered"),
    CheckConstraint("strand IN ('+', '-', '.', '?')", name="strand_known"),
    CheckConstraint("phase IS NULL OR phase IN (0, 1, 2)", name="phase_known"),
    CheckConstraint("segment_index >= 0", name="segment_index_non_negative"),
    CheckConstraint(f"linkage_method IS NULL OR linkage_method IN ({_LINKAGE_IN})",
                    name="linkage_method_known"),
)

# `verify` groups by source_type; `query` filters on feature_type. Two indexes
# because they are two different columns, and adding the second after the
# large-genome benchmark means re-running it.
Index("ix_feature_assembly_source_type", feature.c.assembly_id, feature.c.source_type)
Index("ix_feature_assembly_feature_type", feature.c.assembly_id, feature.c.feature_type)
Index("ix_feature_assembly_source_id", feature.c.assembly_id, feature.c.source_id)
Index("ix_feature_parent", feature.c.parent_id)
# Region queries: "what is on scaffold 12 between X and Y". A plain B-tree on
# (region, start, end) is deliberate -- SQLite's RTREE is a compile-time
# option with no MySQL/PostgreSQL equivalent, and using it would void the
# single-code-path claim in D-002.
Index("ix_feature_region_span", feature.c.sequence_region_id, feature.c.start, feature.c.end)

# Arbitrary column-9 attributes and GenBank qualifiers, kept so `explain` can
# show what the input actually said. Not a substitute for a modelled column.
feature_attribute = Table(
    "feature_attribute",
    metadata,
    Column("feature_attribute_id", Integer, primary_key=True),
    Column("feature_id", Integer, ForeignKey("feature.feature_id"), nullable=False),
    Column("key", String(128), nullable=False),
    Column("value", Text, nullable=False),
)
Index("ix_feature_attribute_feature_key", feature_attribute.c.feature_id,
      feature_attribute.c.key)
Index("ix_feature_attribute_key_value", feature_attribute.c.key, feature_attribute.c.value)


# --------------------------------------------------------------------------
# Functional annotation (milestone 7).
#
# Deliberately outside the feature-graph engine (D-013): eggNOG and
# InterProScan produce per-protein tables, not feature graphs, and forcing them
# into the profile abstraction would bend it. They attach to features rather
# than becoming them, so `verify`'s feature counts are untouched by loading or
# not loading them.
# --------------------------------------------------------------------------

annotation_run = Table(
    "annotation_run",
    metadata,
    Column("annotation_run_id", Integer, primary_key=True),
    Column("assembly_id", Integer, ForeignKey("assembly.assembly_id"), nullable=False),
    Column("tool", String(64), nullable=False),  # eggnog | interproscan
    Column("tool_version", String(64), nullable=True),
    Column("path", Text, nullable=False),
    Column("sha256", String(64), nullable=False),
    # How many query proteins the file held, how many matched a feature, and how
    # many did not. Stored rather than derived: an unmatched query leaves no row
    # behind, so without this the fact that 40% of a file went nowhere would be
    # unrecoverable after the load (D-039).
    Column("query_count", Integer, nullable=False),
    Column("matched_count", Integer, nullable=False),
    Column("unmatched_count", Integer, nullable=False),
    Column("loaded_at", DateTime(timezone=True), nullable=False, server_default=func.now()),
    UniqueConstraint("assembly_id", "path", name="uq_annotation_run_assembly_path"),
)

#: One row per hit. Uniform across both tools: a GO term from eggNOG and a Pfam
#: domain from InterProScan are both "this protein has this accession", which is
#: why one table serves rather than a table per tool plus a term table.
protein_annotation = Table(
    "protein_annotation",
    metadata,
    Column("protein_annotation_id", Integer, primary_key=True),
    Column(
        "annotation_run_id",
        Integer,
        ForeignKey("annotation_run.annotation_run_id"),
        nullable=False,
    ),
    Column("feature_id", Integer, ForeignKey("feature.feature_id"), nullable=False),
    # The identifier the tool was given, kept verbatim. The join to `feature` is
    # by lookup and may have gone through protein_id or transcript_id, so this is
    # the only record of what the functional file actually said.
    Column("query_id", String(512), nullable=False),
    # Which analysis produced it: 'Pfam', 'GO', 'KEGG_ko', 'CAZy', 'PANTHER'...
    Column("analysis", String(64), nullable=False),
    Column("accession", String(128), nullable=False),
    Column("description", Text, nullable=True),
    # Protein coordinates, 1-based inclusive. NULL for whole-protein assignments
    # such as a GO term; set for a located domain.
    Column("protein_start", Integer, nullable=True),
    Column("protein_end", Integer, nullable=True),
    Column("score", Float, nullable=True),
    Column("evalue", Float, nullable=True),
    UniqueConstraint(
        "annotation_run_id",
        "feature_id",
        "analysis",
        "accession",
        "protein_start",
        name="uq_protein_annotation_hit",
    ),
    CheckConstraint(
        "protein_start IS NULL OR protein_start >= 1", name="protein_start_positive"
    ),
    CheckConstraint(
        "protein_end IS NULL OR protein_start IS NULL OR protein_end >= protein_start",
        name="protein_span_ordered",
    ),
)
Index(
    "ix_protein_annotation_analysis_accession",
    protein_annotation.c.analysis,
    protein_annotation.c.accession,
)
Index("ix_protein_annotation_feature", protein_annotation.c.feature_id)


#: Tables that hold no user data and are populated by ``init`` itself.
SEED_TABLES = frozenset({"vocabulary", "vocabulary_meta"})
