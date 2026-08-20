"""Milestone 9b: exporting from the database, and provenance surviving it.

Milestone 9 delivered the serialiser and tested it against parser output. What it
did not deliver was the path a user actually takes: reading the *stored* features
and writing them out, after the input file is gone. That is what `cgload export`
is for, and these tests cover it.

What they assert, stated before they were written:

* A loaded assembly exports and re-imports to the same feature counts.
* The export declares itself and is detected as cgload output, not as the
  dialect it was derived from.
* A computed parent edge comes back as computed, not as one the file stated --
  the D-032 failure arriving by the export path.
* A mixed-dialect assembly exports as one file that still says, per feature,
  which program produced it.
* Nothing is invented: an AUGUSTUS export has no exons and says so.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest
from click.testing import CliRunner
from cgload.cli.main import cli
from cgload.config import load_config
from cgload.db.engine import build_engine
from cgload.db.export_source import load_context, load_features
from cgload.db.init import create_schema
from cgload.db.schema import feature
from cgload.db.verification import count_database
from cgload.export import EXPORT_PRAGMA, is_cgload_export
from cgload.loader import load_organism
from cgload.parsers import profiles as profile_table
from cgload.parsers.normalise import counts_by_source_type, normalise
from cgload.parsers.tokenizer import tokenize
from cgload.verify import count_source_naively
from sqlalchemy import select

FIXTURES = Path(__file__).parent / "fixtures"
PAIR_GFF = FIXTURES / "pair_fgraminearum_NC_026477.gff3"
PAIR_GBFF = FIXTURES / "pair_fgraminearum_NC_026477.gbff"
AUGUSTUS = FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3"
MYCMAY = FIXTURES / "refseq_mycmay_GCF_000328475.gff3"
EGGNOG = FIXTURES / "eggnog_mycmay_emapper214.annotations"


def _write_fasta(annotations: list[Path], target: Path) -> None:
    import re

    lengths: dict[str, int] = {}
    for path in annotations:
        text = path.read_text()
        if path.suffix == ".gbff":
            name = re.search(r"^VERSION\s+(\S+)", text, re.M).group(1)
            lengths[name] = max(
                int(value)
                for line in text.splitlines()
                if line.startswith("     ") and ".." in line
                for value in re.findall(r"\b(\d{3,})\b", line)
            )
        else:
            for line in text.splitlines():
                if line.startswith("#"):
                    continue
                columns = line.split("\t")
                if len(columns) >= 9:
                    lengths[columns[0]] = max(
                        lengths.get(columns[0], 0), int(columns[4])
                    )
    with target.open("w") as handle:
        for name, longest in sorted(lengths.items()):
            handle.write(f">{name}\n")
            handle.write((("ACGT" * 15) + "\n") * ((longest + 1000) // 60 + 1))


def _load(engine, tmp_path: Path, annotations: list[tuple[Path, str | None]],
          functional: list[Path] | None = None):
    fasta = tmp_path / "assembly.fa"
    _write_fasta([path for path, _ in annotations], fasta)
    lines = [
        "[organism]", 'genus = "Testus"', 'species = "exampleii"', 'strain = "S"',
        "", "[assembly]", 'name = "asm"', f'fasta = "{fasta.name}"', "",
    ]
    for path, prefix in annotations:
        (tmp_path / path.name).write_text(path.read_text())
        lines += ["[[files.annotation]]", f'path = "{path.name}"']
        if prefix:
            lines.append(f'id_prefix = "{prefix}"')
        lines.append("")
    for path in functional or []:
        (tmp_path / path.name).write_text(path.read_text())
        lines += ["[[files.functional]]", f'path = "{path.name}"', ""]
    config = tmp_path / "organism.toml"
    config.write_text("\n".join(lines))
    return load_organism(engine, load_config(config), tmp_path / "data")


@pytest.fixture
def engine(tmp_path):
    built = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(built)
    return built


def _export(engine, assembly_id: int, target: Path) -> Path:
    result = CliRunner().invoke(
        cli,
        ["export", "--assembly-id", str(assembly_id), "--db-url", str(engine.url),
         "--output", str(target)],
    )
    assert result.exit_code == 0, result.output
    return target


# -- the round trip through the database ----------------------------------


@pytest.mark.parametrize("fixture", [PAIR_GFF, PAIR_GBFF, AUGUSTUS, MYCMAY])
def test_a_loaded_assembly_exports_and_reimports_to_the_same_counts(
    engine, tmp_path, fixture
):
    """The path a user takes. Milestone 9's round trip went parser -> exporter;
    this one goes parser -> database -> exporter -> parser, which is the only
    version that tests `cgload export`."""
    report = _load(engine, tmp_path, [(fixture, None)])
    with engine.connect() as connection:
        stored = count_database(connection, report.assembly_id)

    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")

    assert dict(count_source_naively(written)) == dict(stored), fixture.name

    tokenized = tokenize(written)
    detection = profile_table.detect(tokenized)
    reimported = normalise(
        tokenized.features, detection.profile, comments=tokenized.comments
    )
    assert dict(counts_by_source_type(reimported)) == dict(stored)


def test_the_export_is_detected_as_cgload_not_as_its_source_dialect(engine, tmp_path):
    """Column 2 still says AUGUSTUS, correctly -- that is where the genes came
    from. But the file no longer follows AUGUSTUS's conventions, so the header
    pragma has to outrank column 2 (D-052)."""
    report = _load(engine, tmp_path, [(AUGUSTUS, None)])
    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")

    assert is_cgload_export(written.read_text())
    tokenized = tokenize(written)
    assert profile_table.detect(tokenized).profile.name == "cgload"
    assert "AUGUSTUS" in {item.source for item in tokenized.features}


# -- provenance across the round trip -------------------------------------


def test_a_computed_parent_edge_comes_back_computed(engine, tmp_path):
    """The bug this milestone nearly shipped. A GenBank load records 43 CDS
    features linked by containment; the export writes Parent attributes for them,
    so a self-stating profile would relabel every one `explicit_parent` on
    re-import -- D-032's false-provenance failure arriving by the export path.
    """
    report = _load(engine, tmp_path, [(PAIR_GBFF, None)])
    with engine.connect() as connection:
        before = Counter(
            connection.execute(
                select(feature.c.linkage_method).where(
                    feature.c.assembly_id == report.assembly_id,
                    feature.c.segment_index == 0,
                )
            ).scalars()
        )
    assert before["structural"] > 0, "fixture no longer exercises containment"

    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")
    tokenized = tokenize(written)
    reimported = normalise(
        tokenized.features,
        profile_table.detect(tokenized).profile,
        comments=tokenized.comments,
    )
    after = Counter(item.linkage_method for item in reimported.features)

    assert after["structural"] == before["structural"]
    assert after["locus_tag"] == before["locus_tag"]
    assert "explicit_parent" not in after


def test_the_exported_linkage_attribute_matches_what_the_parser_reads():
    """Two literals, two modules, deliberately not imported across the boundary.
    They have to agree, so a test says so rather than a comment."""
    from cgload.export import LINKAGE_ATTRIBUTE as written
    from cgload.parsers.profiles import LINKAGE_ATTRIBUTE as read

    assert written == read


def test_the_export_pragma_matches_what_detection_looks_for():
    from cgload.parsers.profiles import EXPORT_PRAGMA as expected

    assert expected == EXPORT_PRAGMA


# -- mixed sources ---------------------------------------------------------


def test_a_mixed_dialect_assembly_exports_as_one_file_keeping_each_source(
    engine, tmp_path
):
    """Provenance is per feature (D-014), so an assembly built from a RefSeq
    download and an AUGUSTUS run exports as one file where each gene still says
    where it came from."""
    report = _load(engine, tmp_path, [(PAIR_GFF, None), (AUGUSTUS, "AUG_")])
    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")

    text = written.read_text()
    assert "# source programs present: AUGUSTUS, RefSeq" in text

    tokenized = tokenize(written)
    assert {item.source for item in tokenized.features} == {"AUGUSTUS", "RefSeq"}

    with engine.connect() as connection:
        stored = count_database(connection, report.assembly_id)
    assert dict(count_source_naively(written)) == dict(stored)


def test_nothing_is_invented_and_the_header_says_what_is_absent(engine, tmp_path):
    """AUGUSTUS emits no exons, so neither does its export. The header records
    that, so 'the input had none' stays distinguishable from 'they were lost'."""
    report = _load(engine, tmp_path, [(AUGUSTUS, None)])
    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")
    text = written.read_text()

    assert "exon" not in dict(count_source_naively(written))
    assert "does not emit these types" in text
    assert "exon" in text.split("does not emit these types")[1][:200]


def test_anonymous_features_stay_anonymous(engine, tmp_path):
    """175 of AUGUSTUS's rows have no ID. Inventing one on export would breach
    D-023, so the round-trip assertion is on counts by type, never on identifier
    sets -- and this test pins that the count survives."""
    report = _load(engine, tmp_path, [(AUGUSTUS, None)])
    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")

    tokenized = tokenize(written)
    anonymous = [item for item in tokenized.features if item.id is None]
    assert len(anonymous) == 175
    assert {item.type for item in anonymous} == {"intron", "start_codon", "stop_codon"}


# -- the reading half in isolation ----------------------------------------


def test_stored_segments_regroup_into_features(engine, tmp_path):
    """A six-segment CDS is one feature on the way out, as it was on the way in."""
    report = _load(engine, tmp_path, [(PAIR_GFF, None)])
    with engine.connect() as connection:
        features = load_features(connection, report.assembly_id)
        stored = count_database(connection, report.assembly_id)

    assert len(features) == sum(stored.values())
    assert any(len(item.segments) > 1 for item in features)
    for item in features:
        starts = [segment.start for segment in item.segments]
        assert starts == sorted(starts), "segments must come out in coordinate order"


def test_the_context_carries_the_vocabulary_hash(engine, tmp_path):
    """An export names the vocabulary it was written under, so a re-import can
    tell whether the type mapping has changed since (D-017)."""
    report = _load(engine, tmp_path, [(PAIR_GFF, None)])
    with engine.connect() as connection:
        context = load_context(connection, report.assembly_id)
    assert len(context.vocabulary_hash) == 64
    assert context.sequence_regions

    written = _export(engine, report.assembly_id, tmp_path / "out.gff3")
    assert context.vocabulary_hash in written.read_text()


def test_exporting_a_missing_assembly_exits_two(engine, tmp_path):
    result = CliRunner().invoke(
        cli, ["export", "--assembly-id", "999", "--db-url", str(engine.url)]
    )
    assert result.exit_code == 2
    assert "no assembly with id" in result.output


# -- the TSV format --------------------------------------------------------


def test_functional_annotation_exports_as_tsv_with_the_go_caveat(engine, tmp_path):
    """D-045. Anyone summing the GO rows needs to know they include inherited
    ancestors before they do it."""
    report = _load(engine, tmp_path, [(MYCMAY, None)], functional=[EGGNOG])
    target = tmp_path / "out.tsv"
    result = CliRunner().invoke(
        cli,
        ["export", "--assembly-id", str(report.assembly_id), "--format", "tsv",
         "--db-url", str(engine.url), "--output", str(target)],
    )
    assert result.exit_code == 0, result.output

    text = target.read_text()
    assert "inherited ancestors" in text
    rows = [line for line in text.splitlines() if not line.startswith("#")]
    assert rows[0].split("\t")[0] == "feature_id"
    assert len(rows) - 1 == report.functional[0].hit_count


def test_export_is_no_longer_pending():
    from cgload.cli.main import PENDING

    assert "export" not in PENDING
    assert set(PENDING) == {"query", "remove"}
