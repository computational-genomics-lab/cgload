"""Milestone 4 and 5 tests: the loader and the `load` command.

What these assert, stated before they were written:

* All three real dialects load end to end and produce the counts the parser
  tests already pin, so the loader adds storage without changing meaning.
* A failure anywhere — bad dialect, colliding IDs, wrong assembly, a rejected
  row — leaves the database *and* the data directory exactly as they were.
* Segments become rows; features stay features. The number `stats` reports is
  the number the parser counted.
* Parents resolve into real row ids, across files as well as within one.
* Declared proteins reach disk, and a bad one stops the load.
* Every input error exits 2 with a message naming the file, never a traceback.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import pytest
from click.testing import CliRunner
from cgload.cli.main import cli
from cgload.config import load_config
from cgload.db.counts import feature_count_expression
from cgload.db.engine import build_engine
from cgload.db.init import create_schema
from cgload.db.schema import assembly, feature, feature_attribute, source_file
from cgload.loader import LoadError, load_organism
from sqlalchemy import func, select

FIXTURES = Path(__file__).parent / "fixtures"
REFSEQ_GFF = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
FUNANNOTATE_GFF = FIXTURES / "funannotate_aphaeospermum_pheo_arth1.gff3"
AUGUSTUS_GFF = FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3"


def write_fasta(path: Path, name: str, length: int = 60000) -> Path:
    random.seed(11)
    sequence = "".join(random.choice("ACGT") for _ in range(length))
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for index in range(0, length, 60):
            handle.write(sequence[index : index + 60] + "\n")
    return path


def config_for(tmp_path: Path, seqid: str, annotations: list[tuple[Path, str | None]],
               *, name: str = "asm1", strain: str = "S1") -> Path:
    fasta = write_fasta(tmp_path / f"{seqid}.fa", seqid)
    lines = [
        "[organism]",
        'genus = "Testus"',
        'species = "exampleii"',
        f'strain = "{strain}"',
        "ncbi_taxon_id = 12345",
        "",
        "[assembly]",
        f'name = "{name}"',
        f'fasta = "{fasta.name}"',
        "",
    ]
    for source, prefix in annotations:
        target = tmp_path / source.name
        if not target.exists():
            target.write_text(source.read_text())
        lines += ["[[files.annotation]]", f'path = "{target.name}"']
        if prefix:
            lines.append(f'id_prefix = "{prefix}"')
        lines.append("")
    path = tmp_path / "organism.toml"
    path.write_text("\n".join(lines))
    return path


@pytest.fixture
def engine(tmp_path):
    eng = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(eng)
    return eng


# -- the three dialects, end to end ---------------------------------------


@pytest.mark.parametrize(
    ("gff", "seqid", "expected"),
    [
        (
            REFSEQ_GFF,
            "NC_026477.1",
            {"region": 1, "gene": 95, "pseudogene": 3, "mRNA": 11, "tRNA": 61,
             "rRNA": 26, "exon": 147, "CDS": 10},
        ),
        (
            FUNANNOTATE_GFF,
            "pheo_arth1",
            {"gene": 46, "mRNA": 21, "tRNA": 25, "exon": 165, "CDS": 21},
        ),
        (
            AUGUSTUS_GFF,
            "pheo_arth1",
            {"gene": 49, "transcript": 49, "CDS": 49, "intron": 77,
             "start_codon": 49, "stop_codon": 49},
        ),
    ],
)
def test_each_dialect_loads_with_the_counts_the_parser_reported(
    engine, tmp_path, gff, seqid, expected
):
    """Milestone 5's exit criterion, per profile. The loader must not change the
    feature count the parser arrived at."""
    config = load_config(config_for(tmp_path, seqid, [(gff, None)]))
    report = load_organism(engine, config, tmp_path / "data")

    assert report.files[0].counts_by_type == expected
    assert report.feature_count == sum(expected.values())

    with engine.connect() as connection:
        stored = dict(
            connection.execute(
                select(feature.c.source_type, feature_count_expression())
                .group_by(feature.c.source_type)
            ).all()
        )
    assert stored == expected


def test_segments_become_rows_while_features_stay_features(engine, tmp_path):
    """56 CDS rows in the RefSeq file, 10 CDS features. Both numbers must be
    recoverable, and it is the feature count that `verify` will compare."""
    config = load_config(config_for(tmp_path, "NC_026477.1", [(REFSEQ_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        rows = connection.execute(
            select(func.count()).select_from(feature).where(feature.c.source_type == "CDS")
        ).scalar_one()
        features = connection.execute(
            select(feature_count_expression()).where(feature.c.source_type == "CDS")
        ).scalar_one()
    assert (rows, features) == (56, 10)


def test_anonymous_rows_each_get_a_row_and_a_null_id(engine, tmp_path):
    """D-024, through the database this time."""
    config = load_config(config_for(tmp_path, "pheo_arth1", [(AUGUSTUS_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        anonymous = connection.execute(
            select(func.count()).select_from(feature).where(feature.c.source_id.is_(None))
        ).scalar_one()
    assert anonymous == 175


def test_segment_index_is_coordinate_order_in_the_database(engine, tmp_path):
    config = load_config(config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.source_id, feature.c.segment_index, feature.c.start)
            .where(feature.c.source_type == "CDS")
            .order_by(feature.c.source_id, feature.c.segment_index)
        ).all()

    grouped: dict[str, list[int]] = {}
    for row in rows:
        grouped.setdefault(row.source_id, []).append(row.start)
    assert any(len(starts) > 1 for starts in grouped.values())
    for starts in grouped.values():
        assert starts == sorted(starts)


def test_parents_resolve_to_real_rows(engine, tmp_path):
    config = load_config(config_for(tmp_path, "NC_026477.1", [(REFSEQ_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")

    child = feature.alias("child")
    parent = feature.alias("parent")
    with engine.connect() as connection:
        pairs = connection.execute(
            select(child.c.source_type, parent.c.source_type)
            .select_from(child.join(parent, child.c.parent_id == parent.c.feature_id))
            .distinct()
        ).all()
        orphaned = connection.execute(
            select(func.count())
            .select_from(feature)
            .where(feature.c.linkage_method.is_not(None), feature.c.parent_id.is_(None))
        ).scalar_one()

    assert ("exon", "mRNA") in pairs
    assert ("CDS", "mRNA") in pairs
    assert ("mRNA", "pseudogene") in pairs
    assert orphaned == 0


def test_attributes_are_stored_once_per_feature_not_per_segment(engine, tmp_path):
    """RefSeq repeats every attribute on every segment row; storing them per
    segment would multiply the table `explain` reads for no gain."""
    config = load_config(config_for(tmp_path, "NC_026477.1", [(REFSEQ_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.source_id, feature_attribute.c.key)
            .select_from(
                feature_attribute.join(
                    feature, feature.c.feature_id == feature_attribute.c.feature_id
                )
            )
            .where(feature.c.source_type == "CDS", feature_attribute.c.key == "protein_id")
        ).all()
    assert len(rows) == len({row.source_id for row in rows})


def test_decoded_attribute_value_survives_to_the_database(engine, tmp_path):
    config = load_config(config_for(tmp_path, "NC_026477.1", [(REFSEQ_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")

    with engine.connect() as connection:
        value = connection.execute(
            select(feature_attribute.c.value).where(feature_attribute.c.key == "strain")
        ).scalar_one()
    assert value == "PH-1; NRRL 31084"


# -- two files, one assembly (D-014) --------------------------------------


def test_two_dialects_load_into_one_assembly(engine, tmp_path):
    config = load_config(
        config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None), (AUGUSTUS_GFF, "AUG_")])
    )
    report = load_organism(engine, config, tmp_path / "data")

    assert {entry.dialect for entry in report.files} == {"funannotate", "augustus"}
    with engine.connect() as connection:
        by_dialect = dict(
            connection.execute(
                select(source_file.c.detected_dialect, source_file.c.id_prefix).where(
                    source_file.c.file_format == "gff3"
                )
            ).all()
        )
        programs = set(
            connection.execute(select(feature.c.source_program).distinct()).scalars()
        )
    assert by_dialect == {"funannotate": None, "augustus": "AUG_"}
    assert programs == {"funannotate", "AUGUSTUS"}


def test_id_prefix_keeps_two_files_distinct_and_reaches_the_proteins(engine, tmp_path):
    """The prefix must apply to the comment-derived protein keys too, or the
    oracle looks up `g1` while the features are stored as `AUG_g1` and a correct
    file fails to load."""
    config = load_config(
        config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None), (AUGUSTUS_GFF, "AUG_")])
    )
    report = load_organism(engine, config, tmp_path / "data")

    proteins = (report.data_directory / "annotation" / "proteins.faa").read_text()
    assert ">AUG_g1\n" in proteins
    assert ">g1\n" not in proteins


def test_colliding_ids_across_files_are_refused_naming_both(engine, tmp_path):
    copy = tmp_path / "second_copy.gff3"
    copy.write_text(AUGUSTUS_GFF.read_text())
    config = load_config(
        config_for(tmp_path, "pheo_arth1", [(AUGUSTUS_GFF, None), (copy, None)])
    )
    with pytest.raises(LoadError, match="defined in both"):
        load_organism(engine, config, tmp_path / "data")


def test_the_same_file_listed_twice_is_refused(engine, tmp_path):
    """Would double every count. Caught before the database, not by a
    constraint."""
    config_path = config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None)])
    text = config_path.read_text()
    config_path.write_text(
        f'{text}\n[[files.annotation]]\npath = "{FUNANNOTATE_GFF.name}"\n'
    )
    with pytest.raises(LoadError, match="listed twice"):
        load_organism(engine, load_config(config_path), tmp_path / "data")


# -- refusals --------------------------------------------------------------


def test_annotation_for_a_different_assembly_is_refused(engine, tmp_path):
    """The seqid in the annotation is absent from the FASTA: the two files
    describe different assemblies, and loading would store unreachable rows."""
    config = load_config(config_for(tmp_path, "some_other_scaffold", [(FUNANNOTATE_GFF, None)]))
    with pytest.raises(LoadError, match="different assemblies"):
        load_organism(engine, config, tmp_path / "data")


def test_an_unrecognised_dialect_is_refused(engine, tmp_path):
    unknown = tmp_path / "mystery.gff3"
    unknown.write_text(
        "##gff-version 3\nchr1\tMysteryTool\tgene\t1\t99\t.\t+\t.\tID=g1\n"
    )
    config = load_config(config_for(tmp_path, "chr1", [(unknown, None)]))
    with pytest.raises(LoadError, match="no confident match"):
        load_organism(engine, config, tmp_path / "data")


def test_a_bad_protein_declaration_stops_the_load(engine, tmp_path, monkeypatch):
    """The oracle wired into the load path, not just the parser tests."""
    import cgload.parsers.normalise as normalise_module

    real = normalise_module.check_declared_proteins
    monkeypatch.setattr(
        "cgload.loader.check_declared_proteins",
        lambda result: ["g1: CDS totals 851 bp, not a whole number of codons"]
        if real(result) == [] and result.proteins
        else real(result),
    )
    config = load_config(config_for(tmp_path, "pheo_arth1", [(AUGUSTUS_GFF, None)]))
    with pytest.raises(LoadError, match="do not follow from"):
        load_organism(engine, config, tmp_path / "data")


def test_reloading_the_same_assembly_is_refused(engine, tmp_path):
    config = load_config(config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None)]))
    load_organism(engine, config, tmp_path / "data")
    with pytest.raises(Exception, match="already"):
        load_organism(engine, config, tmp_path / "data")


# -- all or nothing, with features this time ------------------------------


def test_a_failure_after_features_are_inserted_rolls_back_everything(engine, tmp_path):
    """The milestone-4 exit criterion. Fail late — after the assembly, the
    regions and thousands of feature rows are in — and both the database and the
    data directory must come out untouched."""
    data_root = tmp_path / "data"
    config = load_config(
        config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None), (AUGUSTUS_GFF, "AUG_")])
    )

    import cgload.loader as loader_module

    original = loader_module._insert_attributes
    calls = {"n": 0}

    def fail_on_second(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise RuntimeError("simulated failure on the second file")
        return original(*args, **kwargs)

    monkey = pytest.MonkeyPatch()
    monkey.setattr(loader_module, "_insert_attributes", fail_on_second)
    try:
        with pytest.raises(RuntimeError, match="simulated"):
            load_organism(engine, config, data_root)
    finally:
        monkey.undo()

    with engine.connect() as connection:
        for table in (assembly, feature, feature_attribute, source_file):
            assert connection.execute(select(func.count()).select_from(table)).scalar_one() == 0

    assert [p for p in data_root.rglob("*") if p.is_file()] == []


# -- the command ----------------------------------------------------------


def test_load_command_reports_what_it_did(engine, tmp_path):
    config_path = config_for(tmp_path, "pheo_arth1", [(AUGUSTUS_GFF, None)])
    result = CliRunner().invoke(
        cli,
        ["load", "--config", str(config_path), "--data-dir", str(tmp_path / "data"),
         "--db-url", str(engine.url)],
    )
    assert result.exit_code == 0, result.output
    assert "augustus" in result.output
    assert "322 features" in result.output
    assert "49 declared proteins" in result.output


def test_load_command_exits_two_without_a_traceback(engine, tmp_path):
    unknown = tmp_path / "mystery.gff3"
    unknown.write_text("##gff-version 3\nchr1\tNope\tgene\t1\t9\t.\t+\t.\tID=g\n")
    config_path = config_for(tmp_path, "chr1", [(unknown, None)])
    result = CliRunner().invoke(
        cli,
        ["load", "--config", str(config_path), "--data-dir", str(tmp_path / "data"),
         "--db-url", str(engine.url)],
    )
    assert result.exit_code == 2
    assert "Traceback" not in result.output
    assert "no confident match" in result.output


def test_load_command_refuses_an_uninitialised_database(tmp_path):
    config_path = config_for(tmp_path, "pheo_arth1", [(FUNANNOTATE_GFF, None)])
    result = CliRunner().invoke(
        cli,
        ["load", "--config", str(config_path), "--data-dir", str(tmp_path / "data"),
         "--db-url", f"sqlite:///{tmp_path / 'absent.db'}"],
    )
    assert result.exit_code == 2
    assert "cgload init" in result.output


def test_load_and_stats_agree(engine, tmp_path):
    config_path = config_for(tmp_path, "NC_026477.1", [(REFSEQ_GFF, None)])
    runner = CliRunner()
    loaded = runner.invoke(
        cli,
        ["load", "--config", str(config_path), "--data-dir", str(tmp_path / "data"),
         "--db-url", str(engine.url)],
    )
    assert loaded.exit_code == 0, loaded.output
    stats = runner.invoke(cli, ["stats", "--db-url", str(engine.url)])
    assert stats.exit_code == 0
    assert "exon 147" in stats.output
    assert "CDS 10" in stats.output


def test_manifest_records_the_dialect_and_the_evidence(engine, tmp_path):
    config = load_config(config_for(tmp_path, "pheo_arth1", [(AUGUSTUS_GFF, None)]))
    report = load_organism(engine, config, tmp_path / "data")
    manifest = json.loads((report.data_directory / "manifest.json").read_text())

    annotation = [entry for entry in manifest["files"] if entry["role"] == "annotation"]
    assert len(annotation) == 1
    assert annotation[0]["dialect"] == "augustus"
    assert annotation[0]["evidence"]
    assert annotation[0]["features"] == 322
    assert manifest["proteins"] == 49


def test_load_is_no_longer_pending():
    from cgload.cli.main import PENDING

    assert "load" not in PENDING
    assert set(PENDING) == {"query", "remove"}


# -- the bundled example (milestone 4 exit criterion) ---------------------


def test_the_bundled_example_loads(engine, tmp_path):
    """`examples/minimal/` is promised by the README and the paper's
    availability section, so it must actually load — and be checked here, not by
    a human remembering to try it."""
    source = Path(__file__).resolve().parents[1] / "examples" / "minimal"
    assert source.is_dir(), "examples/minimal is missing"

    working = tmp_path / "example"
    working.mkdir()
    for item in source.iterdir():
        (working / item.name).write_text(item.read_text())

    config = load_config(working / "organism.toml")
    report = load_organism(engine, config, tmp_path / "data")
    assert report.files[0].counts_by_type == {
        "gene": 3, "mRNA": 2, "tRNA": 1, "exon": 4, "CDS": 2
    }
