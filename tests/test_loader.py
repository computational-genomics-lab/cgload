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
from cgload.cli.main import cli
from cgload.config import load_config
from cgload.db.counts import feature_count_expression
from cgload.db.engine import build_engine
from cgload.db.init import create_schema
from cgload.db.schema import assembly, feature, feature_attribute, source_file
from cgload.loader import LoadError, load_organism
from click.testing import CliRunner
from sqlalchemy import func, select

FIXTURES = Path(__file__).parent / "fixtures"
REFSEQ_GFF = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
FUNANNOTATE_GFF = FIXTURES / "funannotate_aphaeospermum_pheo_arth1.gff3"
AUGUSTUS_GFF = FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3"


def write_fasta(path: Path, name: str, length: int = 60000) -> Path:
    # A seeded random block, repeated: the RefSeq fixture's region row declares
    # the chromosome's full 8 Mb, and generating that one base at a time
    # dominated the run.
    random.seed(11)
    block = "".join(random.choice("ACGT") for _ in range(6000))
    sequence = (block * (length // len(block) + 1))[:length]
    with path.open("w") as handle:
        handle.write(f">{name}\n")
        for index in range(0, length, 60):
            handle.write(sequence[index : index + 60] + "\n")
    return path


def length_needed(annotations: list[Path], seqid: str) -> int:
    """The longest coordinate any annotation gives for ``seqid``.

    The loader refuses a row past the end of a linear sequence (AUDIT-6), so a
    test FASTA must be at least as long as the annotation says the sequence is.
    """
    longest = 0
    for path in annotations:
        for line in path.read_text().splitlines():
            if line.startswith("##sequence-region") and line.split()[1] == seqid:
                longest = max(longest, int(line.split()[3]))
            elif line and not line.startswith("#"):
                columns = line.split("\t")
                if len(columns) >= 9 and columns[0] == seqid:
                    longest = max(longest, int(columns[3]), int(columns[4]))
    return longest


def config_for(tmp_path: Path, seqid: str, annotations: list[tuple[Path, str | None]],
               *, name: str = "asm1", strain: str = "S1") -> Path:
    needed = length_needed([source for source, _prefix in annotations], seqid)
    fasta = write_fasta(tmp_path / f"{seqid}.fa", seqid, max(60000, needed + 1000))
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


# -- --strict-vocabulary (D-087) --------------------------------------------
# Two Sequence Ontology terms no file of the 30-genome RefSeq panel uses, so the
# RefSeq profile does not declare them (D-107).

UNDECLARED = (
    "##gff-version 3\n"
    "NC_000001.1\tRefSeq\tregion\t1\t20000\t.\t+\t.\tID=NC_000001.1:1..20000\n"
    "NC_000001.1\tGnomon\tgene\t1000\t3000\t.\t+\t.\tID=gene-A\n"
    "NC_000001.1\tGnomon\tmRNA\t1000\t3000\t.\t+\t.\tID=rna-A;Parent=gene-A\n"
    "NC_000001.1\tGnomon\texon\t1000\t3000\t.\t+\t.\tID=exon-A;Parent=rna-A\n"
    "NC_000001.1\tRefSeq\triboswitch\t5000\t5100\t.\t+\t.\tID=rs-1\n"
    "NC_000001.1\tRefSeq\tattenuator\t6000\t6100\t.\t+\t.\tID=att-1\n"
)


def _undeclared_project(tmp_path: Path) -> Path:
    gff = tmp_path / "undeclared.gff3"
    gff.write_text(UNDECLARED)
    return config_for(tmp_path, "NC_000001.1", [(gff, None)])


def test_undeclared_types_load_with_a_warning_by_default(engine, tmp_path):
    report = load_organism(engine, load_config(_undeclared_project(tmp_path)), tmp_path / "data")
    assert report.files[0].counts_by_type["riboswitch"] == 1
    assert any("not declared by the refseq profile: attenuator, riboswitch" in warning
               for warning in report.warnings)


def test_strict_vocabulary_refuses_listing_every_undeclared_type(engine, tmp_path):
    config = load_config(_undeclared_project(tmp_path))
    with pytest.raises(LoadError, match="2 feature type.*: attenuator, riboswitch"):
        load_organism(engine, config, tmp_path / "data", strict_vocabulary=True)
    with engine.connect() as connection:
        assert connection.execute(select(func.count()).select_from(feature)).scalar_one() == 0
    assert not (tmp_path / "data").exists() or not any((tmp_path / "data").iterdir())


def test_strict_vocabulary_on_the_command_line(tmp_path):
    config = _undeclared_project(tmp_path)
    url = f"sqlite:///{tmp_path / 'cli.db'}"
    assert CliRunner().invoke(cli, ["init", "--db-url", url]).exit_code == 0
    result = CliRunner().invoke(cli, ["load", "--db-url", url, "--config", str(config),
                                      "--data-dir", str(tmp_path / "d"), "--strict-vocabulary"])
    assert result.exit_code == 2
    assert "attenuator, riboswitch" in result.output


# -- ncRNA classes from the INSDC vocabulary (D-088) ------------------------

VERTEBRATE_NCRNA = (
    "##gff-version 3\n"
    "NC_000007.14\tRefSeq\tregion\t1\t60000\t.\t+\t.\tID=NC_000007.14:1..60000;genome=chromosome\n"
    "NC_000007.14\tBestRefSeq\tgene\t1000\t1200\t.\t-\t.\tID=gene-RNY1;gene_biotype=Y_RNA\n"
    "NC_000007.14\tBestRefSeq\tY_RNA\t1000\t1200\t.\t-\t.\tID=rna-NR_004391.1;Parent=gene-RNY1\n"
    "NC_000007.14\tBestRefSeq\texon\t1000\t1200\t.\t-\t.\tID=exon-NR_004391.1-1;Parent=rna-NR_004391.1\n"
    "NC_000007.14\tBestRefSeq\tgene\t3000\t3100\t.\t+\t.\tID=gene-VTRNA1-1;gene_biotype=vault_RNA\n"
    "NC_000007.14\tBestRefSeq\tvault_RNA\t3000\t3100\t.\t+\t.\tID=rna-NR_026703.1;Parent=gene-VTRNA1-1\n"
    "NC_000007.14\tBestRefSeq\texon\t3000\t3100\t.\t+\t.\tID=exon-NR_026703.1-1;Parent=rna-NR_026703.1\n"
    "NC_000007.14\tBestRefSeq\tgene\t5000\t5300\t.\t-\t.\tID=gene-RMRP;gene_biotype=RNase_MRP_RNA\n"
    "NC_000007.14\tBestRefSeq\tRNase_MRP_RNA\t5000\t5300\t.\t-\t.\tID=rna-NR_003051.3;Parent=gene-RMRP\n"
    "NC_000007.14\tBestRefSeq\texon\t5000\t5300\t.\t-\t.\tID=exon-NR_003051.3-1;Parent=rna-NR_003051.3\n"
    "NC_000007.14\tBestRefSeq\tgene\t8000\t9000\t.\t+\t.\tID=gene-A;gene_biotype=protein_coding\n"
    "NC_000007.14\tBestRefSeq\tmRNA\t8000\t9000\t.\t+\t.\tID=rna-A;Parent=gene-A\n"
    "NC_000007.14\tBestRefSeq\texon\t8000\t9000\t.\t+\t.\tID=exon-A-1;Parent=rna-A\n"
    "NC_000007.14\tBestRefSeq\tCDS\t8100\t8900\t.\t+\t0\tID=cds-A;Parent=rna-A\n"
)


def test_vertebrate_ncrna_classes_load_with_their_exons(engine, tmp_path):
    """Y_RNA, vault_RNA and RNase_MRP_RNA loaded inert, so their exons were
    refused ("a exon may only attach to ..."), and with them every vertebrate
    RefSeq file. They are ncRNA transcripts, as snRNA is."""
    gff = tmp_path / "vertebrate.gff3"
    gff.write_text(VERTEBRATE_NCRNA)
    config = load_config(config_for(tmp_path, "NC_000007.14", [(gff, None)]))
    report = load_organism(engine, config, tmp_path / "data")
    counts = report.files[0].counts_by_type
    assert counts["Y_RNA"] == counts["vault_RNA"] == counts["RNase_MRP_RNA"] == 1
    assert counts["exon"] == 4
    with engine.connect() as connection:
        canonical = connection.execute(
            select(feature.c.feature_type).where(feature.c.source_type == "Y_RNA")
        ).scalar_one()
        unlinked = connection.execute(
            select(func.count()).select_from(feature).where(
                feature.c.source_type == "exon", feature.c.parent_id.is_(None))
        ).scalar_one()
    assert canonical == "ncRNA"
    assert unlinked == 0


# -- programmed frameshifts in RefSeq GFF3 (D-089) --------------------------

FRAMESHIFT = (
    "##gff-version 3\n"
    "NC_000007.14\tRefSeq\tregion\t1\t60000\t.\t+\t.\tID=NC_000007.14:1..60000;genome=chromosome\n"
    "NC_000007.14\tBestRefSeq\tgene\t1000\t5000\t.\t+\t.\tID=gene-FS;gene_biotype=protein_coding\n"
    "NC_000007.14\tBestRefSeq\tmRNA\t1000\t5000\t.\t+\t.\tID=rna-FS;Parent=gene-FS\n"
    "NC_000007.14\tBestRefSeq\texon\t1000\t5000\t.\t+\t.\tID=exon-FS-1;Parent=rna-FS\n"
    "NC_000007.14\tBestRefSeq\tCDS\t1100\t2279\t.\t+\t0\tID=cds-FS;Parent=rna-FS{extra}\n"
    "NC_000007.14\tBestRefSeq\tCDS\t2279\t3448\t.\t+\t0\tID=cds-FS;Parent=rna-FS{extra}\n"
)


def _frameshift_config(tmp_path: Path, extra: str) -> Path:
    gff = tmp_path / "frameshift.gff3"
    gff.write_text(FRAMESHIFT.format(extra=extra))
    return config_for(tmp_path, "NC_000007.14", [(gff, None)])


def test_a_ribosomal_slippage_cds_may_share_a_base(engine, tmp_path):
    """The two rows share base 2279, read twice by the ribosome. NCBI writes
    GenBank's /ribosomal_slippage as exception=ribosomal slippage."""
    config = load_config(_frameshift_config(tmp_path, ";exception=ribosomal slippage"))
    report = load_organism(engine, config, tmp_path / "data")
    assert report.files[0].counts_by_type["CDS"] == 1
    with engine.connect() as connection:
        rows = connection.execute(
            select(feature.c.start, feature.c.end)
            .where(feature.c.source_id == "cds-FS")
            .order_by(feature.c.segment_index)
        ).all()
    assert [tuple(row) for row in rows] == [(1100, 2279), (2279, 3448)]


def test_an_unexplained_overlap_is_still_refused(engine, tmp_path):
    for extra in ("", ";exception=trans-splicing"):
        config = load_config(_frameshift_config(tmp_path, extra))
        with pytest.raises(LoadError, match="cannot overlap"):
            load_organism(engine, config, tmp_path / f"data{len(extra)}")
