"""Milestone 10: parent envelopes that contradict their children.

Every fixture here is reconstructed from *Zea mays* B73 GCF_902167145.1, with
the true coordinates taken from the GenBank file for the same assembly.
"""

from __future__ import annotations

import pytest
from cgload.envelope import (
    apply_policy,
    find_violations,
    format_report,
    repair_envelopes,
)
from cgload.parsers.normalise import NormalisedFeature, Segment


def feature(source_id, source_type, start, end, parent=None, strand="-", segments=None):
    spans = segments or [(start, end)]
    return NormalisedFeature(
        source_id=source_id,
        source_type=source_type,
        feature_type=source_type,
        seqid="NC_007982.1",
        strand=strand,
        source_program="RefSeq",
        segments=[
            Segment(line_number=0, start=s, end=e, phase=None, score=None)
            for s, e in spans
        ],
        parent_source_id=parent,
    )


@pytest.fixture
def zeam_mp186():
    """`nad1`, trans-spliced across four loci, one on the opposite strand.

    Published GFF3 envelopes: gene 266,974..267,232 (the last join segment
    only), mRNA 50,490..267,232. GenBank join for the same gene spans
    50,490..548,772. The exon rows are correct in both files.
    """
    return [
        feature("gene-ZeamMp186", "gene", 266974, 267232),
        feature("rna-ZeamMp186", "mRNA", 50490, 267232, "gene-ZeamMp186", strand="?"),
        feature("exon-1", "exon", 50490, 50874, "rna-ZeamMp186"),
        feature("exon-2", "exon", 320928, 321010, "rna-ZeamMp186"),
        feature("exon-3", "exon", 322404, 322595, "rna-ZeamMp186"),
        feature("exon-4", "exon", 548714, 548772, "rna-ZeamMp186"),
        feature("exon-5", "exon", 266974, 267232, "rna-ZeamMp186"),
    ]


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def test_violations_are_found(zeam_mp186):
    """Three exons escape the mRNA, and the mRNA escapes the gene."""
    violations = find_violations(zeam_mp186)
    assert len(violations) == 4
    assert {v.parent_id for v in violations} == {"gene-ZeamMp186", "rna-ZeamMp186"}


def test_every_row_passes_the_per_row_checks(zeam_mp186):
    """This is why the check is needed. Measured on the full maize annotation:
    36 structural problems, 34 of which pass every per-row rule and load
    cleanly with envelopes wrong by up to 281 kb."""
    for item in zeam_mp186:
        assert item.start <= item.end
        assert item.strand in ("+", "-", ".", "?")
        assert item.end <= 569630  # NC_007982.1 length


def test_strand_is_not_the_marker():
    """An early hypothesis held that `?` identified the affected rows. On the
    real annotation the affected parents are six `-`, two `?` and one `+`, so a
    fix scoped to `?` would have missed two thirds of the cases."""
    features = [
        feature("gene-plus", "gene", 111614, 115826, strand="+"),
        feature("rna-plus", "mRNA", 111614, 511384, "gene-plus", strand="+"),
    ]
    assert len(find_violations(features)) == 1


def test_a_correct_hierarchy_reports_nothing():
    features = [
        feature("gene-ok", "gene", 100, 900),
        feature("rna-ok", "mRNA", 100, 900, "gene-ok"),
        feature("exon-ok", "exon", 200, 400, "rna-ok"),
    ]
    assert find_violations(features) == []


# ---------------------------------------------------------------------------
# Repair
# ---------------------------------------------------------------------------


def test_repair_recovers_the_genbank_truth(zeam_mp186):
    """Bottom-up repair reaches 50,490..548,772, which is exactly the extent of
    the GenBank join for this gene."""
    repaired, repairs = repair_envelopes(zeam_mp186)
    envelopes = {f.source_id: (f.start, f.end) for f in repaired}
    assert envelopes["rna-ZeamMp186"] == (50490, 548772)
    assert envelopes["gene-ZeamMp186"] == (50490, 548772)
    assert len(repairs) == 2
    assert find_violations(repaired) == []


def test_repair_must_be_bottom_up(zeam_mp186):
    """The trap. Both the mRNA and the gene envelope are wrong, so repairing
    the gene from its immediate child -- the still-broken mRNA -- gives
    50,490..267,232: plausible, and wrong by 281 kb. A top-down pass produces a
    number that looks fixed and is not, which is worse than refusing.
    """
    repaired, _repairs = repair_envelopes(zeam_mp186)
    gene = next(f for f in repaired if f.source_id == "gene-ZeamMp186")
    assert gene.end != 267232, "gene repaired from the broken mRNA"
    assert gene.end == 548772


def test_repair_only_widens(zeam_mp186):
    """A repair may never shrink an envelope. Narrowing would discard a claim
    the file makes; widening only admits children it already claims."""
    original = {f.source_id: (f.start, f.end) for f in zeam_mp186}
    repaired, _repairs = repair_envelopes(zeam_mp186)
    for item in repaired:
        old_start, old_end = original[item.source_id]
        assert item.start <= old_start
        assert item.end >= old_end


def test_leaves_and_correct_parents_are_untouched(zeam_mp186):
    repaired, repairs = repair_envelopes(zeam_mp186)
    assert {r.feature_id for r in repairs} == {"rna-ZeamMp186", "gene-ZeamMp186"}
    exons_before = {f.source_id: (f.start, f.end) for f in zeam_mp186 if f.source_type == "exon"}
    exons_after = {f.source_id: (f.start, f.end) for f in repaired if f.source_type == "exon"}
    assert exons_before == exons_after


def test_multi_segment_parents_are_never_repaired():
    """Widening one segment of several means choosing which, and the file gives
    no basis for that choice. It stays a violation and refuses the load."""
    features = [
        feature("gene-multi", "gene", 0, 0, segments=[(100, 200), (400, 500)]),
        feature("rna-multi", "mRNA", 100, 900, "gene-multi"),
    ]
    repaired, repairs = repair_envelopes(features)
    assert repairs == []
    assert len(find_violations(repaired)) == 1


# ---------------------------------------------------------------------------
# Policy
# ---------------------------------------------------------------------------


def test_default_policy_refuses_without_repairing(zeam_mp186):
    features, repairs, remaining = apply_policy(zeam_mp186, repair=False)
    assert repairs == []
    assert len(remaining) == 4
    assert [f.start for f in features] == [f.start for f in zeam_mp186]


def test_repair_policy_clears_the_violations(zeam_mp186):
    _features, repairs, remaining = apply_policy(zeam_mp186, repair=True)
    assert len(repairs) == 2
    assert remaining == []


def test_a_repair_that_does_not_clear_still_refuses():
    """The re-check is the point. If repair leaves a violation standing, the
    cause is not a truncated envelope, and proceeding would be the silent
    partial success this project exists to prevent."""
    features = [
        feature("gene-multi", "gene", 0, 0, segments=[(100, 200), (400, 500)]),
        feature("rna-multi", "mRNA", 100, 900, "gene-multi"),
    ]
    _features, repairs, remaining = apply_policy(features, repair=True)
    assert repairs == []
    assert len(remaining) == 1


def test_report_names_every_change(zeam_mp186):
    _features, repairs, remaining = apply_policy(zeam_mp186, repair=True)
    text = format_report(repairs, remaining)
    assert "50,490..548,772" in text
    assert "281,540" in text
    assert "recomputed from children" in text


def test_clean_input_needs_no_policy():
    features = [
        feature("gene-ok", "gene", 100, 900),
        feature("rna-ok", "mRNA", 200, 800, "gene-ok"),
    ]
    result, repairs, remaining = apply_policy(features, repair=False)
    assert repairs == [] and remaining == []
    assert len(result) == 2


# ---------------------------------------------------------------------------
# Wired into `load`: the path a user actually takes
# ---------------------------------------------------------------------------

import json  # noqa: E402
from pathlib import Path as _Path  # noqa: E402

from click.testing import CliRunner  # noqa: E402
from cgload.cli.main import cli  # noqa: E402
from cgload.config import load_config  # noqa: E402
from cgload.db.engine import build_engine  # noqa: E402
from cgload.db.init import create_schema  # noqa: E402
from cgload.db.schema import assembly as assembly_table  # noqa: E402
from cgload.db.schema import feature as feature_table  # noqa: E402
from cgload.loader import LoadError, load_organism  # noqa: E402
from cgload.parsers import profiles as profile_table  # noqa: E402
from cgload.parsers.normalise import normalise  # noqa: E402
from cgload.parsers.tokenizer import tokenize  # noqa: E402
from sqlalchemy import select  # noqa: E402

FIXTURES = _Path(__file__).parent / "fixtures"

BROKEN = FIXTURES / "broken" / "trans_spliced_envelope.gff3"

#: The correct span for nad1, taken from the GenBank release of the same
#: assembly -- which is the only reason we know the pieces are right and the
#: envelope is wrong.
GENBANK_TRUTH = (50490, 548772)


def _prepare(tmp_path):
    (tmp_path / BROKEN.name).write_text(BROKEN.read_text())
    fasta = tmp_path / "a.fa"
    with fasta.open("w") as handle:
        handle.write(">NC_007982.1\n")
        handle.write((("ACGT" * 15) + "\n") * (600000 // 60 + 1))
    config = tmp_path / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Zea"\nspecies = "mays"\nstrain = "B73"\n\n'
        '[assembly]\nname = "mito"\nfasta = "a.fa"\n\n'
        f'[[files.annotation]]\npath = "{BROKEN.name}"\n'
    )
    return config


def _engine(tmp_path):
    built = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(built)
    return built


def test_a_broken_envelope_stops_the_load_by_default(tmp_path):
    """Refusing is right by default: loading known-wrong coordinates silently is
    worse than failing."""
    engine = _engine(tmp_path)
    with pytest.raises(LoadError, match="lie outside the parent"):
        load_organism(engine, load_config(_prepare(tmp_path)), tmp_path / "data")

    with engine.connect() as connection:
        assert (
            connection.execute(select(assembly_table)).first() is None
        ), "the load must roll back, not partially store"


def test_the_refusal_names_the_features_and_the_flag(tmp_path):
    engine = _engine(tmp_path)
    with pytest.raises(LoadError) as raised:
        load_organism(engine, load_config(_prepare(tmp_path)), tmp_path / "data")
    message = str(raised.value)
    assert "rna-ZeamMp186" in message
    assert "--repair-envelopes" in message
    assert "recomputed from their children" not in message, (
        "the refusal path must not claim repairs it did not make"
    )


def test_repair_stores_the_span_genbank_agrees_with(tmp_path):
    """The whole point, end to end. The repaired gene must match the GenBank
    release of the same assembly exactly -- if it matched the transcript's own
    wrong number instead, this passes 50,490..267,232 and looks plausible."""
    engine = _engine(tmp_path)
    report = load_organism(
        engine,
        load_config(_prepare(tmp_path)),
        tmp_path / "data",
        repair_envelopes=True,
    )

    with engine.connect() as connection:
        rows = dict(
            connection.execute(
                select(feature_table.c.source_id, feature_table.c.start).where(
                    feature_table.c.source_type == "gene"
                )
            ).all()
        )
        spans = connection.execute(
            select(feature_table.c.start, feature_table.c.end).where(
                feature_table.c.source_id == "gene-ZeamMp186"
            )
        ).one()
    assert (spans.start, spans.end) == GENBANK_TRUTH
    assert rows["gene-ok"] == 1000, "a correct gene must not be touched"
    assert len(report.envelope_repairs) == 2


def test_every_repair_is_declared_in_the_manifest(tmp_path):
    """A repaired load is declared, not merely successful: someone reading the
    database later has to be able to find out that cgload changed coordinates
    the provider published."""
    engine = _engine(tmp_path)
    report = load_organism(
        engine,
        load_config(_prepare(tmp_path)),
        tmp_path / "data",
        repair_envelopes=True,
    )
    manifest = json.loads((report.data_directory / "manifest.json").read_text())
    repairs = {entry["feature"]: entry for entry in manifest["envelope_repairs"]}

    assert set(repairs) == {"gene-ZeamMp186", "rna-ZeamMp186"}
    assert repairs["gene-ZeamMp186"]["from"] == [266974, 267232]
    assert repairs["gene-ZeamMp186"]["to"] == list(GENBANK_TRUTH)


def test_a_clean_load_records_no_repairs(tmp_path):
    engine = _engine(tmp_path)
    """The manifest key exists either way, so 'no repairs' is a recorded fact
    rather than a missing key."""
    fixture = FIXTURES / "pair_fgraminearum_NC_026477.gff3"
    (tmp_path / fixture.name).write_text(fixture.read_text())
    lengths: dict[str, int] = {}
    for line in fixture.read_text().splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) >= 9:
            lengths[columns[0]] = max(lengths.get(columns[0], 0), int(columns[4]))
    fasta = tmp_path / "b.fa"
    with fasta.open("w") as handle:
        for name, longest in lengths.items():
            handle.write(f">{name}\n")
            handle.write((("ACGT" * 15) + "\n") * ((longest + 1000) // 60 + 1))
    config = tmp_path / "clean.toml"
    config.write_text(
        '[organism]\ngenus = "F"\nspecies = "g"\nstrain = "S"\n\n[assembly]\n'
        f'name = "c"\nfasta = "{fasta.name}"\n\n[[files.annotation]]\n'
        f'path = "{fixture.name}"\n'
    )
    report = load_organism(
        engine, load_config(config), tmp_path / "data", repair_envelopes=True
    )
    assert report.envelope_repairs == []
    manifest = json.loads((report.data_directory / "manifest.json").read_text())
    assert manifest["envelope_repairs"] == []


def test_the_repair_flag_is_reachable_from_the_command(tmp_path):
    engine = _engine(tmp_path)
    config = _prepare(tmp_path)
    runner = CliRunner()

    refused = runner.invoke(
        cli, ["load", "--config", str(config), "--data-dir", str(tmp_path / "d1"),
              "--db-url", str(engine.url)]
    )
    assert refused.exit_code == 2
    assert "--repair-envelopes" in refused.output

    repaired = runner.invoke(
        cli, ["load", "--config", str(config), "--data-dir", str(tmp_path / "d2"),
              "--db-url", str(engine.url), "--repair-envelopes"]
    )
    assert repaired.exit_code == 0, repaired.output
    assert "repaired envelope" in repaired.output


def test_no_correct_fixture_reports_a_violation():
    """A check that fires on correct data is worse than no check. Every real
    fixture in the repository must come back clean."""
    for path in sorted(FIXTURES.iterdir()):
        if path.suffix not in (".gff3", ".gbff"):
            continue
        if path.suffix == ".gbff":
            from cgload.parsers.genbank import tokenize_genbank

            tokenized = tokenize_genbank(path)
        else:
            tokenized = tokenize(path)
        detection = profile_table.detect(tokenized)
        result = normalise(
            tokenized.features, detection.profile, comments=tokenized.comments
        )
        assert find_violations(result.features) == [], path.name
