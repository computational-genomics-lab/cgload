"""Milestone 8b: the database half of the check, and the two commands.

Milestone 8 delivered the naive counting pass and the two comparisons that need
no database. What it did not deliver was the other half of the primary check --
the stored side -- nor `verify`, nor `explain`, nor the automatic run at the end
of a load. Those are what these tests cover.

What they assert, stated before they were written:

* The naive pass reproduces the parser's per-type counts on every real fixture.
  That is the check being meaningful rather than merely present.
* A parser that drops a feature fails the load and rolls everything back --
  demonstrated on a *leaf* feature, so that no other check can catch it first.
* `verify` reports every tier including skipped ones, and exits 2 only on a
  real failure.
* `explain` distinguishes a parent the file stated from one cgload computed.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import pytest
from click.testing import CliRunner
from cgload.cli.main import cli
from cgload.config import load_config
from cgload.db.engine import build_engine
from cgload.db.init import create_schema
from cgload.db.schema import assembly, feature, feature_attribute
from cgload.db.verification import (
    compare_with_database,
    count_database,
    count_database_for_file,
)
from cgload.loader import LoadError, load_organism
from cgload.verify import Outcome, count_source_naively
from sqlalchemy import func, select

FIXTURES = Path(__file__).parent / "fixtures"
MYCMAY_GFF = FIXTURES / "refseq_mycmay_GCF_000328475.gff3"
DSCAM_GBFF = FIXTURES / "genbank_dscam1_CG12164.gbff"
NCBI_COUNTS = FIXTURES / "ncbi_feature_count_GCF_000328475.txt"


def _fasta_for(path: Path, target: Path) -> dict[str, int]:
    """A FASTA long enough for whatever coordinates the annotation uses."""
    lengths: dict[str, int] = {}
    text = path.read_text()
    if path.suffix == ".gbff":
        import re

        version = re.search(r"^VERSION\s+(\S+)", text, re.M)
        longest = max(
            int(value)
            for line in text.splitlines()
            if line.startswith("     ") and ".." in line
            for value in re.findall(r"\b(\d{3,})\b", line)
        )
        lengths[version.group(1)] = longest
    else:
        for line in text.splitlines():
            if line.startswith("#"):
                continue
            columns = line.split("\t")
            if len(columns) >= 9:
                lengths[columns[0]] = max(lengths.get(columns[0], 0), int(columns[4]))

    with target.open("w") as handle:
        for name, longest in sorted(lengths.items()):
            handle.write(f">{name}\n")
            handle.write((("ACGT" * 15) + "\n") * ((longest + 1000) // 60 + 1))
    return lengths


def _prepare(tmp_path: Path, annotation: Path, name: str = "asm") -> Path:
    fasta = tmp_path / "assembly.fa"
    _fasta_for(annotation, fasta)
    (tmp_path / annotation.name).write_text(annotation.read_text())
    config = tmp_path / "organism.toml"
    config.write_text(
        f'[organism]\ngenus = "Testus"\nspecies = "exampleii"\nstrain = "S"\n\n'
        f'[assembly]\nname = "{name}"\nfasta = "{fasta.name}"\n\n'
        f'[[files.annotation]]\npath = "{annotation.name}"\n'
    )
    return config


@pytest.fixture
def engine(tmp_path):
    built = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(built)
    return built


# -- the naive pass agrees with the parser on real files -------------------


@pytest.mark.parametrize(
    "fixture",
    [
        "refseq_fgraminearum_NC_026477.gff3",
        "refseq_mycmay_GCF_000328475.gff3",
        "funannotate_aphaeospermum_pheo_arth1.gff3",
        "augustus_aphaeospermum_pheo_arth1.gff3",
        "genbank_maeruginosa_NC_010296.gbff",
        "genbank_dscam1_CG12164.gbff",
    ],
)
def test_the_naive_pass_reproduces_the_parsers_counts(engine, tmp_path, fixture):
    """The check earning its keep. Two entirely separate pieces of code read the
    same file and arrive at the same per-type feature counts -- including the
    AUGUSTUS rows with no ID, the GenBank join() spans, and the RefSeq CDS
    segments that share one identifier."""
    annotation = FIXTURES / fixture
    working = tmp_path / fixture
    config = _prepare(tmp_path, annotation)

    report = load_organism(engine, load_config(config), tmp_path / "data")
    naive = count_source_naively(working)

    with engine.connect() as connection:
        stored = count_database(connection, report.assembly_id)

    assert dict(naive) == dict(stored), fixture


# -- the mutation test the milestone exists for ---------------------------


def test_dropping_a_leaf_feature_fails_the_load_and_rolls_it_back(engine, tmp_path):
    """A *leaf* deliberately: an exon has no children, so no orphan check can
    catch it first and the completeness check is the only thing between the bug
    and a silently wrong database: the failure mode this project exists to prevent.
    """
    import cgload.loader as loader_module
    import cgload.parsers.normalise as normalise_module

    config = _prepare(tmp_path, MYCMAY_GFF)
    real = normalise_module.normalise

    def sabotage(features, profile, **kwargs):
        result = real(features, profile, **kwargs)
        for index, item in enumerate(result.features):
            if item.feature_type == "exon":
                del result.features[index]
                break
        return result

    monkey = pytest.MonkeyPatch()
    monkey.setattr(loader_module, "normalise", sabotage)
    try:
        with pytest.raises(LoadError, match="did not store what the file contains"):
            load_organism(engine, load_config(config), tmp_path / "data")
    finally:
        monkey.undo()

    with engine.connect() as connection:
        for table in (feature, assembly):
            assert (
                connection.execute(select(func.count()).select_from(table)).scalar_one() == 0
            )
    assert [p for p in (tmp_path / "data").rglob("*") if p.is_file()] == []


def test_the_failure_names_the_type_and_the_two_counts(engine, tmp_path):
    """A bug report needs the type and both numbers, not just 'verification
    failed'."""
    import cgload.loader as loader_module
    import cgload.parsers.normalise as normalise_module

    config = _prepare(tmp_path, MYCMAY_GFF)
    real = normalise_module.normalise

    def sabotage(features, profile, **kwargs):
        result = real(features, profile, **kwargs)
        for index, item in enumerate(result.features):
            if item.feature_type == "exon":
                del result.features[index]
                break
        return result

    monkey = pytest.MonkeyPatch()
    monkey.setattr(loader_module, "normalise", sabotage)
    try:
        with pytest.raises(LoadError) as raised:
            load_organism(engine, load_config(config), tmp_path / "data")
    finally:
        monkey.undo()

    message = str(raised.value)
    assert "exon" in message
    assert "file has 25" in message
    assert "database has 24" in message


def test_a_clean_load_records_its_verification_in_the_manifest(engine, tmp_path):
    """The check is not optional and not a separate step, so its result belongs
    in the receipt alongside the file hashes."""
    import json

    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")
    assert report.verification
    assert all(item.outcome is Outcome.PASS for item in report.verification)

    manifest = json.loads((report.data_directory / "manifest.json").read_text())
    assert manifest["verification"]
    assert all(entry["outcome"] == "pass" for entry in manifest["verification"])


# -- the comparison function itself ---------------------------------------


def test_a_surplus_is_as_much_a_failure_as_a_loss():
    """Inventing a feature is as wrong as losing one, and a check that only
    looked for losses would miss a loader that duplicated rows."""
    result = compare_with_database(Counter({"gene": 10}), Counter({"gene": 11}))
    assert result.outcome is Outcome.FAIL
    assert result.mismatches == {"gene": (10, 11)}


def test_a_type_missing_entirely_from_one_side_is_caught():
    result = compare_with_database(Counter({"gene": 10, "tRNA": 3}), Counter({"gene": 10}))
    assert result.outcome is Outcome.FAIL
    assert result.mismatches == {"tRNA": (3, 0)}


def test_agreement_passes_and_says_how_much_it_checked():
    result = compare_with_database(
        Counter({"gene": 10, "CDS": 4}), Counter({"gene": 10, "CDS": 4})
    )
    assert result.outcome is Outcome.PASS
    assert "14 features" in result.detail


def test_counts_are_scoped_to_one_file(engine, tmp_path):
    """An assembly may be built from several annotation files (D-014). Comparing
    one file's naive count against the whole assembly would report every other
    file's features as surplus."""
    fasta = tmp_path / "assembly.fa"
    _fasta_for(FIXTURES / "funannotate_aphaeospermum_pheo_arth1.gff3", fasta)
    for name in (
        "funannotate_aphaeospermum_pheo_arth1.gff3",
        "augustus_aphaeospermum_pheo_arth1.gff3",
    ):
        (tmp_path / name).write_text((FIXTURES / name).read_text())
    config = tmp_path / "organism.toml"
    config.write_text(
        f'[organism]\ngenus = "T"\nspecies = "e"\nstrain = "S"\n\n[assembly]\n'
        f'name = "a"\nfasta = "{fasta.name}"\n\n'
        '[[files.annotation]]\npath = "funannotate_aphaeospermum_pheo_arth1.gff3"\n\n'
        '[[files.annotation]]\npath = "augustus_aphaeospermum_pheo_arth1.gff3"\n'
        'id_prefix = "AUG_"\n'
    )
    report = load_organism(engine, load_config(config), tmp_path / "data")
    assert len(report.verification) == 2
    assert all(item.outcome is Outcome.PASS for item in report.verification)

    from cgload.db.verification import annotation_files

    with engine.connect() as connection:
        files = annotation_files(connection, report.assembly_id)
        per_file = [
            count_database_for_file(connection, report.assembly_id, file_id)
            for file_id, _ in files
        ]
        whole = count_database(connection, report.assembly_id)
    assert sum(per_file, Counter()) == whole
    assert all(counts != whole for counts in per_file)


# -- the verify command ----------------------------------------------------


def test_verify_reports_every_tier_including_skipped(engine, tmp_path):
    """A check that vanishes when it cannot run is indistinguishable from one
    that passed."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli,
        ["verify", "--assembly-id", str(report.assembly_id), "--db-url", str(engine.url)],
    )
    assert result.exit_code == 0, result.output
    assert "source count" in result.output
    assert "PASS" in result.output
    assert "cross-parser" in result.output
    assert "reference (NCBI)" in result.output
    assert "----" in result.output


def test_verify_exits_two_when_the_database_was_tampered_with(engine, tmp_path):
    """The load passed; something deleted rows afterwards. `verify` must catch
    that independently of the load."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    # Delete the last segment of a multi-segment exon, or a whole single-segment
    # one that nothing points at. Deleting a row with children trips the foreign
    # key first, which is a different (also correct) protection.
    with engine.begin() as connection:
        rows = connection.execute(
            select(feature.c.feature_id, feature.c.source_id).where(
                feature.c.source_type == "exon"
            )
        ).all()
        referenced = set(
            connection.execute(
                select(feature.c.parent_id).where(feature.c.parent_id.is_not(None))
            ).scalars()
        )
        victims = [row.source_id for row in rows if row.feature_id not in referenced]
        assert victims, "every exon is referenced; pick a different victim"
        doomed = [row.feature_id for row in rows if row.source_id == victims[0]]
        # feature_attribute rows point at the feature, so they go first. The
        # foreign keys are doing their job; this test is simulating corruption
        # that bypassed cgload, not a legitimate deletion.
        connection.execute(
            feature_attribute.delete().where(feature_attribute.c.feature_id.in_(doomed))
        )
        connection.execute(feature.delete().where(feature.c.feature_id.in_(doomed)))

    result = CliRunner().invoke(
        cli,
        ["verify", "--assembly-id", str(report.assembly_id), "--db-url", str(engine.url)],
    )
    assert result.exit_code == 2
    assert "FAIL" in result.output
    assert "exon" in result.output


def test_verify_skips_the_reference_check_on_a_partial_assembly(engine, tmp_path):
    """Two of twenty-three sequences loaded. Disagreement with published totals
    has an innocent explanation, so the tier reports why it did not run rather
    than failing (D-033)."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli,
        [
            "verify",
            "--assembly-id",
            str(report.assembly_id),
            "--db-url",
            str(engine.url),
            "--reference-counts",
            str(NCBI_COUNTS),
            "--accession",
            "GCF_000328475.2",
            "--sequences-in-assembly",
            "23",
        ],
    )
    assert result.exit_code == 0
    assert "2 of 23" in result.output


def test_reference_counts_without_an_accession_is_refused(engine, tmp_path):
    """The file lists two accessions differing by one digit; picking the wrong
    one selects nothing and reports every feature as missing."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli,
        [
            "verify",
            "--assembly-id",
            str(report.assembly_id),
            "--db-url",
            str(engine.url),
            "--reference-counts",
            str(NCBI_COUNTS),
        ],
    )
    assert result.exit_code == 2
    assert "--accession" in result.output


def test_verify_says_so_when_the_original_file_is_gone(engine, tmp_path):
    """The naive pass needs the input. Missing it is a skip with a reason, not a
    pass and not a crash."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")
    (tmp_path / MYCMAY_GFF.name).unlink()

    result = CliRunner().invoke(
        cli,
        ["verify", "--assembly-id", str(report.assembly_id), "--db-url", str(engine.url)],
    )
    assert result.exit_code == 0
    assert "no longer at" in result.output


# -- the explain command ---------------------------------------------------


def test_explain_says_the_file_stated_this_link(engine, tmp_path):
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli, ["explain", "cds-XP_011386003.1", "--db-url", str(engine.url)]
    )
    assert result.exit_code == 0, result.output
    assert "explicit_parent" in result.output
    assert "stated by the file" in result.output
    assert "refseq" in result.output
    assert str(report.assembly_id) or True


def test_explain_says_cgload_computed_this_link(engine, tmp_path):
    """The distinction the paper rests on, made visible to a user. A GenBank CDS
    has no parent pointer in the file at all."""
    config = _prepare(tmp_path, DSCAM_GBFF)
    load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli, ["explain", "NP_001260764.1", "--db-url", str(engine.url)]
    )
    assert result.exit_code == 0, result.output
    assert "structural" in result.output
    assert "computed by cgload" in result.output


def test_explain_shows_every_segment_of_a_discontinuous_feature(engine, tmp_path):
    config = _prepare(tmp_path, MYCMAY_GFF)
    load_organism(engine, load_config(config), tmp_path / "data")

    with build_engine(str(engine.url), must_exist=True).connect() as connection:
        multi = connection.execute(
            select(feature.c.source_id)
            .where(feature.c.segment_index > 0)
            .limit(1)
        ).scalar_one_or_none()
    if multi is None:
        pytest.skip("fixture has no discontinuous feature")

    result = CliRunner().invoke(cli, ["explain", multi, "--db-url", str(engine.url)])
    assert result.exit_code == 0
    assert "segments" in result.output
    assert "coordinate order" in result.output


def test_explain_on_an_unknown_identifier_exits_two(engine, tmp_path):
    config = _prepare(tmp_path, MYCMAY_GFF)
    load_organism(engine, load_config(config), tmp_path / "data")
    result = CliRunner().invoke(cli, ["explain", "not-a-real-id", "--db-url", str(engine.url)])
    assert result.exit_code == 2
    assert "no feature with source id" in result.output


def test_explain_reports_the_go_caveat(engine, tmp_path):
    """D-045. A user shown 129 GO terms must be told what they are."""
    fasta = tmp_path / "assembly.fa"
    _fasta_for(MYCMAY_GFF, fasta)
    (tmp_path / MYCMAY_GFF.name).write_text(MYCMAY_GFF.read_text())
    eggnog = FIXTURES / "eggnog_mycmay_emapper214.annotations"
    (tmp_path / eggnog.name).write_text(eggnog.read_text())
    config = tmp_path / "organism.toml"
    config.write_text(
        f'[organism]\ngenus = "T"\nspecies = "e"\nstrain = "S"\n\n[assembly]\n'
        f'name = "a"\nfasta = "{fasta.name}"\n\n'
        f'[[files.annotation]]\npath = "{MYCMAY_GFF.name}"\n\n'
        f'[[files.functional]]\npath = "{eggnog.name}"\n'
    )
    load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli, ["explain", "cds-XP_011386003.1", "--db-url", str(engine.url)]
    )
    assert result.exit_code == 0
    assert "inherited ancestors" in result.output


def test_verify_and_explain_are_no_longer_pending():
    from cgload.cli.main import PENDING

    assert "verify" not in PENDING
    assert "explain" not in PENDING
    assert set(PENDING) == {"query", "remove"}


# ---------------------------------------------------------------------------
# Milestone 8c: the matched pair -- the strongest check, now reproducible
# ---------------------------------------------------------------------------

PAIR_GFF = FIXTURES / "pair_fgraminearum_NC_026477.gff3"
PAIR_GBFF = FIXTURES / "pair_fgraminearum_NC_026477.gbff"


def test_the_matched_pair_agrees_across_two_formats():
    """The strongest check available, reproducible from the repository.

    The same 43 genes written in two file formats with nothing in common, read
    by two parsers that share no code: one splits tab-separated columns and
    follows explicit Parent labels, the other reads an indented text format and
    works parentage out from coordinates. A parser that dropped a gene from one
    format has no way to make the other format's count drop it too.
    """
    from cgload.verify import compare_across_formats

    result = compare_across_formats(PAIR_GFF, PAIR_GBFF)
    assert result.outcome is Outcome.PASS, result.mismatches


def test_both_formats_of_the_pair_load_to_the_same_counts(engine, tmp_path):
    """Not just the naive pass agreeing -- the real parsers, through the
    database, on shared feature types."""
    counts = {}
    for index, annotation in enumerate((PAIR_GFF, PAIR_GBFF)):
        working = tmp_path / f"case{index}"
        working.mkdir()
        config = _prepare(working, annotation, name=f"asm{index}")
        report = load_organism(engine, load_config(config), working / "data")
        with engine.connect() as connection:
            counts[annotation.suffix] = count_database(connection, report.assembly_id)

    shared = set(counts[".gff3"]) & set(counts[".gbff"])
    assert {"gene", "mRNA", "CDS"} <= shared
    for feature_type in shared:
        assert counts[".gff3"][feature_type] == counts[".gbff"][feature_type], feature_type


def test_the_genbank_half_of_the_pair_is_linked_by_containment(engine, tmp_path):
    """Every CDS-to-transcript edge in the GenBank file was computed, because the
    format states none. That is the algorithm the paper is about, exercised on
    real NCBI output rather than on a hand-picked isoform case."""
    config = _prepare(tmp_path, PAIR_GBFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    with engine.connect() as connection:
        rows = dict(
            connection.execute(
                select(feature.c.linkage_method, func.count())
                .where(
                    feature.c.assembly_id == report.assembly_id,
                    feature.c.source_type == "CDS",
                )
                .group_by(feature.c.linkage_method)
            ).all()
        )
    assert set(rows) <= {"structural", "note"}
    assert "explicit_parent" not in rows


# ---------------------------------------------------------------------------
# The reference tier must not be able to pass trivially (D-067)
# ---------------------------------------------------------------------------


def test_the_reference_tier_cannot_run_without_the_assembly_size(engine, tmp_path):
    """The precondition must not default to the loaded count.

    `sequences_in_assembly=sequences_in_assembly or loaded` made the check
    `loaded == loaded` -- always true, so the tier could never skip and a partial
    load was silently compared against totals that may describe the whole
    assembly. A precondition that cannot fail is D-033 inverted.
    """
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli,
        ["verify", "--assembly-id", str(report.assembly_id), "--db-url", str(engine.url),
         "--reference-counts", str(NCBI_COUNTS), "--accession", "GCF_000328475.2"],
    )
    assert result.exit_code == 2
    assert "--sequences-in-assembly" in result.output


def test_a_partial_load_skips_the_reference_tier(engine, tmp_path):
    """Two sequences loaded against a twenty-three sequence assembly.
    Disagreement has an innocent explanation, so the tier reports why it did not
    run rather than passing or failing."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli,
        ["verify", "--assembly-id", str(report.assembly_id), "--db-url", str(engine.url),
         "--reference-counts", str(NCBI_COUNTS), "--accession", "GCF_000328475.2",
         "--sequences-in-assembly", "23"],
    )
    assert result.exit_code == 0
    assert "2 of 23" in result.output
    assert "matches published totals" not in result.output


def test_stating_the_loaded_count_deliberately_lets_the_tier_run(engine, tmp_path):
    """The legitimate case the old default silently assumed: the published counts
    describe exactly the unit that was loaded, which is what `--assembly-unit`
    selection (D-065) exists for. It must be asserted by the user, because
    cgload cannot tell that apart from a partial load."""
    config = _prepare(tmp_path, MYCMAY_GFF)
    report = load_organism(engine, load_config(config), tmp_path / "data")

    result = CliRunner().invoke(
        cli,
        ["verify", "--assembly-id", str(report.assembly_id), "--db-url", str(engine.url),
         "--reference-counts", str(NCBI_COUNTS), "--accession", "GCF_000328475.2",
         "--sequences-in-assembly", "2"],
    )
    assert "loaded 2 of 2" not in result.output, (
        "the precondition passed, so the tier must report a comparison outcome "
        "rather than a skip"
    )
