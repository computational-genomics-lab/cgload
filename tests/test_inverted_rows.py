"""Inverted coordinate rows: refused by default, repaired under the flag.

*Zea mays* B73 GCF_902167145.1 contains one row where the file writes start
after end:

    NC_007982.1  RefSeq  mRNA  691776  267232  .  ?  .  ID=rna-ZeamMp017;...

This is the same NCBI trans-splicing converter bug as the 34 too-small
envelopes, arriving in the one form a per-row check happens to catch. Neither
number is an extent: 267,232 is the end of one segment, 691,776 is not a
coordinate on this molecule at all (NC_007982.1 is 569,630 bp). The feature's
children span 122,146..548,772, which is the correct answer.

Before this change, `--repair-envelopes` could not help: the tokenizer refused
the row long before any envelope logic ran, so cgload repaired one form of the
bug and refused the other while both were recoverable from the same evidence.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from cgload.envelope import apply_policy, find_violations
from cgload.parsers.normalise import NormalisedFeature, Segment
from cgload.parsers.tokenizer import TokenizeError, tokenize

SEQUENCE_LENGTH = 569630  # NC_007982.1


def gff3(tmp_path: Path, body: str) -> Path:
    path = tmp_path / "sample.gff3"
    path.write_text("##gff-version 3\n" + body)
    return path


def feature(source_id, source_type, start, end, parent=None, strand="-", inverted=False):
    return NormalisedFeature(
        source_id=source_id,
        source_type=source_type,
        feature_type=source_type,
        seqid="NC_007982.1",
        strand=strand,
        source_program="RefSeq",
        segments=[
            Segment(
                line_number=0,
                start=start,
                end=end,
                phase=None,
                score=None,
                inverted=inverted,
            )
        ],
        parent_source_id=parent,
    )


@pytest.fixture
def zeam_mp017():
    """`ZeamMp017` as maize publishes it, with its five real exon spans."""
    return [
        feature("gene-ZeamMp017", "gene", 266974, 267232),
        feature(
            "rna-ZeamMp017", "mRNA", 691776, 267232, "gene-ZeamMp017",
            strand="?", inverted=True,
        ),
        feature("exon-1", "exon", 122146, 122530, "rna-ZeamMp017"),
        feature("exon-2", "exon", 320928, 321010, "rna-ZeamMp017"),
        feature("exon-3", "exon", 322404, 322595, "rna-ZeamMp017"),
        feature("exon-4", "exon", 548714, 548772, "rna-ZeamMp017"),
        feature("exon-5", "exon", 266974, 267232, "rna-ZeamMp017"),
    ]


# ---------------------------------------------------------------------------
# The tokenizer gate
# ---------------------------------------------------------------------------


def test_inverted_row_is_refused_by_default(tmp_path):
    """Unchanged behaviour: without the flag, an inverted row stops the load."""
    body = "chr1\tRefSeq\tCDS\t500\t100\t.\t+\t.\tID=c1\n"
    with pytest.raises(TokenizeError, match="precedes start"):
        tokenize(gff3(tmp_path, body))


def test_tolerated_row_keeps_its_coordinates_exactly_as_written(tmp_path):
    """Nothing is swapped. A swap would invent an extent the file never
    claimed, and 691,776..267,232 reversed is 267,232..691,776 -- which runs
    122,146 bases past the end of the molecule and is not the answer either."""
    body = "chr1\tRefSeq\tCDS\t500\t100\t.\t+\t.\tID=c1\n"
    result = tokenize(gff3(tmp_path, body), tolerate_inverted=True)
    row = result.features[0]
    assert (row.start, row.end) == (500, 100)
    assert row.inverted is True


def test_ordinary_rows_are_not_flagged(tmp_path):
    body = "chr1\tRefSeq\tCDS\t100\t500\t.\t+\t.\tID=c1\n"
    result = tokenize(gff3(tmp_path, body), tolerate_inverted=True)
    assert result.features[0].inverted is False


# ---------------------------------------------------------------------------
# Detection and repair
# ---------------------------------------------------------------------------


def test_inverted_feature_is_a_violation_on_its_own(zeam_mp017):
    """It contradicts itself before any parent is consulted, so it is reported
    even when no containment rule is broken."""
    lone = [feature("solo", "mRNA", 691776, 267232, inverted=True)]
    assert len(find_violations(lone)) == 1


def test_repair_recomputes_purely_from_children(zeam_mp017):
    """The written pair is not an extent -- both numbers are segment
    coordinates -- so widening from them is meaningless. The children alone
    give 122,146..548,772."""
    _features, repairs, remaining = apply_policy(zeam_mp017, repair=True)
    assert remaining == []
    by_id = {r.feature_id: r for r in repairs}
    assert (by_id["rna-ZeamMp017"].new_start, by_id["rna-ZeamMp017"].new_end) == (
        122146,
        548772,
    )


def test_repaired_extent_lies_inside_the_molecule(zeam_mp017):
    """691,776 exceeds NC_007982.1's 569,630 bases, so the published row cannot
    be a coordinate pair on this sequence at all. The repair must land inside."""
    features, _repairs, _remaining = apply_policy(zeam_mp017, repair=True)
    for item in features:
        assert item.end <= SEQUENCE_LENGTH
        assert item.start >= 1


def test_repair_propagates_upward(zeam_mp017):
    """The gene envelope is also wrong, and is repaired from the repaired mRNA
    rather than the broken one -- the same bottom-up requirement as D-046."""
    features, _repairs, _remaining = apply_policy(zeam_mp017, repair=True)
    gene = next(f for f in features if f.source_id == "gene-ZeamMp017")
    assert (gene.start, gene.end) == (122146, 548772)


def test_the_inverted_flag_is_cleared_after_repair(zeam_mp017):
    """A repaired feature must not still be marked as contradicting itself, or
    the post-repair re-check would refuse a row it just fixed."""
    features, _repairs, remaining = apply_policy(zeam_mp017, repair=True)
    assert remaining == []
    mrna = next(f for f in features if f.source_id == "rna-ZeamMp017")
    assert all(not segment.inverted for segment in mrna.segments)


def test_default_policy_still_refuses(zeam_mp017):
    _features, repairs, remaining = apply_policy(zeam_mp017, repair=False)
    assert repairs == []
    assert remaining


def test_an_inverted_row_with_no_children_cannot_be_repaired():
    """There is nothing to recompute from, so it stays a violation and refuses
    the load rather than being silently accepted or guessed at."""
    orphan = [feature("solo", "mRNA", 691776, 267232, inverted=True)]
    _features, repairs, remaining = apply_policy(orphan, repair=True)
    assert repairs == []
    assert len(remaining) == 1


# ---------------------------------------------------------------------------
# Wired into `load`: --repair-envelopes has to reach the tokenizer
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
from sqlalchemy import select  # noqa: E402

WIRED_FIXTURES = _Path(__file__).parent / "fixtures"
INVERTED = WIRED_FIXTURES / "broken" / "inverted_row_envelope.gff3"

#: What ZeamMp017's children actually span, and therefore the only defensible
#: extent. Reversing the written pair gives 267,232..691,776, which runs 122,146
#: bases past the end of a 569,630 bp molecule -- a swap is not a guess, it is a
#: demonstrably wrong one (D-056).
CHILD_SPAN = (122146, 548772)
MOLECULE_LENGTH = 569630


def _setup(tmp_path):
    (tmp_path / INVERTED.name).write_text(INVERTED.read_text())
    fasta = tmp_path / "a.fa"
    with fasta.open("w") as handle:
        handle.write(">NC_007982.1\n")
        handle.write((("ACGT" * 15) + "\n") * (MOLECULE_LENGTH // 60 + 1))
    config = tmp_path / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Zea"\nspecies = "mays"\nstrain = "B73"\n\n'
        '[assembly]\nname = "mito"\nfasta = "a.fa"\n\n'
        f'[[files.annotation]]\npath = "{INVERTED.name}"\n'
    )
    engine = build_engine(f"sqlite:///{tmp_path / 'g.db'}")
    create_schema(engine)
    return engine, config


def test_an_inverted_row_is_refused_without_the_flag(tmp_path):
    """Unchanged default behaviour: the tokenizer rejects it at read time."""
    engine, config = _setup(tmp_path)
    with pytest.raises(Exception, match="precedes start"):
        load_organism(engine, load_config(config), tmp_path / "data")
    with engine.connect() as connection:
        assert connection.execute(select(assembly_table)).first() is None


def test_the_flag_reaches_the_tokenizer_not_only_the_envelope_stage(tmp_path):
    """The wiring bug this test exists for: an inverted row is refused at read
    time, long before any envelope logic runs, so a flag that only reached the
    envelope stage could not help the one maize row that needs it."""
    engine, config = _setup(tmp_path)
    report = load_organism(
        engine, load_config(config), tmp_path / "data", repair_envelopes=True
    )
    assert len(report.envelope_repairs) == 2


def test_the_repaired_extent_comes_from_the_children_not_from_a_swap(tmp_path):
    """A swap would give 267,232..691,776, which runs past the end of the
    molecule. The children's span is the only defensible answer."""
    engine, config = _setup(tmp_path)
    load_organism(engine, load_config(config), tmp_path / "data", repair_envelopes=True)

    with engine.connect() as connection:
        spans = dict(
            connection.execute(
                select(feature_table.c.source_id, feature_table.c.start).where(
                    feature_table.c.source_type.in_(("gene", "mRNA"))
                )
            ).all()
        )
        gene = connection.execute(
            select(feature_table.c.start, feature_table.c.end).where(
                feature_table.c.source_id == "gene-ZeamMp017"
            )
        ).one()
    assert (gene.start, gene.end) == CHILD_SPAN
    assert gene.end <= MOLECULE_LENGTH, "a swapped pair would exceed the molecule"
    assert spans["gene-ok"] == 1000, "a correct gene must not be touched"


def test_the_gene_is_repaired_from_the_repaired_transcript(tmp_path):
    """Bottom-up. Repairing the gene from the *broken* mRNA would give
    691,776..267,232 again, which is not even an interval."""
    engine, config = _setup(tmp_path)
    report = load_organism(
        engine, load_config(config), tmp_path / "data", repair_envelopes=True
    )
    by_feature = {item.feature_id: item for item in report.envelope_repairs}
    for identifier in ("rna-ZeamMp017", "gene-ZeamMp017"):
        repair = by_feature[identifier]
        assert (repair.new_start, repair.new_end) == CHILD_SPAN, identifier


def test_the_written_pair_is_preserved_in_the_manifest(tmp_path):
    """The provider's numbers are recorded verbatim, so an audit can see what
    the file actually said rather than only what cgload stored."""
    engine, config = _setup(tmp_path)
    report = load_organism(
        engine, load_config(config), tmp_path / "data", repair_envelopes=True
    )
    manifest = json.loads((report.data_directory / "manifest.json").read_text())
    entry = next(
        item for item in manifest["envelope_repairs"] if item["feature"] == "gene-ZeamMp017"
    )
    assert entry["from"] == [691776, 267232]
    assert entry["to"] == list(CHILD_SPAN)


def test_an_inverted_row_with_no_children_still_refuses(tmp_path):
    """Nothing to recompute from, so the flag cannot help and must not pretend
    to. The message says why rather than reusing the parent/child wording."""
    (tmp_path / "orphan.gff3").write_text(
        "##gff-version 3\n"
        "NC_007982.1\tRefSeq\tgene\t691776\t267232\t.\t?\t.\tID=gene-lonely;locus_tag=X\n"
        "NC_007982.1\tRefSeq\tgene\t1000\t2000\t.\t+\t.\tID=gene-ok2;locus_tag=OK2\n"
    )
    fasta = tmp_path / "a.fa"
    with fasta.open("w") as handle:
        handle.write(">NC_007982.1\n")
        handle.write((("ACGT" * 15) + "\n") * (MOLECULE_LENGTH // 60 + 1))
    config = tmp_path / "orphan.toml"
    config.write_text(
        '[organism]\ngenus = "Z"\nspecies = "m"\nstrain = "B"\n\n[assembly]\n'
        'name = "m2"\nfasta = "a.fa"\n\n[[files.annotation]]\npath = "orphan.gff3"\n'
    )
    engine = build_engine(f"sqlite:///{tmp_path / 'h.db'}")
    create_schema(engine)

    with pytest.raises(LoadError) as raised:
        load_organism(
            engine, load_config(config), tmp_path / "data", repair_envelopes=True
        )
    message = str(raised.value)
    assert "runs backwards" in message
    assert "no children to recompute" in message
    assert "lies outside" not in message, "self-contradiction must not use parent wording"


def test_the_command_reports_each_repair(tmp_path):
    engine, config = _setup(tmp_path)
    result = CliRunner().invoke(
        cli,
        ["load", "--config", str(config), "--data-dir", str(tmp_path / "d"),
         "--db-url", str(engine.url), "--repair-envelopes"],
    )
    assert result.exit_code == 0, result.output
    assert "691,776-267,232 -> 122,146-548,772" in result.output
