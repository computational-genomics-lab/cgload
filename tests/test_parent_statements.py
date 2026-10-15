"""D-092: the parent-edge tier counts the file's Parent statements itself.

Before this, the tier read only the database. Whether a feature "had a parent in
the file" was whatever the parser had recorded, so a parser that dropped a
``Parent`` attribute stored the feature as top-level and the tier agreed with
it. The naive pass now counts, per type, the features whose rows name a parent,
and the tier compares that with the links stored.
"""

from __future__ import annotations

import random
import sqlite3
from pathlib import Path

import cgload.parsers.tokenizer as tokenizer_module
from cgload.cli.main import cli
from cgload.verify import naive_counts
from click.testing import CliRunner

FIXTURES = Path(__file__).parent / "fixtures"
FUSARIUM = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
GENBANK = FIXTURES / "pair_fgraminearum_NC_026477.gbff"


def _fasta(path: Path, lengths: dict[str, int]) -> None:
    rnd = random.Random(7)
    with path.open("w") as handle:
        for name, length in lengths.items():
            handle.write(f">{name}\n")
            sequence = "".join(rnd.choice("ACGT") for _ in range(length))
            for index in range(0, length, 60):
                handle.write(sequence[index : index + 60] + "\n")


def _project(tmp_path: Path, source: Path) -> Path:
    gff = tmp_path / source.name
    gff.write_text(source.read_text())
    lengths = {}
    for line in gff.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            lengths[parts[1]] = int(parts[3])
    _fasta(tmp_path / "assembly.fa", lengths)
    config = tmp_path / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Fusarium"\nspecies = "graminearum"\nstrain = "PH-1"\n\n'
        '[assembly]\nname = "asm"\nfasta = "assembly.fa"\n\n'
        f'[[files.annotation]]\npath = "{gff.name}"\n'
    )
    return config


def _cli(tmp_path: Path, *args: str):
    url = f"sqlite:///{tmp_path / 'g.db'}"
    return CliRunner().invoke(cli, [*args[:1], "--db-url", url, *args[1:]])


def _load(tmp_path: Path, config: Path):
    assert _cli(tmp_path, "init").exit_code == 0
    return _cli(tmp_path, "load", "--config", str(config), "--data-dir", str(tmp_path / "data"))


def test_the_naive_pass_counts_features_that_state_a_parent(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t900\t.\t+\t.\tID=g1\n"
        "chr1\t.\tmRNA\t1\t900\t.\t+\t.\tID=t1;Parent=g1\n"
        # One CDS written on three rows is one feature with a parent, not three.
        "chr1\t.\tCDS\t1\t100\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\t.\tCDS\t200\t300\t.\t+\t2\tID=c1;Parent=t1\n"
        "chr1\t.\tCDS\t400\t500\t.\t+\t1\tID=c1;Parent=t1\n"
        # Rows with no ID count one by one, as they do in the feature count.
        "chr1\t.\tintron\t101\t199\t.\t+\t.\tParent=t1\n"
        "chr1\t.\tintron\t301\t399\t.\t+\t.\tParent=t1\n"
        # A space after the semicolon still names a parent.
        "chr1\t.\texon\t1\t100\t.\t+\t.\tID=e1; Parent=t1\n"
        # A parent mentioned only inside another attribute's value is not one.
        "chr1\t.\tregion\t1\t900\t.\t+\t.\tID=r1;Note=no Parent=here\n"
    )
    counted = naive_counts(gff)
    assert counted.features == {"gene": 1, "mRNA": 1, "CDS": 1, "intron": 2, "exon": 1,
                                "region": 1}
    assert counted.with_parent == {"mRNA": 1, "CDS": 1, "intron": 2, "exon": 1}


def test_genbank_states_no_parents():
    assert naive_counts(GENBANK).with_parent is None


def test_load_and_verify_report_the_statements_they_compared(tmp_path):
    config = _project(tmp_path, FUSARIUM)
    result = _load(tmp_path, config)
    assert result.exit_code == 0, result.output
    verified = _cli(tmp_path, "verify", "--assembly-id", "1")
    assert verified.exit_code == 0, verified.output
    line = next(row for row in verified.output.splitlines() if "parent edges" in row)
    assert "PASS" in line and "the file states" in line


def test_verify_fails_when_a_link_and_its_record_are_both_gone(tmp_path):
    """The gap D-092 closes. With the link and the record of how it was made both
    removed, the feature looks top-level, and the tier used to pass."""
    config = _project(tmp_path, FUSARIUM)
    assert _load(tmp_path, config).exit_code == 0
    connection = sqlite3.connect(tmp_path / "g.db")
    changed = connection.execute(
        "UPDATE feature SET parent_id = NULL, linkage_method = NULL "
        "WHERE source_id = 'exon-XM_011327722.1-2'"
    ).rowcount
    connection.commit()
    connection.close()
    assert changed == 1
    verified = _cli(tmp_path, "verify", "--assembly-id", "1")
    assert verified.exit_code != 0, verified.output
    assert "exon: the file gives" in verified.output


def test_a_parser_that_drops_a_parent_attribute_is_rolled_back(tmp_path, monkeypatch):
    real = tokenizer_module.parse_attributes

    def lose_one_parent(raw, path, line_number):
        attributes = real(raw, path, line_number)
        if attributes.get("ID") == ["exon-XM_011327722.1-2"]:
            attributes.pop("Parent", None)
        return attributes

    monkeypatch.setattr(tokenizer_module, "parse_attributes", lose_one_parent)
    config = _project(tmp_path, FUSARIUM)
    result = _load(tmp_path, config)
    assert result.exit_code == 2, result.output
    assert "did not store the parent links" in result.output
    assert "exon: the file gives" in result.output
    connection = sqlite3.connect(tmp_path / "g.db")
    stored = connection.execute("SELECT count(*) FROM feature").fetchone()[0]
    connection.close()
    assert stored == 0


def test_verify_says_when_it_could_check_links_only_inside_the_database(tmp_path):
    config = _project(tmp_path, FUSARIUM)
    assert _load(tmp_path, config).exit_code == 0
    (tmp_path / FUSARIUM.name).write_text(FUSARIUM.read_text() + "# edited afterwards\n")
    verified = _cli(tmp_path, "verify", "--assembly-id", "1")
    line = next(row for row in verified.output.splitlines() if "parent edges" in row)
    assert "within the database only" in line
