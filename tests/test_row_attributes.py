"""D-096: a feature written on several rows keeps each row's attributes and strand.

What these assert, stated before they were written:

* An alignment written on several rows (`cDNA_match`) keeps each row's own
  `Target` and `Gap`, and a partial CDS keeps `start_range` on the one row that
  had it. Before D-096 every row was exported with the first row's attributes;
  the benchmark's round trip found 3,249 such features in the human RefSeq file.
* A trans-spliced gene whose rows lie on both strands exports each row on the
  strand it was written with, while the feature as a whole still reads `?`
  (D-068) in `query`.
* A row that repeats the first row's attributes, as RefSeq's CDS rows usually
  do, stores nothing extra.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from cgload.cli.main import cli
from cgload.db.engine import build_engine
from cgload.db.schema import feature, feature_attribute
from cgload.export import export_gff3
from cgload.parsers import profiles as dialects
from cgload.parsers.normalise import normalise
from cgload.parsers.tokenizer import tokenize
from click.testing import CliRunner
from sqlalchemy import func, select

GFF = (
    "##gff-version 3\n"
    "##sequence-region chr1 1 10000\n"
    "##sequence-region chrM 1 5000\n"
    "chr1\tRefSeq\tregion\t1\t10000\t.\t+\t.\tID=chr1:1..10000;genome=chromosome\n"
    "chrM\tRefSeq\tregion\t1\t5000\t.\t+\t.\tID=chrM:1..5000;genome=mitochondrion\n"
    "chr1\tRefSeq\tgene\t100\t2000\t.\t+\t.\tID=gene-a;Name=a\n"
    "chr1\tRefSeq\tmRNA\t100\t2000\t.\t+\t.\tID=rna-a;Parent=gene-a;transcript_id=NM_1.1\n"
    "chr1\tRefSeq\texon\t100\t300\t.\t+\t.\tID=exon-a-1;Parent=rna-a\n"
    "chr1\tRefSeq\texon\t1000\t2000\t.\t+\t.\tID=exon-a-2;Parent=rna-a\n"
    "chr1\tRefSeq\tCDS\t100\t300\t.\t+\t0\tID=cds-a;Parent=rna-a;protein_id=NP_1.1;"
    "partial=true;start_range=.,100\n"
    "chr1\tRefSeq\tCDS\t1000\t1900\t.\t+\t2\tID=cds-a;Parent=rna-a;protein_id=NP_1.1;"
    "partial=true\n"
    "chr1\tRefSeq\tgene\t3000\t5000\t.\t-\t.\tID=gene-b;Name=b\n"
    "chr1\tRefSeq\tmRNA\t3000\t5000\t.\t-\t.\tID=rna-b;Parent=gene-b\n"
    "chr1\tRefSeq\texon\t3000\t3500\t.\t-\t.\tID=exon-b-1;Parent=rna-b\n"
    "chr1\tRefSeq\texon\t4500\t5000\t.\t-\t.\tID=exon-b-2;Parent=rna-b\n"
    "chr1\tRefSeq\tCDS\t3200\t3500\t.\t-\t0\tID=cds-b;Parent=rna-b;protein_id=NP_2.1\n"
    "chr1\tRefSeq\tCDS\t4500\t4800\t.\t-\t0\tID=cds-b;Parent=rna-b;protein_id=NP_2.1\n"
    "chr1\tRefSeq\tcDNA_match\t100\t300\t.\t+\t.\tID=aln-1;Target=NM_1.1 1 201 +;"
    "Gap=M201;pct_identity_gap=100\n"
    "chr1\tRefSeq\tcDNA_match\t1000\t2000\t.\t+\t.\tID=aln-1;Target=NM_1.1 202 1202 +;"
    "Gap=M1001;pct_identity_gap=100\n"
    "chrM\tRefSeq\tgene\t500\t800\t.\t+\t.\tID=gene-t;gene=nad1;exception=trans-splicing\n"
    "chrM\tRefSeq\tgene\t3000\t3300\t.\t-\t.\tID=gene-t;gene=nad1;exception=trans-splicing\n"
    "chrM\tRefSeq\tmRNA\t500\t3300\t.\t?\t.\tID=rna-t;Parent=gene-t;exception=trans-splicing\n"
    "chrM\tRefSeq\texon\t500\t800\t.\t+\t.\tID=exon-t-1;Parent=rna-t\n"
    "chrM\tRefSeq\texon\t3000\t3300\t.\t-\t.\tID=exon-t-2;Parent=rna-t\n"
)


def rows_of(text: str) -> Counter:
    """Every data row as (columns 1, 3-8, attributes as a sorted tuple), with
    cgload's own linkage stamp left out: what the file says, row by row."""
    found = Counter()
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        attributes = tuple(sorted(
            pair for pair in columns[8].split(";")
            if pair and not pair.startswith("cgload_linkage_method=")
        ))
        found[(columns[0], *columns[2:8], attributes)] += 1
    return found


def write_inputs(directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "a.gff3").write_text(GFF)
    with (directory / "a.fa").open("w") as handle:
        for name, length in (("chr1", 10000), ("chrM", 5000)):
            sequence = ("ACGT" * (length // 4 + 1))[:length]
            handle.write(f">{name}\n")
            for index in range(0, length, 60):
                handle.write(sequence[index : index + 60] + "\n")
    config = directory / "organism.toml"
    config.write_text(
        '[organism]\ngenus = "Testus"\nspecies = "rowwise"\nstrain = "S1"\n\n'
        '[assembly]\nname = "asm1"\nfasta = "a.fa"\n\n'
        '[[files.annotation]]\npath = "a.gff3"\n'
    )
    return config


def run(*args: str):
    return CliRunner().invoke(cli, list(args))


def loaded(tmp_path: Path) -> str:
    config = write_inputs(tmp_path / "in")
    url = f"sqlite:///{tmp_path / 'g.db'}"
    assert run("init", "--db-url", url).exit_code == 0
    result = run("load", "--db-url", url, "--config", str(config), "--data-dir",
                 str(tmp_path / "data"))
    assert result.exit_code == 0, result.output
    return url


def test_the_round_trip_through_the_database_returns_every_row_as_written(tmp_path):
    url = loaded(tmp_path)
    out = tmp_path / "out.gff3"
    exported = run("export", "--db-url", url, "--assembly-id", "1", "--format", "gff3",
                   "--output", str(out))
    assert exported.exit_code == 0, exported.output
    assert rows_of(out.read_text()) == rows_of(GFF)


def test_the_indexed_annotation_in_the_data_directory_is_row_for_row_too(tmp_path):
    import gzip

    loaded(tmp_path)
    annotation = next((tmp_path / "data").glob("**/annotation.gff3.gz"))
    with gzip.open(annotation, "rt") as handle:
        assert rows_of(handle.read()) == rows_of(GFF)


def test_the_parser_keeps_what_differs_and_only_that(tmp_path):
    path = tmp_path / "a.gff3"
    path.write_text(GFF)
    tokenized = tokenize(path)
    result = normalise(tokenized.features, dialects.REFSEQ, comments=tokenized.comments)
    by_id = {item.source_id: item for item in result.features}
    alignment = by_id["aln-1"]
    assert alignment.segments[1].attributes["Target"] == ["NM_1.1 202 1202 +"]
    # RefSeq's usual case: the second CDS row repeats the first exactly.
    assert all(segment.attributes is None for segment in by_id["cds-b"].segments)
    gene = by_id["gene-t"]
    assert gene.strand == "?"
    assert [segment.strand for segment in gene.segments] == ["+", "-"]
    text, _summary = export_gff3(result.features)
    assert rows_of(text) == rows_of(GFF)


def test_only_rows_that_differ_add_attribute_rows(tmp_path):
    url = loaded(tmp_path)
    with build_engine(url).connect() as connection:
        owners = dict(connection.execute(
            select(feature.c.source_id, func.count(func.distinct(feature_attribute.c.feature_id)))
            .join(feature_attribute, feature_attribute.c.feature_id == feature.c.feature_id)
            .group_by(feature.c.source_id)).all())
    assert owners["cds-b"] == 1  # identical rows: the first segment only
    assert owners["cds-a"] == 2  # start_range on one row only
    assert owners["aln-1"] == 2  # each row its own Target and Gap


def test_query_reports_a_trans_spliced_gene_as_both_strands(tmp_path):
    url = loaded(tmp_path)
    result = run("query", "--db-url", url, "--assembly-id", "1", "--id", "gene-t")
    assert result.exit_code == 0, result.output
    row = next(line for line in result.output.splitlines() if line.startswith("gene-t\t"))
    assert row.split("\t")[6] == "?"


# NCBI writes the rows of a minus-strand feature 5' to 3', so from the highest
# coordinate down. cgload orders segments by coordinate (finding 7), so the
# file's first row is not the first segment. Found by the benchmark's preflight,
# 22 October 2026: about half the alignments and partial CDS whose rows differ
# lost their lowest row's own attributes (D-104).
DESCENDING = (
    "##gff-version 3\n"
    "##sequence-region chr1 1 10000\n"
    "chr1\tRefSeq\tregion\t1\t10000\t.\t+\t.\tID=chr1:1..10000;genome=chromosome\n"
    "chr1\tRefSeq\tgene\t6000\t9000\t.\t-\t.\tID=gene-c;Name=c\n"
    "chr1\tRefSeq\tmRNA\t6000\t9000\t.\t-\t.\tID=rna-c;Parent=gene-c;transcript_id=NM_3.1\n"
    "chr1\tRefSeq\texon\t8000\t9000\t.\t-\t.\tID=exon-c-1;Parent=rna-c\n"
    "chr1\tRefSeq\texon\t6000\t6500\t.\t-\t.\tID=exon-c-2;Parent=rna-c\n"
    "chr1\tRefSeq\tCDS\t8000\t8800\t.\t-\t0\tID=cds-c;Parent=rna-c;protein_id=NP_3.1;"
    "partial=true;end_range=8800,.\n"
    "chr1\tRefSeq\tCDS\t6200\t6500\t.\t-\t2\tID=cds-c;Parent=rna-c;protein_id=NP_3.1;"
    "partial=true\n"
    "chr1\tRefSeq\tcDNA_match\t8000\t9000\t.\t-\t.\tID=aln-2;Target=NM_3.1 1 1001 +;"
    "Gap=M1001\n"
    "chr1\tRefSeq\tcDNA_match\t7000\t7400\t.\t-\t.\tID=aln-2;Target=NM_3.1 1002 1402 +;"
    "Gap=M401\n"
    "chr1\tRefSeq\tcDNA_match\t6000\t6500\t.\t-\t.\tID=aln-2;Target=NM_3.1 1403 1903 +;"
    "Gap=M501\n"
)


def loaded_from(tmp_path: Path, text: str) -> str:
    directory = tmp_path / "in"
    config = write_inputs(directory)
    (directory / "a.gff3").write_text(text)
    url = f"sqlite:///{tmp_path / 'g.db'}"
    assert run("init", "--db-url", url).exit_code == 0
    result = run("load", "--db-url", url, "--config", str(config), "--data-dir",
                 str(tmp_path / "data"))
    assert result.exit_code == 0, result.output
    return url


def test_rows_written_highest_first_keep_their_own_attributes(tmp_path):
    import gzip

    url = loaded_from(tmp_path, DESCENDING)
    out = tmp_path / "out.gff3"
    exported = run("export", "--db-url", url, "--assembly-id", "1", "--format", "gff3",
                   "--output", str(out))
    assert exported.exit_code == 0, exported.output
    assert rows_of(out.read_text()) == rows_of(DESCENDING)
    annotation = next((tmp_path / "data").glob("**/annotation.gff3.gz"))
    with gzip.open(annotation, "rt") as handle:
        assert rows_of(handle.read()) == rows_of(DESCENDING)


def test_the_first_segment_carries_the_feature_attributes(tmp_path):
    """The database reads a feature's attributes from its first segment, so the
    parser hands over the lowest row's as the feature's, and a later segment
    keeps its own only where it differs from those."""
    path = tmp_path / "a.gff3"
    path.write_text(DESCENDING)
    tokenized = tokenize(path)
    result = normalise(tokenized.features, dialects.REFSEQ, comments=tokenized.comments)
    by_id = {item.source_id: item for item in result.features}
    alignment = by_id["aln-2"]
    assert alignment.attributes["Target"] == ["NM_3.1 1403 1903 +"]
    assert [segment.attributes and segment.attributes["Target"]
            for segment in alignment.segments] == [
        None, ["NM_3.1 1002 1402 +"], ["NM_3.1 1 1001 +"]]
    cds = by_id["cds-c"]
    assert "end_range" not in cds.attributes
    assert cds.segments[0].attributes is None
    assert cds.segments[1].attributes["end_range"] == ["8800", "."]
