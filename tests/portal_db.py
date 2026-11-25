"""A small database for the portal's tests, built from the fixtures.

Three assemblies:

1. *Drosophila melanogaster*, the Dscam1 region of chromosome 2R, from GenBank,
   with the record's own sequence: real bases, so translations can be compared
   with the proteins NCBI wrote in the record.
2. *Fusarium graminearum*, a slice of chromosome 4 in RefSeq GFF3, with
   eggNOG-mapper and InterProScan results; the sequence is made up.
3. A file written to attack the portal: attribute values holding HTML and
   script, an ID with characters a URL must escape.

The search index is built for the first two only, so the third shows what the
portal does without one.
"""

from __future__ import annotations

import random
import re
from pathlib import Path

from cgload.cli.main import cli
from click.testing import CliRunner

FIXTURES = Path(__file__).parent / "fixtures"
DSCAM = FIXTURES / "genbank_dscam1_CG12164.gbff"
FUSARIUM = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"
EGGNOG = FIXTURES / "eggnog_fgraminearum.emapper.annotations"
INTERPRO = FIXTURES / "interproscan_fgraminearum.tsv"

HOSTILE_NAME = '<script>alert("name")</script>'
HOSTILE_PRODUCT = '"><img src=x onerror=alert(1)>'
HOSTILE_ID = "gene-<b>&evil/../1"
_GENE_ID = HOSTILE_ID.replace("&", "%26")
_ROWS = [
    ("gene", 100, 900, ".", f"ID={_GENE_ID};Name={HOSTILE_NAME};Dbxref=GeneID:javascript:alert(1)"),
    ("mRNA", 100, 900, ".",
     f"ID=rna-1;Parent={_GENE_ID};product={HOSTILE_PRODUCT.replace('=', '%3D')}"),
    ("exon", 100, 300, ".", "ID=exon-1;Parent=rna-1"),
    ("exon", 500, 900, ".", "ID=exon-2;Parent=rna-1"),
    ("CDS", 150, 300, "0", "ID=cds-1;Parent=rna-1"),
    ("CDS", 500, 800, "0", "ID=cds-1;Parent=rna-1"),
]
HOSTILE_GFF = "##gff-version 3\n##sequence-region chrX 1 5000\n" + "".join(
    f"chrX\tfunannotate\t{kind}\t{start}\t{end}\t.\t+\t{phase}\t{attributes}\n"
    for kind, start, end, phase, attributes in _ROWS)


def _run(*args: str) -> None:
    result = CliRunner().invoke(cli, list(args), catch_exceptions=False)
    assert result.exit_code == 0, result.output


def _origin_fasta(gbff: Path, target: Path) -> str:
    text = gbff.read_text()
    seqid = re.search(r"^VERSION\s+(\S+)", text, re.M).group(1)
    bases = "".join(re.sub(r"[^a-z]", "", line)
                    for line in text.split("\nORIGIN", 1)[1].splitlines())
    _write_fasta(target, {seqid: bases.upper()})
    return seqid


def _write_fasta(target: Path, sequences: dict[str, str]) -> None:
    with target.open("w") as handle:
        for name, bases in sequences.items():
            handle.write(f">{name}\n")
            for index in range(0, len(bases), 60):
                handle.write(bases[index : index + 60] + "\n")


def _made_up(lengths: dict[str, int], seed: int = 7) -> dict[str, str]:
    block = "".join(random.Random(seed).choice("ACGT") for _ in range(6000))
    return {name: (block * (length // len(block) + 1))[:length]
            for name, length in lengths.items()}


def _declared_lengths(gff: Path) -> dict[str, int]:
    lengths = {}
    for line in gff.read_text().splitlines():
        if line.startswith("##sequence-region"):
            parts = line.split()
            lengths[parts[1]] = int(parts[3])
    return lengths


def build(directory: Path, url: str | None = None) -> tuple[Path, Path]:
    """Load the three assemblies into ``url`` (a new SQLite file by default);
    return the SQLite file (or None for a server) and the data directory."""
    directory.mkdir(parents=True, exist_ok=True)
    db = directory / "portal.db" if url is None else None
    data = directory / "data"
    url = url or f"sqlite:///{db}"
    _run("init", "--db-url", url)

    fly = directory / "fly"
    fly.mkdir()
    (fly / DSCAM.name).write_text(DSCAM.read_text())
    _origin_fasta(DSCAM, fly / "region.fa")
    (fly / "organism.toml").write_text(
        '[organism]\ngenus = "Drosophila"\nspecies = "melanogaster"\n'
        'strain = "reference"\nncbi_taxon_id = 7227\n\n'
        '[assembly]\nname = "Dscam1 region"\nfasta = "region.fa"\n\n'
        f'[[files.annotation]]\npath = "{DSCAM.name}"\n')
    _run("load", "--config", str(fly / "organism.toml"), "--db-url", url,
         "--data-dir", str(data))

    fusarium = directory / "fusarium"
    fusarium.mkdir()
    for path in (FUSARIUM, EGGNOG, INTERPRO):
        (fusarium / path.name).write_text(path.read_text())
    _write_fasta(fusarium / "assembly.fa", _made_up(_declared_lengths(FUSARIUM)))
    (fusarium / "organism.toml").write_text(
        '[organism]\ngenus = "Fusarium"\nspecies = "graminearum"\nstrain = "PH-1"\n\n'
        '[assembly]\nname = "chromosome 4 slice"\nfasta = "assembly.fa"\n\n'
        f'[[files.annotation]]\npath = "{FUSARIUM.name}"\n\n'
        f'[[files.functional]]\npath = "{EGGNOG.name}"\ntool = "eggnog"\n\n'
        f'[[files.functional]]\npath = "{INTERPRO.name}"\ntool = "interproscan"\n')
    _run("load", "--config", str(fusarium / "organism.toml"), "--db-url", url,
         "--data-dir", str(data))

    hostile = directory / "hostile"
    hostile.mkdir()
    (hostile / "hostile.gff3").write_text(HOSTILE_GFF)
    _write_fasta(hostile / "assembly.fa", _made_up({"chrX": 5000}, seed=3))
    (hostile / "organism.toml").write_text(
        '[organism]\ngenus = "Hostilia"\nspecies = "testis"\nstrain = "S"\n\n'
        '[assembly]\nname = "attack"\nfasta = "assembly.fa"\n\n'
        '[[files.annotation]]\npath = "hostile.gff3"\n')
    _run("load", "--config", str(hostile / "organism.toml"), "--db-url", url,
         "--data-dir", str(data))

    _run("index", "--db-url", url, "--assembly-id", "1", "--assembly-id", "2")
    return db, data
