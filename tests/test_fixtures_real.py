"""Assertions about the real fixtures.

These are not parser tests -- the parser does not exist yet. They pin the
structural facts the milestone-3 profiles are being built against, so that if a
fixture is ever replaced with one that does not exhibit the trap, the test that
depends on the trap fails loudly rather than passing vacuously.
"""

from __future__ import annotations

import collections
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"

AUGUSTUS = FIXTURES / "augustus_aphaeospermum_pheo_arth1.gff3"
FUNANNOTATE = FIXTURES / "funannotate_aphaeospermum_pheo_arth1.gff3"
REFSEQ = FIXTURES / "refseq_fgraminearum_NC_026477.gff3"


def rows(path: Path):
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.split("\t")
        attributes = dict(
            field.split("=", 1)
            for field in columns[8].split(";")
            if field.strip() and "=" in field
        )
        yield columns, attributes


@pytest.mark.parametrize("path", [AUGUSTUS, FUNANNOTATE, REFSEQ])
def test_fixture_is_present_and_non_trivial(path):
    assert path.exists(), f"{path.name} missing; see tests/fixtures/PROVENANCE.md"
    assert sum(1 for _ in rows(path)) > 100


def test_augustus_emits_no_exons():
    """The dialect declares no exon type; zero exons is expected, not a drop."""
    types = {columns[2] for columns, _ in rows(AUGUSTUS)}
    assert "exon" not in types
    assert {"gene", "transcript", "CDS", "intron", "start_codon", "stop_codon"} <= types


def test_augustus_anonymous_rows_exist_and_are_numerous():
    """D-024. If this ever returns zero, the anonymous-feature path is untested."""
    anonymous = collections.Counter(
        columns[2] for columns, attributes in rows(AUGUSTUS) if "ID" not in attributes
    )
    assert set(anonymous) == {"intron", "start_codon", "stop_codon"}
    assert sum(anonymous.values()) > 100


def test_augustus_codons_lie_inside_a_cds_of_the_same_transcript():
    """D-022's evidence. Mapping codon rows to CDS would double-count bases."""
    cds = collections.defaultdict(list)
    for columns, attributes in rows(AUGUSTUS):
        if columns[2] == "CDS":
            cds[attributes["Parent"]].append((int(columns[3]), int(columns[4])))

    checked = 0
    for columns, attributes in rows(AUGUSTUS):
        if columns[2] not in ("start_codon", "stop_codon"):
            continue
        start, end = int(columns[3]), int(columns[4])
        spans = cds[attributes["Parent"]]
        assert any(s <= start and end <= e for s, e in spans), (
            f"{columns[2]} {start}-{end} is not inside a CDS of {attributes['Parent']}"
        )
        checked += 1
    assert checked > 50


def test_funannotate_transcript_suffix_would_collide_with_its_gene():
    """D-023's evidence: stripping -T1 turns an mRNA ID into its parent's ID."""
    collisions = 0
    for columns, attributes in rows(FUNANNOTATE):
        if columns[2] in ("mRNA", "tRNA") and attributes["ID"].endswith("-T1"):
            assert attributes["ID"].removesuffix("-T1") == attributes["Parent"]
            collisions += 1
    assert collisions > 20


def test_refseq_has_a_url_encoded_attribute_value():
    """Percent-decoding must happen after splitting on ';'."""
    region = next(columns for columns, _ in rows(REFSEQ) if columns[2] == "region")
    assert "%3B" in region[8]


def test_dbxref_is_a_comma_separated_list():
    multi = [
        attributes["Dbxref"]
        for _, attributes in rows(REFSEQ)
        if "," in attributes.get("Dbxref", "")
    ]
    assert multi, "no multi-valued Dbxref in the fixture"


def test_funannotate_spells_it_dbxref_with_a_capital_b():
    keys = {key for _, attributes in rows(FUNANNOTATE) for key in attributes}
    assert "DBxref" in keys
    assert "Dbxref" not in keys


@pytest.mark.parametrize("path", [AUGUSTUS, FUNANNOTATE, REFSEQ])
def test_only_cds_rows_share_an_id(path):
    """CDS segments repeat one ID; exons do not. An implementation assuming
    either rule held universally would be wrong half the time."""
    by_id = collections.defaultdict(list)
    for columns, attributes in rows(path):
        if "ID" in attributes:
            by_id[attributes["ID"]].append(columns[2])
    repeated = {i: set(t) for i, t in by_id.items() if len(t) > 1}
    assert repeated, f"{path.name} has no discontinuous feature"
    assert all(types == {"CDS"} for types in repeated.values())


@pytest.mark.parametrize("path", [AUGUSTUS, FUNANNOTATE, REFSEQ])
def test_discontinuous_groups_satisfy_the_d014_exemption(path):
    """Same type, strand and parent, non-overlapping. If a real file ever
    violates this, the exemption is wrong, not the file."""
    groups = collections.defaultdict(list)
    for columns, attributes in rows(path):
        if "ID" in attributes:
            groups[attributes["ID"]].append((columns, attributes))

    for identifier, group in groups.items():
        if len(group) < 2:
            continue
        assert len({c[0] for c, _ in group}) == 1, identifier
        assert len({c[2] for c, _ in group}) == 1, identifier
        assert len({c[6] for c, _ in group}) == 1, identifier
        assert len({a.get("Parent") for _, a in group}) == 1, identifier
        spans = sorted((int(c[3]), int(c[4])) for c, _ in group)
        for (_, end), (next_start, _) in zip(spans, spans[1:], strict=False):
            assert end < next_start, f"{identifier} has overlapping segments"


def test_funannotate_emits_minus_strand_segments_in_descending_file_order():
    """Evidence that segment_index must be coordinate order, not file order."""
    groups = collections.defaultdict(list)
    for columns, attributes in rows(FUNANNOTATE):
        if "ID" in attributes:
            groups[attributes["ID"]].append(int(columns[3]))
    descending = [
        i for i, starts in groups.items()
        if len(starts) > 1 and starts == sorted(starts, reverse=True)
    ]
    assert descending, "fixture no longer exercises reverse-order segments"


@pytest.mark.parametrize("path", [AUGUSTUS, FUNANNOTATE, REFSEQ])
def test_every_parent_reference_resolves(path):
    """The fixtures are graph-closed slices. A dangling Parent here would mean
    the slice was cut badly, not that the loader should tolerate one."""
    known = {a["ID"] for _, a in rows(path) if "ID" in a}
    for _columns, attributes in rows(path):
        for parent in attributes.get("Parent", "").split(","):
            if parent:
                assert parent in known, f"{path.name}: dangling Parent {parent!r}"


def test_refseq_child_parent_type_pairs_are_all_permitted():
    """Checks the real hierarchy against ALLOWED_PARENTS, including the
    surprising ones: mRNA and tRNA hanging off a pseudogene."""
    from cgload.db.vocabulary import ALLOWED_PARENTS, SYNONYMS

    def canonical(term: str) -> str:
        return SYNONYMS.get(term, term)

    types = {a["ID"]: columns[2] for columns, a in rows(REFSEQ) if "ID" in a}
    seen = set()
    for columns, attributes in rows(REFSEQ):
        parent = attributes.get("Parent")
        if not parent:
            continue
        child_type = canonical(columns[2])
        parent_type = canonical(types[parent])
        seen.add((child_type, parent_type))
        if child_type in ALLOWED_PARENTS:
            assert parent_type in ALLOWED_PARENTS[child_type], (
                f"{child_type} <- {parent_type} occurs in RefSeq but is not permitted"
            )
    assert ("mRNA", "pseudogene") in seen
    assert ("tRNA", "pseudogene") in seen


# -- AUGUSTUS comment blocks (D-025) ---------------------------------------


def comment_lines(path: Path):
    return [line for line in path.read_text().splitlines() if line.startswith("#")]


def data_lines(path: Path):
    return [
        line
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ]


def test_augustus_is_mostly_comments():
    """788 comment lines to 399 data lines. Any check that counted lines rather
    than features would be wrong by a factor of three, so the naive counting
    pass of D-008 must apply the same skip rule as the tokenizer."""
    assert len(comment_lines(AUGUSTUS)) > 1.5 * len(data_lines(AUGUSTUS))


def test_no_comment_line_contains_a_tab():
    """The skip rule is 'first character is #'. If a comment carried tabs, a
    field-count heuristic could mistake it for data; it does not, so the simple
    rule is safe."""
    assert all("\t" not in line for line in comment_lines(AUGUSTUS))


def extract_proteins(path: Path) -> dict[str, str]:
    """Reconstruct AUGUSTUS's declared translations from its comment blocks.

    This is the only place AUGUSTUS emits protein sequence. The block opens with
    ``# protein sequence = [`` and continues over ``#``-prefixed lines until a
    closing ``]``.
    """
    proteins: dict[str, str] = {}
    gene: str | None = None
    buffer: list[str] | None = None
    for line in path.read_text().splitlines():
        if line.startswith("# start gene "):
            gene = line.split()[-1]
            continue
        if line.startswith("# protein sequence = ["):
            buffer = [line.split("[", 1)[1]]
        elif buffer is not None:
            buffer.append(line[1:].strip())
        else:
            continue
        if buffer[-1].endswith("]"):
            buffer[-1] = buffer[-1][:-1]
            proteins[str(gene)] = "".join(buffer)
            buffer = None
    return proteins


def test_every_augustus_gene_declares_a_protein():
    proteins = extract_proteins(AUGUSTUS)
    genes = {
        attributes["ID"] for columns, attributes in rows(AUGUSTUS) if columns[2] == "gene"
    }
    assert set(proteins) == genes
    assert all(protein.isalpha() for protein in proteins.values())


def test_declared_protein_length_matches_the_cds_spans():
    """The free oracle. AUGUSTUS's declared protein length must equal
    sum(CDS spans) / 3 - 1 (the stop codon is inside the CDS, D-022). This
    checks our coordinate arithmetic against the tool's own output, using
    nothing but the file -- an independent cross-check of exactly the kind
    D-008 asks for, available at no cost on every AUGUSTUS load."""
    proteins = extract_proteins(AUGUSTUS)

    cds_bases: dict[str, int] = collections.defaultdict(int)
    for columns, attributes in rows(AUGUSTUS):
        if columns[2] == "CDS":
            cds_bases[attributes["Parent"]] += int(columns[4]) - int(columns[3]) + 1

    checked = 0
    for gene, protein in proteins.items():
        total = cds_bases[f"{gene}.t1"]
        assert total % 3 == 0, f"{gene}: CDS total {total} is not a whole number of codons"
        assert total // 3 - 1 == len(protein), (
            f"{gene}: CDS {total} bp implies {total // 3 - 1} aa, declared {len(protein)}"
        )
        checked += 1
    assert checked == 49
