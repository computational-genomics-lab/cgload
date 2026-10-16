#!/usr/bin/env python3
"""Fault injection: which damage does cgload notice, and which does it not?

    python3 benchmark/faults.py ACCESSION --inputs DIR --work DIR --results DIR

Needs the SQLite load that run_cgload.py made for ACCESSION. Two kinds of fault:

* Damage after a load (F01-F18). Each fault is made on a fresh copy of the
  loaded database (and, where the fault is in the data directory, of that
  directory), and ``cgload verify`` is run on the copy.
* Damage inside the parser (P01-P08). The file is loaded again into a new
  database with cgload's tokenizer wrapped so that it hands on altered rows.
  The question is whether the load refuses and leaves the database empty, and,
  if it loads, whether ``verify`` then notices.

Each fault has an expected outcome written down here, before any run, from what
cgload claims: it checks completeness (every feature and every stated parent
link stored, the stored assembly unchanged), not fidelity (D-075). A fault that
keeps every count is expected to pass unnoticed. One row per fault goes to
RESULTS/faults.tsv, with what was expected and what happened.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import random
import shutil
import sqlite3
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import benchcommon as bc  # noqa: E402

#: id, name, layer, expected outcome, what it does.
DB_FAULTS = (
    ("F01", "delete_one_feature", "database", "detected", "delete one exon and its attributes"),
    ("F02", "delete_one_percent", "database", "detected", "delete 1% of exons"),
    ("F03", "delete_gene_model", "database", "detected", "delete a gene and everything below it"),
    ("F04", "unlink_one", "database", "detected", "set one exon's parent to none"),
    ("F05", "unlink_and_unrecord", "database", "detected",
     "set one exon's parent and its linkage record to none"),
    ("F06", "repoint_parent", "database", "not detected",
     "move one exon to a different mRNA"),
    ("F07", "flip_strand", "database", "not detected", "flip the strand of one gene"),
    ("F08", "shift_start", "database", "not detected", "move one exon's start by one base"),
    ("F09", "change_type", "database", "detected", "turn one exon into a CDS"),
    ("F10", "rename_feature", "database", "not detected", "give one exon a new ID"),
    ("F11", "duplicate_feature", "database", "detected", "add a copy of one exon under a new ID"),
    ("F12", "delete_attribute", "database", "not detected", "delete one attribute value"),
    ("F13", "truncate_fasta", "data", "detected", "cut the last line off the stored FASTA"),
    ("F14", "delete_fai", "data", "detected", "delete the FASTA index"),
    ("F15", "change_one_base", "data", "detected", "change one base of the stored FASTA"),
    ("F16", "input_edited", "input", "detected",
     "add a comment line to the input GFF3 after the load"),
    ("F17", "data_dir_moved", "data", "detected",
     "move the data directory; verify without --data-dir"),
    ("F18", "unlink_all", "database", "detected", "set every parent to none"),
    ("F19", "delete_one_row", "database", "detected",
     "delete one row of a CDS written on several rows"),
)

PARSER_FAULTS = (
    ("P01", "drop_rows", "detected", "the tokenizer drops every 1000th row"),
    ("P02", "drop_one_segment", "detected",
     "the tokenizer drops one row of a feature written on several rows"),
    ("P03", "drop_parent", "detected", "the tokenizer drops the Parent of one exon"),
    ("P04", "retype_row", "detected", "the tokenizer reads one exon row as a CDS"),
    ("P05", "flip_strand", "not detected", "the tokenizer flips the strand of one exon row"),
    ("P06", "invent_feature", "detected", "the tokenizer adds a copy of one exon under a new ID"),
    ("P07", "shift_start", "not detected", "the tokenizer moves one exon's start by one base"),
    ("P08", "swap_ids", "not detected", "the tokenizer swaps the IDs of two exons"),
)

COLUMNS = ("accession", "fault_id", "fault", "layer", "description", "target", "expected",
           "observed", "as_expected", "exit_code", "tier_source_count", "tier_parent_edges",
           "tier_stored_assembly", "load_refused", "database_empty", "message")


def _rng(accession: str, fault: str) -> random.Random:
    seed = int(hashlib.sha256(f"{accession}:{fault}".encode()).hexdigest()[:12], 16)
    return random.Random(seed)


def leaf_exons(connection: sqlite3.Connection) -> list[str]:
    """Exons written on one row, with a parent, and with no children."""
    return [row[0] for row in connection.execute("""
        SELECT f.source_id FROM feature f
        WHERE f.source_type = 'exon' AND f.source_id IS NOT NULL AND f.parent_id IS NOT NULL
          AND NOT EXISTS (SELECT 1 FROM feature c WHERE c.parent_id = f.feature_id)
        GROUP BY f.source_id HAVING count(*) = 1 ORDER BY f.source_id""")]


def apply_db_fault(fault: str, db: Path, rng: random.Random, inputs: Path, copy: Path) -> str:
    connection = sqlite3.connect(db)
    exons = leaf_exons(connection)
    target = rng.choice(exons) if exons else None

    def ids_of(source_id):
        return [row[0] for row in connection.execute(
            "SELECT feature_id FROM feature WHERE source_id = ?", (source_id,))]

    def delete(feature_ids):
        for chunk in range(0, len(feature_ids), 500):
            part = feature_ids[chunk:chunk + 500]
            marks = ",".join("?" * len(part))
            connection.execute(
                f"DELETE FROM feature_attribute WHERE feature_id IN ({marks})", part)  # noqa: S608
            connection.execute(f"DELETE FROM feature WHERE feature_id IN ({marks})", part)  # noqa: S608

    if fault == "delete_one_feature":
        delete(ids_of(target))
    elif fault == "delete_one_percent":
        chosen = rng.sample(exons, max(1, len(exons) // 100))
        target = f"{len(chosen)} exons"
        delete([i for name in chosen for i in ids_of(name)])
    elif fault == "delete_gene_model":
        genes = [row[0] for row in connection.execute("""
            SELECT g.source_id FROM feature g JOIN feature m ON m.parent_id = g.feature_id
            WHERE g.source_type = 'gene' AND m.source_type = 'mRNA'
            GROUP BY g.source_id ORDER BY g.source_id""")]
        target = rng.choice(genes)
        rows = [row[0] for row in connection.execute("""
            WITH RECURSIVE below(id) AS (
                SELECT feature_id FROM feature WHERE source_id = ?
                UNION SELECT f.feature_id FROM feature f JOIN below b ON f.parent_id = b.id)
            SELECT id FROM below""", (target,))]
        delete(rows)
    elif fault == "unlink_one":
        connection.execute("UPDATE feature SET parent_id = NULL WHERE source_id = ?", (target,))
    elif fault == "unlink_and_unrecord":
        connection.execute("UPDATE feature SET parent_id = NULL, linkage_method = NULL "
                           "WHERE source_id = ?", (target,))
    elif fault == "repoint_parent":
        current = connection.execute(
            "SELECT parent_id FROM feature WHERE source_id = ?", (target,)).fetchone()[0]
        others = [row[0] for row in connection.execute(
            "SELECT feature_id FROM feature WHERE source_type = 'mRNA' AND segment_index = 0 "
            "AND feature_id != ? ORDER BY feature_id LIMIT 1000", (current,))]
        connection.execute("UPDATE feature SET parent_id = ? WHERE source_id = ?",
                           (rng.choice(others), target))
    elif fault == "flip_strand":
        genes = [row[0] for row in connection.execute(
            "SELECT source_id FROM feature WHERE source_type = 'gene' AND strand IN ('+', '-') "
            "ORDER BY source_id")]
        target = rng.choice(genes)
        connection.execute("UPDATE feature SET strand = CASE strand WHEN '+' THEN '-' "
                           "ELSE '+' END WHERE source_id = ?", (target,))
    elif fault == "shift_start":
        connection.execute('UPDATE feature SET start = start + 1 WHERE source_id = ? '
                           'AND start < "end"', (target,))
    elif fault == "change_type":
        connection.execute("UPDATE feature SET source_type = 'CDS', feature_type = 'CDS' "
                           "WHERE source_id = ?", (target,))
    elif fault == "rename_feature":
        connection.execute("UPDATE feature SET source_id = source_id || '-renamed' "
                           "WHERE source_id = ?", (target,))
    elif fault == "duplicate_feature":
        columns = [row[1] for row in connection.execute("PRAGMA table_info(feature)")
                   if row[1] != "feature_id"]
        names = ", ".join(f'"{c}"' for c in columns)
        chosen = ", ".join("source_id || '-copy'" if c == "source_id" else f'"{c}"'
                           for c in columns)
        connection.execute(f"INSERT INTO feature ({names}) SELECT {chosen} FROM feature "  # noqa: S608
                           "WHERE source_id = ?", (target,))
    elif fault == "delete_attribute":
        attribute = connection.execute(
            "SELECT a.feature_attribute_id FROM feature_attribute a JOIN feature f "
            "ON f.feature_id = a.feature_id WHERE f.source_id = ? LIMIT 1", (target,)).fetchone()
        connection.execute("DELETE FROM feature_attribute WHERE feature_attribute_id = ?",
                           (attribute[0],))
    elif fault == "unlink_all":
        target = "every feature"
        connection.execute("UPDATE feature SET parent_id = NULL")
    elif fault == "delete_one_row":
        multi = [row[0] for row in connection.execute(
            "SELECT source_id FROM feature WHERE source_type = 'CDS' AND source_id IS NOT NULL "
            "GROUP BY source_id HAVING count(*) > 1 ORDER BY source_id")]
        target = rng.choice(multi)
        (last,) = connection.execute(
            "SELECT feature_id FROM feature WHERE source_id = ? "
            "ORDER BY segment_index DESC LIMIT 1", (target,)).fetchone()
        delete([last])
    elif fault == "truncate_fasta":
        fasta = next(copy.glob("**/assembly.fa"))
        lines = fasta.read_bytes().splitlines(keepends=True)
        fasta.write_bytes(b"".join(lines[:-1]))
        target = fasta.name
    elif fault == "delete_fai":
        index = next(copy.glob("**/assembly.fa.fai"))
        index.unlink()
        target = index.name
    elif fault == "change_one_base":
        fasta = next(copy.glob("**/assembly.fa"))
        with open(fasta, "r+b") as handle:
            handle.readline()
            position = handle.tell() + 10
            handle.seek(position)
            base = handle.read(1)
            handle.seek(position)
            handle.write(b"C" if base != b"C" else b"G")
        target = f"{fasta.name} byte {position}"
    elif fault == "input_edited":
        source = inputs / "genomic.gff"
        edited = copy / "genomic.gff"
        shutil.copyfile(source, edited)
        with edited.open("a") as handle:
            handle.write("# edited after the load\n")
        connection.execute("UPDATE source_file SET path = ? WHERE path = ?",
                           (str(edited), str(source)))
        target = "genomic.gff"
    elif fault == "data_dir_moved":
        connection.execute("UPDATE assembly SET data_root = ?", (str(copy / "moved-away"),))
        target = "assembly.data_root"
    connection.commit()
    connection.close()
    return target or ""


def run_verify(name: str, db: Path, data_dir: Path | None, steps: Path) -> tuple[dict, dict]:
    command = ["cgload", "verify", "--assembly-id", "1"]
    if data_dir is not None:
        command += ["--data-dir", str(data_dir)]
    measured = bc.measured(name, command, steps, env={"CGLOAD_DB_URL": f"sqlite:///{db}"})
    return measured, bc.tiers(bc.text_of(steps, name))


def assembly_dir(data: Path) -> Path:
    return next(data.glob("**/manifest.json")).parent


def db_faults(args, source_db: Path, source_data: Path, inputs: Path) -> None:
    table = args.results / "faults.tsv"
    for fault_id, name, layer, expected, description in DB_FAULTS:
        if bc.done(table, accession=args.accession, fault_id=fault_id):
            continue
        work = args.work / "faults" / args.accession / fault_id
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
        db = work / "cgload.db"
        shutil.copyfile(source_db, db)
        data_dir = source_data
        if layer == "data":
            data_dir = work / "data"
            shutil.copytree(source_data, data_dir)
        target = apply_db_fault(name, db, _rng(args.accession, fault_id), inputs, work)
        steps = work / "steps"
        verify_dir = None if name == "data_dir_moved" else data_dir
        measured, found = run_verify("verify", db, verify_dir, steps)
        text = bc.text_of(steps, "verify")
        if name == "input_edited":
            observed = ("detected" if "has changed since it was loaded" in text
                        else "not detected")
        else:
            observed = "detected" if measured["exit_code"] != 0 else "not detected"
        row = dict.fromkeys(COLUMNS)
        row.update(accession=args.accession, fault_id=fault_id, fault=name, layer=layer,
                   description=description, target=target, expected=expected,
                   observed=observed, as_expected=observed == expected,
                   exit_code=measured["exit_code"],
                   tier_source_count=found.get("source count"),
                   tier_parent_edges=found.get("parent edges"),
                   tier_stored_assembly=found.get("stored assembly"),
                   message=bc.first_lines("\n".join(
                       line for line in text.splitlines() if "FAIL" in line or "changed" in line
                   ), 2))
        bc.append_row(table, row)
        print(f"{args.accession} {fault_id} {name}: expected {expected}, {observed}")
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)


# --------------------------------------------------------------------------- parser faults

def apply_parser_fault(fault: str, features: list) -> str:
    """Alter the tokenizer's rows in place; return what was changed."""
    exon_rows = [i for i, f in enumerate(features) if f.type == "exon" and f.id
                 and f.parents]
    counts: dict[str, int] = {}
    for f in features:
        if f.id:
            counts[f.id] = counts.get(f.id, 0) + 1
    single_exons = [i for i in exon_rows if counts.get(features[i].id) == 1]
    if fault == "drop_rows":
        doomed = set(range(999, len(features), 1000)) or {len(features) - 1}
        features[:] = [f for i, f in enumerate(features) if i not in doomed]
        return f"{len(doomed)} rows"
    if fault == "drop_one_segment":
        multi = next(f.id for f in features if f.id and counts[f.id] > 1)
        last = max(i for i, f in enumerate(features) if f.id == multi)
        del features[last]
        return multi
    index = single_exons[len(single_exons) // 2]
    row = features[index]
    if fault == "drop_parent":
        attributes = {k: v for k, v in row.attributes.items() if k.lower() != "parent"}
        features[index] = replace(row, attributes=attributes)
    elif fault == "retype_row":
        features[index] = replace(row, type="CDS", phase=0)
    elif fault == "flip_strand":
        features[index] = replace(row, strand="-" if row.strand == "+" else "+")
    elif fault == "invent_feature":
        attributes = dict(row.attributes)
        attributes["ID"] = [row.id + "-invented"]
        features.insert(index + 1, replace(row, attributes=attributes))
    elif fault == "shift_start":
        features[index] = replace(row, start=row.start + 1 if row.start < row.end else row.start)
    elif fault == "swap_ids":
        other_index = single_exons[len(single_exons) // 2 + 1]
        other = features[other_index]
        first = dict(row.attributes)
        second = dict(other.attributes)
        first["ID"], second["ID"] = second["ID"], first["ID"]
        features[index] = replace(row, attributes=first)
        features[other_index] = replace(other, attributes=second)
        return f"{row.id} <-> {other.id}"
    return row.id


def parser_run(argv: list[str]) -> int:
    """Run inside a subprocess: patch cgload's tokenizer, then run `cgload load`."""
    fault, target_file, *load_args = argv
    import cgload.loader as loader_module
    from cgload.cli.main import cli

    real = loader_module.tokenize

    def faulty(path, *, tolerate_inverted=False):
        result = real(path, tolerate_inverted=tolerate_inverted)
        Path(target_file).write_text(apply_parser_fault(fault, result.features))
        return result

    loader_module.tokenize = faulty
    return cli.main(args=["load", *load_args], prog_name="cgload", standalone_mode=False) or 0


def parser_faults(args, repair: bool) -> None:
    table = args.results / "faults.tsv"
    config_source = args.work / "cgload" / "sqlite" / args.accession / "organism.toml"
    for fault_id, name, expected, description in PARSER_FAULTS:
        if bc.done(table, accession=args.accession, fault_id=fault_id):
            continue
        work = args.work / "faults" / args.accession / fault_id
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
        db = work / "cgload.db"
        env = {"CGLOAD_DB_URL": f"sqlite:///{db}"}
        steps = work / "steps"
        bc.measured("init", ["cgload", "init"], steps, env=env)
        target_file = work / "target.txt"
        load_args = ["--config", str(config_source), "--data-dir", str(work / "data")]
        if repair:
            load_args.append("--repair-envelopes")
        load = bc.measured("load", [sys.executable, str(Path(__file__).resolve()),
                                    "_parser-run", name, str(target_file), *load_args],
                           steps, env=env)
        with sqlite3.connect(db) as connection:
            empty = connection.execute("SELECT count(*) FROM assembly").fetchone()[0] == 0
        found = {}
        text = bc.text_of(steps, "load")
        if load["exit_code"] == 0:
            measured, found = run_verify("verify", db, None, steps)
            observed = "detected" if measured["exit_code"] != 0 else "not detected"
            text = bc.text_of(steps, "verify")
            exit_code = measured["exit_code"]
        else:
            observed = "detected"
            exit_code = load["exit_code"]
        row = dict.fromkeys(COLUMNS)
        row.update(accession=args.accession, fault_id=fault_id, fault=name, layer="parser",
                   description=description,
                   target=target_file.read_text() if target_file.exists() else "",
                   expected=expected, observed=observed, as_expected=observed == expected,
                   exit_code=exit_code, tier_source_count=found.get("source count"),
                   tier_parent_edges=found.get("parent edges"),
                   tier_stored_assembly=found.get("stored assembly"),
                   load_refused=load["exit_code"] != 0, database_empty=empty,
                   message=bc.first_lines(text, 2))
        bc.append_row(table, row)
        print(f"{args.accession} {fault_id} {name}: expected {expected}, {observed}")
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["_parser-run"]:
        return parser_run(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--keep", action="store_true", help="keep each fault's copy")
    args = parser.parse_args(argv)
    loaded = args.work / "cgload" / "sqlite" / args.accession
    source_db = loaded / "cgload.db"
    if not source_db.exists():
        raise SystemExit(f"{source_db} is missing: run run_cgload.py --backend sqlite first")
    repair = bc.done(args.results / "cgload_scale.tsv", accession=args.accession,
                     backend="sqlite", outcome="loaded_with_repair")
    os.environ.pop("CGLOAD_DB_URL", None)
    db_faults(args, source_db, loaded / "data", args.inputs / args.accession)
    parser_faults(args, repair)
    return 0


if __name__ == "__main__":
    sys.exit(main())
