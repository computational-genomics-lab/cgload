#!/usr/bin/env python3
"""Do cgload's structural checks and the survey's checker agree on the same file?

    python3 benchmark/agreement.py ACCESSION --inputs DIR --results DIR

The survey's checker (survey/cgsurvey.py) flagged rows; cgload judges features.
For one genome this script lists, side by side:

* every row the survey flagged (survey/results/hits.tsv), with the survey's
  class (survey/classify_hits.py's rule, imported, not copied), and what cgload
  does with the feature on that row: accept it, or refuse the file because of
  it; and, under --repair-envelopes, repair it, store it as written with a
  warning, or still refuse;
* every feature cgload objects to that the survey did not flag.

cgload's view comes from its own code, run the way ``cgload load
--repair-envelopes`` runs it, without a database: tokenise with inverted rows
tolerated, detect the dialect, normalise, find envelope violations, repair,
set aside disagreements that exist only across a circular origin (D-091), and
check each row against its sequence's length (D-081). Lengths are the
``##sequence-region`` lengths the survey used.

The comparison is valid only on the file the survey read; the file's MD5 is
checked against survey/results/survey.tsv and a mismatch is recorded instead.
Rows go to RESULTS/agreement.tsv; one summary row per genome to
RESULTS/agreement_summary.tsv.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "survey"))

import benchcommon as bc  # noqa: E402
from classify_hits import row_class  # noqa: E402

SURVEY = HERE.parent / "survey" / "results"
COLUMNS = ("accession", "side", "line", "feature", "type", "survey_kind", "survey_class",
           "cgload_plain", "cgload_repair", "agreement", "detail")
SUMMARY_COLUMNS = ("accession", "file_md5", "survey_md5", "same_file", "survey_rows",
                   "cgload_features", "both_accept", "both_object",
                   "survey_documented_cgload_objects", "survey_error_cgload_accepts",
                   "cgload_only", "note")


def md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def survey_rows(accession: str, results: Path = SURVEY) -> list[dict]:
    with open(results / "hits.tsv", newline="") as handle:
        return [row for row in csv.DictReader(handle, delimiter="\t")
                if row["accession"] == accession]


def survey_md5(accession: str, results: Path = SURVEY) -> str | None:
    with open(results / "survey.tsv", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["accession"] == accession:
                return row["file_md5"]
    return None


def feature_of(gff_row: str) -> tuple[str, str]:
    """(ID, type) of a hit's row, as hits.tsv writes it: columns joined by ' | '."""
    columns = [part.strip() for part in gff_row.split(" | ")]
    ftype = columns[2] if len(columns) > 2 else ""
    for chunk in (columns[8] if len(columns) > 8 else "").split(";"):
        key, _, value = chunk.partition("=")
        if key.strip() == "ID":
            return unquote(value), ftype
    return "", ftype


def cgload_view(gff: Path) -> dict:
    """What cgload makes of the file, per feature ID."""
    from cgload.envelope import apply_policy, find_violations
    from cgload.loader import _circular_lengths, _split_origin_artefacts
    from cgload.parsers import profiles as dialects
    from cgload.parsers.normalise import NormaliseError, normalise
    from cgload.parsers.tokenizer import TokenizeError, tokenize

    try:
        tokenized = tokenize(gff, tolerate_inverted=True)
        detection = dialects.detect(tokenized)
        normalised = normalise(
            tokenized.features, detection.profile, comments=tokenized.comments,
            unexpected=dialects.unexpected_types(detection.profile, tokenized.features),
            missing=dialects.missing_types(detection.profile, tokenized.features))
    except (TokenizeError, NormaliseError, dialects.DetectionError) as error:
        return {"refused_outright": str(error)[:500], "findings": {}}

    features = normalised.features
    findings: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for violation in find_violations(features):
        if violation.is_self_contradiction:
            findings[violation.parent_id].append(("refused", "written backwards"))
        else:
            findings[violation.child_id].append(("refused", f"outside {violation.parent_id}"))
            findings[violation.parent_id].append(("refused", f"excludes {violation.child_id}"))

    kept, repairs, remaining = apply_policy(features, repair=True)
    across, remaining = _split_origin_artefacts(kept, remaining)
    repair_view: dict[str, str] = {}
    for repair in repairs:
        repair_view[repair.feature_id] = "repaired"
    for item in across:
        for name in (item["parent"], item["child"]):
            repair_view.setdefault(name, "stored as written (D-091)")
    for violation in remaining:
        for name in (violation.parent_id, violation.child_id):
            if name:
                repair_view[name] = "still refused"

    declared = {name: end - start + 1 for name, (start, end) in tokenized.declared_regions.items()}
    circular = set(_circular_lengths(kept))
    for item in kept:
        length = declared.get(item.seqid)
        if length is None:
            continue
        for segment in item.segments:
            if segment.start <= length and segment.end <= length:
                continue
            name = item.source_id or f"line {segment.line_number}"
            if item.seqid not in circular:
                findings[name].append(("refused", "past the end of a linear sequence"))
                repair_view[name] = "still refused"
            elif segment.start > length:
                findings[name].append(("accepted with a warning",
                                       "wholly past the end of a circular sequence"))
                repair_view.setdefault(name, "stored as written, with a warning")
    return {"findings": findings, "repair": repair_view,
            "types": {item.source_id: item.source_type for item in features if item.source_id}}


def after_repair(name: str, plain: str, repair_view: dict[str, str]) -> str:
    """What --repair-envelopes does about this feature. A child refused only for
    lying outside its parent is accepted once the parent is recomputed."""
    if name in repair_view:
        return repair_view[name]
    return "accepted" if plain == "accepted" else "accepted once its parent is repaired"


def compare(args) -> None:
    gff = args.inputs / args.accession / "genomic.gff"
    table = args.results / "agreement.tsv"
    summary_table = args.results / "agreement_summary.tsv"
    if bc.done(summary_table, accession=args.accession):
        print(f"{args.accession}: already in {summary_table}, skipped")
        return
    expected, actual = survey_md5(args.accession, args.survey_results), md5(gff)
    summary = dict.fromkeys(SUMMARY_COLUMNS)
    summary.update({"accession": args.accession, "file_md5": actual, "survey_md5": expected,
                    "same_file": expected == actual})
    if expected != actual:
        summary["note"] = "the file differs from the one the survey read; not compared"
        bc.append_row(summary_table, summary)
        print(f"{args.accession}: {summary['note']}")
        return
    view = cgload_view(gff)
    if "refused_outright" in view:
        summary["note"] = "cgload refused the file before the checks: " + view["refused_outright"]
        bc.append_row(summary_table, summary)
        return
    findings, repair_view = view["findings"], view["repair"]

    worst: dict[str, dict] = {}
    rank = {"A": 1, "B": 2, "C": 3, "D": 4}
    for hit in survey_rows(args.accession, args.survey_results):
        hit["class"] = row_class(hit)
        if hit["line"] not in worst or rank[hit["class"]] > rank[worst[hit["line"]]["class"]]:
            worst[hit["line"]] = hit
    agreement = Counter()
    seen_by_survey = set()
    for line, hit in sorted(worst.items(), key=lambda item: int(item[0])):
        name, ftype = feature_of(hit["gff_row"])
        seen_by_survey.add(name)
        plain = "refused" if any(v == "refused" for v, _ in findings.get(name, [])) else (
            findings[name][0][0] if findings.get(name) else "accepted")
        repaired = after_repair(name, plain, repair_view)
        survey_says_error = hit["class"] in ("C", "D")
        if plain == "accepted":
            verdict = ("both accept as valid" if not survey_says_error
                       else "survey: error; cgload accepts")
        else:
            verdict = ("both object" if survey_says_error
                       else "survey: documented; cgload objects")
        agreement[verdict] += 1
        bc.append_row(table, {
            "accession": args.accession, "side": "survey", "line": line, "feature": name,
            "type": ftype, "survey_kind": hit["kind"], "survey_class": hit["class"],
            "cgload_plain": plain, "cgload_repair": repaired, "agreement": verdict,
            "detail": "; ".join(sorted({d for _v, d in findings.get(name, [])}))})
    only_cgload = sorted(name for name in findings if name not in seen_by_survey)
    for name in only_cgload:
        verdicts = findings[name]
        plain = "refused" if any(v == "refused" for v, _ in verdicts) else verdicts[0][0]
        agreement["cgload only"] += 1
        bc.append_row(table, {
            "accession": args.accession, "side": "cgload", "line": "", "feature": name,
            "type": view["types"].get(name, ""), "survey_kind": "", "survey_class": "",
            "cgload_plain": plain, "cgload_repair": after_repair(name, plain, repair_view),
            "agreement": "cgload only",
            "detail": "; ".join(sorted({d for _v, d in verdicts}))})
    summary.update({
        "survey_rows": len(worst),
        "cgload_features": len(findings),
        "both_accept": agreement["both accept as valid"],
        "both_object": agreement["both object"],
        "survey_documented_cgload_objects": agreement["survey: documented; cgload objects"],
        "survey_error_cgload_accepts": agreement["survey: error; cgload accepts"],
        "cgload_only": agreement["cgload only"],
        "note": "",
    })
    bc.append_row(summary_table, summary)
    print(f"{args.accession}: {summary}")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accession")
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--survey-results", type=Path, default=SURVEY,
                        help="the survey's hits.tsv and survey.tsv (default %(default)s)")
    args = parser.parse_args(argv)
    compare(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
