#!/usr/bin/env python3
"""The query-timing tasks, done with gffutils, as one command each.

    python3 benchmark/gffutils_query.py DB gene ID
    python3 benchmark/gffutils_query.py DB region SEQID START END
    python3 benchmark/gffutils_query.py DB count
    python3 benchmark/gffutils_query.py DB export SEQID OUT.gff3

Written the way the gffutils documentation shows: ``db[id]`` with
``db.children``, ``db.region``, ``db.count_features_of_type``, and
``str(feature)`` for a GFF3 line. ``gene`` and ``region`` print what they
found as GFF3 lines, as cgload's query prints its rows; ``count`` prints each
type's count and then the total; ``export`` writes the file and prints how
many lines it wrote. gffutils' own command, ``gffutils-cli``, is not used: in
gffutils 0.14 its ``region`` raises NotImplementedError and its ``children``
prints grandchildren twice.
"""

from __future__ import annotations

import sys


def main(argv=None) -> int:
    import gffutils

    argv = sys.argv[1:] if argv is None else argv
    db = gffutils.FeatureDB(argv[0])
    task = argv[1]
    if task == "gene":
        gene = db[argv[2]]
        found = [gene, *db.children(gene)]
    elif task == "region":
        found = list(db.region(seqid=argv[2], start=int(argv[3]), end=int(argv[4])))
    elif task == "count":
        counts = {name: db.count_features_of_type(name) for name in db.featuretypes()}
        for name, count in sorted(counts.items()):
            print(f"{name}\t{count}")
        print(sum(counts.values()))
        return 0
    elif task == "export":
        written = 0
        with open(argv[3], "w") as handle:
            for feature in db.region(seqid=argv[2]):
                handle.write(str(feature) + "\n")
                written += 1
        print(written)
        return 0
    else:
        raise SystemExit(f"unknown task {task}")
    out = sys.stdout
    for feature in found:
        out.write(str(feature) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
