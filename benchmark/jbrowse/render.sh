#!/usr/bin/env bash
# Draw a region of cgload's stored assembly and annotation with JBrowse 2.
#
#   bash benchmark/jbrowse/render.sh ASSEMBLY.fa ANNOTATION.gff3.gz REGION OUT.svg
#
# ASSEMBLY.fa and ANNOTATION.gff3.gz are the files cgload writes into its data
# directory, with their .fai and .tbi beside them; nothing else is needed. This
# is the claim the paper makes for "browsable": JBrowse 2 reads cgload's files
# directly, with no server code of ours. Needs Node.js 18 or later and npm.
#
# @jbrowse/img 3.7.0 is the newest release that runs on Node 22 (4.x needs 23).
# Two of its dependencies have since moved on, so package.json pins yargs 17,
# and swaps node-fetch v2, which cannot read the data: URL that the BGZF
# decoder loads its WebAssembly from, for Node's own fetch (fetch-shim/).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
[ $# -eq 4 ] || { echo "usage: render.sh ASSEMBLY.fa ANNOTATION.gff3.gz REGION OUT.svg" >&2; exit 2; }
if [ ! -x "$here/node_modules/.bin/jb2export" ]; then
  (cd "$here" && npm install --no-audit --no-fund --loglevel=error)
fi
# npm links a file: override relative to the package that asked for it, which
# leaves node-fetch pointing nowhere; put the shim in place itself.
rm -rf "$here/node_modules/node-fetch"
cp -r "$here/fetch-shim" "$here/node_modules/node-fetch"
"$here/node_modules/.bin/jb2export" --fasta "$1" --gffgz "$2" --loc "$3" --out "$4"
echo "wrote $4"
