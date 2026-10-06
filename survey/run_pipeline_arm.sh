#!/usr/bin/env bash
# Survey the pipeline arm: same assemblies as the RefSeq arm, annotated by
# funannotate and BRAKER. This script does NOT run the pipelines; it
#   (1) downloads the RefSeq FASTA + GFF for each accession (the FASTA is what
#       you feed the pipelines, so seqids match), and
#   (2) surveys whatever pipeline GFF3s you have placed at
#       pipeline/<accession>/<dialect>.gff3
# Missing pipeline outputs are logged as MISSING, never skipped silently.
#
#   bash run_pipeline_arm.sh fetch      # step 1 only
#   bash run_pipeline_arm.sh survey     # step 2, after the pipelines have run
set -uo pipefail
MODE=${1:-survey}
ARM=${2:-pipeline_arm.tsv}
OUT=${3:-results}
STAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)
if [ "$MODE" = "survey" ]; then

  # --- provenance guard -------------------------------------------------------
  # Results must trace back to exactly one committed version of the survey
  # inputs. This directory is a live checkout of main, so an uncommitted edit (or
  # a pull halfway through a run) would otherwise mix two versions of the code or
  # the panel into one results table without any record of it.
  here=$(cd "$(dirname "$0")" && pwd)
  INPUTS="cgsurvey.py run_panel.sh run_pipeline_arm.sh panel.tsv pipeline_arm.tsv"
  if git -C "$here" rev-parse --git-dir >/dev/null 2>&1; then
    commit=$(git -C "$here" rev-parse HEAD)
    dirty=$(git -C "$here" status --porcelain -- $INPUTS)
    if [ -n "$dirty" ] && [ "${ALLOW_DIRTY:-0}" != "1" ]; then
      echo "refusing to run: survey inputs have uncommitted changes:" >&2
      echo "$dirty" >&2
      echo "commit them first, or set ALLOW_DIRTY=1 for a trial run you will not report." >&2
      exit 1
    fi
  else
    commit="not-a-git-checkout"; dirty=""
  fi
  mkdir -p "$OUT"
  [ -f "$OUT/provenance.tsv" ] || printf 'started\tscript\tcommit\tuncommitted_inputs\tinput_md5\n' > "$OUT/provenance.tsv"
  md5s=$(cd "$here" && python3 -c 'import hashlib,sys,os; print(" ".join(f + ":" + hashlib.md5(open(f, "rb").read()).hexdigest()[:12] for f in sys.argv[1:] if os.path.exists(f)))' $INPUTS)
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(basename "$0")" "$commit" \
    "$( [ -n "$dirty" ] && echo yes || echo no )" "$md5s" >> "$OUT/provenance.tsv"
  echo "# survey inputs at commit $commit" >&2
  # ----------------------------------------------------------------------------
fi
mkdir -p "$OUT" downloads pipeline

tail -n +2 "$ARM" | while IFS=$'\t' read -r acc species dialects note; do
  [ -z "${acc:-}" ] && continue
  d="downloads/$acc"
  if [ "$MODE" = "fetch" ]; then
    if [ ! -f "$d/genome.fna" ]; then
      echo "=== fetching $species ($acc)"
      datasets download genome accession "$acc" --include genome,gff3 \
        --filename "downloads/$acc.full.zip" || { echo "  FETCH FAILED"; continue; }
      mkdir -p "$d"
      unzip -o -j "downloads/$acc.full.zip" "ncbi_dataset/data/$acc/*" -d "$d" >/dev/null
      mv "$d"/*_genomic.fna "$d/genome.fna" 2>/dev/null
    fi
    mkdir -p "pipeline/$acc"
    echo "  ready: $d/genome.fna  ->  put outputs in pipeline/$acc/{funannotate,braker}.gff3"
    continue
  fi
  for dialect in ${dialects//,/ }; do
    g="pipeline/$acc/$dialect.gff3"
    if [ ! -f "$g" ]; then
      echo "MISSING  $species  $dialect  ($g)" | tee -a "$OUT/pipeline_missing.log"
      continue
    fi
    python3 cgsurvey.py --gff "$g" --fasta "$d/genome.fna" \
      --molecules-from "$d/genomic.gff" \
      --arm pipeline --dialect "$dialect" \
      --accession "$acc" --species "$species" --clade fungus \
      --downloaded "$STAMP" --out-dir "$OUT"
  done
done
