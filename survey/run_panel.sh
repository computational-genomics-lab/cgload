#!/usr/bin/env bash
# Run the pre-registered panel end to end.
#
# Requires the NCBI datasets CLI:
#   conda install -c conda-forge ncbi-datasets-cli
#
# Every genome in panel.tsv is downloaded, checksummed, surveyed and recorded --
# including the clean ones. Nothing is dropped for being boring. If a download
# or a resolution fails, the row is written with status=FAILED rather than
# omitted, so the denominator stays honest.

set -uo pipefail
# CHECK_ONLY=1 bash run_panel.sh   -> resolve accessions only, no downloads.
PANEL=${1:-panel.tsv}
CHECK_ONLY=${CHECK_ONLY:-0}
OUT=${2:-results}
STAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)
here=$(cd "$(dirname "$0")" && pwd)
if [ "$CHECK_ONLY" != "1" ]; then

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
  if [ -e "$OUT/survey.tsv" ] || [ -e "$OUT/failures.tsv" ]; then
    mv "$OUT" "$OUT.prev-$(date -u +%Y%m%dT%H%M%SZ)"
    echo "# previous results kept in $OUT.prev-*" >&2
  fi
  mkdir -p "$OUT"
  [ -f "$OUT/provenance.tsv" ] || printf 'started\tscript\tcommit\tuncommitted_inputs\tinput_md5\n' > "$OUT/provenance.tsv"
  md5s=$(cd "$here" && python3 -c 'import hashlib,sys,os; print(" ".join(f + ":" + hashlib.md5(open(f, "rb").read()).hexdigest()[:12] for f in sys.argv[1:] if os.path.exists(f)))' $INPUTS)
  printf '%s\t%s\t%s\t%s\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$(basename "$0")" "$commit" \
    "$( [ -n "$dirty" ] && echo yes || echo no )" "$md5s" >> "$OUT/provenance.tsv"
  echo "# survey inputs at commit $commit" >&2
  # ----------------------------------------------------------------------------
fi
mkdir -p "$OUT" downloads

echo "# panel run started $STAMP" | tee "$OUT/run.log"

n_panel=$(tail -n +2 "$PANEL" | awk -F'\t' 'NF >= 4 && $4 != ""' | wc -l | tr -d ' ')

tail -n +2 "$PANEL" | while IFS=$'\t' read -r clade species asm acc conf orgs why; do
  [ -z "${acc:-}" ] && continue
  echo "=== $species ($acc)" | tee -a "$OUT/run.log"

  # 1. Resolve through NCBI's own formatter (no guessing at JSON field names).
  info=$(datasets summary genome accession "$acc" --as-json-lines 2>>"$OUT/run.log" \
         | dataformat tsv genome --elide-header \
             --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name \
             2>>"$OUT/run.log" | head -n 1)
  if [ -z "$info" ]; then
    echo "  !! RESOLVE FAILED: accession not found at NCBI" | tee -a "$OUT/run.log"
    [ "$CHECK_ONLY" = "1" ] || printf '%s\t%s\t%s\tRESOLVE_FAILED\n' "$clade" "$species" "$acc" >> "$OUT/failures.tsv"
    continue
  fi
  got_asm=$(printf '%s' "$info" | cut -f2)
  ann_name=$(printf '%s' "$info" | cut -f3)
  ann_date=$(printf '%s' "$info" | cut -f4)
  got_org=$(printf '%s' "$info" | cut -f5)
  echo "  declared=$asm resolved=$got_asm | organism=$got_org | annotation=${ann_name:-NONE} $ann_date" \
    | tee -a "$OUT/run.log"
  problem=""
  if [ "$asm" != "$got_asm" ]; then
    echo "  !! ASSEMBLY MISMATCH: fix this panel row" | tee -a "$OUT/run.log"; problem=ASSEMBLY_MISMATCH
  fi
  if [ -z "$ann_name" ]; then
    echo "  !! NO ANNOTATION on this assembly: replace it in the panel" | tee -a "$OUT/run.log"; problem=NO_ANNOTATION
  fi
  [ "$CHECK_ONLY" = "1" ] && continue
  if [ -n "$problem" ]; then
    printf '%s\t%s\t%s\t%s\n' "$clade" "$species" "$acc" "$problem" >> "$OUT/failures.tsv"
    continue
  fi

  # 2. Download annotation only (cached across runs).
  if [ ! -s "downloads/$acc/genomic.gff" ]; then
    if ! datasets download genome accession "$acc" --include gff3 --no-progressbar \
         --filename "downloads/$acc.zip" >>"$OUT/run.log" 2>&1; then
      printf '%s\t%s\t%s\tDOWNLOAD_FAILED\n' "$clade" "$species" "$acc" >> "$OUT/failures.tsv"
      continue
    fi
    mkdir -p "downloads/$acc"
    unzip -o -j "downloads/$acc.zip" "ncbi_dataset/data/$acc/genomic.gff" \
      -d "downloads/$acc" >>"$OUT/run.log" 2>&1
  fi
  if [ ! -s "downloads/$acc/genomic.gff" ]; then
    echo "  !! NO GFF in the download" | tee -a "$OUT/run.log"
    printf '%s\t%s\t%s\tNO_GFF\n' "$clade" "$species" "$acc" >> "$OUT/failures.tsv"
    continue
  fi

  # 3. Survey. A crash is recorded as a failure, never silently dropped.
  if ! python3 "$here/cgsurvey.py" --gff "downloads/$acc/genomic.gff" \
       --accession "$acc" --species "$species" --clade "$clade" \
       --assembly-name "$got_asm" --annotation-release "$ann_name $ann_date" \
       --arm refseq --dialect refseq \
       --downloaded "$STAMP" --out-dir "$OUT" 2>&1 | tee -a "$OUT/run.log"; then
    printf '%s\t%s\t%s\tSURVEY_FAILED\n' "$clade" "$species" "$acc" >> "$OUT/failures.tsv"
  fi
done

if [ "$CHECK_ONLY" = "1" ]; then
  n_bad=$(grep -c '  !! ' "$OUT/run.log" || true)
  echo "# check done: $n_panel panel rows, $n_bad problem(s)."
  [ "$n_bad" -gt 0 ] && grep -B1 '  !! ' "$OUT/run.log" | grep -v '^--$'
  [ "$n_bad" -eq 0 ] && echo "# panel is clean: commit and push panel.tsv now."
else
  surveyed=0; failed=0
  [ -f "$OUT/survey.tsv" ] && surveyed=$(( $(wc -l < "$OUT/survey.tsv") - 1 ))
  [ -f "$OUT/failures.tsv" ] && failed=$(wc -l < "$OUT/failures.tsv" | tr -d ' ')
  echo "# done: $n_panel panel rows = $surveyed surveyed + $failed failed"
  [ "$failed" -gt 0 ] && { echo "# failures:"; cat "$OUT/failures.tsv"; }
  [ $((surveyed + failed)) -ne "$n_panel" ] && echo "!! counts do not add up: read $OUT/run.log"
fi
