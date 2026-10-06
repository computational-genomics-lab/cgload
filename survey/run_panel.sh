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

tail -n +2 "$PANEL" | while IFS=$'\t' read -r clade species asm acc conf orgs why; do
  [ -z "${acc:-}" ] && continue
  echo "=== $species ($acc)" | tee -a "$OUT/run.log"

  # 1. Resolve: confirm the accession really is the assembly we named.
  #    A silent substitution here is the one error that would poison the panel.
  real=$(datasets summary genome accession "$acc" --as-json-lines 2>/dev/null \
         | python3 -c 'import sys,json
for l in sys.stdin:
    d=json.loads(l)
    print(d.get("assembly_info",{}).get("assembly_name",""),
          d.get("organism",{}).get("organism_name",""), sep="\t")
    break' )
  if [ -z "$real" ]; then
    echo "  RESOLVE FAILED" | tee -a "$OUT/run.log"
    printf '%s\t%s\t%s\tRESOLVE_FAILED\n' "$clade" "$species" "$acc" >> "$OUT/failures.tsv"
    continue
  fi
  got_asm=$(printf '%s' "$real" | cut -f1)
  echo "  declared=$asm resolved=$got_asm" | tee -a "$OUT/run.log"
  if [ "$asm" != "$got_asm" ]; then
    echo "  !! assembly name mismatch -- panel row needs correcting BEFORE the run" \
      | tee -a "$OUT/run.log"
  fi

  [ "$CHECK_ONLY" = "1" ] && continue

  # 2. Download annotation only.
  if [ ! -f "downloads/$acc/genomic.gff" ]; then
    datasets download genome accession "$acc" --include gff3 \
      --filename "downloads/$acc.zip" >>"$OUT/run.log" 2>&1 || {
        printf '%s\t%s\t%s\tDOWNLOAD_FAILED\n' "$clade" "$species" "$acc" >> "$OUT/failures.tsv"
        continue; }
    mkdir -p "downloads/$acc"
    unzip -o -j "downloads/$acc.zip" "ncbi_dataset/data/$acc/genomic.gff" \
      -d "downloads/$acc" >>"$OUT/run.log" 2>&1
  fi

  # 3. Pin the annotation release from the GFF header, not from the website.
  rel=$(grep -m1 -i 'annotation-source\|#!annotation-date\|Annotation Release' \
        "downloads/$acc/genomic.gff" | tr -d '\r' | sed 's/^#*!* *//')

  # 4. Survey.
  python3 cgsurvey.py --gff "downloads/$acc/genomic.gff" \
    --accession "$acc" --species "$species" --clade "$clade" \
    --assembly-name "$got_asm" --annotation-release "$rel" \
    --arm refseq --dialect refseq \
    --downloaded "$STAMP" --out-dir "$OUT" 2>&1 | tee -a "$OUT/run.log"
done

if [ "$CHECK_ONLY" = "1" ]; then
  echo "# check done. Fix every 'mismatch' and every RESOLVE FAILED in panel.tsv, then commit it."
  echo "# mismatches: $(grep -c 'mismatch' "$OUT/run.log" || true)"
else
  echo "# done. rows in survey.tsv: $(( $(wc -l < "$OUT/survey.tsv") - 1 ))"
fi
