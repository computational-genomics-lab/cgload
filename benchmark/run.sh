#!/usr/bin/env bash
# The cgload benchmark, start to finish, as PROTOCOL.md fixes it.
#
#   bash benchmark/run.sh [WORK]          (WORK defaults to ./benchmark-run)
#
# Run from the root of a clean checkout of the tagged commit. Every step logs
# to WORK/logs/, writes its tables to WORK/results/, and skips what an earlier
# run already finished, so after a stop the same command carries on. At the
# end WORK/cgload-benchmark-results.tar.gz holds everything to send back.
#
# Optional settings, as environment variables:
#   GENOMES            accessions to run (default: all 30 in benchmark/genomes.tsv)
#   JOBS               parallel jobs for the steps that are not timed (default 8)
#   R_JOBS             parallel txdbmaker runs, each a large R process (default 3)
#   CHADO_JOBS         parallel python-chado loads (default 6)
#   SURVEY_DOWNLOADS   the survey's downloads/ folder, holding ACCESSION/genomic.gff
#                      (default: survey/downloads in this checkout, if present)
#   ZENODO_RECORD      the survey's Zenodo record, a second source for the GFF3
#   CASE_STUDY         a cgload organism config for the case study (Phase 9)
#   PG_CONNINFO        use this PostgreSQL ("host=... port=... user=...") instead
#                      of starting one under WORK
#   USE_CURRENT_ENV=1  use the tools already on PATH instead of making the
#                      conda environments (needs BENCH_CHADO_PYTHON, Rscript, tabix)
#   SKIP               comma-separated steps to leave out, e.g. "chado,kills"
#   ALLOW_DIRTY=1      run from a checkout with uncommitted changes (not for results)
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
BENCH="$REPO/benchmark"
WORK=$(realpath -m "${1:-$PWD/benchmark-run}")
RESULTS="$WORK/results"
INPUTS="$WORK/inputs"
LOGS="$WORK/logs"
API=${NCBI_DATASETS_API:-https://api.ncbi.nlm.nih.gov/datasets}
JOBS=${JOBS:-8}
R_JOBS=${R_JOBS:-3}
CHADO_JOBS=${CHADO_JOBS:-6}
SKIP=",${SKIP:-},"
mkdir -p "$RESULTS" "$INPUTS" "$LOGS" "$WORK/tools" || { echo "cannot create $WORK"; exit 1; }
# One run per folder. Two copies writing the same tables and databases spoil
# both, so a second copy stops here. The lock goes when this script ends.
exec 9> "$WORK/.run.lock"
if ! flock -n 9; then
  echo "Another run.sh is already running in $WORK. Let it finish, or stop it"
  echo "(pkill -f benchmark/run.sh), before starting this one."
  exit 1
fi
exec > >(tee -a "$LOGS/run.log") 2>&1

say() { echo "$(date -u +%FT%TZ)  $*"; }
skipped() { [[ "$SKIP" == *",$1,"* ]]; }
declare -a SUMMARY=()
step() {  # step NAME command... : run, log, record the exit code, carry on
  local name=$1; shift
  if skipped "$name"; then SUMMARY+=("$(printf '%-28s %s' "$name" "skipped (SKIP)")"); return 0; fi
  say "=== $name"
  "$@" >> "$LOGS/$name.log" 2>&1
  local code=$?
  SUMMARY+=("$(printf '%-28s %s' "$name" "$([ $code -eq 0 ] && echo ok || echo "FAILED ($code), see logs/$name.log")")")
  say "=== $name: exit $code"
  return 0
}

cd "$REPO" || exit 1
if [ -n "$(git status --porcelain --untracked-files=no)" ] && [ "${ALLOW_DIRTY:-0}" != 1 ]; then
  echo "This checkout has uncommitted changes. The benchmark runs from a committed,"
  echo "tagged version only (PROTOCOL.md, section 3). Commit or stash, or set ALLOW_DIRTY=1"
  echo "for a trial run whose results will not be used."
  exit 1
fi
say "# cgload benchmark, $(git describe --always --dirty --tags) ($(git rev-parse HEAD))"
say "# host: $(uname -srm); cpus: $(nproc); $(free -g | awk '/^Mem:/ {print "RAM " $2 " GB"}')"
say "# disk under $WORK: $(df -h "$WORK" | awk 'NR == 2 {print $4 " free of " $2}')"

# BENCH_GENOMES points at another genome table, for the benchmark's own tests.
GENOME_TABLE=${BENCH_GENOMES:-$BENCH/genomes.tsv}
if [ -n "${GENOMES:-}" ]; then
  read -r -a ACCESSIONS <<<"$GENOMES"
else
  mapfile -t ACCESSIONS < <(grep -v '^#' "$GENOME_TABLE" | tail -n +2 | cut -f1)
fi
role() {  # role NAME : the accessions that take part in that part of the benchmark
  local wanted=" ${ACCESSIONS[*]} "
  grep -v '^#' "$GENOME_TABLE" | tail -n +2 | awk -F'\t' -v r="$1" '$8 ~ "(^|,)"r"(,|$)" {print $1}' \
    | while read -r acc; do [[ "$wanted" == *" $acc "* ]] && echo "$acc"; done
}
say "# genomes: ${#ACCESSIONS[@]}"

# --- 1. environments -----------------------------------------------------------------
MAMBA=""
find_mamba() {
  local candidate
  for candidate in micromamba mamba conda; do
    if command -v "$candidate" >/dev/null 2>&1; then MAMBA=$(command -v "$candidate"); return 0; fi
  done
  if [ -x "$WORK/tools/micromamba" ]; then MAMBA="$WORK/tools/micromamba"; return 0; fi
  say "no conda, mamba or micromamba on PATH; downloading micromamba"
  curl -fsSL --retry 3 -o "$WORK/tools/micromamba" \
    https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64 \
    && chmod +x "$WORK/tools/micromamba" && MAMBA="$WORK/tools/micromamba"
}
make_env() {  # make_env PREFIX FILE
  local prefix=$1 file=$2
  if [ -x "$prefix/bin/python" ] || [ -x "$prefix/bin/R" ]; then echo "$prefix exists"; return 0; fi
  export MAMBA_ROOT_PREFIX="$WORK/tools/mamba-root"
  case "$(basename "$MAMBA")" in
    micromamba) "$MAMBA" create -y -p "$prefix" -f "$file" ;;
    *) "$MAMBA" env create -p "$prefix" -f "$file" ;;
  esac
}
environments() {
  find_mamba || { echo "could not get micromamba"; return 1; }
  echo "using $MAMBA"
  local failed=0
  make_env "$WORK/env" "$BENCH/environment.yml" || failed=1
  make_env "$WORK/env-r" "$BENCH/environment-r.yml" || echo "the R environment failed; txdbmaker will not run"
  make_env "$WORK/env-chado" "$BENCH/environment-chado.yml" || echo "the Chado environment failed; Chado will not run"
  "$WORK/env/bin/pip" install -q -e "${REPO}[postgres]" || failed=1
  mkdir -p "$RESULTS/environment"
  for prefix in env env-r env-chado; do
    [ -d "$WORK/$prefix" ] && "$MAMBA" list -p "$WORK/$prefix" --explicit > "$RESULTS/environment/$prefix.txt" 2>/dev/null
  done
  "$WORK/env/bin/pip" freeze > "$RESULTS/environment/env-pip.txt"
  [ -x "$WORK/env-chado/bin/pip" ] && "$WORK/env-chado/bin/pip" freeze > "$RESULTS/environment/env-chado-pip.txt"
  return $failed
}
if [ "${USE_CURRENT_ENV:-0}" = 1 ]; then
  PY=$(command -v python3); SUMMARY+=("$(printf '%-28s %s' environments "tools on PATH (USE_CURRENT_ENV=1)")")
else
  step environments environments
  PY="$WORK/env/bin/python"
  export PATH="$WORK/env/bin:$PATH"
  export BENCH_RSCRIPT=${BENCH_RSCRIPT:-$WORK/env-r/bin/Rscript}
  export BENCH_CHADO_PYTHON=${BENCH_CHADO_PYTHON:-$WORK/env-chado/bin/python}
  export BENCH_TABIX=${BENCH_TABIX:-$WORK/env/bin/tabix}
  export BENCH_PSQL=${BENCH_PSQL:-$WORK/env/bin/psql}
fi
[ -x "$PY" ] || { echo "no Python for the benchmark at $PY"; exit 1; }
"$PY" -c "import cgload, gffutils" || { echo "cgload or gffutils is missing from $PY"; exit 1; }
{
  echo "cgload $(cgload --version)"; echo "commit $(git rev-parse HEAD)"
  "$PY" -c "import gffutils, sys; print('gffutils', gffutils.__version__); print('python', sys.version.split()[0])"
} > "$RESULTS/versions.txt" 2>&1

# --- 2. PostgreSQL -------------------------------------------------------------------
PG_STARTED=0
stop_postgres() { [ "$PG_STARTED" = 1 ] && pg_ctl -D "$WORK/pgdata" -m fast stop >/dev/null 2>&1; }
trap stop_postgres EXIT
postgres() {
  if [ -n "${PG_CONNINFO:-}" ]; then export BENCH_PG_CONNINFO="$PG_CONNINFO"; return 0; fi
  # The socket lives under /tmp: a socket path longer than 107 bytes is refused.
  local sock port=${PG_PORT:-55433}
  sock="${TMPDIR:-/tmp}/cgb-$(id -u)"
  mkdir -p "$sock"
  if [ ! -d "$WORK/pgdata" ]; then
    initdb -D "$WORK/pgdata" -A trust -U "$(id -un)" >/dev/null || return 1
  fi
  # 9>&-: the server must not inherit the run's lock, or a server left running
  # would block the next run.
  pg_ctl -D "$WORK/pgdata" -o "-p $port -k $sock -c listen_addresses=''" -l "$LOGS/postgres-server.log" -w start 9>&- || return 1
  PG_STARTED=1
  BENCH_PG_CONNINFO="host=$sock port=$port user=$(id -un)"
  export BENCH_PG_CONNINFO
  echo "PostgreSQL $(postgres --version 2>/dev/null) at $BENCH_PG_CONNINFO"
}
step postgres postgres
if [ -z "${BENCH_PG_CONNINFO:-}" ]; then
  say "PostgreSQL did not start (see logs/postgres.log); half the benchmark needs it, so the run stops."
  exit 1
fi
echo "postgres: $BENCH_PG_CONNINFO" >> "$RESULTS/versions.txt"

# --- 3. inputs -----------------------------------------------------------------------
SURVEY_DOWNLOADS=${SURVEY_DOWNLOADS:-$REPO/survey/downloads}
survey_md5() { awk -F'\t' -v a="$1" 'NR > 1 && $7 == a {print $11}' "$REPO/survey/results/survey.tsv"; }
fetch_ncbi() {  # fetch_ncbi ACC DIR TYPES...  (GENOME_FASTA GENOME_GFF GENOME_GBFF)
  local acc=$1 dir=$2; shift 2
  local query="" type
  for type in "$@"; do query+="&include_annotation_type=$type"; done
  curl -fsSL --retry 5 --retry-delay 30 -o "$dir/download.zip" \
    "$API/v2/genome/accession/$acc/download?${query#&}" || return 1
  unzip -tq "$dir/download.zip" >/dev/null || return 1
  unzip -o -q -j "$dir/download.zip" "ncbi_dataset/data/$acc/*" -d "$dir" || return 1
  rm -f "$dir/download.zip"
}
inputs_one() {
  local acc=$1 dir="$INPUTS/$1" want source="" got
  mkdir -p "$dir"
  want=$(survey_md5 "$acc")
  if [ ! -s "$dir/genomic.gff" ]; then
    if [ -s "$SURVEY_DOWNLOADS/$acc/genomic.gff" ]; then
      cp "$SURVEY_DOWNLOADS/$acc/genomic.gff" "$dir/genomic.gff"; source="survey downloads"
    elif [ -n "${ZENODO_RECORD:-}" ] && curl -fsSL --retry 3 -o "$dir/genomic.gff.gz" \
        "https://zenodo.org/records/$ZENODO_RECORD/files/${acc}_genomic.gff.gz?download=1"; then
      gunzip -f "$dir/genomic.gff.gz"; source="zenodo $ZENODO_RECORD"
    else
      fetch_ncbi "$acc" "$dir" GENOME_GFF || { echo "$acc: no GFF3"; return 1; }
      source="NCBI $(date -u +%F)"
    fi
    echo "$source" > "$dir/gff_source.txt"
  fi
  if ! ls "$dir"/*_genomic.fna >/dev/null 2>&1 && [ ! -s "$dir/assembly.fa" ]; then
    fetch_ncbi "$acc" "$dir" GENOME_FASTA || { echo "$acc: no FASTA"; return 1; }
  fi
  if [ ! -s "$dir/assembly.fa" ]; then mv "$dir"/*_genomic.fna "$dir/assembly.fa" || return 1; fi
  if [[ " $(role crossparser | tr '\n' ' ') " == *" $acc "* ]] && [ ! -s "$dir/genomic.gbff" ]; then
    fetch_ncbi "$acc" "$dir" GENOME_GBFF || echo "$acc: no GenBank file; the cross-parser tier will not run"
  fi
  got=$(md5sum "$dir/genomic.gff" | cut -d' ' -f1)
  printf '%s\t%s\t%s\t%s\t%s\n' "$acc" "$(cat "$dir/gff_source.txt" 2>/dev/null)" "$got" \
    "$([ "$got" = "$want" ] && echo yes || echo NO)" "$(md5sum "$dir/assembly.fa" | cut -d' ' -f1)"
}
inputs() {
  local table="$RESULTS/inputs.tsv" acc failed=0
  printf 'accession\tgff_source\tgff_md5\tsame_as_survey\tfasta_md5\n' > "$table"
  for acc in "${ACCESSIONS[@]}"; do
    inputs_one "$acc" >> "$table" || failed=1
  done
  column -t -s $'\t' "$table" 2>/dev/null || cat "$table"
  return $failed
}
step inputs inputs

# --- 4. the reference count ------------------------------------------------------------
truth() {
  local acc todo=()
  mkdir -p "$RESULTS/truth"
  for acc in "${ACCESSIONS[@]}"; do
    [ -s "$RESULTS/truth/$acc.json" ] || todo+=("$acc")
  done
  [ ${#todo[@]} -eq 0 ] && return 0
  printf '%s\n' "${todo[@]}" | xargs -P "$JOBS" -I{} "$PY" "$BENCH/truth.py" \
    "$INPUTS/{}/genomic.gff" "$RESULTS/truth/{}.sqlite" --summary "$RESULTS/truth/{}.json" \
    > /dev/null
}
step truth truth

# --- 5. agreement with the survey --------------------------------------------------------
agreement() {
  role agreement | xargs -P "$JOBS" -I{} "$PY" "$BENCH/agreement.py" {} \
    --inputs "$INPUTS" --results "$RESULTS"
}
step agreement agreement

# --- 6. cgload at scale: one genome at a time, nothing else running ------------------------
scale() {
  local backend=$1 acc
  for acc in $(role scale); do
    "$PY" "$BENCH/run_cgload.py" "$acc" --backend "$backend" --inputs "$INPUTS" --work "$WORK" \
      --results "$RESULTS" > /dev/null || echo "$acc $backend: runner failed"
  done
}
step scale-sqlite scale sqlite
step scale-postgres scale postgres

# --- 7. faults and interrupted loads (cgload) ----------------------------------------------
faults() {
  local acc
  for acc in $(role faults); do
    "$PY" "$BENCH/faults.py" "$acc" --inputs "$INPUTS" --work "$WORK" --results "$RESULTS"
  done
}
kills() {  # kills TOOL TRIALS [ACCESSION...]; default: the genomes with the kills role
  local tool=$1 trials=$2 acc accessions
  shift 2
  if [ $# -gt 0 ]; then accessions=("$@"); else mapfile -t accessions < <(role kills); fi
  for acc in "${accessions[@]}"; do
    "$PY" "$BENCH/kills.py" "$acc" --tool "$tool" --trials "$trials" --inputs "$INPUTS" \
      --work "$WORK" --results "$RESULTS"
  done
}
step faults faults
kills_cgload() { kills cgload 10; kills cgload-force 10; }
step kills-cgload kills_cgload

# --- 8. the other loaders ---------------------------------------------------------------
KEEP=$(role timing | paste -sd, -)
compare() {  # compare TOOL JOBS
  local tool=$1 jobs=$2
  role compare | xargs -P "$jobs" -I{} "$PY" "$BENCH/run_compare.py" {} \
    --tool "$tool" --inputs "$INPUTS" --work "$WORK" --results "$RESULTS" --keep-db "$KEEP" \
    --timeout "${3:-43200}" > /dev/null
}
step gffutils-default compare gffutils-default "$JOBS"
step gffutils compare gffutils "$JOBS"
step kills-gffutils kills gffutils 10
step txdbmaker compare txdbmaker "$R_JOBS"

# --- 9. query timing: alone on the machine -----------------------------------------------
timing() {
  local acc
  for acc in $(role timing); do
    "$PY" "$BENCH/query_timing.py" "$acc" --work "$WORK" --results "$RESULTS"
  done
}
step timing timing

# --- 10. Chado: long, so last -------------------------------------------------------------
if [ -z "${BENCH_CHADO_SCHEMA:-}" ]; then
  export BENCH_CHADO_SCHEMA="$WORK/tools/chado-default_schema-1.31.sql"
  export BENCH_SO_OBO="$WORK/tools/so.obo"
fi
chado_inputs() {
  # The GMOD Chado 1.31 schema and the Sequence Ontology, each at a fixed commit.
  [ -s "$BENCH_CHADO_SCHEMA" ] || curl -fsSL --retry 3 -o "$BENCH_CHADO_SCHEMA" \
    https://raw.githubusercontent.com/GMOD/Chado/ef04588c62be103359252af6ef07337b0c8df500/schemas/1.31/default_schema.sql || return 1
  [ -s "$BENCH_SO_OBO" ] || curl -fsSL --retry 3 -o "$BENCH_SO_OBO" \
    https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/4340c143bac3578bba0c013d02e8c5f0e51ad14a/Ontology_Files/so.obo || return 1
  md5sum "$BENCH_CHADO_SCHEMA" "$BENCH_SO_OBO"
}
step chado-inputs chado_inputs
step chado compare chado "$CHADO_JOBS" "${CHADO_TIMEOUT:-86400}"
kills_chado() {  # F. graminearum only (PROTOCOL.md, section 2): Chado loads are slow
  local acc
  for acc in $(role kills_chado); do kills chado 5 "$acc"; done
}
step kills-chado kills_chado

# --- 11. the case study ---------------------------------------------------------------
if [ -n "${CASE_STUDY:-}" ]; then
  step case-study "$PY" "$BENCH/case_study.py" "$CASE_STUDY" --work "$WORK" --results "$RESULTS"
else
  SUMMARY+=("$(printf '%-28s %s' case-study "not run: CASE_STUDY is not set")")
fi

# --- 12. tables and the bundle to send back -------------------------------------------------
step summarise "$PY" "$BENCH/summarise.py" "$RESULTS"
cp "$LOGS"/*.log "$RESULTS/" 2>/dev/null
tar -czf "$WORK/cgload-benchmark-results.tar.gz" -C "$WORK" --exclude='results/truth/*.sqlite' results
echo
say "=== summary"
printf '%s\n' "${SUMMARY[@]}"
say "send back: $WORK/cgload-benchmark-results.tar.gz ($(du -h "$WORK/cgload-benchmark-results.tar.gz" | cut -f1))"
