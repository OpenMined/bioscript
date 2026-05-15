#!/usr/bin/env bash
# test-vntyper.sh — prove Java Kestrel and BioScript/Rust Kestrel both call
# every shipped real-data VNtyper fixture the way upstream says they should.
#
#   ./test-vntyper.sh --java --bam             # Java: all fixtures, assert correct
#   ./test-vntyper.sh --rust --bam             # BioScript/Rust: same
#   ./test-vntyper.sh --java --rust --bam      # both + correctness parity (default)
#   ./test-vntyper.sh --java --rust --bam --case 66bf   # one fixture
#   ./test-vntyper.sh --small                  # fast small-fixture suites only
#   ./test-vntyper.sh --all                    # small + both engines + vendor
#
# Every fixture upstream ships a kestrel_assertions entry for
# (ports/vntyper/vntyper/tests/test_data_config.json) is run and asserted
# against upstream's expected Confidence and Alt/ActiveRegion/Depth_Score
# tolerances. Positives must detect the variant; negatives must stay
# Negative. A wrong call is a hard FAIL, never a skip.
#
# "Java" = the Java-Kestrel pipeline (java + kestrel.jar; BAM also needs
# samtools + bcftools). "Rust" = the BioScript native pipeline through
# kestrel-rs via python/bioscript/_native.so. For FASTQ, "Java" is the same
# coordinator with the Java engine selected (no separate entry point).
#
# Parity contract: for every fixture both engines must make the upstream-
# correct call and agree on the positive/negative classification. Exact
# REF/ALT can differ (same dup frameshift reported against an equivalent
# motif reference) and the BAM TSV sha differs by the tracked samtools-rs
# FASTQ-extraction gap — neither is a parity failure; a wrong or
# disagreeing call is.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

RUN_JAVA=0
RUN_RUST=0
INPUT_BAM=0
INPUT_FASTQ=0
RUN_SMALL=0
RUN_ALL=0
RUN_VENDOR=0
DO_REBUILD=0
VERBOSE=0
CASE_FILTER=""

if [[ -t 1 ]]; then
  C_RED=$'\033[31m';  C_GRN=$'\033[32m';  C_YLW=$'\033[33m'
  C_BLU=$'\033[34m';  C_BLD=$'\033[1m';   C_DIM=$'\033[2m';  C_RST=$'\033[0m'
else
  C_RED=; C_GRN=; C_YLW=; C_BLU=; C_BLD=; C_DIM=; C_RST=
fi

usage() {
  cat <<EOF
${C_BLD}test-vntyper.sh${C_RST} — Java Kestrel vs BioScript/Rust Kestrel VNtyper parity.

${C_BLD}USAGE${C_RST}
  ./test-vntyper.sh [ENGINE...] [INPUT...] [OPTIONS]

${C_BLD}ENGINE${C_RST}
  -j, --java        Run the Java-Kestrel pipeline and print its output.
                    Needs java + kestrel.jar (BAM also needs samtools+bcftools).
  -r, --rust        Run the BioScript native (kestrel-rs) pipeline and print
                    its output. Needs python/bioscript/_native.so.
                    --java and --rust together runs both, then diffs the
                    normalized output and exits non-zero on any divergence.

${C_BLD}INPUT${C_RST}
  -b, --bam         BAM entry point (all upstream fixtures; the path
                    upstream's kestrel_assertions are defined for).
  -f, --fastq       FASTQ entry point (fixtures that ship a FASTQ pair).

${C_BLD}MODES${C_RST}
      --small       Just the fast small-fixture unittest suites.
      --all         --small + --java --rust --bam --fastq + --vendor.
      --vendor      Also run the kestrel-rs vendor parity gate
                    (KESTREL_RUN_VNTYPER_FASTQ_PARITY=1, cargo test --release).

${C_BLD}OPTIONS${C_RST}
      --case SUB    Restrict to fixtures whose name contains SUB
                    (e.g. 66bf, dfc3, a5c1, b178, 7a61, 40cf).
      --rebuild     Rebuild python/bioscript/_native.so via maturin first.
  -v, --verbose     Stream full engine output (default: tail on failure).
  -h, --help        Show this help.

${C_BLD}OUTPUT${C_RST}
  Per-engine JSON lands under /tmp/vntyper-run-<timestamp>/. The terminal
  shows each fixture's expected vs actual call and OK/FAIL. With both
  engines a per-fixture correctness-parity table is printed
  and the script exits non-zero if any case diverges.

${C_BLD}NOTES${C_RST}
  - Large fixtures live under ports/vntyper/test-data/. Missing data/tools
    are reported as a concrete prerequisite list, not a silent skip.
EOF
}

if [[ $# -eq 0 ]]; then
  usage; exit 0
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    -j|--java)    RUN_JAVA=1 ;;
    -r|--rust)    RUN_RUST=1 ;;
    -b|--bam)     INPUT_BAM=1 ;;
    -f|--fastq)   INPUT_FASTQ=1 ;;
    --small)      RUN_SMALL=1 ;;
    --all)        RUN_ALL=1; RUN_SMALL=1; RUN_JAVA=1; RUN_RUST=1; INPUT_BAM=1; INPUT_FASTQ=1; RUN_VENDOR=1 ;;
    --vendor)     RUN_VENDOR=1 ;;
    --case)       shift; CASE_FILTER="${1:-}" ;;
    --case=*)     CASE_FILTER="${1#*=}" ;;
    --rebuild)    DO_REBUILD=1 ;;
    -v|--verbose) VERBOSE=1 ;;
    -h|--help)    usage; exit 0 ;;
    *) echo "${C_RED}Unknown flag: $1${C_RST}" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

# --case is a fixture-name substring filter (e.g. 66bf, dfc3, a5c1, 7a61).
# No value restriction: an unknown filter just yields "no fixtures match".
if [[ $RUN_SMALL -eq 0 && $RUN_JAVA -eq 0 && $RUN_RUST -eq 0 && $RUN_VENDOR -eq 0 ]]; then
  echo "${C_RED}Pick at least one of --java / --rust / --small / --all / --vendor.${C_RST}" >&2
  exit 2
fi
if [[ ($RUN_JAVA -eq 1 || $RUN_RUST -eq 1) && $INPUT_BAM -eq 0 && $INPUT_FASTQ -eq 0 ]]; then
  echo "${C_RED}Engine selected but no input — add --bam and/or --fastq.${C_RST}" >&2
  exit 2
fi

TS=$(date +%Y%m%d-%H%M%S)
OUT_DIR="/tmp/vntyper-run-$TS"
mkdir -p "$OUT_DIR"

export PYTHONPATH="$ROOT/python:$ROOT/ports/vntyper/bioscript${PYTHONPATH:+:$PYTHONPATH}"
export CC=${CC:-cc}
export AR=${AR:-ar}

HELPER="$ROOT/ports/vntyper/tests/run_parity_pipeline.py"
DIFFER="$ROOT/ports/vntyper/tests/diff_parity_outputs.py"

declare -a STEP_LABELS=()
declare -a STEP_STATUS=()
declare -a STEP_SECS=()

print_header() {
  printf '\n%s%s━━ %s ━━%s\n' "$C_BLU" "$C_BLD" "$*" "$C_RST"
}

# run_step <label> <log-name> <cmd...>
run_step() {
  local label="$1"; shift
  local log="$OUT_DIR/$1"; shift
  local start end secs status rc

  print_header "$label"
  echo "${C_DIM}\$ $*${C_RST}"
  echo "${C_DIM}log: $log${C_RST}"

  start=$(date +%s)
  set +e
  if [[ $VERBOSE -eq 1 ]]; then
    "$@" 2>&1 | tee "$log"
    rc=${PIPESTATUS[0]}
  else
    "$@" >"$log" 2>&1
    rc=$?
  fi
  set -e
  end=$(date +%s)
  secs=$((end - start))

  if [[ $rc -eq 0 ]]; then
    status="PASS"; echo "${C_GRN}✓ PASS${C_RST} (${secs}s)"
  else
    status="FAIL"; echo "${C_RED}✗ FAIL${C_RST} (${secs}s, exit $rc)"
    if [[ $VERBOSE -eq 0 ]]; then
      echo "${C_DIM}── last 30 lines ──${C_RST}"
      tail -n 30 "$log" || true
    fi
  fi

  STEP_LABELS+=("$label")
  STEP_STATUS+=("$status")
  STEP_SECS+=("$secs")
  return 0
}

# Print the human-readable per-case summary that the helper wrote to JSON.
show_engine_output() {
  local json="$1"
  [[ -f "$json" ]] || return 0
  python3 - "$json" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
print(f"  all_correct={data.get('all_correct')}")
for stem, c in sorted(data.get("cases", {}).items()):
    exp = c.get("expected", {})
    called = c.get("called")
    verdict = "OK  " if c.get("correct") else "FAIL"
    if called:
        call = (f"{called.get('CHROM')}:{called.get('POS')} "
                f"{called.get('REF')}>{called.get('ALT')} "
                f"conf={called.get('Confidence')} "
                f"alt={called.get('Estimated_Depth_AlternateVariant')}")
    else:
        call = "no positive call"
    print(f"  [{verdict}] {stem}")
    print(f"         expect {exp.get('confidence')!r} -> {call} "
          f"({c.get('wall_seconds')}s)")
    for reason in c.get("reasons", []):
        print(f"         - {reason}")
PY
}

if [[ $DO_REBUILD -eq 1 ]]; then
  print_header "Rebuilding native extension"
  if ! command -v maturin >/dev/null 2>&1; then
    echo "${C_RED}maturin not on PATH — install with: pipx install maturin${C_RST}"
    exit 3
  fi
  rm -rf /tmp/bioscript-maturin-wheel
  maturin build --release -o /tmp/bioscript-maturin-wheel
  install -m 755 "$ROOT/rust/target/release/lib_native.so" "$ROOT/python/bioscript/_native.so"
  echo "${C_GRN}✓ _native.so updated${C_RST}"
fi

if [[ $RUN_SMALL -eq 1 ]]; then
  run_step "small Python tests" "small-python.log" \
    python -m unittest discover -s python/tests -p 'test_*.py'
  run_step "small VNtyper port tests" "small-vntyper.log" \
    python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'
fi

# Engine runs. For each selected input run each selected engine, then diff.
declare -a INPUTS=()
[[ $INPUT_BAM -eq 1 ]]   && INPUTS+=("bam")
[[ $INPUT_FASTQ -eq 1 ]] && INPUTS+=("fastq")

CASE_ARGS=()
[[ -n "$CASE_FILTER" ]] && CASE_ARGS=(--fixture "$CASE_FILTER")

for input in "${INPUTS[@]}"; do
  java_json="$OUT_DIR/java-$input.json"
  rust_json="$OUT_DIR/rust-$input.json"

  if [[ $RUN_JAVA -eq 1 ]]; then
    run_step "Java $input pipeline" "java-$input.log" \
      python3 "$HELPER" --engine java --input "$input" \
        "${CASE_ARGS[@]}" --json "$java_json" \
        --out-dir "$OUT_DIR/java-$input"
    print_header "Java $input output"
    show_engine_output "$java_json"
  fi

  if [[ $RUN_RUST -eq 1 ]]; then
    run_step "Rust $input pipeline" "rust-$input.log" \
      python3 "$HELPER" --engine rust --input "$input" \
        "${CASE_ARGS[@]}" --json "$rust_json" \
        --out-dir "$OUT_DIR/rust-$input"
    print_header "Rust $input output"
    show_engine_output "$rust_json"
  fi

  if [[ $RUN_JAVA -eq 1 && $RUN_RUST -eq 1 ]]; then
    run_step "Java↔Rust $input parity" "parity-$input.log" \
      python3 "$DIFFER" "$java_json" "$rust_json" \
        --label-left java --label-right rust
    print_header "Java↔Rust $input parity"
    cat "$OUT_DIR/parity-$input.log"
  fi
done

if [[ $RUN_VENDOR -eq 1 ]]; then
  KESTREL_DIR="$ROOT/vendor/rust/kestrel-rs"
  PARITY_OUT="$OUT_DIR/kestrel-vendor-parity"
  mkdir -p "$PARITY_OUT"
  run_step "kestrel-rs vendor parity gate" "vendor-kestrel.log" \
    env KESTREL_RUN_VNTYPER_FASTQ_PARITY=1 \
        KESTREL_VNTYPER_PARITY_OUT="$PARITY_OUT" \
    bash -c "cd '$KESTREL_DIR' && cargo test --release -p kestrel --test vntyper_fastq_parity -- --nocapture"
fi

print_header "Summary"
pass=0; fail=0
total_secs=0
for i in "${!STEP_LABELS[@]}"; do
  s=${STEP_STATUS[$i]}
  total_secs=$(( total_secs + STEP_SECS[$i] ))
  case "$s" in
    PASS) color=$C_GRN; pass=$((pass+1)) ;;
    *)    color=$C_RED; fail=$((fail+1)) ;;
  esac
  printf "  %s%-4s%s  %5ss  %s\n" "$color" "$s" "$C_RST" "${STEP_SECS[$i]}" "${STEP_LABELS[$i]}"
done

printf '\n  %d steps, %s%d pass%s, %s%d fail%s, %ss wall\n' \
  "${#STEP_LABELS[@]}" \
  "$C_GRN" "$pass" "$C_RST" "$C_RED" "$fail" "$C_RST" "$total_secs"
printf '  logs:  %s\n\n' "$OUT_DIR"

[[ $fail -eq 0 ]]
