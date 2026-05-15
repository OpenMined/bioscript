#!/usr/bin/env bash
# test-vntyper.sh — run VNtyper through the BioScript pipeline with either
# the Java Kestrel engine, the Rust kestrel-rs engine (via the native
# extension), or both for side-by-side comparison.
#
# Quick examples:
#   ./test-vntyper.sh --rust --fastq            # Rust native FASTQ gate
#   ./test-vntyper.sh --rust --bam              # Rust native BAM gate
#   ./test-vntyper.sh --java --bam              # external Java BAM gate
#   ./test-vntyper.sh --rust --bam --strict     # + strict TSV/report parity
#   ./test-vntyper.sh --java --rust --bam       # run both, compare
#   ./test-vntyper.sh --small                   # fast small-fixture suite only
#   ./test-vntyper.sh --all                     # everything (heavy, ~minutes)

set -euo pipefail

# Resolve repo root from script location.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

# ---- Defaults ---------------------------------------------------------------
RUN_JAVA=0
RUN_RUST=0
INPUT_BAM=0
INPUT_FASTQ=0
RUN_SMALL=0
RUN_ALL=0
RUN_STRICT=0
RUN_VENDOR=0
DO_REBUILD=0
VERBOSE=0

# ---- ANSI helpers (only if stdout is a terminal) ----------------------------
if [[ -t 1 ]]; then
  C_RED=$'\033[31m';  C_GRN=$'\033[32m';  C_YLW=$'\033[33m'
  C_BLU=$'\033[34m';  C_BLD=$'\033[1m';   C_DIM=$'\033[2m';  C_RST=$'\033[0m'
else
  C_RED=; C_GRN=; C_YLW=; C_BLU=; C_BLD=; C_DIM=; C_RST=
fi

usage() {
  cat <<EOF
${C_BLD}test-vntyper.sh${C_RST} — run VNtyper through BioScript with selectable engines.

${C_BLD}USAGE${C_RST}
  ./test-vntyper.sh [ENGINE] [INPUT] [OPTIONS]

${C_BLD}ENGINE${C_RST} (at least one)
  -j, --java        Run the Java-Kestrel external pipeline
                    (BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1; needs java + samtools + bcftools).
  -r, --rust        Run the Rust-Kestrel native pipeline
                    (BIOSCRIPT_RUN_NATIVE_*_PARITY=1; uses kestrel-rs via _native.so).
                    --java and --rust together runs both back-to-back.

${C_BLD}INPUT${C_RST} (at least one when an engine is selected)
  -b, --bam         BAM input fixture (positive + negative).
  -f, --fastq       FASTQ input fixture (Rust only — there is no Java-only FASTQ gate).

${C_BLD}MODES${C_RST}
      --small       Just the fast small-fixture suite (no opt-in gates).
      --all         Run small + every opt-in parity gate (Rust + Java + strict).
      --strict      Add BIOSCRIPT_RUN_NATIVE_BAM_OUTPUT_PARITY=1
                    (byte-equal TSV/report fingerprint for the Rust BAM path).
      --vendor      Also run the kestrel-rs vendor parity gate
                    (KESTREL_RUN_VNTYPER_FASTQ_PARITY=1, cargo test --release).

${C_BLD}OPTIONS${C_RST}
      --rebuild     Rebuild python/bioscript/_native.so via maturin before running.
  -v, --verbose     Show full unittest output (default: tail on failure only).
  -h, --help        Show this help.

${C_BLD}OUTPUT${C_RST}
  Logs land under /tmp/vntyper-run-<timestamp>/. Each step shows wall time
  and pass/fail. A summary table is printed at the end.

${C_BLD}NOTES${C_RST}
  - The opt-in gates need large fixtures under ports/vntyper/test-data/.
    They will skip with a clear message listing the missing file/tool.
  - --rust requires python/bioscript/_native.so to exist (use --rebuild
    if you changed Rust sources).
EOF
}

# ---- Arg parsing ------------------------------------------------------------
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
    --all)        RUN_ALL=1; RUN_SMALL=1; RUN_JAVA=1; RUN_RUST=1; INPUT_BAM=1; INPUT_FASTQ=1; RUN_STRICT=1; RUN_VENDOR=1 ;;
    --strict)     RUN_STRICT=1 ;;
    --vendor)     RUN_VENDOR=1 ;;
    --rebuild)    DO_REBUILD=1 ;;
    -v|--verbose) VERBOSE=1 ;;
    -h|--help)    usage; exit 0 ;;
    *) echo "${C_RED}Unknown flag: $1${C_RST}" >&2; usage >&2; exit 2 ;;
  esac
  shift
done

# Validate selection. --small alone is fine. An engine without an input is not.
if [[ $RUN_SMALL -eq 0 && $RUN_JAVA -eq 0 && $RUN_RUST -eq 0 && $RUN_VENDOR -eq 0 ]]; then
  echo "${C_RED}Pick at least one of --java / --rust / --small / --all / --vendor.${C_RST}" >&2
  exit 2
fi
if [[ ($RUN_JAVA -eq 1 || $RUN_RUST -eq 1) && $INPUT_BAM -eq 0 && $INPUT_FASTQ -eq 0 ]]; then
  echo "${C_RED}Engine selected but no input — add --bam and/or --fastq.${C_RST}" >&2
  exit 2
fi
if [[ $RUN_JAVA -eq 1 && $INPUT_FASTQ -eq 1 && $INPUT_BAM -eq 0 ]]; then
  echo "${C_YLW}Warning: there is no Java-only FASTQ gate; --java will be skipped for FASTQ.${C_RST}" >&2
fi

# ---- Setup ------------------------------------------------------------------
TS=$(date +%Y%m%d-%H%M%S)
OUT_DIR="/tmp/vntyper-run-$TS"
mkdir -p "$OUT_DIR"

export PYTHONPATH="$ROOT/python:$ROOT/ports/vntyper/bioscript${PYTHONPATH:+:$PYTHONPATH}"
export CC=${CC:-cc}
export AR=${AR:-ar}

# Tracks results. Indexed by step label.
declare -a STEP_LABELS=()
declare -a STEP_STATUS=()  # PASS / FAIL / SKIP
declare -a STEP_SECS=()
declare -a STEP_LOGS=()

print_header() {
  printf '\n%s\n' "${C_BLU}${C_BLD}━━ %s ━━${C_RST}" | sed "s/%s/$*/"
}

# Run one step: label, log-filename, command...
run_step() {
  local label="$1"; shift
  local log_name="$1"; shift
  local log="$OUT_DIR/$log_name"
  local start end secs status

  print_header "$label"
  echo "${C_DIM}\$ $*${C_RST}"
  echo "${C_DIM}log: $log${C_RST}"

  start=$(date +%s)
  set +e
  if [[ $VERBOSE -eq 1 ]]; then
    "$@" 2>&1 | tee "$log"
    local rc=${PIPESTATUS[0]}
  else
    "$@" >"$log" 2>&1
    local rc=$?
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
  STEP_LOGS+=("$log")
  return 0  # never fail the script; summary at end shows status
}

# ---- Optional rebuild -------------------------------------------------------
if [[ $DO_REBUILD -eq 1 ]]; then
  print_header "Rebuilding native extension"
  rm -rf /tmp/bioscript-maturin-wheel
  if ! command -v maturin >/dev/null 2>&1; then
    echo "${C_RED}maturin not on PATH — install with: pipx install maturin${C_RST}"
    exit 3
  fi
  maturin build --release -o /tmp/bioscript-maturin-wheel
  install -m 755 "$ROOT/rust/target/release/lib_native.so" "$ROOT/python/bioscript/_native.so"
  echo "${C_GRN}✓ _native.so updated${C_RST}"
fi

# ---- Steps ------------------------------------------------------------------
# Small / fast tests
if [[ $RUN_SMALL -eq 1 || $RUN_ALL -eq 1 ]]; then
  run_step "small Python tests" "small-python.log" \
    python -m unittest discover -s python/tests -p 'test_*.py'

  run_step "small VNtyper port tests" "small-vntyper.log" \
    python -m unittest discover -s ports/vntyper/tests -p 'test_*.py'
fi

# Rust native pipeline
if [[ $RUN_RUST -eq 1 ]]; then
  if [[ $INPUT_BAM -eq 1 ]]; then
    BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 \
      run_step "Rust BAM parity gate" "rust-bam.log" \
      env BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 \
      python -m unittest -v ports.vntyper.tests.test_native_bam_pipeline_gate
  fi
  if [[ $INPUT_FASTQ -eq 1 ]]; then
    run_step "Rust FASTQ parity gate" "rust-fastq.log" \
      env BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1 \
      python -m unittest -v ports.vntyper.tests.test_native_fastq_pipeline_gate
  fi
fi

# Java external pipeline (BAM only — there is no Java-only FASTQ gate)
if [[ $RUN_JAVA -eq 1 && $INPUT_BAM -eq 1 ]]; then
  run_step "Java external BAM parity gate" "java-bam.log" \
    env BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1 \
    python -m unittest -v ports.vntyper.tests.test_native_bam_pipeline_gate
fi

# Strict BAM TSV/report fingerprint
if [[ $RUN_STRICT -eq 1 && ($RUN_RUST -eq 1 || $RUN_ALL -eq 1) && $INPUT_BAM -eq 1 ]]; then
  run_step "Rust BAM strict output parity" "rust-bam-strict.log" \
    env BIOSCRIPT_RUN_NATIVE_BAM_OUTPUT_PARITY=1 \
        BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1 \
    python -m unittest -v \
      ports.vntyper.tests.test_native_bam_pipeline_gate.VntyperNativeBamPipelineGateTests.test_native_bam_output_fingerprints_match_expected_outputs
fi

# Vendor kestrel-rs gate
if [[ $RUN_VENDOR -eq 1 ]]; then
  KESTREL_DIR="$ROOT/vendor/rust/kestrel-rs"
  PARITY_OUT="$OUT_DIR/kestrel-vendor-parity"
  mkdir -p "$PARITY_OUT"
  run_step "kestrel-rs vendor parity gate" "vendor-kestrel.log" \
    env KESTREL_RUN_VNTYPER_FASTQ_PARITY=1 \
        KESTREL_VNTYPER_PARITY_OUT="$PARITY_OUT" \
    bash -c "cd '$KESTREL_DIR' && cargo test --release -p kestrel --test vntyper_fastq_parity -- --nocapture"
fi

# ---- Summary ----------------------------------------------------------------
print_header "Summary"
pass=0; fail=0; skip=0
total_secs=0
for i in "${!STEP_LABELS[@]}"; do
  s=${STEP_STATUS[$i]}
  total_secs=$(( total_secs + STEP_SECS[$i] ))
  case "$s" in
    PASS) color=$C_GRN; pass=$((pass+1)) ;;
    FAIL) color=$C_RED; fail=$((fail+1)) ;;
    *)    color=$C_YLW; skip=$((skip+1)) ;;
  esac
  printf "  %s%-6s%s  %5ss  %s\n" "$color" "$s" "$C_RST" "${STEP_SECS[$i]}" "${STEP_LABELS[$i]}"
done

printf '\n  total: %d steps, %s%d pass%s, %s%d fail%s, %ss wall\n' \
  "${#STEP_LABELS[@]}" \
  "$C_GRN" "$pass" "$C_RST" "$C_RED" "$fail" "$C_RST" "$total_secs"
printf '  logs:  %s\n\n' "$OUT_DIR"

# If both --java and --rust were run for BAM, suggest a diff command.
if [[ $RUN_JAVA -eq 1 && $RUN_RUST -eq 1 && $INPUT_BAM -eq 1 ]]; then
  echo "${C_DIM}Compare Java vs Rust BAM logs:"
  echo "  diff <(grep -E '^(ok|FAIL|ERROR)' $OUT_DIR/java-bam.log) \\"
  echo "       <(grep -E '^(ok|FAIL|ERROR)' $OUT_DIR/rust-bam.log)${C_RST}"
fi

[[ $fail -eq 0 ]]
