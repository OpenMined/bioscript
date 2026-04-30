#!/usr/bin/env bash
set -euo pipefail

# Coverage runner for the BioScript Rust crates using cargo-llvm-cov.
# - Mirrors test.sh by operating inside ./rust
# - Runs focused BioScript tests by default, with optional all-tests mode
# - Generates HTML report and LCOV file
# - Prints a sorted summary to stdout

FULL_CLEAN_FLAG=0
OPEN_HTML_FLAG=${OPEN_HTML:-0}
LARGE_FLAG=0
ALL_TESTS_FLAG=0
NO_LINT_FLAG=0
FOCUSED_TEST=""

usage() {
  cat <<'EOF'
Usage: ./coverage.sh [--full-clean|-c] [--open] [--large] [--all-tests] [--no-lint] [--focused-test name]

  --full-clean, -c  Run cargo clean and remove coverage dirs before running
  --open            Open HTML report locally (no-op in CI)
  --large           Include tests that require large local fixtures
  --all-tests       Run all tests for the first-party BioScript crates
  --no-lint         Skip cargo fmt and clippy checks
  --focused-test    Run one focused integration test target:
                    file_formats, formats_lib, inspect, prepare, cli, schema, core,
                    runtime_security, or runtime_resources

Environment:
  AUTO_INSTALL_LLVM_COV=0    Do not auto-install cargo-llvm-cov
  AUTO_INSTALL_LLVM_TOOLS=0  Do not auto-install llvm-tools-preview
  FULL_CLEAN=1               Same as --full-clean
  LCOV_OUT=path              Override LCOV output path
  OPEN_HTML=1                Same as --open
  TEST_THREADS=n             Test thread count for focused integration tests
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --full-clean|-c)
      FULL_CLEAN_FLAG=1
      ;;
    --open)
      OPEN_HTML_FLAG=1
      ;;
    --large)
      LARGE_FLAG=1
      ;;
    --all-tests)
      ALL_TESTS_FLAG=1
      ;;
    --no-lint)
      NO_LINT_FLAG=1
      ;;
    --focused-test)
      if [[ $# -lt 2 ]]; then
        echo "--focused-test requires a test target name" >&2
        usage >&2
        exit 2
      fi
      FOCUSED_TEST="$2"
      shift
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
  shift
done

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$ROOT_DIR/rust"

PACKAGES=(
  bioscript-cli
  bioscript-core
  bioscript-formats
  bioscript-runtime
  bioscript-schema
)

PKG_ARGS=()
for package in "${PACKAGES[@]}"; do
  PKG_ARGS+=(-p "$package")
done

TEST_RUSTFLAGS="${RUSTFLAGS:-} -Aunused-assignments -Amissing-docs"

if [[ "$NO_LINT_FLAG" != "1" ]]; then
  echo "==> Formatting and linting"
  cargo fmt --check "${PKG_ARGS[@]}"
  RUSTFLAGS="$TEST_RUSTFLAGS" cargo clippy "${PKG_ARGS[@]}" --all-targets --color=never -- -D warnings
fi

echo "==> Checking cargo-llvm-cov availability"
if ! cargo llvm-cov --version >/dev/null 2>&1; then
  if [[ "${AUTO_INSTALL_LLVM_COV:-1}" == "1" ]]; then
    echo "==> Installing cargo-llvm-cov (first run only)"
    if ! cargo install cargo-llvm-cov; then
      echo "Failed to install cargo-llvm-cov. Install manually with:" >&2
      echo "  cargo install cargo-llvm-cov" >&2
      exit 1
    fi
  else
    echo "cargo-llvm-cov is not installed. Install with:" >&2
    echo "  cargo install cargo-llvm-cov" >&2
    exit 1
  fi
fi

if ! rustup component list --installed | grep -Eq '^(llvm-tools-preview|llvm-tools)'; then
  if [[ "${AUTO_INSTALL_LLVM_TOOLS:-1}" == "1" ]]; then
    echo "==> Installing rustup component: llvm-tools-preview (first run only)"
    if ! rustup component add llvm-tools-preview; then
      echo "Failed to install llvm-tools-preview. Install manually with:" >&2
      echo "  rustup component add llvm-tools-preview" >&2
      exit 1
    fi
  else
    echo "llvm-tools-preview is missing. Enable auto-install or run:" >&2
    echo "  rustup component add llvm-tools-preview" >&2
    exit 1
  fi
fi

if [[ "${FULL_CLEAN:-0}" == "1" || "$FULL_CLEAN_FLAG" == "1" ]]; then
  echo "==> FULL_CLEAN=1: performing cargo clean and removing coverage dirs"
  cargo clean
  rm -rf target/llvm-cov target/coverage target/llvm-cov-target
fi

echo "==> Cleaning previous coverage artifacts"
cargo llvm-cov clean --workspace
mkdir -p target/coverage

LCOV_OUT=${LCOV_OUT:-target/coverage/lcov.info}
TEST_THREADS=${TEST_THREADS:-1}
IGNORE_REGEX=${IGNORE_REGEX:-'(^|/)(monty|noodles|vendor)/'}
OPEN_FLAG=""
if [[ "$OPEN_HTML_FLAG" == "1" && "${CI:-0}" != "true" ]]; then
  OPEN_FLAG="--open"
fi

COV_ENV=(RUSTFLAGS="$TEST_RUSTFLAGS")
if [[ "$LARGE_FLAG" == "1" ]]; then
  COV_ENV+=(BIOSCRIPT_RUN_LARGE_TESTS=1)
fi

echo "==> Running coverage"
if [[ -n "$FOCUSED_TEST" ]]; then
  case "$FOCUSED_TEST" in
    file_formats)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --test file_formats -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    formats_lib)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --lib
      ;;
    inspect)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --test inspect -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    prepare)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --test prepare -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    cli)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-cli --test cli -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    schema)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-schema --test validate_variants -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    core)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-core --lib
      ;;
    runtime_security)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-runtime --test security -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    runtime_resources)
      env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-runtime --test resources_coverage -- --nocapture --test-threads="$TEST_THREADS"
      ;;
    *)
      echo "Unknown focused test target: $FOCUSED_TEST" >&2
      usage >&2
      exit 2
      ;;
  esac
elif [[ "$ALL_TESTS_FLAG" == "1" ]]; then
  env "${COV_ENV[@]}" cargo llvm-cov --no-report "${PKG_ARGS[@]}" --all-targets
else
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --test file_formats -- --nocapture --test-threads="$TEST_THREADS"
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --lib
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --test inspect -- --nocapture --test-threads="$TEST_THREADS"
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-formats --test prepare -- --nocapture --test-threads="$TEST_THREADS"
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-cli --test cli -- --nocapture --test-threads="$TEST_THREADS"
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-schema --test validate_variants -- --nocapture --test-threads="$TEST_THREADS"
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-core --lib
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-runtime --test security -- --nocapture --test-threads="$TEST_THREADS"
  env "${COV_ENV[@]}" cargo llvm-cov --no-report -p bioscript-runtime --test resources_coverage -- --nocapture --test-threads="$TEST_THREADS"
fi

echo "==> Generating HTML report"
cargo llvm-cov report "${PKG_ARGS[@]}" --html --ignore-filename-regex "$IGNORE_REGEX" $OPEN_FLAG

echo "==> Exporting LCOV"
cargo llvm-cov report "${PKG_ARGS[@]}" --ignore-filename-regex "$IGNORE_REGEX" --lcov --output-path "$LCOV_OUT"

echo "==> Coverage summary (sorted by coverage %)"
SUMMARY_OUTPUT=$(cargo llvm-cov report "${PKG_ARGS[@]}" --ignore-filename-regex "$IGNORE_REGEX" --summary-only)
printf '%s\n' "$SUMMARY_OUTPUT" | head -n 3
printf '%s\n' "$SUMMARY_OUTPUT" | tail -n +4 | grep -v "^TOTAL" | sort -t'%' -k3 -n
printf '%s\n' "$SUMMARY_OUTPUT" | grep "^TOTAL"

HTML_DIR="target/llvm-cov/html"
if [[ -d "$HTML_DIR" ]]; then
  echo "HTML report: rust/$HTML_DIR/index.html"
else
  echo "HTML report directory not found. cargo-llvm-cov typically writes to target/llvm-cov/html" >&2
fi

echo "LCOV file: rust/$LCOV_OUT"
