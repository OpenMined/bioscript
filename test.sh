#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPORT=0
LARGE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --security)
      cd "$ROOT_DIR/rust"
      exec cargo test -p bioscript-runtime --test security -- --nocapture
      ;;
    --report)
      REPORT=1
      shift
      ;;
    --large)
      LARGE=1
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

cd "$ROOT_DIR/rust"

TEST_RUSTFLAGS="${RUSTFLAGS:-} -Aunused-assignments -Amissing-docs"

cargo_test() {
  RUSTFLAGS="$TEST_RUSTFLAGS" cargo test "$@"
}

if [[ "$LARGE" == "1" ]]; then
  BIOSCRIPT_RUN_LARGE_TESTS=1 cargo_test -p bioscript-formats --test file_formats --test inspect -- --nocapture
else
  cargo_test -p bioscript-formats --test file_formats --test inspect -- --nocapture
fi
cargo_test -p bioscript-cli --test cli -- --nocapture

if [[ "$REPORT" == "1" ]]; then
  cargo build -p bioscript-cli
  cd "$ROOT_DIR"
  REPORT_PATH="$ROOT_DIR/test-reports/file-formats.html"
  python3 "$ROOT_DIR/tools/generate_file_format_report.py" \
    --bioscript "$ROOT_DIR/rust/target/debug/bioscript" \
    --root "$ROOT_DIR" \
    --output "$REPORT_PATH"
  echo "HTML report: $REPORT_PATH"
fi
