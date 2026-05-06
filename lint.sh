#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

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

cargo fmt --check "${PKG_ARGS[@]}"

filter_vendored() {
  awk 'BEGIN{RS=""; ORS="\n\n"} !/\/noodles\/|\/vendor\/|`noodles-|`lexical-/'
}

CLIPPY_STDERR="$(mktemp)"
FILTERED_STDERR="$(mktemp)"
cleanup() {
  rm -f "$CLIPPY_STDERR" "$FILTERED_STDERR"
}
trap cleanup EXIT

set +e
cargo clippy "${PKG_ARGS[@]}" --all-targets --color=never -- -D warnings 2> "$CLIPPY_STDERR"
CLIPPY_STATUS=$?
set -e

filter_vendored < "$CLIPPY_STDERR" > "$FILTERED_STDERR"
if [[ -s "$FILTERED_STDERR" ]]; then
  cat "$FILTERED_STDERR" >&2
elif [[ "$CLIPPY_STATUS" -ne 0 ]]; then
  cat "$CLIPPY_STDERR" >&2
fi

if [[ "$CLIPPY_STATUS" -ne 0 ]]; then
  exit "$CLIPPY_STATUS"
fi

cargo test -p bioscript-core --test source_size -- --nocapture
