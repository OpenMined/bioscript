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

CLIPPY_OUTPUT="$(mktemp)"
FILTERED_OUTPUT="$(mktemp)"
cleanup() {
  rm -f "$CLIPPY_OUTPUT" "$FILTERED_OUTPUT"
}
trap cleanup EXIT

set +e
cargo clippy "${PKG_ARGS[@]}" --all-targets --color=never -- -D warnings > "$CLIPPY_OUTPUT" 2>&1
CLIPPY_STATUS=$?
set -e

filter_vendored < "$CLIPPY_OUTPUT" > "$FILTERED_OUTPUT"
if [[ "$CLIPPY_STATUS" -ne 0 ]]; then
  cat "$CLIPPY_OUTPUT" >&2
elif [[ -s "$FILTERED_OUTPUT" ]]; then
  cat "$FILTERED_OUTPUT" >&2
fi

if [[ "$CLIPPY_STATUS" -ne 0 ]]; then
  exit "$CLIPPY_STATUS"
fi

cargo test -p bioscript-core --test source_size -- --nocapture
