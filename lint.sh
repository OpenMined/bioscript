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

cargo clippy "${PKG_ARGS[@]}" --all-targets --color=never -- -D warnings 2> >(filter_vendored >&2)
