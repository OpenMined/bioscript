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

FMT_ARGS=()
for package in "${PACKAGES[@]}"; do
  FMT_ARGS+=(-p "$package")
done

cargo fmt --check "${FMT_ARGS[@]}"
cargo clippy --workspace --all-targets -- -D warnings
