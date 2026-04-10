#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

cd "$ROOT_DIR/rust"

if [[ "${1:-}" == "--security" ]]; then
  exec cargo test -p bioscript --test security -- --nocapture
fi

exec cargo test -p bioscript -- --nocapture
