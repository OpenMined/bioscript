#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <version>" >&2
  echo "Example: $0 0.2.0" >&2
  exit 1
fi

VERSION="$1"
if [[ ! "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+(-[0-9A-Za-z.-]+)?$ ]]; then
  echo "Invalid semver version: $VERSION" >&2
  exit 1
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

printf "%s\n" "$VERSION" > VERSION

for manifest in \
  rust/bioscript-cli/Cargo.toml \
  rust/bioscript-core/Cargo.toml \
  rust/bioscript-ffi/Cargo.toml \
  rust/bioscript-formats/Cargo.toml \
  rust/bioscript-runtime/Cargo.toml \
  rust/bioscript-schema/Cargo.toml \
  rust/bioscript-wasm/Cargo.toml
do
  perl -0pi -e 's/^version = ".*?"$/version = "'"$VERSION"'"/m' "$manifest"
done

cargo metadata --manifest-path rust/Cargo.toml --format-version 1 >/dev/null

echo "Bumped BioScript to $VERSION"
