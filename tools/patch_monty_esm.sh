#!/usr/bin/env bash
# Post-build patch for the monty wasi-browser ESM loader.
#
# `pydantic/monty` (vendored as a submodule at bioscript/monty) emits a browser
# loader via napi-rs that calls `instantiateNapiModuleSync`, which triggers
# `new WebAssembly.Module(bytes)` on the main thread. Chromium refuses that for
# WASM binaries >8 MB, and the monty runtime is well over that.
#
# The fix is purely cosmetic to the generated file (swap the sync init import
# for the async one the same runtime already exports). We keep this script in
# bioscript/ — not inside the monty submodule — so we don't have to fork monty.
# Run it every time we rebuild or update monty.
#
# Usage:
#   bash bioscript/tools/patch_monty_esm.sh <path-to-monty.wasi-browser.mjs>
#
# If no path is given, we patch the file next to us in the monty submodule
# (best-effort location produced by `napi build --esm`).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_TARGET="${SCRIPT_DIR}/../monty/crates/monty-js/monty.wasi-browser.mjs"
TARGET="${1:-${DEFAULT_TARGET}}"

if [ ! -f "${TARGET}" ]; then
  echo "[patch_monty_esm] target not found: ${TARGET}" >&2
  echo "[patch_monty_esm] skipping (build monty-js first or pass a path)" >&2
  exit 0
fi

if ! grep -q 'instantiateNapiModuleSync' "${TARGET}"; then
  echo "[patch_monty_esm] already patched; leaving ${TARGET} alone"
  exit 0
fi

# Portable sed -i (BSD vs GNU).
if sed --version >/dev/null 2>&1; then
  SED_INPLACE=(sed -i)
else
  SED_INPLACE=(sed -i '')
fi

"${SED_INPLACE[@]}" \
  -e 's/instantiateNapiModuleSync as __emnapiInstantiateNapiModuleSync/instantiateNapiModule as __emnapiInstantiateNapiModule/' \
  -e 's/__emnapiInstantiateNapiModuleSync(/await __emnapiInstantiateNapiModule(/' \
  "${TARGET}"

# Sanity check: make sure the patch actually applied both substitutions.
if grep -q 'instantiateNapiModuleSync' "${TARGET}"; then
  echo "[patch_monty_esm] patch failed — sync import still present in ${TARGET}" >&2
  exit 1
fi

echo "[patch_monty_esm] rewrote sync → async init in ${TARGET}"
