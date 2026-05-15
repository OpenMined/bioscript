#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

VERSION="${1:-$(tr -d '[:space:]' < VERSION)}"
IMAGE="${BIOSCRIPT_IMAGE:-ghcr.io/openmined/bioscript}"
PLATFORMS="${PLATFORMS:-linux/amd64,linux/arm64}"
PUSH="${PUSH:-0}"
LOAD="${LOAD:-0}"

BUILD_ARGS=(
  docker buildx build
  --platform "$PLATFORMS"
  -f docker/Dockerfile
  -t "${IMAGE}:${VERSION}"
  -t "${IMAGE}:latest"
  .
)

if [[ "$PUSH" == "1" ]]; then
  BUILD_ARGS+=(--push)
elif [[ "$LOAD" == "1" ]]; then
  if [[ "$PLATFORMS" == *,* ]]; then
    echo "LOAD=1 requires a single platform, got: $PLATFORMS" >&2
    exit 1
  fi
  BUILD_ARGS+=(--load)
else
  BUILD_ARGS+=(--output "type=oci,dest=docker/bioscript-${VERSION}.oci.tar")
fi

"${BUILD_ARGS[@]}"
