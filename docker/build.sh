#!/bin/bash
set -e

# Build BioScript Docker image
#
# Usage: ./build.sh [version]
#
# If no version is provided, reads from __init__.py

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/.."

# Get version from argument or __init__.py
if [ -n "$1" ]; then
    VERSION="$1"
else
    VERSION=$(grep '^__version__ = ' python/src/bioscript/__init__.py | cut -d'"' -f2)
fi

echo "Building bioscript:${VERSION}..."

docker build \
    -f docker/Dockerfile \
    -t "bioscript:${VERSION}" \
    -t "bioscript:latest" \
    -t "ghcr.io/openmined/bioscript:${VERSION}" \
    -t "ghcr.io/openmined/bioscript:latest" \
    .

echo "✓ Built bioscript:${VERSION}"
echo "  Tagged: bioscript:${VERSION}, bioscript:latest"
echo "  Tagged: ghcr.io/openmined/bioscript:${VERSION}, ghcr.io/openmined/bioscript:latest"
echo ""
echo "Test with:"
echo "  docker run --rm bioscript:${VERSION}"
echo "  docker run --rm ghcr.io/openmined/bioscript:latest"
