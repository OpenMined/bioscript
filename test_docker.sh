#!/bin/bash
set -e

# Test bioscript using Docker
# Rebuilds the container and runs the classifier

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Building Docker image..." >&2
cd "$REPO_ROOT/docker"
bash build.sh >&2

echo "" >&2
echo "Running classification in Docker..." >&2
echo "" >&2

# Run with Docker
docker run --rm \
    -v "$REPO_ROOT/examples/apol1:/data" \
    -w /data \
    bioscript:latest \
    bash process_samplesheet.sh samplesheet.csv classify_apol1.py
