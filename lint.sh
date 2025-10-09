#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/python"

echo "Running ruff format..."
uv run ruff format .

echo "Running ruff check with fixes..."
uv run ruff check . --fix

echo "Running mypy..."
uv run mypy src/bioscript

echo "Running vulture to detect dead code..."
uv run vulture src tests --min-confidence 80

echo "✓ All linting checks passed!"
