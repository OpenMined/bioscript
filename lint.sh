#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Setup venv (must be sequential)
uv venv --quiet 2>/dev/null || true
source .venv/bin/activate
uv pip install -e ./python --quiet 2>/dev/null
uv pip install pytest ruff mypy vulture --quiet 2>/dev/null

cd "$SCRIPT_DIR/python"

echo "Running ruff format..."
ruff format .

echo "Running ruff check with fixes..."
ruff check . --fix

echo "Running mypy..."
mypy src/bioscript

echo "Running vulture to detect dead code..."
vulture src tests --min-confidence 80

echo "✓ All linting checks passed!"
