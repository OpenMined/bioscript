#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

set -e

export UV_VENV_CLEAR=1
VENV_PATH="$SCRIPT_DIR/.venv"
export UV_PROJECT_ENVIRONMENT="$VENV_PATH"

uv venv "$VENV_PATH"
uv pip install -e ./python
uv pip install pytest ruff mypy vulture

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
