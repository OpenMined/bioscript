"""Test utilities for BioScript classifiers."""

from __future__ import annotations

import importlib.util
import inspect
import sys
from pathlib import Path
from typing import Any, Callable, Iterator

from .data import GenotypeGenerator
from .types import GRCh, VariantRow


class VariantFixture:
    """Helper for creating test variant data with cleaner syntax."""

    def __init__(
        self,
        variant_templates: list[dict[str, str | int | float]],
        assembly: GRCh | str | None = None,
    ):
        """Initialize fixture with variant templates."""
        self.gen = GenotypeGenerator(variant_templates, assembly=assembly)
        self.templates = variant_templates

    def variants(self, genotypes: list[str]) -> Iterator[VariantRow]:
        """Generate variants with given genotypes."""
        # Get the base variants from the generator
        base_variants = self.gen(genotypes)

        # Add raw_line to each variant
        for variant, genotype, template in zip(base_variants, genotypes, self.templates):
            # Construct raw_line in TSV format
            variant.raw_line = (
                f"{template['rsid']}\t{template['chromosome']}\t{template['position']}\t{genotype}"
            )
            yield variant

    def __call__(self, genotypes: list[str]) -> Iterator[VariantRow]:
        """Shorthand for variants()."""
        return self.variants(genotypes)


def discover_tests(module_path: Path | str) -> dict[str, Callable]:
    """
    Discover test functions in a classifier module.

    Looks for functions starting with 'test_' and returns them as a dict.

    Args:
        module_path: Path to the classifier module

    Returns:
        Dictionary mapping test names to test functions
    """
    module_path = Path(module_path)
    if not module_path.exists():
        raise FileNotFoundError(f"Module not found: {module_path}")

    # Load module
    spec = importlib.util.spec_from_file_location(module_path.stem, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {module_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_path.stem] = module
    spec.loader.exec_module(module)

    # Find all test_* functions
    tests = {}
    for name, obj in inspect.getmembers(module):
        if name.startswith("test_") and callable(obj):
            tests[name] = obj

    return tests


def run_tests(module_path: Path | str, verbose: bool = True) -> dict[str, Any]:
    """
    Run all tests in a classifier module.

    First tries to use pytest if available, otherwise runs tests directly.

    Args:
        module_path: Path to the classifier module
        verbose: Print test results (default: True)

    Returns:
        Dictionary with test results:
            - passed: List of test names that passed
            - failed: List of tuples (test_name, error_message)
            - total: Total number of tests
            - success: True if all tests passed
    """
    module_path = Path(module_path)

    # Try pytest first
    try:
        import pytest

        if verbose:
            print(f"Running tests with pytest: {module_path}")

        # Run pytest with capture disabled for better output
        result = pytest.main(["-v", str(module_path)])

        # pytest returns 0 for success, non-zero for failures
        return {
            "success": result == 0,
            "pytest": True,
            "exit_code": result,
        }

    except ImportError:
        # Fallback: run tests directly
        if verbose:
            print(f"Running tests (pytest not available): {module_path}")

        tests = discover_tests(module_path)
        passed = []
        failed = []

        for test_name, test_func in tests.items():
            try:
                if verbose:
                    print(f"  {test_name} ... ", end="", flush=True)
                test_func()
                passed.append(test_name)
                if verbose:
                    print("✓")
            except AssertionError as e:
                failed.append((test_name, str(e)))
                if verbose:
                    print(f"✗\n    {e}")
            except Exception as e:
                failed.append((test_name, f"ERROR: {e}"))
                if verbose:
                    print(f"✗\n    ERROR: {e}")

        if verbose:
            print(f"\n{len(passed)} passed, {len(failed)} failed")

        return {
            "passed": passed,
            "failed": failed,
            "total": len(tests),
            "success": len(failed) == 0,
            "pytest": False,
        }


def _is_definition_or_assignment(code: str) -> bool:
    """
    Check if code is a definition or assignment (not a function call).

    Keeps:
    - import/from statements
    - class definitions
    - function definitions
    - variable assignments (x = ...)
    - __bioscript__ dictionary

    Skips:
    - Standalone function calls
    - Expression statements
    - Export-related code
    - Jupyter magic commands
    """
    import ast

    # Skip export-related imports
    export_keywords = ["export_from_notebook", "run_tests", "discover_tests"]
    if any(keyword in code for keyword in export_keywords):
        return False

    # Check each line for Jupyter magic commands
    for line in code.split("\n"):
        stripped = line.strip()
        if stripped.startswith(("!", "%", "%%")):
            return False

    try:
        tree = ast.parse(code)
    except SyntaxError:
        # If we can't parse it, it's probably invalid Python (like magic commands)
        # Skip it to be safe
        return False

    # Check each statement in the code
    for node in tree.body:
        # Keep imports, class/function definitions, assignments
        if isinstance(
            node,
            (
                ast.Import,
                ast.ImportFrom,
                ast.ClassDef,
                ast.FunctionDef,
                ast.AsyncFunctionDef,
                ast.Assign,
                ast.AnnAssign,
                ast.AugAssign,
            ),
        ):
            continue
        # Skip standalone expressions (function calls, etc.)
        elif isinstance(node, ast.Expr):
            return False
        # Skip other statement types (for, while, if, etc. at top level)
        else:
            return False

    return True


def export_from_notebook(
    notebook_path: Path | str,
    output_path: Path | str | None = None,
    include_tests: bool = True,
) -> Path:
    """
    Export classifier code from a Jupyter notebook to a Python module.

    Extracts only definitions and assignments from code cells:
    - Imports (except export-related)
    - Class definitions
    - Function definitions
    - Variable assignments

    Skips:
    - Function calls (like test_g0_homozygous())
    - Export-related code (export_from_notebook, run_tests)
    - Expression statements

    Args:
        notebook_path: Path to the Jupyter notebook (.ipynb)
        output_path: Output path for Python file (default: same name as notebook)
        include_tests: Include test_* functions in export (default: True)

    Returns:
        Path to the exported Python file

    Example:
        >>> export_from_notebook("apol1_dev.ipynb", "classify_apol1.py")
    """
    import json

    notebook_path = Path(notebook_path)
    if not notebook_path.exists():
        raise FileNotFoundError(f"Notebook not found: {notebook_path}")

    output_path = notebook_path.with_suffix(".py") if output_path is None else Path(output_path)

    # Load notebook
    with open(notebook_path) as f:
        notebook = json.load(f)

    # Extract code cells
    code_lines = []
    for cell in notebook.get("cells", []):
        if cell.get("cell_type") == "code":
            source = cell.get("source", [])
            # Join source lines and add to output
            cell_code = "".join(source)

            # Skip cells that are just comments or empty
            cell_code = cell_code.strip()
            if not cell_code:
                continue

            # Skip Jupyter magic commands (!, %, %%)
            if cell_code.startswith(("!", "%", "%%")):
                continue

            # Skip pure comment cells (but keep code with inline comments)
            lines = cell_code.split("\n")
            has_code = any(line.strip() and not line.strip().startswith("#") for line in lines)
            if not has_code:
                continue

            # Filter out test functions if requested
            if not include_tests and "def test_" in cell_code:
                continue

            # Only keep definitions and assignments (skip function calls)
            if not _is_definition_or_assignment(cell_code):
                continue

            code_lines.append(cell_code)
            code_lines.append("")  # Blank line between cells

    # Write to output file
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write("\n".join(code_lines))

    return output_path


# Convenience exports for Jupyter notebooks
__all__ = [
    "VariantFixture",
    "discover_tests",
    "run_tests",
    "export_from_notebook",
]
