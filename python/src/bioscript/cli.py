"""BioScript CLI for running genetic variant classifiers."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import sys
import traceback
from pathlib import Path

from .reader import load_variants_tsv
from .testing import export_from_notebook, run_tests
from .types import MatchList


def load_classifier_module(script_path: Path):
    """
    Dynamically load a classifier script.

    Expected export: __bioscript__ dictionary with:
        variant_calls: List of VariantCall objects (required for auto mode)
        classifier: Callable that takes matches and returns string (required for auto mode)
        name (optional): Column name for output (defaults to script filename)
        main (optional): Custom function(*args, **kwargs) -> dict for full control

    Args:
        script_path: Path to the classifier script

    Returns:
        Dictionary with 'name' and either 'handler' or 'config'
    """
    spec = importlib.util.spec_from_file_location(script_path.stem, script_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {script_path}")

    module = importlib.util.module_from_spec(spec)
    sys.modules[script_path.stem] = module
    spec.loader.exec_module(module)

    # Load __bioscript__ dict
    if not hasattr(module, "__bioscript__"):
        raise AttributeError(
            f"{script_path} must export '__bioscript__' dict with "
            "'variant_calls' and 'classifier' (or 'main')"
        )

    config = module.__bioscript__

    # Check if custom main function provided
    if "main" in config:
        return {
            "name": config.get("name", script_path.stem),
            "main": config["main"],
        }

    # Auto mode - require variant_calls and classifier
    if "variant_calls" not in config:
        raise AttributeError(f"{script_path}: __bioscript__ must include 'variant_calls'")
    if "classifier" not in config:
        raise AttributeError(f"{script_path}: __bioscript__ must include 'classifier'")

    return {
        "name": config.get("name", script_path.stem),
        "variant_calls": config["variant_calls"],
        "classifier": config["classifier"],
    }


def test_command(args):
    """Run tests in classifier modules."""
    all_passed = True

    for script_path_str in args.classifiers:
        script_path = Path(script_path_str)
        if not script_path.exists():
            print(f"Error: Classifier script not found: {script_path}", file=sys.stderr)
            sys.exit(1)

        print(f"\n{'='*60}")
        print(f"Testing: {script_path}")
        print('='*60)

        result = run_tests(script_path, verbose=True)

        if not result["success"]:
            all_passed = False

    # Exit with error code if any tests failed
    if not all_passed:
        sys.exit(1)


def export_command(args):
    """Export classifier from Jupyter notebook."""
    notebook_path = Path(args.notebook)
    if not notebook_path.exists():
        print(f"Error: Notebook not found: {args.notebook}", file=sys.stderr)
        sys.exit(1)

    output_path = args.output if args.output else None

    try:
        result = export_from_notebook(
            notebook_path,
            output_path=output_path,
            include_tests=not args.no_tests,
        )
        print(f"✓ Exported to: {result}")

        # Run tests if requested
        if args.test and not args.no_tests:
            print(f"\nRunning tests in exported file...")
            test_result = run_tests(result, verbose=True)
            if not test_result["success"]:
                sys.exit(1)

    except Exception as e:
        print(f"Error exporting notebook: {e}", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(1)


def classify_command(args):
    """Run classification on SNP file with multiple classifiers."""
    # Load SNP file
    snp_file_path = Path(args.file)
    if not snp_file_path.exists():
        print(f"Error: File not found: {args.file}", file=sys.stderr)
        sys.exit(1)

    # Results dictionary - only add participant_id if provided
    results = {}
    if args.participant_id:
        results[args.participant_col] = args.participant_id

    # Process each classifier
    for script_path_str in args.classifiers:
        script_path = Path(script_path_str)
        if not script_path.exists():
            print(f"Error: Classifier script not found: {script_path}", file=sys.stderr)
            sys.exit(1)

        try:
            # Load classifier module
            module_config = load_classifier_module(script_path)
            name = module_config["name"]

            # Check if custom main function
            if "main" in module_config:
                # Call main with full control
                main_func = module_config["main"]
                try:
                    # Build kwargs - only include participant_id if provided
                    main_kwargs = {
                        "snp_file": str(snp_file_path),
                        "file": str(snp_file_path),
                    }
                    if args.participant_id:
                        main_kwargs["participant_id"] = args.participant_id

                    result = main_func(**main_kwargs)

                    # Handle different return types
                    if isinstance(result, dict):
                        results.update(result)
                    elif isinstance(result, str):
                        results[name] = result
                    elif isinstance(result, Path):
                        # File output - verify exists
                        if not result.exists():
                            raise FileNotFoundError(
                                f"main() returned path {result} but file does not exist"
                            )
                        results[name] = str(result)
                    else:
                        results[name] = str(result)

                except Exception as e:
                    print(
                        f"Error in {script_path} main(): {e}",
                        file=sys.stderr,
                    )
                    print(traceback.format_exc(), file=sys.stderr)
                    results[name] = "ERROR"

            else:
                # Auto mode - load variants and classify
                variant_calls = module_config["variant_calls"]
                classifier = module_config["classifier"]

                try:
                    # Load and match variants
                    variants = load_variants_tsv(snp_file_path)
                    calls = MatchList(variant_calls=variant_calls)
                    matches = calls.match_rows(variants)

                    # Call classifier (uses __call__ interface)
                    result = classifier(matches)
                    results[name] = result

                except Exception as e:
                    print(
                        f"Error in {script_path} classification: {e}",
                        file=sys.stderr,
                    )
                    print(traceback.format_exc(), file=sys.stderr)
                    results[name] = "ERROR"

        except (ImportError, AttributeError) as e:
            print(f"Error loading {script_path}: {e}", file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            sys.exit(1)

    # Output based on format
    try:
        if args.out == "tsv":
            writer = csv.DictWriter(sys.stdout, fieldnames=results.keys(), delimiter="\t")
            writer.writeheader()
            writer.writerow(results)
        elif args.out == "json":
            import json

            print(json.dumps(results, indent=2))
        else:
            # Simple key=value output
            for key, value in results.items():
                print(f"{key}={value}")
    except Exception as e:
        print(f"Error writing output: {e}", file=sys.stderr)
        print(traceback.format_exc(), file=sys.stderr)
        sys.exit(1)


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="BioScript - Genetic variant classification tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run single classifier
  bioscript classify --file=snps.txt classify_apol1.py

  # With participant ID
  bioscript classify --participant_id=P001 --file=snps.txt classify_apol1.py

  # Chain multiple classifiers with TSV output
  bioscript classify --participant_id=P001 --file=snps.txt \\
    classify_apol1.py classify_apol2.py --out=tsv

  # Custom participant column name
  bioscript classify --participant_id=P001 --file=snps.txt \\
    classify_apol1.py --participant_col=sample_id --out=tsv
        """,
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Test command
    test_parser = subparsers.add_parser(
        "test", help="Run tests in classifier modules"
    )
    test_parser.add_argument(
        "classifiers",
        nargs="+",
        help="Paths to classifier scripts with test_* functions",
    )

    # Export command
    export_parser = subparsers.add_parser(
        "export", help="Export classifier from Jupyter notebook"
    )
    export_parser.add_argument(
        "notebook",
        help="Path to Jupyter notebook (.ipynb)",
    )
    export_parser.add_argument(
        "-o", "--output",
        help="Output path for Python file (default: same name as notebook)",
    )
    export_parser.add_argument(
        "--no-tests",
        action="store_true",
        help="Exclude test functions from export",
    )
    export_parser.add_argument(
        "--test",
        action="store_true",
        help="Run tests after export",
    )

    # Classify command
    classify_parser = subparsers.add_parser(
        "classify", help="Run variant classification on SNP file"
    )
    classify_parser.add_argument(
        "--participant_id",
        help="Optional participant ID for output column",
    )
    classify_parser.add_argument(
        "--file", required=True, help="Path to SNP genotype file (TSV format)"
    )
    classify_parser.add_argument(
        "classifiers",
        nargs="+",
        help="Paths to classifier scripts",
    )
    classify_parser.add_argument(
        "--out",
        choices=["tsv", "json", "simple"],
        default="simple",
        help="Output format (default: simple)",
    )
    classify_parser.add_argument(
        "--participant_col",
        default="participant_id",
        help="Column name for participant ID in output (default: participant_id)",
    )

    args = parser.parse_args()

    # Route to command handler
    if args.command == "test":
        test_command(args)
    elif args.command == "export":
        export_command(args)
    elif args.command == "classify":
        classify_command(args)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
