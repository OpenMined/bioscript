"""Utilities for writing variant match results to files."""

from __future__ import annotations

import csv
from pathlib import Path


def write_csv(
    filename: str | Path,
    matches: list,
    headers: list[str] | None = None,
    delimiter: str = "\t",
):
    """Write matches to a CSV/TSV file.

    Args:
        filename: Path to output file
        matches: List of VariantMatch objects or dicts
        headers: Optional list of column headers. If not provided, infers from first match.
        delimiter: Field delimiter (default: tab for TSV)
    """
    path = Path(filename)

    if not matches:
        with path.open("w", newline="") as f:
            if headers:
                writer = csv.DictWriter(f, fieldnames=headers, delimiter=delimiter)
                writer.writeheader()
        return

    if headers is None:
        if hasattr(matches[0], "as_dict"):
            headers = list(matches[0].as_dict().keys())
        elif isinstance(matches[0], dict):
            headers = list(matches[0].keys())
        else:
            raise ValueError("Cannot infer headers - please provide them explicitly")

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter=delimiter)
        writer.writeheader()
        for match in matches:
            if hasattr(match, "as_dict"):
                writer.writerow(match.as_dict())
            elif isinstance(match, dict):
                writer.writerow(match)
            else:
                raise ValueError(f"Cannot write match of type {type(match)}")


def write_tsv(filename: str | Path, matches: list, headers: list[str] | None = None):
    """Write matches to a TSV file (convenience wrapper for write_csv with tab delimiter).

    Args:
        filename: Path to output TSV file
        matches: List of VariantMatch objects or dicts
        headers: Optional list of column headers. If not provided, infers from first match.
    """
    write_csv(filename, matches, headers=headers, delimiter="\t")
