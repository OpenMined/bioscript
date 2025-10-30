"""Utility helpers shared across BioScript classifiers."""

from __future__ import annotations

import inspect
from pathlib import Path

import pandas as pd


def optional_str(value, upper: bool = False) -> str | None:
    """Convert pandas value to optional string, handling NaN values.

    Args:
        value: Pandas value that may be NaN
        upper: If True, convert result to uppercase

    Returns:
        None if value is NaN, otherwise stripped string (optionally uppercased)
    """
    if pd.isna(value):
        return None
    result = str(value).strip()
    return result.upper() if upper else result


def optional_int(value) -> int | None:
    """Convert pandas value to optional int, handling NaN values.

    Args:
        value: Pandas value that may be NaN

    Returns:
        None if value is NaN, otherwise int
    """
    return None if pd.isna(value) else int(value)


def assets_dir(module_file: str | Path | None = None) -> Path:
    """Return the directory containing the calling classifier module.

    When ``module_file`` is provided (e.g. ``__file__``) the path is resolved
    relative to it. Otherwise the caller's module file is inferred from the
    stack, falling back to the current working directory.
    """

    if module_file is not None:
        return Path(module_file).resolve().parent

    frame_info = inspect.stack()[1]
    module = inspect.getmodule(frame_info.frame)
    if module and hasattr(module, "__file__"):
        return Path(module.__file__).resolve().parent

    return Path.cwd()
