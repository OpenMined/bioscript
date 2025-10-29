"""Utility functions for working with pandas DataFrames and optional values."""

from __future__ import annotations

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
