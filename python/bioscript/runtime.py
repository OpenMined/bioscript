"""Backend selection for Python-side BioScript shims."""

from __future__ import annotations

import os
from enum import Enum


class BackendMode(str, Enum):
    AUTO = "auto"
    PYTHON = "python"
    RUST = "rust"


def selected_backend() -> BackendMode:
    raw = os.environ.get("BIOSCRIPT_BACKEND", BackendMode.AUTO.value).strip().lower()
    try:
        return BackendMode(raw)
    except ValueError as exc:
        allowed = ", ".join(mode.value for mode in BackendMode)
        raise ValueError(f"BIOSCRIPT_BACKEND must be one of: {allowed}") from exc
