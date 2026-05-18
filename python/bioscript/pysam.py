"""BioScript-supported pysam compatibility subset."""

from __future__ import annotations

import importlib
from typing import Any

from .runtime import BackendMode, ModuleBackendPolicy, selected_backend

BACKEND_POLICY = ModuleBackendPolicy(
    auto="use real pysam when installed; otherwise native pysam shim is pending",
    python="requires real pysam",
    rust="native pysam shim is pending",
)


def _real_pysam() -> Any:
    return importlib.import_module("pysam")


class AlignmentFile:
    """Proxy for the supported `pysam.AlignmentFile` subset."""

    def __init__(self, path: str, mode: str = "r", **kwargs: Any) -> None:
        backend = selected_backend()
        if backend in {BackendMode.AUTO, BackendMode.PYTHON}:
            try:
                self._inner = _real_pysam().AlignmentFile(path, mode, **kwargs)
                return
            except ModuleNotFoundError:
                if backend == BackendMode.PYTHON:
                    raise
        raise NotImplementedError("Rust-backed bioscript.pysam is not available yet")

    def __enter__(self) -> "AlignmentFile":
        if hasattr(self._inner, "__enter__"):
            self._inner.__enter__()
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> object:
        return self._inner.__exit__(exc_type, exc, tb)

    def fetch(self, *args: Any, **kwargs: Any) -> Any:
        return self._inner.fetch(*args, **kwargs)

    def close(self) -> None:
        self._inner.close()
