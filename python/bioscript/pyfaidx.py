"""BioScript-supported pyfaidx compatibility subset."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

from .runtime import BackendMode, ModuleBackendPolicy, selected_backend

BACKEND_POLICY = ModuleBackendPolicy(
    auto="use real pyfaidx when installed; otherwise use the pure Python FASTA fallback",
    python="requires real pyfaidx",
    rust="native pyfaidx shim is pending",
)


def _real_pyfaidx() -> Any:
    return importlib.import_module("pyfaidx")


class Fasta:
    """Small `pyfaidx.Fasta` subset with optional real-library delegation."""

    def __init__(self, path: str | Path, **kwargs: Any) -> None:
        backend = selected_backend()
        if backend in {BackendMode.AUTO, BackendMode.PYTHON}:
            try:
                self._inner = _real_pyfaidx().Fasta(path, **kwargs)
                self._simple = None
                return
            except ModuleNotFoundError:
                if backend == BackendMode.PYTHON:
                    raise
        if backend == BackendMode.RUST:
            raise NotImplementedError("Rust-backed bioscript.pyfaidx is not available yet")
        self._inner = None
        self._simple = _SimpleFasta(Path(path))

    def __getitem__(self, contig: str) -> Any:
        if self._inner is not None:
            return self._inner[contig]
        return self._simple[contig]


class _SimpleFasta:
    def __init__(self, path: Path) -> None:
        self.records = _read_fasta(path)

    def __getitem__(self, contig: str) -> "_SimpleRecord":
        try:
            return _SimpleRecord(self.records[contig])
        except KeyError as exc:
            raise KeyError(contig) from exc


class _SimpleRecord:
    def __init__(self, sequence: str) -> None:
        self.seq = sequence

    def __getitem__(self, key: slice) -> "_SimpleSequence":
        if not isinstance(key, slice):
            raise TypeError("BioScript pyfaidx fallback only supports slicing")
        return _SimpleSequence(self.seq[key])


class _SimpleSequence:
    def __init__(self, sequence: str) -> None:
        self.seq = sequence

    def __str__(self) -> str:
        return self.seq

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return self.seq == other
        return NotImplemented


def _read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(chunks)
            name = line[1:].split()[0]
            chunks = []
        elif name is None:
            raise ValueError("FASTA sequence appeared before first header")
        else:
            chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks)
    if not records:
        raise ValueError("FASTA did not contain any records")
    return records
