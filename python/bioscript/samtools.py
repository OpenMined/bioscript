"""BioScript-supported samtools command-builder subset."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .runtime import ModuleBackendPolicy

BACKEND_POLICY = ModuleBackendPolicy(
    auto="command builders are pure Python; native helpers require bioscript._native",
    python="command builders are pure Python; native helpers require bioscript._native",
    rust="native helpers require bioscript._native with the samtools-rs backend",
)


def view_region(bam: str, region: str, output_bam: str, include_unmapped: bool = False) -> list[str]:
    args = ["samtools", "view", "-b", _path_arg(bam), region, "-o", _path_arg(output_bam)]
    if include_unmapped:
        args.extend(["-f", "4"])
    return args


def fastq(bam: str, fastq_1: str, fastq_2: str) -> list[str]:
    return ["samtools", "fastq", "-1", _path_arg(fastq_1), "-2", _path_arg(fastq_2), _path_arg(bam)]


def depth(bam: str, region: str, include_zero: bool = False) -> list[str]:
    args = ["samtools", "depth"]
    if include_zero:
        args.append("-a")
    args.extend(["-r", region, _path_arg(bam)])
    return args


def index(bam: str) -> list[str]:
    return ["samtools", "index", _path_arg(bam)]


def view_region_native(bam: str, region: str, output_bam: str, index: str | None = None) -> int:
    native = _native()
    return int(
        native.samtools_view_region_native(
            _path_arg(bam),
            _optional_path(index),
            region,
            _path_arg(output_bam),
        )
    )


def depth_native(bam: str, region: str, index: str | None = None) -> dict[str, float]:
    native = _native()
    return dict(native.samtools_depth_native(_path_arg(bam), _optional_path(index), region))


def fastq_native(
    bam: str,
    region: str,
    fastq_1: str,
    fastq_2: str,
    index: str | None = None,
) -> dict[str, int]:
    native = _native()
    return {
        key: int(value)
        for key, value in native.samtools_fastq_native(
            _path_arg(bam),
            _optional_path(index),
            region,
            _path_arg(fastq_1),
            _path_arg(fastq_2),
        ).items()
    }


def _path_arg(path: str) -> str:
    value = str(Path(path))
    if "\0" in value:
        raise ValueError("path arguments cannot contain NUL bytes")
    return value


def _optional_path(path: str | None) -> str | None:
    if path is None:
        return None
    return _path_arg(path)


def _native() -> Any:
    try:
        from . import _native as native
    except ImportError as exc:
        raise NotImplementedError("BioScript native samtools backend is not installed") from exc
    return native
