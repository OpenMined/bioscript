"""BioScript-supported bcftools command-planning subset."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .runtime import ModuleBackendPolicy

BACKEND_POLICY = ModuleBackendPolicy(
    auto="command builders are pure Python; native helpers require bioscript._native",
    python="command builders are pure Python; native helpers require bioscript._native",
    rust="native helpers require bioscript._native backed by bcftools-rs",
)


def sort(input_vcf: str, output_vcf_gz: str) -> list[str]:
    return ["bcftools", "sort", "-Oz", "-o", _path_arg(output_vcf_gz), _path_arg(input_vcf)]


def plan_sort(input_vcf: str, output_vcf_gz: str) -> list[str]:
    return sort(input_vcf, output_vcf_gz)


def index(vcf_gz: str) -> list[str]:
    return ["bcftools", "index", "-t", _path_arg(vcf_gz)]


def plan_index(vcf_gz: str) -> list[str]:
    return index(vcf_gz)


def view(input_vcf: str, output_vcf: str, output_type: str = "z") -> list[str]:
    return [
        "bcftools",
        "view",
        "-O",
        output_type,
        "-o",
        _path_arg(output_vcf),
        _path_arg(input_vcf),
    ]


def plan_view(input_vcf: str, output_vcf: str, output_type: str = "z") -> list[str]:
    return view(input_vcf, output_vcf, output_type)


def view_filter(input_vcf: str, output_vcf_gz: str, include_expr: str) -> list[str]:
    return [
        "bcftools",
        "view",
        "-i",
        include_expr,
        "-Oz",
        "-o",
        _path_arg(output_vcf_gz),
        _path_arg(input_vcf),
    ]


def plan_view_filter(input_vcf: str, output_vcf_gz: str, include_expr: str) -> list[str]:
    return view_filter(input_vcf, output_vcf_gz, include_expr)


def norm(input_vcf: str, reference_fasta: str, output_vcf_gz: str) -> list[str]:
    return [
        "bcftools",
        "norm",
        "-f",
        _path_arg(reference_fasta),
        "-Oz",
        "-o",
        _path_arg(output_vcf_gz),
        _path_arg(input_vcf),
    ]


def plan_norm(input_vcf: str, reference_fasta: str, output_vcf_gz: str) -> list[str]:
    return norm(input_vcf, reference_fasta, output_vcf_gz)


def view_header_native(input_vcf: str, output_vcf: str) -> None:
    native = _native()
    native.bcftools_view_header_native(_path_arg(input_vcf), _path_arg(output_vcf))


def view_native(input_vcf: str, output_vcf: str, output_type: str = "v") -> None:
    native = _native()
    native.bcftools_view_native(_path_arg(input_vcf), _path_arg(output_vcf), output_type)


def sort_native(
    input_vcf: str,
    output_vcf: str,
    *,
    output_type: str = "z",
    write_index: bool = True,
) -> None:
    native = _native()
    native.bcftools_sort_native(
        _path_arg(input_vcf),
        _path_arg(output_vcf),
        output_type,
        write_index,
    )


def index_native(
    vcf_gz: str,
    output_index: str | None = None,
    *,
    tbi: bool = True,
    force: bool = True,
) -> None:
    native = _native()
    native.bcftools_index_native(_path_arg(vcf_gz), _optional_path(output_index), tbi, force)


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
        raise NotImplementedError("BioScript native bcftools backend is not installed") from exc
    return native
