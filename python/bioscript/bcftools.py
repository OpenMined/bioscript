"""BioScript-supported bcftools command-builder subset."""

from __future__ import annotations

from pathlib import Path


def sort(input_vcf: str, output_vcf_gz: str) -> list[str]:
    return ["bcftools", "sort", "-Oz", "-o", _path_arg(output_vcf_gz), _path_arg(input_vcf)]


def index(vcf_gz: str) -> list[str]:
    return ["bcftools", "index", "-t", _path_arg(vcf_gz)]


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


def _path_arg(path: str) -> str:
    value = str(Path(path))
    if "\0" in value:
        raise ValueError("path arguments cannot contain NUL bytes")
    return value
