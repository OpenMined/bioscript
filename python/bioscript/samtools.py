"""BioScript-supported samtools command-builder subset."""

from __future__ import annotations

from pathlib import Path


def view_region(bam: str, region: str, output_bam: str, include_unmapped: bool = False) -> list[str]:
    args = ["samtools", "view", "-b", _path_arg(bam), region, "-o", _path_arg(output_bam)]
    if include_unmapped:
        args.extend(["-f", "4"])
    return args


def fastq(bam: str, fastq_1: str, fastq_2: str) -> list[str]:
    return ["samtools", "fastq", "-1", _path_arg(fastq_1), "-2", _path_arg(fastq_2), _path_arg(bam)]


def depth(bam: str, region: str) -> list[str]:
    return ["samtools", "depth", "-r", region, _path_arg(bam)]


def index(bam: str) -> list[str]:
    return ["samtools", "index", _path_arg(bam)]


def _path_arg(path: str) -> str:
    value = str(Path(path))
    if "\0" in value:
        raise ValueError("path arguments cannot contain NUL bytes")
    return value
