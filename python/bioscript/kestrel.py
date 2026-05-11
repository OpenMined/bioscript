"""BioScript-supported Kestrel compatibility subset."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable


def build_command(
    jar_path: str,
    reference_vntr: str,
    output_vcf: str,
    output_sam: str,
    temp_dir: str,
    sample_name: str,
    fastq_1: str,
    fastq_2: str,
    *,
    java_program: str = "java",
    memory: str = "12g",
    kmer_size: int = 20,
    max_align_states: int = 40,
    max_hap_states: int = 40,
    log_level: str = "INFO",
    additional_args: Iterable[str] = (),
) -> list[str]:
    """Build the structured argv list for VNtyper's Kestrel invocation."""

    _validate_program(java_program)
    args = [
        java_program,
        f"-Xmx{memory}",
        "-jar",
        _path_arg(jar_path),
        "-k",
        str(kmer_size),
        "--maxalignstates",
        str(max_align_states),
        "--maxhapstates",
        str(max_hap_states),
        "-r",
        _path_arg(reference_vntr),
        "-o",
        _path_arg(output_vcf),
        f"-s{sample_name}",
        _path_arg(fastq_1),
        _path_arg(fastq_2),
        "--hapfmt",
        "sam",
        "-p",
        _path_arg(output_sam),
        "--logstderr",
        "--logstdout",
        "--loglevel",
        log_level.upper(),
        "--temploc",
        _path_arg(temp_dir),
    ]
    args.extend(str(arg) for arg in additional_args)
    return args


def run(*args: object, **kwargs: object) -> dict[str, object]:
    """Return the planned command for now; tool execution is runtime-owned."""

    argv = build_command(*args, **kwargs)
    return {
        "argv": argv,
        "vcf": kwargs.get("output_vcf") if "output_vcf" in kwargs else None,
        "sam": kwargs.get("output_sam") if "output_sam" in kwargs else None,
    }


def read_vcf(path: str) -> list[dict[str, str]]:
    """Read a small Kestrel VCF into dictionaries."""

    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            values = line.split("\t")
            rows.append({key: values[idx] if idx < len(values) else "" for idx, key in enumerate(header)})
    return rows


def call_sequences_native(
    reference_name: str,
    reference_sequence: str,
    read_sequences: Iterable[str],
    kmer_size: int,
    *,
    sample_name: str = "sample1",
    source_version: str = "native",
    reference_md5: str = ".",
    minimum_difference: int = 5,
    difference_quantile: float = 0.90,
    anchor_both_ends: bool = True,
    decay_min: float = 0.55,
    decay_alpha: float = 0.80,
    peak_scan_length: int = 7,
    scan_limit_factor: float = 7.0,
    recover_right_anchor: bool = True,
    min_kmer_count: int = 1,
    max_haplotypes: int = 40,
    max_bases: int = 500,
    max_repeat_count: int = 0,
    max_saved_states: int = 40,
    locus_depth: int = 1,
) -> str:
    """Run the native synthetic reads-to-VCF Kestrel path."""

    native = _native()
    return str(
        native.kestrel_call_sequences_native(
            reference_name,
            reference_sequence,
            list(read_sequences),
            int(kmer_size),
            sample_name,
            source_version,
            reference_md5,
            int(minimum_difference),
            float(difference_quantile),
            bool(anchor_both_ends),
            float(decay_min),
            float(decay_alpha),
            int(peak_scan_length),
            float(scan_limit_factor),
            bool(recover_right_anchor),
            int(min_kmer_count),
            int(max_haplotypes),
            int(max_bases),
            int(max_repeat_count),
            int(max_saved_states),
            int(locus_depth),
        )
    )


def call_fastq_native(
    reference_name: str,
    reference_sequence: str,
    fastq_paths: Iterable[str],
    kmer_size: int,
    *,
    sample_name: str = "sample1",
    source_version: str = "native",
    reference_md5: str = ".",
    minimum_difference: int = 5,
    difference_quantile: float = 0.90,
    anchor_both_ends: bool = True,
    decay_min: float = 0.55,
    decay_alpha: float = 0.80,
    peak_scan_length: int = 7,
    scan_limit_factor: float = 7.0,
    recover_right_anchor: bool = True,
    min_kmer_count: int = 1,
    max_haplotypes: int = 40,
    max_bases: int = 500,
    max_repeat_count: int = 0,
    max_saved_states: int = 40,
    locus_depth: int = 1,
) -> str:
    """Run the native FASTQ-to-VCF Kestrel path."""

    native = _native()
    return str(
        native.kestrel_call_fastq_native(
            reference_name,
            reference_sequence,
            [_path_arg(path) for path in fastq_paths],
            int(kmer_size),
            sample_name,
            source_version,
            reference_md5,
            int(minimum_difference),
            float(difference_quantile),
            bool(anchor_both_ends),
            float(decay_min),
            float(decay_alpha),
            int(peak_scan_length),
            float(scan_limit_factor),
            bool(recover_right_anchor),
            int(min_kmer_count),
            int(max_haplotypes),
            int(max_bases),
            int(max_repeat_count),
            int(max_saved_states),
            int(locus_depth),
        )
    )


def _path_arg(path: str) -> str:
    value = str(Path(path))
    if "\0" in value:
        raise ValueError("path arguments cannot contain NUL bytes")
    return value


def _validate_program(program: str) -> None:
    if not program.strip():
        raise ValueError("program cannot be empty")
    if "/" in program or any(ch in program for ch in "|&;<>`$\n\r"):
        raise ValueError(f"program must be a simple executable name: {program!r}")


def _native() -> Any:
    try:
        from . import _native as native
    except ImportError as exc:
        raise NotImplementedError("BioScript native Kestrel backend is not installed") from exc
    return native
