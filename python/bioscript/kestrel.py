"""BioScript-supported Kestrel compatibility subset."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Iterable

from .runtime import ModuleBackendPolicy

BACKEND_POLICY = ModuleBackendPolicy(
    auto="command builders and FASTA parsing are pure Python; native calls require bioscript._native",
    python="command builders and FASTA parsing are pure Python; native calls require bioscript._native",
    rust="native calls require bioscript._native backed by kestrel-rs",
)


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


def plan_command(*args: object, **kwargs: object) -> list[str]:
    return build_command(*args, **kwargs)


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


def load_reference_regions(path: str) -> list[tuple[str, str, str]]:
    """Read FASTA records as native Kestrel reference triples."""

    records: list[tuple[str, str, str]] = []
    current_name: str | None = None
    current_parts: list[str] = []
    with open(path, encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    records.append(_reference_region(current_name, current_parts))
                current_name = line[1:].split()[0]
                if not current_name:
                    raise ValueError("FASTA record name cannot be empty")
                current_parts = []
                continue
            if current_name is None:
                raise ValueError("FASTA sequence appeared before a record header")
            current_parts.append(line)
    if current_name is not None:
        records.append(_reference_region(current_name, current_parts))
    if not records:
        raise ValueError(f"FASTA file contains no records: {path}")
    return records


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
    max_gap_size: int | None = None,
    recover_right_anchor: bool = True,
    call_ambiguous_regions: bool = True,
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
            _optional_int(max_gap_size),
            bool(recover_right_anchor),
            bool(call_ambiguous_regions),
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
    max_gap_size: int | None = None,
    recover_right_anchor: bool = True,
    call_ambiguous_regions: bool = True,
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
            _optional_int(max_gap_size),
            bool(recover_right_anchor),
            bool(call_ambiguous_regions),
            int(min_kmer_count),
            int(max_haplotypes),
            int(max_bases),
            int(max_repeat_count),
            int(max_saved_states),
            int(locus_depth),
        )
    )


def call_fastq_references_native(
    references: Iterable[tuple[str, str, str]],
    fastq_paths: Iterable[str],
    kmer_size: int,
    *,
    sample_name: str = "sample1",
    source_version: str = "native",
    minimum_difference: int = 5,
    difference_quantile: float = 0.90,
    anchor_both_ends: bool = True,
    decay_min: float = 0.55,
    decay_alpha: float = 0.80,
    peak_scan_length: int = 7,
    scan_limit_factor: float = 7.0,
    max_gap_size: int | None = None,
    recover_right_anchor: bool = True,
    call_ambiguous_regions: bool = True,
    min_kmer_count: int = 1,
    max_haplotypes: int = 40,
    max_bases: int = 500,
    max_repeat_count: int = 0,
    max_saved_states: int = 40,
    locus_depth: int = 1,
) -> str:
    """Run the native FASTQ-to-VCF Kestrel path over multiple references."""

    native = _native()
    reference_rows = [(str(name), str(sequence), str(md5)) for name, sequence, md5 in references]
    return str(
        native.kestrel_call_fastq_references_native(
            reference_rows,
            [_path_arg(path) for path in fastq_paths],
            int(kmer_size),
            sample_name,
            source_version,
            int(minimum_difference),
            float(difference_quantile),
            bool(anchor_both_ends),
            float(decay_min),
            float(decay_alpha),
            int(peak_scan_length),
            float(scan_limit_factor),
            _optional_int(max_gap_size),
            bool(recover_right_anchor),
            bool(call_ambiguous_regions),
            int(min_kmer_count),
            int(max_haplotypes),
            int(max_bases),
            int(max_repeat_count),
            int(max_saved_states),
            int(locus_depth),
        )
    )


def run_native(
    reference_fasta: str,
    fastq_paths: Iterable[str],
    output_vcf: str,
    *,
    kmer_size: int = 20,
    sample_name: str = "sample1",
    minimum_difference: int = 5,
    difference_quantile: float = 0.90,
    min_kmer_count: int = 5,
    max_haplotypes: int = 40,
    max_bases: int = 500,
    max_saved_states: int = 40,
) -> str:
    """Run native Kestrel over FASTQs and write the resulting VCF."""

    vcf = call_fastq_references_native(
        load_reference_regions(reference_fasta),
        fastq_paths,
        kmer_size,
        sample_name=sample_name,
        minimum_difference=minimum_difference,
        difference_quantile=difference_quantile,
        min_kmer_count=min_kmer_count,
        max_haplotypes=max_haplotypes,
        max_bases=max_bases,
        max_saved_states=max_saved_states,
    )
    output = Path(_path_arg(output_vcf))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(vcf, encoding="utf-8")
    return str(output)


def _path_arg(path: str) -> str:
    value = str(Path(path))
    if "\0" in value:
        raise ValueError("path arguments cannot contain NUL bytes")
    return value


def _optional_int(value: int | None) -> int | None:
    if value is None:
        return None
    return int(value)


def _reference_region(name: str, sequence_parts: list[str]) -> tuple[str, str, str]:
    sequence = "".join(sequence_parts)
    if not sequence:
        raise ValueError(f"FASTA record contains no sequence: {name}")
    return (name, sequence, hashlib.md5(sequence.encode("ascii")).hexdigest())


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
