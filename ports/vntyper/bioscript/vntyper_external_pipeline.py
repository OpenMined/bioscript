"""External-tool-backed VNtyper pipeline runner.

The command builders live in `vntyper_commands`; this module is the narrow
execution layer for the BAM path. It intentionally accepts an injectable runner
so tests can validate command order and output materialization without requiring
samtools, bcftools, or Kestrel.
"""

from __future__ import annotations

import csv
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from bioscript import bcftools, kestrel, samtools

try:
    from . import vntyper_commands, vntyper_config, vntyper_port
except ImportError:
    import vntyper_commands
    import vntyper_config
    import vntyper_port


Runner = Callable[..., object]

KESTREL_TSV_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "Estimated_Depth_AlternateVariant",
    "Estimated_Depth_Variant_ActiveRegion",
    "Depth_Score",
    "Confidence",
    "Flag",
    "is_valid_frameshift",
    "alt_filter_pass",
    "passes_vntyper_filters",
]

NATIVE_KESTREL_MAX_HAPLOTYPES = vntyper_config.NATIVE_KESTREL_MAX_HAPLOTYPES
NATIVE_KESTREL_MAX_SAVED_STATES = vntyper_config.NATIVE_KESTREL_MAX_SAVED_STATES
NATIVE_KESTREL_MAX_BASES = vntyper_config.NATIVE_KESTREL_MAX_BASES
NATIVE_KESTREL_MIN_KMER_COUNT = vntyper_config.NATIVE_KESTREL_MIN_KMER_COUNT


@dataclass(frozen=True)
class ExternalPipelineResult:
    participant_id: str
    output_dir: str
    commands: list[list[str]]
    kestrel_vcf: str
    kestrel_tsv: str
    report_json: str


def run_vntyper(
    bam: str,
    reference_build: str = "hg19",
    output_dir: str = "vntyper-output",
    participant_id: str | None = None,
    **kwargs: object,
) -> ExternalPipelineResult:
    sample = participant_id or Path(bam).stem
    return run_bam_pipeline(
        bam,
        sample,
        output_dir,
        assembly=reference_build,
        **kwargs,
    )


def run_vntyper_fastq(
    r1: str,
    r2: str,
    reference_build: str = "hg19",
    output_dir: str = "vntyper-output",
    participant_id: str | None = None,
    **kwargs: object,
) -> ExternalPipelineResult:
    sample = participant_id or Path(r1).name.split("_")[0]
    return run_fastq_kestrel(
        r1,
        r2,
        sample,
        output_dir,
        assembly=reference_build,
        **kwargs,
    )


def run_bam_pipeline(
    input_bam: str,
    participant_id: str,
    output_dir: str,
    assembly: str = "hg19",
    chromosome_convention: str | None = None,
    kestrel_jar: str = vntyper_commands.DEFAULT_KESTREL_JAR,
    muc1_reference: str = vntyper_commands.DEFAULT_MUC1_REFERENCE,
    dry_run: bool = False,
    runner: Runner | None = None,
    use_native_samtools: bool = False,
    use_native_kestrel: bool = False,
    use_native_bcftools: bool = False,
    native_samtools: object | None = None,
    native_kestrel: object | None = None,
    native_bcftools: object | None = None,
) -> ExternalPipelineResult:
    out_dir = Path(output_dir)
    plan = vntyper_commands.plan_bam_pipeline(
        input_bam,
        participant_id,
        assembly=assembly,
        work_dir=str(out_dir),
        chromosome_convention=chromosome_convention,
        kestrel_jar=kestrel_jar,
        muc1_reference=muc1_reference,
    )
    commands = pipeline_commands(
        input_bam,
        plan,
        muc1_reference,
        use_native_samtools,
        use_native_kestrel,
        use_native_bcftools,
    )

    result = ExternalPipelineResult(
        participant_id=plan.participant_id,
        output_dir=str(out_dir),
        commands=commands,
        kestrel_vcf=plan.kestrel_vcf,
        kestrel_tsv=str(out_dir / "kestrel" / "kestrel_result.tsv"),
        report_json=str(out_dir / "report.json"),
    )
    if dry_run:
        return result

    create_output_dirs(result, plan)
    command_runner = runner or subprocess.run
    if use_native_samtools:
        backend = native_samtools or samtools
        index = default_bam_index(input_bam)
        backend.view_region_native(input_bam, plan.bam_region, plan.sliced_bam, index=index)
        backend.fastq_native(input_bam, plan.bam_region, plan.fastq_1, plan.fastq_2, index=index)
        coverage = backend.depth_native(input_bam, plan.vntr_region, index=index)
        if use_native_kestrel:
            run_native_kestrel(native_kestrel or kestrel, muc1_reference, plan, result.kestrel_vcf)
        else:
            command_runner(plan.kestrel_command, check=True)
        if use_native_bcftools:
            run_native_bcftools(native_bcftools or bcftools, plan)
        materialize_post_kestrel_outputs(
            result,
            input_bam,
            assembly,
            coverage,
            input_files=bam_input_files(input_bam, result.kestrel_vcf, plan, use_native_bcftools),
            alignment_pipeline=alignment_pipeline_label(use_native_samtools, use_native_kestrel),
        )
    else:
        depth_output = ""
        for command in external_commands(
            plan,
            include_kestrel=not use_native_kestrel,
            include_bcftools=not use_native_bcftools,
        ):
            if command == plan.samtools_depth_command:
                completed = command_runner(command, check=True, capture_output=True, text=True)
                depth_output = getattr(completed, "stdout", "") or ""
            else:
                command_runner(command, check=True)
        if use_native_kestrel:
            run_native_kestrel(native_kestrel or kestrel, muc1_reference, plan, result.kestrel_vcf)
        if use_native_bcftools:
            run_native_bcftools(native_bcftools or bcftools, plan)
        materialize_post_kestrel_outputs(
            result,
            input_bam,
            assembly,
            coverage_from_depth(depth_output),
            input_files=bam_input_files(input_bam, result.kestrel_vcf, plan, use_native_bcftools),
            alignment_pipeline=alignment_pipeline_label(use_native_samtools, use_native_kestrel),
        )
    return result


def pipeline_commands(
    input_bam: str,
    plan: vntyper_commands.VntyperCommandPlan,
    muc1_reference: str,
    use_native_samtools: bool,
    use_native_kestrel: bool,
    use_native_bcftools: bool,
) -> list[list[str]]:
    if use_native_samtools:
        commands = native_samtools_commands(input_bam, plan)
        if not use_native_kestrel:
            commands.append(plan.kestrel_command)
    else:
        commands = external_commands(
            plan,
            include_kestrel=not use_native_kestrel,
            include_bcftools=not use_native_bcftools,
        )
    if use_native_kestrel:
        commands.append(native_kestrel_command(plan, muc1_reference))
    if use_native_bcftools:
        commands.append(native_bcftools_sort_command(plan.kestrel_vcf, plan.sorted_vcf))
    return commands


def external_commands(
    plan: vntyper_commands.VntyperCommandPlan,
    include_kestrel: bool = True,
    include_bcftools: bool = True,
) -> list[list[str]]:
    commands = [
        plan.samtools_view_command,
        plan.samtools_index_command,
        plan.samtools_fastq_command,
        plan.samtools_depth_command,
    ]
    if include_kestrel:
        commands.append(plan.kestrel_command)
        if include_bcftools:
            commands.extend([plan.bcftools_sort_command, plan.bcftools_index_command])
    return commands


def native_samtools_commands(
    input_bam: str,
    plan: vntyper_commands.VntyperCommandPlan,
) -> list[list[str]]:
    index = default_bam_index(input_bam)
    return [
        [
            "bioscript.samtools.view_region_native",
            input_bam,
            plan.bam_region,
            plan.sliced_bam,
            "--index",
            index,
        ],
        [
            "bioscript.samtools.fastq_native",
            input_bam,
            plan.bam_region,
            plan.fastq_1,
            plan.fastq_2,
            "--index",
            index,
        ],
        [
            "bioscript.samtools.depth_native",
            input_bam,
            plan.vntr_region,
            "--index",
            index,
        ],
    ]


def native_kestrel_command(
    plan: vntyper_commands.VntyperCommandPlan,
    muc1_reference: str,
) -> list[str]:
    return [
        "bioscript.kestrel.run_native",
        muc1_reference,
        plan.fastq_1,
        plan.fastq_2,
        "-o",
        plan.kestrel_vcf,
    ]


def run_native_kestrel(
    backend: object,
    muc1_reference: str,
    plan: vntyper_commands.VntyperCommandPlan,
    output_vcf: str,
) -> None:
    backend.run_native(
        muc1_reference,
        [plan.fastq_1, plan.fastq_2],
        output_vcf,
        kmer_size=20,
        sample_name=plan.participant_id,
        min_kmer_count=NATIVE_KESTREL_MIN_KMER_COUNT,
        max_haplotypes=NATIVE_KESTREL_MAX_HAPLOTYPES,
        max_saved_states=NATIVE_KESTREL_MAX_SAVED_STATES,
        max_bases=NATIVE_KESTREL_MAX_BASES,
    )


def run_native_bcftools(
    backend: object,
    plan: vntyper_commands.VntyperCommandPlan,
) -> None:
    backend.sort_native(
        plan.kestrel_vcf,
        plan.sorted_vcf,
        output_type="z",
        write_index=True,
    )


def bam_input_files(
    input_bam: str,
    kestrel_vcf: str,
    plan: vntyper_commands.VntyperCommandPlan,
    include_sorted_vcf: bool,
) -> dict[str, str]:
    files = {"bam": input_bam, "vcf": kestrel_vcf}
    if include_sorted_vcf:
        files["sorted_vcf"] = plan.sorted_vcf
    return files


def alignment_pipeline_label(use_native_samtools: bool, use_native_kestrel: bool) -> str:
    if use_native_samtools and use_native_kestrel:
        return "native bioscript samtools/kestrel"
    if use_native_samtools:
        return "native bioscript samtools/kestrel"
    if use_native_kestrel:
        return "external samtools/native bioscript kestrel"
    return "external samtools/kestrel"


def default_bam_index(input_bam: str) -> str:
    return f"{input_bam}.bai"


def run_fastq_kestrel(
    fastq_1: str,
    fastq_2: str,
    participant_id: str,
    output_dir: str,
    assembly: str = "unknown",
    kestrel_jar: str = vntyper_commands.DEFAULT_KESTREL_JAR,
    muc1_reference: str = vntyper_commands.DEFAULT_MUC1_REFERENCE,
    dry_run: bool = False,
    runner: Runner | None = None,
    use_native_kestrel: bool = False,
    use_native_bcftools: bool = False,
    native_kestrel: object | None = None,
    native_bcftools: object | None = None,
) -> ExternalPipelineResult:
    out_dir = Path(output_dir)
    sample = vntyper_commands._safe_sample_name(participant_id)
    kestrel_dir = out_dir / "kestrel"
    kestrel_vcf = str(kestrel_dir / "output.vcf")
    kestrel_sam = str(kestrel_dir / "output.sam")
    sorted_vcf = str(kestrel_dir / "output.sorted.vcf.gz")
    if use_native_kestrel:
        command = native_kestrel_fastq_command(muc1_reference, fastq_1, fastq_2, kestrel_vcf)
    else:
        command = kestrel.build_command(
            kestrel_jar,
            muc1_reference,
            kestrel_vcf,
            kestrel_sam,
            str(kestrel_dir / "tmp"),
            sample,
            fastq_1,
            fastq_2,
        )
    result = ExternalPipelineResult(
        participant_id=sample,
        output_dir=str(out_dir),
        commands=(
            [command, native_bcftools_sort_command(kestrel_vcf, sorted_vcf)]
            if use_native_bcftools
            else [command]
        ),
        kestrel_vcf=kestrel_vcf,
        kestrel_tsv=str(kestrel_dir / "kestrel_result.tsv"),
        report_json=str(out_dir / "report.json"),
    )
    if dry_run:
        return result

    Path(result.kestrel_vcf).parent.mkdir(parents=True, exist_ok=True)
    Path(kestrel_dir / "tmp").mkdir(parents=True, exist_ok=True)
    if use_native_kestrel:
        plan = SimpleFastqKestrelPlan(sample, muc1_reference, fastq_1, fastq_2)
        run_native_kestrel(native_kestrel or kestrel, muc1_reference, plan, result.kestrel_vcf)
    else:
        command_runner = runner or subprocess.run
        command_runner(command, check=True)
    if use_native_bcftools:
        (native_bcftools or bcftools).sort_native(
            result.kestrel_vcf,
            sorted_vcf,
            output_type="z",
            write_index=True,
        )
    materialize_post_kestrel_outputs(
        result,
        f"{fastq_1},{fastq_2}",
        assembly,
        {},
        input_files={
            "fastq_1": fastq_1,
            "fastq_2": fastq_2,
            "vcf": result.kestrel_vcf,
            "sorted_vcf": sorted_vcf,
        }
        if use_native_bcftools
        else {"fastq_1": fastq_1, "fastq_2": fastq_2, "vcf": result.kestrel_vcf},
        alignment_pipeline=(
            "native bioscript kestrel from FASTQ"
            if use_native_kestrel
            else "external kestrel from FASTQ"
        ),
    )
    return result


@dataclass(frozen=True)
class SimpleFastqKestrelPlan:
    participant_id: str
    muc1_reference: str
    fastq_1: str
    fastq_2: str


def native_kestrel_fastq_command(
    muc1_reference: str,
    fastq_1: str,
    fastq_2: str,
    output_vcf: str,
) -> list[str]:
    return [
        "bioscript.kestrel.run_native",
        muc1_reference,
        fastq_1,
        fastq_2,
        "-o",
        output_vcf,
    ]


def native_bcftools_sort_command(input_vcf: str, output_vcf: str) -> list[str]:
    return [
        "bioscript.bcftools.sort_native",
        input_vcf,
        output_vcf,
        "--output-type",
        "z",
        "--write-index",
    ]


def create_output_dirs(result: ExternalPipelineResult, plan: vntyper_commands.VntyperCommandPlan) -> None:
    Path(result.output_dir).mkdir(parents=True, exist_ok=True)
    Path(plan.sliced_bam).parent.mkdir(parents=True, exist_ok=True)
    Path(plan.fastq_1).parent.mkdir(parents=True, exist_ok=True)
    Path(plan.kestrel_vcf).parent.mkdir(parents=True, exist_ok=True)
    Path(plan.kestrel_vcf).parent.joinpath("tmp").mkdir(parents=True, exist_ok=True)


def materialize_post_kestrel_outputs(
    result: ExternalPipelineResult,
    input_bam: str,
    assembly: str,
    coverage: dict[str, float | int] | None = None,
    input_files: dict[str, str] | None = None,
    alignment_pipeline: str = "external samtools/kestrel",
) -> None:
    if not Path(result.kestrel_vcf).exists():
        raise FileNotFoundError(f"Kestrel VCF was not produced: {result.kestrel_vcf}")
    rows = vntyper_port.process_kestrel_vcf(result.kestrel_vcf)
    write_kestrel_result_tsv(result.kestrel_tsv, rows)
    report = vntyper_port.build_report_json(
        sample_name=result.participant_id,
        input_files=input_files or {"bam": input_bam, "vcf": result.kestrel_vcf},
        kestrel_rows=rows,
        coverage=coverage or {},
        metadata={
            "alignment_pipeline": alignment_pipeline,
            "detected_assembly": assembly,
        },
        pipeline_log=[{"command": command} for command in result.commands],
    )
    vntyper_port.write_report_json(result.report_json, report)


def coverage_from_depth(depth_output: str) -> dict[str, float | int]:
    depths = []
    for raw_line in depth_output.splitlines():
        fields = raw_line.split("\t")
        if len(fields) < 3:
            continue
        try:
            depths.append(int(fields[2]))
        except ValueError:
            continue
    if not depths:
        return {}
    zero_count = sum(1 for depth in depths if depth == 0)
    return {
        "mean": statistics.fmean(depths),
        "median": statistics.median(depths),
        "stdev": statistics.pstdev(depths),
        "min": min(depths),
        "max": max(depths),
        "region_length": len(depths),
        "uncovered_bases": zero_count,
        "percent_uncovered": zero_count / len(depths) * 100,
    }


def write_kestrel_result_tsv(path: str, rows: list[dict[str, object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=KESTREL_TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
