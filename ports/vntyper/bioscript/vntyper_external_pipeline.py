"""External-tool-backed VNtyper pipeline runner.

The command builders live in `vntyper_commands`; this module is the narrow
execution layer for the BAM path. It intentionally accepts an injectable runner
so tests can validate command order and output materialization without requiring
samtools, bcftools, or Kestrel.
"""

from __future__ import annotations

import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

try:
    from . import vntyper_commands, vntyper_port
except ImportError:
    import vntyper_commands
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


@dataclass(frozen=True)
class ExternalPipelineResult:
    participant_id: str
    output_dir: str
    commands: list[list[str]]
    kestrel_vcf: str
    kestrel_tsv: str
    report_json: str


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
    commands = [
        plan.samtools_view_command,
        plan.samtools_index_command,
        plan.samtools_fastq_command,
        plan.samtools_depth_command,
        plan.kestrel_command,
        plan.bcftools_sort_command,
        plan.bcftools_index_command,
    ]

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
    for command in commands:
        command_runner(command, check=True)
    materialize_post_kestrel_outputs(result, input_bam, assembly)
    return result


def create_output_dirs(result: ExternalPipelineResult, plan: vntyper_commands.VntyperCommandPlan) -> None:
    Path(result.output_dir).mkdir(parents=True, exist_ok=True)
    Path(plan.sliced_bam).parent.mkdir(parents=True, exist_ok=True)
    Path(plan.fastq_1).parent.mkdir(parents=True, exist_ok=True)
    Path(plan.kestrel_vcf).parent.mkdir(parents=True, exist_ok=True)


def materialize_post_kestrel_outputs(result: ExternalPipelineResult, input_bam: str, assembly: str) -> None:
    if not Path(result.kestrel_vcf).exists():
        raise FileNotFoundError(f"Kestrel VCF was not produced: {result.kestrel_vcf}")
    rows = vntyper_port.process_kestrel_vcf(result.kestrel_vcf)
    write_kestrel_result_tsv(result.kestrel_tsv, rows)
    report = vntyper_port.build_report_json(
        sample_name=result.participant_id,
        input_files={"bam": input_bam, "vcf": result.kestrel_vcf},
        kestrel_rows=rows,
        coverage={},
        metadata={
            "alignment_pipeline": "external samtools/kestrel",
            "detected_assembly": assembly,
        },
        pipeline_log=[{"command": command} for command in result.commands],
    )
    vntyper_port.write_report_json(result.report_json, report)


def write_kestrel_result_tsv(path: str, rows: list[dict[str, object]]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=KESTREL_TSV_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
