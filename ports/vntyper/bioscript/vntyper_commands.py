"""Command planning helpers for the minimal VNtyper BioScript port."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from bioscript import bcftools, kestrel, samtools

try:
    from . import vntyper_config, vntyper_regions
except ImportError:
    import vntyper_config
    import vntyper_regions


DEFAULT_KESTREL_JAR = vntyper_config.DEFAULT_KESTREL_JAR
DEFAULT_MUC1_REFERENCE = vntyper_config.DEFAULT_MUC1_REFERENCE


@dataclass(frozen=True)
class VntyperCommandPlan:
    participant_id: str
    assembly: str
    bam_region: str
    vntr_region: str
    sliced_bam: str
    fastq_1: str
    fastq_2: str
    kestrel_vcf: str
    kestrel_sam: str
    sorted_vcf: str
    samtools_view_command: list[str]
    samtools_index_command: list[str]
    samtools_fastq_command: list[str]
    samtools_depth_command: list[str]
    kestrel_command: list[str]
    bcftools_sort_command: list[str]
    bcftools_index_command: list[str]

    def as_report_row(self) -> dict[str, Any]:
        return {
            "participant_id": self.participant_id,
            "assembly": self.assembly,
            "bam_region": self.bam_region,
            "vntr_region": self.vntr_region,
            "samtools_view_command": self.samtools_view_command,
            "samtools_index_command": self.samtools_index_command,
            "samtools_fastq_command": self.samtools_fastq_command,
            "samtools_depth_command": self.samtools_depth_command,
            "kestrel_command": self.kestrel_command,
            "bcftools_sort_command": self.bcftools_sort_command,
            "bcftools_index_command": self.bcftools_index_command,
        }


def plan_bam_pipeline(
    input_bam: str,
    participant_id: str,
    assembly: str = "hg19",
    work_dir: str = "vntyper",
    chromosome_convention: str | None = None,
    kestrel_jar: str = DEFAULT_KESTREL_JAR,
    muc1_reference: str = DEFAULT_MUC1_REFERENCE,
) -> VntyperCommandPlan:
    bam_region = vntyper_regions.region_string(
        assembly,
        "bam_region_coords",
        convention=chromosome_convention,
    )
    vntr_region = vntyper_regions.region_string(
        assembly,
        "vntr_region_coords",
        convention=chromosome_convention,
    )

    root = Path(work_dir)
    sample = _safe_sample_name(participant_id)
    sliced_bam = str(root / f"{sample}_sliced.bam")
    fastq_1 = str(root / f"{sample}_R1.fastq.gz")
    fastq_2 = str(root / f"{sample}_R2.fastq.gz")
    kestrel_dir = root / "kestrel"
    kestrel_vcf = str(kestrel_dir / "output.vcf")
    kestrel_sam = str(kestrel_dir / "output.sam")
    sorted_vcf = str(kestrel_dir / "output.sorted.vcf.gz")

    return VntyperCommandPlan(
        participant_id=sample,
        assembly=assembly,
        bam_region=bam_region,
        vntr_region=vntr_region,
        sliced_bam=sliced_bam,
        fastq_1=fastq_1,
        fastq_2=fastq_2,
        kestrel_vcf=kestrel_vcf,
        kestrel_sam=kestrel_sam,
        sorted_vcf=sorted_vcf,
        samtools_view_command=samtools.view_region(input_bam, bam_region, sliced_bam),
        samtools_index_command=samtools.index(sliced_bam),
        samtools_fastq_command=samtools.fastq(sliced_bam, fastq_1, fastq_2),
        samtools_depth_command=samtools.depth(sliced_bam, vntr_region, include_zero=True),
        kestrel_command=kestrel.build_command(
            kestrel_jar,
            muc1_reference,
            kestrel_vcf,
            kestrel_sam,
            str(kestrel_dir / "tmp"),
            sample,
            fastq_1,
            fastq_2,
        ),
        bcftools_sort_command=bcftools.sort(kestrel_vcf, sorted_vcf),
        bcftools_index_command=bcftools.index(sorted_vcf),
    )


def _safe_sample_name(participant_id: str) -> str:
    if not participant_id:
        raise ValueError("participant_id is required")
    if "/" in participant_id or "\\" in participant_id or "\0" in participant_id:
        raise ValueError("participant_id must be a simple sample name")
    return participant_id
