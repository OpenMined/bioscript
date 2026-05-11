"""Plan generation of large VNtyper expected outputs.

This script is intentionally not part of normal unit-test discovery. It is a
maintainer helper for files under ignored `ports/vntyper/test-data`.

Dry-run mode does not require samtools, Java, Kestrel, or the BAM files. Use it
to review the exact sample labels, command plans, and expected-output layout
before running an external-tool-backed pipeline.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = ROOT / "ports" / "vntyper" / "test-data"
EXPECTED_ROOT = DATA_ROOT / "expected"
DEFAULT_KESTREL_JAR = DATA_ROOT / "tools" / "kestrel" / "kestrel.jar"
VNTYPER_BIOSCRIPT = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper.bs.py"
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

import vntyper_commands  # noqa: E402
import vntyper_external_pipeline  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-sample", required=True, help="Sample basename without .bam")
    parser.add_argument("--negative-sample", required=True, help="Sample basename without .bam")
    parser.add_argument("--kestrel-jar", default=str(DEFAULT_KESTREL_JAR))
    parser.add_argument("--assembly", default="hg19")
    parser.add_argument(
        "--fastq-only",
        action="store_true",
        help="Generate Kestrel VCF/TSV/report outputs from existing FASTQ pairs without samtools.",
    )
    parser.add_argument(
        "--write-manifest",
        action="store_true",
        help="Write expected/manifest.json even in dry-run mode.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    payload = build_payload(args.positive_sample, args.negative_sample, args.assembly, args.kestrel_jar)
    if args.dry_run:
        print(json.dumps(payload, indent=2))
        if args.write_manifest:
            write_manifest(payload["manifest"])
        return 0

    missing = prerequisites(args.kestrel_jar, payload, fastq_only=args.fastq_only)
    if missing:
        raise SystemExit("Missing prerequisites: " + ", ".join(missing))

    for sample in payload["samples"]:
        if args.fastq_only:
            vntyper_external_pipeline.run_fastq_kestrel(
                sample["input_fastq_1"],
                sample["input_fastq_2"],
                sample["sample"],
                str(EXPECTED_ROOT / sample["label"]),
                kestrel_jar=args.kestrel_jar,
            )
        else:
            vntyper_external_pipeline.run_bam_pipeline(
                sample["input_bam"],
                sample["sample"],
                str(EXPECTED_ROOT / sample["label"]),
                assembly=args.assembly,
                kestrel_jar=args.kestrel_jar,
            )
    write_manifest(payload["manifest"])
    return 0


def build_payload(positive_sample: str, negative_sample: str, assembly: str, kestrel_jar: str) -> dict[str, object]:
    samples = [
        sample_payload("positive", positive_sample, assembly, kestrel_jar),
        sample_payload("negative", negative_sample, assembly, kestrel_jar),
    ]
    return {
        "note": (
            "This harness records the expected-output layout and command plans. "
            "Without --dry-run it executes the external-tool-backed runner and "
            "materializes ignored VCF/TSV/report outputs under test-data/expected."
        ),
        "bioscript_command_plan_commands": [sample["bioscript_command_plan_command"] for sample in samples],
        "samples": samples,
        "manifest": {
            "positive_sample": positive_sample,
            "negative_sample": negative_sample,
            "assembly": assembly,
            "expected_outputs": [
                "positive/kestrel/output.vcf",
                "positive/kestrel/kestrel_result.tsv",
                "negative/kestrel/output.vcf",
                "negative/kestrel/kestrel_result.tsv",
            ],
        },
    }


def sample_payload(label: str, sample: str, assembly: str, kestrel_jar: str) -> dict[str, object]:
    bam = DATA_ROOT / f"{sample}.bam"
    output_root = EXPECTED_ROOT / label
    work_dir = output_root / "work"
    plan = vntyper_commands.plan_bam_pipeline(
        str(bam),
        sample,
        assembly=assembly,
        work_dir=str(work_dir),
        kestrel_jar=kestrel_jar,
    )
    return {
        "label": label,
        "sample": sample,
        "input_bam": str(bam),
        "input_bai": str(DATA_ROOT / f"{sample}.bam.bai"),
        "input_fastq_1": str(DATA_ROOT / f"{sample}_R1.fastq.gz"),
        "input_fastq_2": str(DATA_ROOT / f"{sample}_R2.fastq.gz"),
        "expected_kestrel_vcf": str(output_root / "kestrel" / "output.vcf"),
        "expected_kestrel_tsv": str(output_root / "kestrel" / "kestrel_result.tsv"),
        "bioscript_command_plan_command": [
            "cargo",
            "run",
            "-p",
            "bioscript-cli",
            "--",
            str(VNTYPER_BIOSCRIPT),
            "--root",
            str(ROOT),
            "--input-file",
            str(bam),
            "--output-file",
            str(output_root / "command_plan.tsv"),
            "--participant-id",
            sample,
        ],
        "pipeline_command_plan": plan.as_report_row(),
    }


def write_manifest(manifest: dict[str, object]) -> None:
    EXPECTED_ROOT.mkdir(parents=True, exist_ok=True)
    (EXPECTED_ROOT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def prerequisites(kestrel_jar: str, payload: dict[str, object], fastq_only: bool = False) -> list[str]:
    missing = []
    if not fastq_only and shutil.which("samtools") is None:
        missing.append("samtools")
    if not fastq_only and shutil.which("bcftools") is None:
        missing.append("bcftools")
    if shutil.which("java") is None:
        missing.append("java")
    if not Path(kestrel_jar).exists():
        missing.append(kestrel_jar)
    muc1_reference = ROOT / vntyper_commands.DEFAULT_MUC1_REFERENCE
    if not muc1_reference.exists():
        missing.append(str(muc1_reference))
    for sample in payload["samples"]:
        keys = ["input_fastq_1", "input_fastq_2"] if fastq_only else ["input_bam", "input_bai"]
        for key in keys:
            if not Path(sample[key]).exists():
                missing.append(sample[key])
    return missing


if __name__ == "__main__":
    raise SystemExit(main())
