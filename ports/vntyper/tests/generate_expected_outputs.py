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
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA_ROOT = ROOT / "ports" / "vntyper" / "test-data"
EXPECTED_ROOT = DATA_ROOT / "expected"
VNTYPER_BIOSCRIPT = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper.bs.py"
RUST_ROOT = ROOT / "rust"
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

import vntyper_commands  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-sample", required=True, help="Sample basename without .bam")
    parser.add_argument("--negative-sample", required=True, help="Sample basename without .bam")
    parser.add_argument("--kestrel-jar", default=str(ROOT / "ports" / "vntyper" / "kestrel" / "kestrel.jar"))
    parser.add_argument("--assembly", default="hg19")
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

    missing = prerequisites(args.kestrel_jar, payload)
    if missing:
        raise SystemExit("Missing prerequisites: " + ", ".join(missing))

    for command in payload["bioscript_command_plan_commands"]:
        subprocess.run(command, cwd=RUST_ROOT, check=True)
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
            "The current BioScript entrypoint writes command plans; full VCF/TSV "
            "materialization is enabled once the external pipeline runner lands."
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


def prerequisites(kestrel_jar: str, payload: dict[str, object]) -> list[str]:
    missing = []
    if shutil.which("samtools") is None:
        missing.append("samtools")
    if shutil.which("java") is None:
        missing.append("java")
    if not Path(kestrel_jar).exists():
        missing.append(kestrel_jar)
    for sample in payload["samples"]:
        for key in ["input_bam", "input_bai"]:
            if not Path(sample[key]).exists():
                missing.append(sample[key])
    return missing


if __name__ == "__main__":
    raise SystemExit(main())
