#!/usr/bin/env python3
"""Run one VNtyper engine+input combo and emit a normalized output.

Test plumbing for `test-vntyper.sh`:

    run_parity_pipeline.py --engine java --input fastq --json /tmp/java.json
    run_parity_pipeline.py --engine rust --input fastq --json /tmp/rust.json
    diff -u /tmp/java.json /tmp/rust.json   # empty diff = parity

The pipeline is the same `run_bam_pipeline` / `run_fastq_kestrel` the opt-in
gate tests use; this script just exposes the normalized output (classification,
TSV fingerprint over stable columns, filtered report summary) without paths,
timestamps, or tool-version metadata.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import sys
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
TESTS_ROOT = ROOT / "ports" / "vntyper" / "tests"
MANIFEST_PATH = TESTS_ROOT / "data_manifest.py"
PIPELINE_PATH = BIOSCRIPT_PORT / "vntyper_external_pipeline.py"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))
sys.path.insert(0, str(TESTS_ROOT))

_manifest_spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(_manifest_spec)
_manifest_spec.loader.exec_module(data_manifest)

_pipeline_spec = importlib.util.spec_from_file_location(
    "vntyper_external_pipeline", PIPELINE_PATH
)
vntyper_external_pipeline = importlib.util.module_from_spec(_pipeline_spec)
sys.modules["vntyper_external_pipeline"] = vntyper_external_pipeline
_pipeline_spec.loader.exec_module(vntyper_external_pipeline)

from parity_helpers import normalized_report_summary, normalized_tsv_fingerprint


STABLE_FIELDS = [
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


def _row_is_passing(row):
    return row.get("passes_vntyper_filters") in ("True", True)


def _passing_rows(rows, limit):
    passing = [row for row in rows if _row_is_passing(row)]
    passing.sort(
        key=lambda row: float(row.get("Depth_Score") or 0),
        reverse=True,
    )
    return [
        {field: row.get(field, "") for field in STABLE_FIELDS}
        for row in passing[:limit]
    ]


def _run_one_case(engine, input_kind, label, fixture, output_dir, kestrel_jar):
    use_native = engine == "rust"
    case_dir = output_dir / engine / label
    case_dir.parent.mkdir(parents=True, exist_ok=True)

    started = time.monotonic()
    if input_kind == "bam":
        result = vntyper_external_pipeline.run_bam_pipeline(
            fixture,
            label,
            str(case_dir),
            kestrel_jar=kestrel_jar,
            muc1_reference=str(data_manifest.MUC1_REFERENCE),
            use_native_samtools=use_native,
            use_native_kestrel=use_native,
            use_native_bcftools=use_native,
        )
    else:
        fastq_1, fastq_2 = fixture
        result = vntyper_external_pipeline.run_fastq_kestrel(
            fastq_1,
            fastq_2,
            label,
            str(case_dir),
            assembly="hg19",
            kestrel_jar=kestrel_jar,
            muc1_reference=str(data_manifest.MUC1_REFERENCE),
            use_native_kestrel=use_native,
            use_native_bcftools=use_native,
        )
    elapsed = time.monotonic() - started

    with open(result.report_json, "r", encoding="utf-8") as handle:
        report = json.load(handle)
    with open(result.kestrel_tsv, "r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    return {
        "engine": engine,
        "input": input_kind,
        "case": label,
        "wall_seconds": round(elapsed, 3),
        "classification": report.get("algorithm_results", {}).get("kestrel"),
        "screening_summary": report.get("screening_summary"),
        "tsv_fingerprint": normalized_tsv_fingerprint(rows),
        "report_summary": normalized_report_summary(report),
        "top_passing_rows": _passing_rows(rows, limit=10),
    }


def _resolve_cases(input_kind, case_filter):
    if input_kind == "bam":
        cases = data_manifest.REPRESENTATIVE_BAM_CASES
        fixtures = {label: str(path) for label, path in cases.items()}
    else:
        cases = data_manifest.REPRESENTATIVE_FASTQ_CASES
        fixtures = {label: (str(pair[0]), str(pair[1])) for label, pair in cases.items()}
    if case_filter:
        if case_filter not in fixtures:
            raise SystemExit(
                f"unknown case '{case_filter}'; choices: {sorted(fixtures)}"
            )
        return {case_filter: fixtures[case_filter]}
    return fixtures


def _check_prerequisites(engine, input_kind, fixtures):
    """Return a list of missing prerequisites. Empty list = ready to run."""
    missing = []
    if not data_manifest.MUC1_REFERENCE.exists():
        missing.append(str(data_manifest.MUC1_REFERENCE))
    if engine == "java":
        import shutil

        if shutil.which("java") is None:
            missing.append("java on PATH")
        if not data_manifest.KESTREL_JAR.exists():
            missing.append(str(data_manifest.KESTREL_JAR))
        if input_kind == "bam":
            if shutil.which("samtools") is None:
                missing.append("samtools on PATH")
            if shutil.which("bcftools") is None:
                missing.append("bcftools on PATH")
    else:
        try:
            data_manifest.import_native_module()
        except Exception as exc:
            missing.append(f"bioscript._native importable ({exc})")
    for label, fixture in fixtures.items():
        if input_kind == "bam":
            paths = [Path(fixture), Path(f"{fixture}.bai")]
        else:
            paths = [Path(fixture[0]), Path(fixture[1])]
        for path in paths:
            if not path.exists():
                missing.append(f"{label}: {path}")
    return missing


def _prepend_tool_path():
    """If the vendored samtools/bcftools tree exists, put it first on PATH."""
    tool_bin = data_manifest.LOCAL_TOOL_BIN
    if tool_bin.exists():
        old = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{tool_bin}{os.pathsep}{old}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--engine", choices=["java", "rust"], required=True)
    parser.add_argument("--input", choices=["bam", "fastq"], required=True)
    parser.add_argument(
        "--case",
        choices=["positive", "negative"],
        default=None,
        help="restrict to one fixture (default: both positive and negative)",
    )
    parser.add_argument(
        "--out-dir",
        help="scratch directory for pipeline artifacts (default: a tempdir)",
    )
    parser.add_argument(
        "--json",
        help="path to write the normalized JSON output for diffing",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress per-case progress on stderr (still writes JSON)",
    )
    args = parser.parse_args()

    fixtures = _resolve_cases(args.input, args.case)

    _prepend_tool_path()

    missing = _check_prerequisites(args.engine, args.input, fixtures)
    if missing:
        print(
            f"prerequisites missing for engine={args.engine} input={args.input}:",
            file=sys.stderr,
        )
        for item in missing:
            print(f"  - {item}", file=sys.stderr)
        return 3

    if args.out_dir:
        out_dir = Path(args.out_dir)
    else:
        out_dir = Path(tempfile.mkdtemp(prefix=f"vntyper-parity-{args.engine}-"))
    out_dir.mkdir(parents=True, exist_ok=True)

    if not args.quiet:
        print(
            f"engine={args.engine} input={args.input} cases={sorted(fixtures)}",
            file=sys.stderr,
        )
        print(f"out_dir={out_dir}", file=sys.stderr)

    cases = {}
    for label in sorted(fixtures):
        if not args.quiet:
            print(f"  running case={label} ...", file=sys.stderr)
        cases[label] = _run_one_case(
            args.engine,
            args.input,
            label,
            fixtures[label],
            out_dir,
            str(data_manifest.KESTREL_JAR),
        )
        if not args.quiet:
            summary = cases[label]
            print(
                f"    classification={summary['classification']!r} "
                f"rows={summary['tsv_fingerprint']['row_count']} "
                f"passing={summary['tsv_fingerprint']['passing_count']} "
                f"sha={summary['tsv_fingerprint']['sha256'][:12]} "
                f"wall={summary['wall_seconds']}s",
                file=sys.stderr,
            )

    out = {
        "engine": args.engine,
        "input": args.input,
        "cases": cases,
    }

    payload = json.dumps(out, indent=2, sort_keys=True)
    print(payload)
    if args.json:
        Path(args.json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json).write_text(payload + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
