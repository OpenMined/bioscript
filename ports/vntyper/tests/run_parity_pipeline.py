#!/usr/bin/env python3
"""Run every upstream-asserted VNtyper fixture through one engine.

Test plumbing for `test-vntyper.sh`:

    run_parity_pipeline.py --engine java --input bam --json /tmp/java.json
    run_parity_pipeline.py --engine rust --input bam --json /tmp/rust.json
    diff_parity_outputs.py /tmp/java.json /tmp/rust.json   # parity check

For every fixture upstream ships a `kestrel_assertions` entry for
(`ports/vntyper/vntyper/tests/test_data_config.json`), this runs the same
`run_bam_pipeline` / `run_fastq_kestrel` the opt-in gate tests use, finds
the called variant (top passing row), and asserts it against upstream's
expected Confidence and depth tolerances. Exit is non-zero if any fixture
is mis-called, so a wrong positive/negative is a hard failure, not a skip.
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
import upstream_expectations


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


def _to_float(value):
    try:
        if value is None or value == "" or value == "None":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _row_is_passing(row):
    return row.get("passes_vntyper_filters") in ("True", True)


def _called_variant(rows):
    """The variant the pipeline calls: top passing row by Depth_Score.

    Mirrors upstream taking `rows[0]` of the finalized kestrel_result.tsv
    (passing variants, highest depth-score first). Returns the stable-field
    view of that row, or None when nothing passes (a negative call).
    """
    passing = [r for r in rows if _row_is_passing(r)]
    if not passing:
        return None
    passing.sort(
        key=lambda r: (
            _to_float(r.get("Depth_Score")) or 0.0,
            _to_float(r.get("Estimated_Depth_AlternateVariant")) or 0.0,
        ),
        reverse=True,
    )
    top = passing[0]
    return {field: top.get(field, "") for field in STABLE_FIELDS}


def _passing_rows(rows, limit):
    passing = [r for r in rows if _row_is_passing(r)]
    passing.sort(key=lambda r: _to_float(r.get("Depth_Score")) or 0.0, reverse=True)
    return [
        {field: r.get(field, "") for field in STABLE_FIELDS}
        for r in passing[:limit]
    ]


def _fastq_pair(stem):
    r1 = data_manifest.DATA_ROOT / f"{stem}_R1.fastq.gz"
    r2 = data_manifest.DATA_ROOT / f"{stem}_R2.fastq.gz"
    return (r1, r2) if r1.exists() and r2.exists() else None


def _run_fixture(engine, stem, expectation, input_kind, out_dir, kestrel_jar):
    use_native = engine == "rust"
    assembly = expectation.get("reference_assembly", "hg19")
    case_dir = out_dir / engine / stem
    case_dir.parent.mkdir(parents=True, exist_ok=True)

    started = time.monotonic()
    if input_kind == "bam":
        result = vntyper_external_pipeline.run_bam_pipeline(
            expectation["bam"],
            stem,
            str(case_dir),
            assembly=assembly,
            kestrel_jar=kestrel_jar,
            muc1_reference=str(data_manifest.MUC1_REFERENCE),
            use_native_samtools=use_native,
            use_native_kestrel=use_native,
            use_native_bcftools=use_native,
        )
    else:
        pair = _fastq_pair(stem)
        fastq_1, fastq_2 = str(pair[0]), str(pair[1])
        result = vntyper_external_pipeline.run_fastq_kestrel(
            fastq_1,
            fastq_2,
            stem,
            str(case_dir),
            assembly=assembly,
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

    called = _called_variant(rows)
    correct, reasons = upstream_expectations.evaluate(expectation, called)

    return {
        "engine": engine,
        "input": input_kind,
        "fixture": stem,
        "assembly": assembly,
        "wall_seconds": round(elapsed, 3),
        "expected": {
            "confidence": expectation["confidence"],
            "is_negative": expectation["is_negative"],
            "alt_depth": expectation.get("alt_depth"),
            "region_depth": expectation.get("region_depth"),
            "depth_score": expectation.get("depth_score"),
        },
        "called": called,
        "correct": correct,
        "reasons": reasons,
        "classification": report.get("algorithm_results", {}).get("kestrel"),
        "screening_summary": report.get("screening_summary"),
        "tsv_fingerprint": normalized_tsv_fingerprint(rows),
        "report_summary": normalized_report_summary(report),
        "top_passing_rows": _passing_rows(rows, limit=5),
    }


def _resolve_fixtures(input_kind, fixture_filter):
    expectations = upstream_expectations.load_expectations()
    if not expectations:
        raise SystemExit(
            "no upstream-asserted fixtures found under test-data; "
            "check ports/vntyper/test-data and test_data_config.json"
        )
    chosen = {}
    for stem, exp in expectations.items():
        if fixture_filter and fixture_filter not in stem:
            continue
        if input_kind == "fastq" and _fastq_pair(stem) is None:
            continue
        chosen[stem] = exp
    if not chosen:
        raise SystemExit(
            f"no fixtures match filter={fixture_filter!r} for input={input_kind} "
            f"(known: {sorted(expectations)})"
        )
    return chosen


def _check_prerequisites(engine, input_kind, fixtures):
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
    for stem, exp in fixtures.items():
        if input_kind == "bam":
            paths = [Path(exp["bam"]), Path(f"{exp['bam']}.bai")]
        else:
            pair = _fastq_pair(stem)
            paths = [pair[0], pair[1]] if pair else [Path(f"{stem}_R1.fastq.gz")]
        for path in paths:
            if not Path(path).exists():
                missing.append(f"{stem}: {path}")
    return missing


def _prepend_tool_path():
    tool_bin = data_manifest.LOCAL_TOOL_BIN
    if tool_bin.exists():
        old = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{tool_bin}{os.pathsep}{old}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--engine", choices=["java", "rust"], required=True)
    parser.add_argument("--input", choices=["bam", "fastq"], default="bam")
    parser.add_argument(
        "--fixture",
        default=None,
        help="substring filter on fixture name (e.g. 66bf); default: all",
    )
    parser.add_argument("--out-dir", help="scratch dir (default: a tempdir)")
    parser.add_argument("--json", help="path to write the JSON output")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    fixtures = _resolve_fixtures(args.input, args.fixture)
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
            f"engine={args.engine} input={args.input} "
            f"fixtures={sorted(fixtures)}",
            file=sys.stderr,
        )
        print(f"out_dir={out_dir}", file=sys.stderr)

    cases = {}
    all_correct = True
    for stem in sorted(fixtures):
        if not args.quiet:
            print(f"  running {stem} ...", file=sys.stderr)
        res = _run_fixture(
            args.engine,
            stem,
            fixtures[stem],
            args.input,
            out_dir,
            str(data_manifest.KESTREL_JAR),
        )
        cases[stem] = res
        all_correct = all_correct and res["correct"]
        if not args.quiet:
            verdict = "OK  " if res["correct"] else "FAIL"
            called = res["called"]
            call_str = (
                f"{called['CHROM']}:{called['POS']} {called['REF']}>"
                f"{called['ALT']} conf={called['Confidence']} "
                f"alt={called['Estimated_Depth_AlternateVariant']}"
                if called
                else "no positive call"
            )
            print(
                f"  [{verdict}] {stem} "
                f"expect={res['expected']['confidence']!r} -> {call_str} "
                f"({res['wall_seconds']}s)",
                file=sys.stderr,
            )
            for reason in res["reasons"]:
                print(f"         - {reason}", file=sys.stderr)

    out = {
        "engine": args.engine,
        "input": args.input,
        "all_correct": all_correct,
        "cases": cases,
    }
    payload = json.dumps(out, indent=2, sort_keys=True)
    print(payload)
    if args.json:
        Path(args.json).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json).write_text(payload + "\n", encoding="utf-8")

    return 0 if all_correct else 1


if __name__ == "__main__":
    sys.exit(main())
