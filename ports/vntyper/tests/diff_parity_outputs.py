#!/usr/bin/env python3
"""Diff two `run_parity_pipeline.py` JSON outputs and exit non-zero on divergence.

Engine-identifying fields (`engine`, `alignment_pipeline`, `wall_seconds`) are
stripped from both sides before comparison — those legitimately differ for the
same input. Anything else that differs is a real parity gap and is printed
case by case.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ENGINE_LABEL_FIELDS = {"engine"}
PIPELINE_LABEL_FIELDS = {"alignment_pipeline"}
NOISY_FIELDS = {"wall_seconds"}


def _scrub(obj):
    if isinstance(obj, dict):
        scrubbed = {}
        for key, value in obj.items():
            if key in ENGINE_LABEL_FIELDS or key in NOISY_FIELDS:
                continue
            if key in PIPELINE_LABEL_FIELDS:
                continue
            scrubbed[key] = _scrub(value)
        return scrubbed
    if isinstance(obj, list):
        return [_scrub(item) for item in obj]
    return obj


def _diff_keys(prefix, left, right, hits):
    if type(left) is not type(right):
        hits.append(f"{prefix}: type {type(left).__name__} != {type(right).__name__}")
        return
    if isinstance(left, dict):
        all_keys = sorted(set(left) | set(right))
        for key in all_keys:
            if key not in left:
                hits.append(f"{prefix}.{key}: missing on left ({right[key]!r})")
            elif key not in right:
                hits.append(f"{prefix}.{key}: missing on right ({left[key]!r})")
            else:
                _diff_keys(f"{prefix}.{key}", left[key], right[key], hits)
        return
    if isinstance(left, list):
        if len(left) != len(right):
            hits.append(f"{prefix}: list length {len(left)} != {len(right)}")
            return
        for index, (lhs, rhs) in enumerate(zip(left, right)):
            _diff_keys(f"{prefix}[{index}]", lhs, rhs, hits)
        return
    if left != right:
        hits.append(f"{prefix}: {left!r} != {right!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("left", type=Path, help="first run_parity_pipeline.py JSON")
    parser.add_argument("right", type=Path, help="second run_parity_pipeline.py JSON")
    parser.add_argument(
        "--label-left",
        default=None,
        help="label for left side in output (default: engine field)",
    )
    parser.add_argument(
        "--label-right",
        default=None,
        help="label for right side in output (default: engine field)",
    )
    args = parser.parse_args()

    left_raw = json.loads(args.left.read_text(encoding="utf-8"))
    right_raw = json.loads(args.right.read_text(encoding="utf-8"))

    label_left = args.label_left or left_raw.get("engine", str(args.left))
    label_right = args.label_right or right_raw.get("engine", str(args.right))

    if left_raw.get("input") != right_raw.get("input"):
        print(
            f"input mismatch: {label_left} ran {left_raw.get('input')!r}, "
            f"{label_right} ran {right_raw.get('input')!r}",
            file=sys.stderr,
        )
        return 2

    left = _scrub(left_raw)
    right = _scrub(right_raw)

    cases = sorted(set(left.get("cases", {})) | set(right.get("cases", {})))
    if not cases:
        print("no cases to compare", file=sys.stderr)
        return 2

    any_diff = False
    for case in cases:
        case_left = left.get("cases", {}).get(case)
        case_right = right.get("cases", {}).get(case)
        if case_left is None:
            print(f"[{case}] missing on {label_left}")
            any_diff = True
            continue
        if case_right is None:
            print(f"[{case}] missing on {label_right}")
            any_diff = True
            continue
        hits = []
        _diff_keys("case", case_left, case_right, hits)
        sha_left = case_left.get("tsv_fingerprint", {}).get("sha256", "?")[:12]
        sha_right = case_right.get("tsv_fingerprint", {}).get("sha256", "?")[:12]
        cls_left = case_left.get("classification")
        cls_right = case_right.get("classification")
        rows_left = case_left.get("tsv_fingerprint", {}).get("row_count")
        rows_right = case_right.get("tsv_fingerprint", {}).get("row_count")
        if hits:
            print(
                f"[{case}] DIFF  "
                f"{label_left}: cls={cls_left} rows={rows_left} sha={sha_left}  "
                f"{label_right}: cls={cls_right} rows={rows_right} sha={sha_right}"
            )
            for hit in hits:
                print(f"    {hit}")
            any_diff = True
        else:
            print(
                f"[{case}] MATCH cls={cls_left} rows={rows_left} sha={sha_left}"
            )

    if any_diff:
        print(
            f"\nparity FAIL: {label_left} and {label_right} differ on at least one case",
            file=sys.stderr,
        )
        return 1
    print(
        f"\nparity OK: {label_left} and {label_right} agree on all cases"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
