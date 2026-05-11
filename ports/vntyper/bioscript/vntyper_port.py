"""Dependency-light VNtyper logic for the BioScript port.

This module is written as plain Python-compatible BioScript-style code: lists
and dictionaries instead of pandas DataFrames, and functions instead of classes.
It mirrors the upstream VNtyper post-processing surface that can be tested
without running samtools or Kestrel.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path


NEGATIVE_LABEL = "Negative"

DEFAULT_KESTREL_CONFIG = {
    "confidence_assignment": {
        "depth_score_thresholds": {
            "low": 0.00469,
            "high": 0.00515,
        },
        "alt_depth_thresholds": {
            "low": 20,
            "mid_low": 21,
            "mid_high": 100,
        },
        "var_active_region_threshold": 200,
        "confidence_levels": {
            "low_precision": "Low_Precision",
            "high_precision": "High_Precision",
            "high_precision_star": "High_Precision*",
        },
    },
    "alt_filtering": {
        "gg_alt_value": "GG",
        "gg_depth_score_threshold": 0.00469,
        "exclude_alts": [],
    },
}

DEFAULT_REPORT_CONFIG = {
    "mean_vntr_coverage_threshold": 100,
}


def read_vcf_without_comments(vcf_file):
    rows = []
    header = None
    with open(vcf_file, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            values = line.split("\t")
            row = {}
            for idx, key in enumerate(header):
                row[key] = values[idx] if idx < len(values) else ""
            if "SAMPLE" in row and "Sample" not in row:
                row["Sample"] = row["SAMPLE"]
            rows.append(row)
    return rows


def split_depth_and_calculate_frame_score(rows):
    out = []
    for row in rows:
        next_row = dict(row)
        sample = str(next_row.get("Sample", ""))
        parts = sample.split(":")
        next_row["Del"] = parts[0] if len(parts) > 0 else ""
        next_row["Estimated_Depth_AlternateVariant"] = parts[1] if len(parts) > 1 else "0"
        next_row["Estimated_Depth_Variant_ActiveRegion"] = parts[2] if len(parts) > 2 else "0"
        ref_len = len(str(next_row.get("REF", "")))
        alt_len = len(str(next_row.get("ALT", "")))
        delta = alt_len - ref_len
        next_row["ref_len"] = ref_len
        next_row["alt_len"] = alt_len
        next_row["Frame_Score"] = delta / 3
        next_row["is_frameshift"] = delta % 3 != 0
        out.append(next_row)
    return out


def split_frame_score(rows):
    out = []
    for row in rows:
        next_row = dict(row)
        delta = int(next_row.get("alt_len", 0)) - int(next_row.get("ref_len", 0))
        if delta > 0:
            direction = 1
        elif delta < 0:
            direction = -1
        else:
            direction = 0
        next_row["direction"] = direction
        next_row["frameshift_amount"] = abs(delta) % 3
        out.append(next_row)
    return out


def extract_frameshifts(rows):
    out = []
    for row in rows:
        next_row = dict(row)
        direction = int(next_row.get("direction", 0))
        amount = int(next_row.get("frameshift_amount", 0))
        insertion = direction > 0 and amount == 1
        deletion = direction < 0 and amount == 2
        next_row["is_valid_frameshift"] = insertion or deletion
        out.append(next_row)
    return out


def calculate_depth_score_and_assign_confidence(rows, kestrel_config=None):
    config = kestrel_config or DEFAULT_KESTREL_CONFIG
    assignment = config.get("confidence_assignment", {})
    score_thresholds = assignment.get("depth_score_thresholds", {})
    alt_thresholds = assignment.get("alt_depth_thresholds", {})
    levels = assignment.get("confidence_levels", {})

    low_threshold = float(score_thresholds.get("low", 0.2))
    high_threshold = float(score_thresholds.get("high", 0.4))
    var_region_threshold = float(assignment.get("var_active_region_threshold", 0))
    alt_low = float(alt_thresholds.get("low", 5))
    alt_mid_low = float(alt_thresholds.get("mid_low", 10))
    alt_mid_high = float(alt_thresholds.get("mid_high", 20))

    low_precision = levels.get("low_precision", "Low_Precision")
    high_precision = levels.get("high_precision", "High_Precision")
    high_precision_star = levels.get("high_precision_star", "High_Precision*")

    out = []
    for row in rows:
        next_row = dict(row)
        alt_depth = _float(next_row.get("Estimated_Depth_AlternateVariant", 0))
        region_depth = _float(next_row.get("Estimated_Depth_Variant_ActiveRegion", 0))
        depth_score = alt_depth / region_depth if region_depth != 0 else None
        next_row["Estimated_Depth_AlternateVariant"] = alt_depth
        next_row["Estimated_Depth_Variant_ActiveRegion"] = region_depth
        next_row["Depth_Score"] = depth_score

        confidence = NEGATIVE_LABEL
        if depth_score is not None and depth_score >= low_threshold:
            if region_depth <= var_region_threshold or depth_score == low_threshold:
                confidence = low_precision
            if alt_depth >= alt_mid_high and depth_score >= high_threshold:
                confidence = high_precision_star
            if alt_mid_low <= alt_depth < alt_mid_high and low_threshold <= depth_score <= high_threshold:
                confidence = low_precision
            if alt_depth <= alt_low:
                confidence = low_precision
            if alt_mid_low <= alt_depth < alt_mid_high and depth_score >= high_threshold:
                confidence = high_precision
            if low_threshold < depth_score < high_threshold:
                confidence = low_precision

        next_row["Confidence"] = confidence
        next_row["depth_confidence_pass"] = confidence != NEGATIVE_LABEL
        out.append(next_row)
    return out


def filter_by_alt_values_and_finalize(rows, kestrel_config=None):
    config = kestrel_config or DEFAULT_KESTREL_CONFIG
    alt_filter = config.get("alt_filtering", {})
    gg_alt_value = alt_filter.get("gg_alt_value", "GG")
    gg_depth_threshold = float(alt_filter.get("gg_depth_score_threshold", 0.0))
    exclude_alts = alt_filter.get("exclude_alts", [])

    out = []
    for row in rows:
        if "ALT" not in row or "Depth_Score" not in row:
            raise KeyError("Missing required columns: {'ALT', 'Depth_Score'}")
        next_row = dict(row)
        alt = next_row.get("ALT")
        depth_score = _float(next_row.get("Depth_Score", 0))
        is_gg = alt == gg_alt_value
        next_row["alt_filter_pass"] = (not is_gg or depth_score >= gg_depth_threshold) and alt not in exclude_alts
        out.append(next_row)
    return out


def process_kestrel_vcf(vcf_file, kestrel_config=None):
    rows = read_vcf_without_comments(vcf_file)
    rows = split_depth_and_calculate_frame_score(rows)
    rows = split_frame_score(rows)
    rows = extract_frameshifts(rows)
    rows = calculate_depth_score_and_assign_confidence(rows, kestrel_config)
    rows = filter_by_alt_values_and_finalize(rows, kestrel_config)
    for row in rows:
        row["passes_vntyper_filters"] = (
            bool(row.get("is_valid_frameshift"))
            and bool(row.get("depth_confidence_pass"))
            and bool(row.get("alt_filter_pass"))
        )
        row.setdefault("Flag", "Not flagged")
    return rows


def build_report_json(
    sample_name,
    input_files,
    kestrel_rows,
    coverage=None,
    fastp=None,
    report_config=None,
    pipeline_version="bioscript-vntyper-port",
):
    config = report_config or DEFAULT_REPORT_CONFIG
    coverage = coverage or {}
    fastp = fastp or {}
    mean_cov = coverage.get("mean")
    threshold = config.get("mean_vntr_coverage_threshold", 100)
    quality_pass = mean_cov is None or float(mean_cov) >= float(threshold)
    screening = screening_summary(kestrel_rows, quality_pass)
    return {
        "sample_name": sample_name,
        "version": pipeline_version,
        "report_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "input_files": input_files,
        "coverage": {
            "mean": mean_cov,
            "median": coverage.get("median"),
            "stdev": coverage.get("stdev"),
            "min": coverage.get("min"),
            "max": coverage.get("max"),
            "region_length": coverage.get("region_length"),
            "uncovered_bases": coverage.get("uncovered_bases"),
            "percent_uncovered": coverage.get("percent_uncovered"),
            "quality_pass": quality_pass,
        },
        "fastp": fastp,
        "screening_summary": screening,
        "kestrel_variants": kestrel_rows,
    }


def write_report_json(path, report):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)


def screening_summary(kestrel_rows, quality_pass):
    candidates = [row for row in kestrel_rows if row.get("passes_vntyper_filters")]
    if not candidates:
        if quality_pass:
            return "No variant detected. Quality metrics are acceptable."
        return "No variant detected; however, quality metrics are below threshold."
    best = best_kestrel_call(candidates)
    confidence = best.get("Confidence", NEGATIVE_LABEL)
    flagged = best.get("Flag", "Not flagged") != "Not flagged"
    if confidence in ["High_Precision", "High_Precision*"]:
        if flagged:
            return "Kestrel detected a high-precision pathogenic variant with a flagged result."
        if quality_pass:
            return "Kestrel detected a high-precision pathogenic variant."
        return "Kestrel detected a high-precision pathogenic variant with quality metrics below threshold."
    if confidence == "Low_Precision":
        if flagged:
            return "Kestrel detected a pathogenic variant with low precision and a flagged result."
        return "Kestrel detected a pathogenic variant with low precision."
    return "No variant detected."


def best_kestrel_call(rows):
    if not rows:
        return None
    return sorted(rows, key=lambda row: _float(row.get("Depth_Score", 0)), reverse=True)[0]


def _float(value):
    if value is None or value == "":
        return 0.0
    return float(value)


def main():
    # Placeholder CLI for local smoke checks. The BioScript runtime entry point
    # can call these same functions once local module imports are available.
    fixture = Path(__file__).parents[1] / "tests" / "fixtures" / "kestrel_minimal.vcf"
    rows = process_kestrel_vcf(str(fixture))
    report = build_report_json("fixture", {"vcf": str(fixture)}, rows)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
