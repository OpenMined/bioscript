"""Dependency-light VNtyper logic for the BioScript port.

This module is written as plain Python-compatible BioScript-style code: lists
and dictionaries instead of pandas DataFrames, and functions instead of classes.
It mirrors the upstream VNtyper post-processing surface that can be tested
without running samtools or Kestrel.
"""

from __future__ import annotations

import json
import re
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
    "flagging_rules": {
        "False_Positive_4bp_Insertion": "(REF == 'C') and (ALT == 'CGGCA')",
        "Low_Depth_Conserved_Motifs": "(Depth_Score < 0.4) and (Motif in ['1', '2', '3', '4', '6', '7', '8', '9'])",
    },
    "duplicate_flagging": {
        "enabled": False,
        "flag_name": "Potential_Duplicate",
        "group_by": ["REF", "ALT"],
        "sort_by": [
            {"column": "Depth_Score", "ascending": False},
            {"column": "Motifs", "ascending": True},
            {"column": "POS", "ascending": True},
        ],
    },
}

DEFAULT_REPORT_CONFIG = {
    "mean_vntr_coverage_threshold": 100,
    "algorithm_logic": {
        "kestrel": {
            "rules": [
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                        "Flag": {"operator": "==", "value": "Not flagged"},
                    },
                    "result": "High_Precision",
                },
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["Low_Precision"]},
                        "Flag": {"operator": "==", "value": "Not flagged"},
                    },
                    "result": "Low_Precision",
                },
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["High_Precision", "High_Precision*"]},
                        "Flag": {"operator": "!=", "value": "Not flagged"},
                    },
                    "result": "High_Precision_flagged",
                },
                {
                    "conditions": {
                        "Confidence": {"operator": "in", "value": ["Low_Precision"]},
                        "Flag": {"operator": "!=", "value": "Not flagged"},
                    },
                    "result": "Low_Precision_flagged",
                },
            ],
            "default": "negative",
        },
        "advntr": {
            "rules": [
                {
                    "conditions": {
                        "VID": {"operator": "!=", "value": "Negative"},
                        "Flag": {"operator": "==", "value": "Not flagged"},
                    },
                    "result": "positive",
                },
                {
                    "conditions": {
                        "Flag": {"operator": "not in", "value": ["Not flagged", "Not applicable", "None"]},
                    },
                    "result": "positive flagged",
                },
            ],
            "default": "negative",
        },
    },
    "screening_summary_default": "The screening was negative (no valid Kestrel or adVNTR data).",
    "screening_summary_rules": [
        {
            "conditions": {
                "kestrel_result": "High_Precision",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "Kestrel detected a high-precision pathogenic variant.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the result using orthogonal methods (e.g., SNaPshot, long-read sequencing).",
        },
        {
            "conditions": {
                "kestrel_result": "High_Precision",
                "advntr_result": "none",
                "quality_metrics_pass": False,
            },
            "message": "Kestrel detected a high-precision pathogenic variant with quality metrics below threshold, and adVNTR genotyping was not performed.<br>Further validation using alternative methods (e.g., SNaPshot, long-read sequencing) is strongly recommended.",
        },
        {
            "conditions": {
                "kestrel_result": "High_Precision_flagged",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "Kestrel detected a high-precision pathogenic variant with a flagged result.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the finding using orthogonal methods (e.g., SNaPshot, long-read sequencing).",
        },
        {
            "conditions": {
                "kestrel_result": "Low_Precision",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "Kestrel detected a pathogenic variant with low precision.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the result using alternative methods (e.g., SNaPshot, long-read sequencing).",
        },
        {
            "conditions": {
                "kestrel_result": "negative",
                "advntr_result": "none",
                "quality_metrics_pass": True,
            },
            "message": "No variant detected.<br>Note: adVNTR genotyping was not performed.",
        },
    ],
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
    config = kestrel_config or DEFAULT_KESTREL_CONFIG
    rows = read_vcf_without_comments(vcf_file)
    rows = split_depth_and_calculate_frame_score(rows)
    rows = split_frame_score(rows)
    rows = extract_frameshifts(rows)
    rows = calculate_depth_score_and_assign_confidence(rows, config)
    rows = filter_by_alt_values_and_finalize(rows, config)
    rows = add_flags(
        rows,
        config.get("flagging_rules", {}),
        duplicates_config=config.get("duplicate_flagging", {}),
    )
    for row in rows:
        row["passes_vntyper_filters"] = (
            bool(row.get("is_valid_frameshift"))
            and bool(row.get("depth_confidence_pass"))
            and bool(row.get("alt_filter_pass"))
        )
    return rows


def regex_match(pattern, value):
    try:
        return re.search(pattern, str(value)) is not None
    except re.error:
        return False


def evaluate_condition(row, condition):
    env = {key: _condition_value(value) for key, value in row.items()}
    env["regex_match"] = regex_match
    try:
        return bool(eval(condition, {"__builtins__": {}}, env))
    except Exception:
        return False


def add_flags(rows, flagging_rules, duplicates_config=None):
    out = []
    for row in rows:
        next_row = dict(row)
        flags = []
        for flag_name, condition in flagging_rules.items():
            if evaluate_condition(next_row, condition):
                flags.append(flag_name)
        next_row["Flag"] = ", ".join(flags) if flags else "Not flagged"
        out.append(next_row)
    return mark_potential_duplicates(out, duplicates_config or {})


def mark_potential_duplicates(rows, duplicates_config):
    if not duplicates_config.get("enabled"):
        return rows
    flag_name = duplicates_config.get("flag_name", "Potential_Duplicate")
    group_by = duplicates_config.get("group_by", [])
    sort_by = duplicates_config.get("sort_by", [])
    groups = {}
    for idx, row in enumerate(rows):
        key = tuple(row.get(column) for column in group_by)
        groups.setdefault(key, []).append(idx)
    out = [dict(row) for row in rows]
    for indexes in groups.values():
        if len(indexes) <= 1:
            continue
        ranked = sorted(indexes, key=lambda idx: _duplicate_sort_key(out[idx], sort_by))
        for duplicate_idx in ranked[1:]:
            existing = out[duplicate_idx].get("Flag", "Not flagged")
            out[duplicate_idx]["Flag"] = flag_name if existing == "Not flagged" else f"{existing}, {flag_name}"
    return out


def apply_uniform_filtering_right_motif(
    rows,
    exclude_motifs_right,
    alt_for_motif_right_gg,
    motifs_for_alt_gg,
):
    if not rows:
        return []
    filtered = [dict(row) for row in rows if row.get("Motif") not in exclude_motifs_right]
    if not filtered:
        return []
    filtered = sorted(
        filtered,
        key=lambda row: (_float(row.get("Depth_Score", 0)), _float(row.get("POS", 0))),
        reverse=True,
    )
    deduped = []
    seen = set()
    for row in filtered:
        key = (row.get("POS"), row.get("REF"), row.get("ALT"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(row)
    if any(row.get("ALT") == alt_for_motif_right_gg for row in deduped):
        gg_allowed = [
            row
            for row in deduped
            if row.get("ALT") == alt_for_motif_right_gg and row.get("Motif") in motifs_for_alt_gg
        ]
        non_gg = [row for row in deduped if row.get("ALT") != alt_for_motif_right_gg]
        return gg_allowed + non_gg
    return deduped


def build_report_json(
    sample_name,
    input_files,
    kestrel_rows,
    coverage=None,
    fastp=None,
    report_config=None,
    pipeline_version="bioscript-vntyper-port",
    metadata=None,
    advntr_rows=None,
    pipeline_log=None,
):
    config = report_config or DEFAULT_REPORT_CONFIG
    coverage_qc = build_coverage_qc(coverage or {}, config)
    fastp_qc = build_fastp_qc(fastp or {})
    advntr_rows = advntr_rows or []
    kestrel_result = compute_algorithm_result(kestrel_rows, config, "kestrel")
    advntr_result = "none" if not advntr_rows else compute_algorithm_result(advntr_rows, config, "advntr")
    screening = screening_summary_from_config(
        kestrel_result,
        advntr_result,
        coverage_qc["quality_pass"],
        config,
    )
    report_metadata = build_run_metadata(
        sample_name=sample_name,
        input_files=input_files,
        pipeline_version=pipeline_version,
        metadata=metadata or {},
    )
    return {
        "sample_name": sample_name,
        "version": pipeline_version,
        "report_date": report_metadata["report_date"],
        "metadata": report_metadata,
        "input_files": input_files,
        "coverage": coverage_qc,
        "fastp": fastp_qc,
        "algorithm_results": {
            "kestrel": kestrel_result,
            "advntr": advntr_result,
            "quality_metrics_pass": coverage_qc["quality_pass"],
        },
        "screening_summary": screening,
        "kestrel_variants": kestrel_rows,
        "advntr_variants": advntr_rows,
        "cross_match_summary": build_cross_match_summary(kestrel_result, advntr_result),
        "pipeline_log": pipeline_log or [],
    }


def write_report_json(path, report):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)


def screening_summary(kestrel_rows, quality_pass):
    config = DEFAULT_REPORT_CONFIG
    return screening_summary_from_config(
        compute_algorithm_result(kestrel_rows, config, "kestrel"),
        "none",
        quality_pass,
        config,
    )


def build_cross_match_summary(kestrel_result, advntr_result):
    if advntr_result == "none":
        return {
            "available": False,
            "status": "not_performed",
            "message": "adVNTR genotyping was not performed.",
        }
    kestrel_positive = kestrel_result not in ("negative", "none")
    advntr_positive = advntr_result in ("positive", "positive flagged")
    if kestrel_positive and advntr_positive:
        status = "concordant_positive"
        message = "Kestrel and adVNTR both detected a pathogenic signal."
    elif not kestrel_positive and not advntr_positive:
        status = "concordant_negative"
        message = "Kestrel and adVNTR were both negative."
    elif kestrel_positive:
        status = "kestrel_only"
        message = "Kestrel detected a pathogenic signal that adVNTR did not confirm."
    else:
        status = "advntr_only"
        message = "adVNTR detected a pathogenic signal that Kestrel did not detect."
    return {
        "available": True,
        "status": status,
        "message": message,
    }


def build_run_metadata(sample_name, input_files, pipeline_version, metadata=None):
    metadata = metadata or {}
    return {
        "sample_name": sample_name,
        "vntyper_version": metadata.get("vntyper_version", pipeline_version),
        "report_date": metadata.get("report_date", datetime.now().strftime("%Y-%m-%d %H:%M:%S")),
        "input_files": input_files,
        "alignment_pipeline": metadata.get("alignment_pipeline"),
        "detected_assembly": metadata.get("detected_assembly"),
        "detected_contig": metadata.get("detected_contig"),
        "bam_header_warnings": metadata.get("bam_header_warnings", []),
    }


def build_coverage_qc(coverage, report_config=None):
    config = report_config or DEFAULT_REPORT_CONFIG
    mean_cov = coverage.get("mean")
    threshold = config.get("mean_vntr_coverage_threshold", 100)
    quality_pass = mean_cov is None or float(mean_cov) >= float(threshold)
    return {
        "mean": mean_cov,
        "median": coverage.get("median"),
        "stdev": coverage.get("stdev"),
        "min": coverage.get("min"),
        "max": coverage.get("max"),
        "region_length": coverage.get("region_length"),
        "uncovered_bases": coverage.get("uncovered_bases"),
        "percent_uncovered": coverage.get("percent_uncovered"),
        "threshold": threshold,
        "quality_pass": quality_pass,
        "status": "pass" if quality_pass else "warning",
    }


def build_fastp_qc(fastp):
    if not fastp:
        return {"available": False}
    return {
        "available": True,
        "sequencing_setup": fastp.get("sequencing_setup"),
        "duplication_rate": fastp.get("duplication_rate"),
        "q20_rate": fastp.get("q20_rate"),
        "q30_rate": fastp.get("q30_rate"),
        "passed_filter_read_rate": fastp.get("passed_filter_read_rate"),
        "quality_pass": fastp.get("quality_pass"),
        "status": "pass" if fastp.get("quality_pass", True) else "warning",
    }


def compute_algorithm_result(rows, report_config=None, algorithm="kestrel"):
    config = report_config or DEFAULT_REPORT_CONFIG
    logic = config.get("algorithm_logic", {}).get(algorithm, {})
    default = logic.get("default", "negative")
    for row in rows:
        for rule in logic.get("rules", []):
            if all(_condition_matches(row, field, condition) for field, condition in rule.get("conditions", {}).items()):
                return rule.get("result", default)
    return default


def screening_summary_from_config(kestrel_result, advntr_result, quality_metrics_pass, report_config=None):
    config = report_config or DEFAULT_REPORT_CONFIG
    context = {
        "kestrel_result": kestrel_result,
        "advntr_result": advntr_result,
        "quality_metrics_pass": quality_metrics_pass,
    }
    for rule in config.get("screening_summary_rules", []):
        if rule.get("conditions", {}) == context:
            return rule.get("message", config.get("screening_summary_default", ""))
    return config.get("screening_summary_default", "")


def _condition_matches(row, field, condition):
    if not isinstance(condition, dict):
        return row.get(field) == condition
    operator = condition.get("operator", "==")
    expected = condition.get("value")
    actual = row.get(field)
    if operator == "==":
        return actual == expected
    if operator == "!=":
        return actual != expected
    if operator == "in":
        return actual in expected
    if operator == "not in":
        return actual not in expected
    raise ValueError(f"Unsupported condition operator: {operator}")


def best_kestrel_call(rows):
    if not rows:
        return None
    return sorted(rows, key=lambda row: _float(row.get("Depth_Score", 0)), reverse=True)[0]


def _condition_value(value):
    if value is None or value == "":
        return None
    return value


def _duplicate_sort_key(row, sort_by):
    key = []
    for spec in sort_by:
        value = row.get(spec.get("column"))
        if spec.get("ascending", True):
            key.append(value)
        else:
            key.append(_reverse_sort_value(value))
    return tuple(key)


def _reverse_sort_value(value):
    try:
        return -float(value)
    except (TypeError, ValueError):
        return "".join(chr(255 - ord(char)) for char in str(value))


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
