"""Shared normalized parity helpers for VNtyper large-data gates."""

from __future__ import annotations

import hashlib


def normalized_tsv_fingerprint(rows):
    stable_fields = [
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
    digest = hashlib.sha256()
    normalized_rows = [
        tuple(str(row.get(field, "")) for field in stable_fields)
        for row in rows
    ]
    for row in sorted(normalized_rows):
        digest.update(
            "\t".join(row).encode("utf-8")
        )
        digest.update(b"\n")
    return {
        "row_count": len(rows),
        "passing_count": len(
            [row for row in rows if row.get("passes_vntyper_filters") in ("True", True)]
        ),
        "non_negative_confidence_count": len(
            [row for row in rows if row.get("Confidence") != "Negative"]
        ),
        "sha256": digest.hexdigest(),
    }


def normalized_report_summary(report):
    return {
        "algorithm_results": report.get("algorithm_results"),
        "screening_summary": report.get("screening_summary"),
        "kestrel_variant_count": len(report.get("kestrel_variants", [])),
        "coverage_status": report.get("coverage", {}).get("status"),
        "quality_pass": report.get("coverage", {}).get("quality_pass"),
        "alignment_pipeline": report.get("metadata", {}).get("alignment_pipeline"),
        "detected_assembly": report.get("metadata", {}).get("detected_assembly"),
    }


def parity_context(actual_rows, expected_rows, actual_report, expected_report):
    passing_rows = [
        row for row in actual_rows if row.get("passes_vntyper_filters") in ("True", True)
    ]
    expected_passing_rows = [
        row for row in expected_rows if row.get("passes_vntyper_filters") in ("True", True)
    ]
    top_passing = sorted(
        passing_rows,
        key=lambda row: float(row.get("Depth_Score") or 0),
        reverse=True,
    )[:5]
    return {
        "actual_row_count": len(actual_rows),
        "expected_row_count": len(expected_rows),
        "actual_passing_count": len(passing_rows),
        "expected_passing_count": len(expected_passing_rows),
        "top_passing": [
            {
                "CHROM": row.get("CHROM"),
                "POS": row.get("POS"),
                "REF": row.get("REF"),
                "ALT": row.get("ALT"),
                "Depth_Score": row.get("Depth_Score"),
                "Confidence": row.get("Confidence"),
            }
            for row in top_passing
        ],
        "actual_tsv_fingerprint": normalized_tsv_fingerprint(actual_rows),
        "expected_tsv_fingerprint": normalized_tsv_fingerprint(expected_rows),
        "actual_report_summary": normalized_report_summary(actual_report),
        "expected_report_summary": normalized_report_summary(expected_report),
    }
