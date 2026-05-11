import importlib.util
import csv
import json
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PORT_PATH = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper_port.py"
FIXTURE = ROOT / "ports" / "vntyper" / "tests" / "fixtures" / "kestrel_minimal.vcf"
EXPECTED_TSV = ROOT / "ports" / "vntyper" / "tests" / "fixtures" / "kestrel_minimal_expected.tsv"
EXPECTED_REPORT = ROOT / "ports" / "vntyper" / "tests" / "fixtures" / "kestrel_minimal_expected_report.json"


spec = importlib.util.spec_from_file_location("vntyper_port", PORT_PATH)
vntyper_port = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vntyper_port)


class VntyperPortTests(unittest.TestCase):
    def test_process_kestrel_vcf_marks_expected_filters(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))

        self.assertEqual(len(rows), 3)
        self.assertTrue(rows[0]["is_valid_frameshift"])
        self.assertEqual(rows[0]["Confidence"], "High_Precision*")
        self.assertTrue(rows[0]["passes_vntyper_filters"])

        self.assertTrue(rows[1]["is_valid_frameshift"])
        self.assertEqual(rows[1]["Confidence"], "Low_Precision")
        self.assertTrue(rows[1]["passes_vntyper_filters"])

        self.assertFalse(rows[2]["is_valid_frameshift"])
        self.assertEqual(rows[2]["Confidence"], "Negative")
        self.assertFalse(rows[2]["passes_vntyper_filters"])

    def test_best_kestrel_call_uses_depth_score(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))
        passing = [row for row in rows if row["passes_vntyper_filters"]]
        best = vntyper_port.best_kestrel_call(passing)
        self.assertEqual(best["POS"], "100")
        self.assertEqual(best["Depth_Score"], 0.012)

    def test_report_json_contains_core_ui_fields(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))
        report = vntyper_port.build_report_json(
            sample_name="fixture",
            input_files={"vcf": str(FIXTURE)},
            kestrel_rows=rows,
            coverage={
                "mean": 250,
                "median": 240,
                "stdev": 12,
                "min": 210,
                "max": 280,
                "region_length": 1500,
                "uncovered_bases": 0,
                "percent_uncovered": 0,
            },
        )

        self.assertEqual(report["sample_name"], "fixture")
        self.assertTrue(report["coverage"]["quality_pass"])
        self.assertIn("high-precision pathogenic variant", report["screening_summary"])
        self.assertEqual(len(report["kestrel_variants"]), 3)

    def test_kestrel_fixture_matches_expected_tsv_rows(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))
        columns = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "Estimated_Depth_AlternateVariant",
            "Estimated_Depth_Variant_ActiveRegion",
            "Depth_Score",
            "Confidence",
            "is_valid_frameshift",
            "alt_filter_pass",
            "passes_vntyper_filters",
        ]
        actual = [{column: str(row[column]) for column in columns} for row in rows]
        with EXPECTED_TSV.open("r", encoding="utf-8", newline="") as handle:
            expected = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(actual, expected)

    def test_kestrel_fixture_matches_expected_report_summary(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))
        report = vntyper_port.build_report_json(
            sample_name="fixture",
            input_files={"vcf": str(FIXTURE)},
            kestrel_rows=rows,
            coverage={"mean": 250},
        )
        best = vntyper_port.best_kestrel_call(
            [row for row in rows if row["passes_vntyper_filters"]]
        )
        actual = {
            "screening_summary": report["screening_summary"],
            "coverage_quality_pass": report["coverage"]["quality_pass"],
            "kestrel_variant_count": len(report["kestrel_variants"]),
            "best_call": {
                "CHROM": best["CHROM"],
                "POS": best["POS"],
                "REF": best["REF"],
                "ALT": best["ALT"],
                "Estimated_Depth_AlternateVariant": best["Estimated_Depth_AlternateVariant"],
                "Estimated_Depth_Variant_ActiveRegion": best["Estimated_Depth_Variant_ActiveRegion"],
                "Depth_Score": best["Depth_Score"],
                "Confidence": best["Confidence"],
                "passes_vntyper_filters": best["passes_vntyper_filters"],
            },
        }
        with EXPECTED_REPORT.open("r", encoding="utf-8") as handle:
            expected = json.load(handle)
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
