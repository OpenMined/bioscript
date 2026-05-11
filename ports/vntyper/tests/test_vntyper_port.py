import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PORT_PATH = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper_port.py"
FIXTURE = ROOT / "ports" / "vntyper" / "tests" / "fixtures" / "kestrel_minimal.vcf"


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


if __name__ == "__main__":
    unittest.main()
