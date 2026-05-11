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
        self.assertEqual(report["coverage"]["status"], "pass")
        self.assertEqual(report["algorithm_results"]["kestrel"], "High_Precision_flagged")
        self.assertEqual(report["algorithm_results"]["advntr"], "none")
        self.assertFalse(report["cross_match_summary"]["available"])
        self.assertIn("adVNTR genotyping was not performed", report["screening_summary"])
        self.assertEqual(len(report["kestrel_variants"]), 3)

    def test_report_json_includes_optional_advntr_table_and_cross_match(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))
        advntr_rows = [
            {
                "VID": "MUC1-dupC",
                "Variant": "dupC",
                "SupportingReads": 42,
                "MeanCoverage": 80,
                "Pvalue": 0.001,
                "RU": "MUC1",
                "POS": "100",
                "REF": "C",
                "ALT": "CC",
                "Flag": "Not flagged",
            }
        ]
        report = vntyper_port.build_report_json(
            sample_name="fixture",
            input_files={"vcf": str(FIXTURE)},
            kestrel_rows=rows,
            coverage={"mean": 250},
            advntr_rows=advntr_rows,
        )
        self.assertEqual(report["algorithm_results"]["advntr"], "positive")
        self.assertEqual(report["advntr_variants"], advntr_rows)
        self.assertEqual(report["cross_match_summary"]["status"], "concordant_positive")

    def test_report_json_contains_metadata_and_fastp_qc(self):
        rows = vntyper_port.process_kestrel_vcf(str(FIXTURE))
        report = vntyper_port.build_report_json(
            sample_name="fixture",
            input_files={"bam": "fixture.bam"},
            kestrel_rows=rows,
            coverage={"mean": 10},
            fastp={
                "sequencing_setup": "paired-end",
                "duplication_rate": 0.01,
                "q20_rate": 0.99,
                "q30_rate": 0.95,
                "passed_filter_read_rate": 0.98,
                "quality_pass": True,
            },
            metadata={
                "alignment_pipeline": "external samtools/kestrel",
                "detected_assembly": "hg19",
                "detected_contig": "chr1",
                "bam_header_warnings": ["missing PG"],
                "report_date": "2026-05-11 00:00:00",
            },
        )
        self.assertEqual(report["metadata"]["detected_assembly"], "hg19")
        self.assertEqual(report["metadata"]["detected_contig"], "chr1")
        self.assertEqual(report["metadata"]["bam_header_warnings"], ["missing PG"])
        self.assertEqual(report["coverage"]["status"], "warning")
        self.assertTrue(report["fastp"]["available"])
        self.assertEqual(report["fastp"]["sequencing_setup"], "paired-end")

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
            metadata={
                "alignment_pipeline": "external samtools/kestrel",
                "detected_assembly": "hg19",
                "detected_contig": "chr1",
                "bam_header_warnings": [],
                "report_date": "2026-05-11 00:00:00",
            },
        )
        best = vntyper_port.best_kestrel_call(
            [row for row in rows if row["passes_vntyper_filters"]]
        )
        actual = {
            "screening_summary": report["screening_summary"],
            "coverage": {
                "quality_pass": report["coverage"]["quality_pass"],
                "status": report["coverage"]["status"],
                "threshold": report["coverage"]["threshold"],
            },
            "algorithm_results": report["algorithm_results"],
            "kestrel_variant_count": len(report["kestrel_variants"]),
            "metadata": {
                "vntyper_version": report["metadata"]["vntyper_version"],
                "alignment_pipeline": report["metadata"]["alignment_pipeline"],
                "detected_assembly": report["metadata"]["detected_assembly"],
                "detected_contig": report["metadata"]["detected_contig"],
                "bam_header_warnings": report["metadata"]["bam_header_warnings"],
            },
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
