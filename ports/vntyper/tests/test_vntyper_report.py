import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
PORT_PATH = BIOSCRIPT_PORT / "vntyper_port.py"
REPORT_PATH = BIOSCRIPT_PORT / "vntyper_report.py"
FIXTURE = ROOT / "ports" / "vntyper" / "tests" / "fixtures" / "kestrel_minimal.vcf"

sys.path.insert(0, str(BIOSCRIPT_PORT))

port_spec = importlib.util.spec_from_file_location("vntyper_port", PORT_PATH)
vntyper_port = importlib.util.module_from_spec(port_spec)
port_spec.loader.exec_module(vntyper_port)

report_spec = importlib.util.spec_from_file_location("vntyper_report", REPORT_PATH)
vntyper_report = importlib.util.module_from_spec(report_spec)
report_spec.loader.exec_module(vntyper_report)


class VntyperReportTests(unittest.TestCase):
    def test_html_report_contains_core_sections_without_igv(self):
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
            metadata={
                "alignment_pipeline": "external samtools/kestrel",
                "detected_assembly": "hg19",
                "detected_contig": "chr1",
                "bam_header_warnings": [],
                "report_date": "2026-05-11 00:00:00",
            },
            pipeline_log=["planned samtools view", "planned kestrel"],
        )
        html = vntyper_report.render_html_report(report)
        self.assertIn("<h2>Screening Summary</h2>", html)
        self.assertIn("<h2>Run Metadata</h2>", html)
        self.assertIn("<summary>VNTR Coverage QC</summary>", html)
        self.assertIn("<h2>Kestrel Identified Variants</h2>", html)
        self.assertIn("<summary>Pipeline Log</summary>", html)
        self.assertIn("external samtools/kestrel", html)
        self.assertIn("High_Precision*", html)
        self.assertIn("planned samtools view", html)
        self.assertIn('id="variant-search"', html)
        self.assertIn('id="show-flagged"', html)
        self.assertIn("sortVariants", html)
        self.assertIn("filterVariants", html)
        self.assertIn("confidence-high", html)
        self.assertIn('title="Not flagged"', html)
        self.assertIn("<details open>", html)


if __name__ == "__main__":
    unittest.main()
