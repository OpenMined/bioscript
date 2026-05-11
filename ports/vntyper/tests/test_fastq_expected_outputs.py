import csv
import importlib.util
import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"

sys.path.insert(0, str(BIOSCRIPT_PORT))


spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_manifest)

import vntyper_port


class VntyperFastqExpectedOutputsTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_fastq_kestrel_expected_outputs()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_fastq_kestrel_outputs_are_parseable_for_representative_samples(self):
        for label in ["positive", "negative"]:
            with self.subTest(label=label):
                root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                vcf = root / "kestrel" / "output.vcf"
                tsv = root / "kestrel" / "kestrel_result.tsv"
                report_json = root / "report.json"

                self.assertGreater(vcf.stat().st_size, 0)
                with tsv.open("r", encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                with report_json.open("r", encoding="utf-8") as handle:
                    report = json.load(handle)

                self.assertGreater(len(rows), 0)
                self.assertEqual(len(report["kestrel_variants"]), len(rows))
                self.assertIn(
                    report["algorithm_results"]["kestrel"],
                    ["negative", "Low_Precision", "High_Precision", "High_Precision_flagged"],
                )
                self.assertEqual(report["metadata"]["alignment_pipeline"], "external kestrel from FASTQ")

    def test_reprocessed_java_kestrel_vcf_matches_expected_classification(self):
        for label in ["positive", "negative"]:
            with self.subTest(label=label):
                root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                rows = vntyper_port.process_kestrel_vcf(str(root / "kestrel" / "output.vcf"))
                with (root / "report.json").open("r", encoding="utf-8") as handle:
                    report = json.load(handle)
                rebuilt = vntyper_port.build_report_json(
                    sample_name=report["sample_name"],
                    input_files=report["input_files"],
                    kestrel_rows=rows,
                    metadata=report["metadata"],
                )

                self.assertEqual(
                    rebuilt["algorithm_results"]["kestrel"],
                    report["algorithm_results"]["kestrel"],
                )


if __name__ == "__main__":
    unittest.main()
