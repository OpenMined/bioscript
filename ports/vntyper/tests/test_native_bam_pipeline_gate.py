import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"
PIPELINE_PATH = BIOSCRIPT_PORT / "vntyper_external_pipeline.py"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

manifest_spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(manifest_spec)
manifest_spec.loader.exec_module(data_manifest)

pipeline_spec = importlib.util.spec_from_file_location(
    "vntyper_external_pipeline",
    PIPELINE_PATH,
)
vntyper_external_pipeline = importlib.util.module_from_spec(pipeline_spec)
sys.modules["vntyper_external_pipeline"] = vntyper_external_pipeline
pipeline_spec.loader.exec_module(vntyper_external_pipeline)


class VntyperNativeBamPipelineGateTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_native_bam_pipeline_prerequisites()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_native_bam_pipeline_matches_expected_sample_classification(self):
        for label, bam in self.prereqs["bam_cases"].items():
            with self.subTest(label=label):
                expected_root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                with (expected_root / "report.json").open("r", encoding="utf-8") as handle:
                    expected_report = json.load(handle)

                with tempfile.TemporaryDirectory() as tmp:
                    result = vntyper_external_pipeline.run_bam_pipeline(
                        bam,
                        label,
                        str(Path(tmp) / label),
                        kestrel_jar=self.prereqs["kestrel_jar"],
                        muc1_reference=self.prereqs["muc1_reference"],
                        use_native_samtools=True,
                    )

                    with open(result.report_json, "r", encoding="utf-8") as handle:
                        actual_report = json.load(handle)
                    with open(result.kestrel_tsv, "r", encoding="utf-8", newline="") as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))

                self.assertGreater(len(rows), 0)
                self.assertEqual(
                    actual_report["algorithm_results"]["kestrel"],
                    expected_report["algorithm_results"]["kestrel"],
                )
                self.assertEqual(set(actual_report), set(expected_report))
                self.assertEqual(len(actual_report["kestrel_variants"]), len(rows))
                self.assertEqual(actual_report["screening_summary"], expected_report["screening_summary"])
                self.assertEqual(actual_report["coverage"]["status"], "pass")
                self.assertTrue(actual_report["coverage"]["quality_pass"])
                for key in [
                    "mean",
                    "median",
                    "stdev",
                    "min",
                    "max",
                    "region_length",
                    "uncovered_bases",
                    "percent_uncovered",
                ]:
                    self.assertIsNotNone(actual_report["coverage"][key])
                self.assertGreater(actual_report["coverage"]["region_length"], 0)
                self.assertIn("bam", actual_report["input_files"])
                self.assertIn("vcf", actual_report["input_files"])
                self.assertEqual(
                    actual_report["metadata"]["alignment_pipeline"],
                    "native bioscript samtools/kestrel",
                )
                self.assertEqual(actual_report["metadata"]["detected_assembly"], "hg19")

    def test_native_bam_pipeline_with_native_kestrel_matches_expected_classification(self):
        for label, bam in self.prereqs["bam_cases"].items():
            with self.subTest(label=label):
                expected_root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                with (expected_root / "report.json").open("r", encoding="utf-8") as handle:
                    expected_report = json.load(handle)

                with tempfile.TemporaryDirectory() as tmp:
                    result = vntyper_external_pipeline.run_bam_pipeline(
                        bam,
                        label,
                        str(Path(tmp) / label),
                        kestrel_jar=self.prereqs["kestrel_jar"],
                        muc1_reference=self.prereqs["muc1_reference"],
                        use_native_samtools=True,
                        use_native_kestrel=True,
                    )

                    with open(result.report_json, "r", encoding="utf-8") as handle:
                        actual_report = json.load(handle)
                    with open(result.kestrel_tsv, "r", encoding="utf-8", newline="") as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))

                self.assertGreater(len(rows), 0)
                self.assertEqual(
                    actual_report["algorithm_results"]["kestrel"],
                    expected_report["algorithm_results"]["kestrel"],
                )
                self.assertEqual(set(actual_report), set(expected_report))
                self.assertEqual(len(actual_report["kestrel_variants"]), len(rows))
                self.assertEqual(actual_report["screening_summary"], expected_report["screening_summary"])
                self.assertEqual(
                    actual_report["metadata"]["alignment_pipeline"],
                    "native bioscript samtools/kestrel",
                )
                self.assertEqual(actual_report["metadata"]["detected_assembly"], "hg19")


if __name__ == "__main__":
    unittest.main()
