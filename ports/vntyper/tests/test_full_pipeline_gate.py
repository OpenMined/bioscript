import importlib.util
import csv
import json
import os
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


spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_manifest)

pipeline_spec = importlib.util.spec_from_file_location(
    "vntyper_external_pipeline",
    PIPELINE_PATH,
)
vntyper_external_pipeline = importlib.util.module_from_spec(pipeline_spec)
sys.modules["vntyper_external_pipeline"] = vntyper_external_pipeline
pipeline_spec.loader.exec_module(vntyper_external_pipeline)


class VntyperFullPipelineGateTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_full_pipeline_prerequisites()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_full_pipeline_prerequisites_are_available(self):
        self.assertGreater(self.prereqs["manifest"]["present"], 0)
        self.assertTrue(self.prereqs["samtools"])
        self.assertTrue(self.prereqs["bcftools"])
        self.assertTrue(self.prereqs["java"])
        self.assertTrue(self.prereqs["kestrel_jar"].endswith("kestrel.jar"))
        self.assertTrue(
            self.prereqs["muc1_reference"].endswith(
                "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"
            )
        )
        self.assertGreaterEqual(len(self.prereqs["expected_outputs"]), 6)


class VntyperExternalBamPipelineGateTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_external_bam_pipeline_prerequisites()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_external_bam_pipeline_matches_expected_sample_classification(self):
        old_path = os.environ.get("PATH", "")
        os.environ["PATH"] = f"{self.prereqs['tool_path']}{os.pathsep}{old_path}"
        try:
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
                        )

                        with open(result.report_json, "r", encoding="utf-8") as handle:
                            actual_report = json.load(handle)
                        with open(result.kestrel_tsv, "r", encoding="utf-8", newline="") as handle:
                            rows = list(csv.DictReader(handle, delimiter="\t"))

                    self.assertGreater(len(rows), 0)
                    self.assertEqual(set(actual_report), set(expected_report))
                    self.assertEqual(len(actual_report["kestrel_variants"]), len(rows))
                    self.assertEqual(
                        actual_report["algorithm_results"]["kestrel"],
                        expected_report["algorithm_results"]["kestrel"],
                    )
                    self.assertEqual(
                        actual_report["metadata"]["alignment_pipeline"],
                        "external samtools/kestrel",
                    )
        finally:
            os.environ["PATH"] = old_path


if __name__ == "__main__":
    unittest.main()
