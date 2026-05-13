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


class VntyperNativeFastqPipelineGateTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_native_fastq_pipeline_prerequisites()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_native_fastq_pipeline_with_native_kestrel_and_bcftools_matches_expected_classification(self):
        for label, (fastq_1, fastq_2) in self.prereqs["fastq_cases"].items():
            with self.subTest(label=label):
                expected_root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                with (expected_root / "report.json").open("r", encoding="utf-8") as handle:
                    expected_report = json.load(handle)

                with tempfile.TemporaryDirectory() as tmp:
                    result = vntyper_external_pipeline.run_fastq_kestrel(
                        fastq_1,
                        fastq_2,
                        label,
                        str(Path(tmp) / label),
                        assembly="hg19",
                        muc1_reference=self.prereqs["muc1_reference"],
                        use_native_kestrel=True,
                        use_native_bcftools=True,
                    )

                    with open(result.report_json, "r", encoding="utf-8") as handle:
                        actual_report = json.load(handle)
                    with open(result.kestrel_tsv, "r", encoding="utf-8", newline="") as handle:
                        rows = list(csv.DictReader(handle, delimiter="\t"))

                    sorted_vcf = Path(actual_report["input_files"]["sorted_vcf"])
                    sorted_vcf_index = Path(f"{sorted_vcf}.csi")

                    self.assertTrue(sorted_vcf.exists())
                    self.assertTrue(sorted_vcf_index.exists())

                self.assertGreater(len(rows), 0)
                self.assertEqual(
                    actual_report["algorithm_results"]["kestrel"],
                    expected_report["algorithm_results"]["kestrel"],
                )
                self.assertEqual(set(actual_report), set(expected_report))
                self.assertEqual(len(actual_report["kestrel_variants"]), len(rows))
                self.assertEqual(
                    actual_report["screening_summary"],
                    expected_report["screening_summary"],
                )
                self.assertEqual(
                    actual_report["metadata"]["alignment_pipeline"],
                    "native bioscript kestrel from FASTQ",
                )
                self.assertEqual(actual_report["metadata"]["detected_assembly"], "hg19")


if __name__ == "__main__":
    unittest.main()
