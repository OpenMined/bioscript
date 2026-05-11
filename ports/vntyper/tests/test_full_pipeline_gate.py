import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"


spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_manifest)


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


if __name__ == "__main__":
    unittest.main()
