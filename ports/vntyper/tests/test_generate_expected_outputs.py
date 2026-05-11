import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
GENERATOR_PATH = ROOT / "ports" / "vntyper" / "tests" / "generate_expected_outputs.py"


spec = importlib.util.spec_from_file_location("generate_expected_outputs", GENERATOR_PATH)
generate_expected_outputs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generate_expected_outputs)


class GenerateExpectedOutputsTests(unittest.TestCase):
    def test_dry_run_payload_plans_expected_layout_without_external_tools(self):
        payload = generate_expected_outputs.build_payload(
            "example_6449_hg19_subset",
            "example_66bf_hg19_subset",
            "hg19",
            "ports/vntyper/kestrel/kestrel.jar",
        )

        self.assertEqual(payload["manifest"]["positive_sample"], "example_6449_hg19_subset")
        self.assertEqual(payload["manifest"]["negative_sample"], "example_66bf_hg19_subset")
        self.assertEqual(
            payload["manifest"]["expected_outputs"],
            [
                "positive/kestrel/output.vcf",
                "positive/kestrel/kestrel_result.tsv",
                "negative/kestrel/output.vcf",
                "negative/kestrel/kestrel_result.tsv",
            ],
        )
        self.assertEqual(len(payload["samples"]), 2)
        positive = payload["samples"][0]
        self.assertEqual(positive["label"], "positive")
        self.assertTrue(positive["expected_kestrel_vcf"].endswith("positive/kestrel/output.vcf"))
        self.assertEqual(
            positive["pipeline_command_plan"]["bam_region"],
            "chr1:155158000-155163000",
        )
        self.assertIn("-sexample_6449_hg19_subset", positive["pipeline_command_plan"]["kestrel_command"])


if __name__ == "__main__":
    unittest.main()
