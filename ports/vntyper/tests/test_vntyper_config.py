import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"

sys.path.insert(0, str(BIOSCRIPT_PORT))

config_spec = importlib.util.spec_from_file_location("vntyper_config", BIOSCRIPT_PORT / "vntyper_config.py")
vntyper_config = importlib.util.module_from_spec(config_spec)
config_spec.loader.exec_module(vntyper_config)

port_spec = importlib.util.spec_from_file_location("vntyper_port", BIOSCRIPT_PORT / "vntyper_port.py")
vntyper_port = importlib.util.module_from_spec(port_spec)
sys.modules["vntyper_config"] = vntyper_config
port_spec.loader.exec_module(vntyper_port)


class VntyperConfigTests(unittest.TestCase):
    def test_muc1_regions_and_reference_paths_are_explicit(self):
        self.assertEqual(
            vntyper_config.COORDINATE_SYSTEMS["GRCh37"]["bam_region_coords"],
            "155158000-155163000",
        )
        self.assertEqual(
            vntyper_config.COORDINATE_SYSTEMS["GRCh38"]["vntr_region_coords"],
            "155188000-155192500",
        )
        self.assertIn("All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa", vntyper_config.DEFAULT_MUC1_REFERENCE)

    def test_confidence_thresholds_and_optional_validation_toggles_are_explicit(self):
        assignment = vntyper_config.DEFAULT_KESTREL_CONFIG["confidence_assignment"]

        self.assertEqual(assignment["depth_score_thresholds"]["low"], 0.00469)
        self.assertEqual(assignment["depth_score_thresholds"]["high"], 0.00515)
        self.assertEqual(assignment["alt_depth_thresholds"]["mid_high"], 100)
        self.assertFalse(vntyper_config.OPTIONAL_VALIDATION_DEFAULTS["advntr_enabled"])
        self.assertEqual(vntyper_config.OPTIONAL_VALIDATION_DEFAULTS["advntr_result_when_disabled"], "none")

    def test_report_schema_keys_match_generated_report_surface(self):
        report = vntyper_port.build_report_json("sample1", {"vcf": "output.vcf"}, [])

        self.assertEqual(set(vntyper_config.REPORT_SCHEMA_KEYS), set(report))


if __name__ == "__main__":
    unittest.main()
