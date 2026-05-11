import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PORT_PATH = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper_port.py"

spec = importlib.util.spec_from_file_location("vntyper_port", PORT_PATH)
vntyper_port = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vntyper_port)


class PortedUpstreamUnitTests(unittest.TestCase):
    def test_scoring_marks_non_frameshift_but_retains_row(self):
        rows = [
            {
                "Sample": "Del:10:100",
                "REF": "ATG",
                "ALT": "ATGATG",
                "POS": 123,
            }
        ]
        out = vntyper_port.split_depth_and_calculate_frame_score(rows)
        self.assertEqual(len(out), 1)
        self.assertFalse(out[0]["is_frameshift"])
        self.assertEqual(out[0]["Frame_Score"], 1.0)

    def test_scoring_splits_frame_direction_and_amount(self):
        rows = [
            {"Frame_Score": 1.0, "ref_len": 3, "alt_len": 4},
            {"Frame_Score": -2.0, "ref_len": 6, "alt_len": 4},
        ]
        out = vntyper_port.split_frame_score(rows)
        self.assertEqual(out[0]["direction"], 1)
        self.assertEqual(out[0]["frameshift_amount"], 1)
        self.assertEqual(out[1]["direction"], -1)
        self.assertEqual(out[1]["frameshift_amount"], 2)

    def test_extract_frameshifts_marks_upstream_patterns(self):
        rows = [
            {"direction": 1, "frameshift_amount": 1, "Variant": "ins_ok"},
            {"direction": 1, "frameshift_amount": 2, "Variant": "ins_wrong"},
            {"direction": -1, "frameshift_amount": 2, "Variant": "del_ok"},
            {"direction": -1, "frameshift_amount": 1, "Variant": "del_wrong"},
        ]
        out = vntyper_port.extract_frameshifts(rows)
        self.assertEqual([row["is_valid_frameshift"] for row in out], [True, False, True, False])

    def test_confidence_threshold_boundaries(self):
        low = vntyper_port.DEFAULT_KESTREL_CONFIG["confidence_assignment"]["depth_score_thresholds"]["low"]
        below = [{"Estimated_Depth_AlternateVariant": low * 10000 * 0.5, "Estimated_Depth_Variant_ActiveRegion": 10000}]
        at_threshold = [{"Estimated_Depth_AlternateVariant": low * 10000, "Estimated_Depth_Variant_ActiveRegion": 10000}]

        below_out = vntyper_port.calculate_depth_score_and_assign_confidence(below)
        self.assertEqual(below_out[0]["Confidence"], vntyper_port.NEGATIVE_LABEL)
        self.assertFalse(below_out[0]["depth_confidence_pass"])

        threshold_out = vntyper_port.calculate_depth_score_and_assign_confidence(at_threshold)
        self.assertNotEqual(threshold_out[0]["Confidence"], vntyper_port.NEGATIVE_LABEL)
        self.assertTrue(threshold_out[0]["depth_confidence_pass"])

    def test_confidence_high_precision_star(self):
        conf = vntyper_port.DEFAULT_KESTREL_CONFIG["confidence_assignment"]
        high = conf["depth_score_thresholds"]["high"]
        alt_mid_high = conf["alt_depth_thresholds"]["mid_high"]
        rows = [
            {
                "Estimated_Depth_AlternateVariant": alt_mid_high,
                "Estimated_Depth_Variant_ActiveRegion": int(alt_mid_high / high),
            }
        ]
        out = vntyper_port.calculate_depth_score_and_assign_confidence(rows)
        self.assertEqual(out[0]["Confidence"], "High_Precision*")

    def test_alt_filtering_matches_upstream_gg_and_exclude_rules(self):
        config = {
            "alt_filtering": {
                "gg_alt_value": "GG",
                "gg_depth_score_threshold": 0.02,
                "exclude_alts": ["BAD_ALT", "ZZZ"],
            }
        }
        rows = [
            {"ALT": "GG", "Depth_Score": 0.019},
            {"ALT": "GG", "Depth_Score": 0.02},
            {"ALT": "XYZ", "Depth_Score": 0.5},
            {"ALT": "BAD_ALT", "Depth_Score": 0.5},
        ]
        out = vntyper_port.filter_by_alt_values_and_finalize(rows, config)
        self.assertEqual([row["alt_filter_pass"] for row in out], [False, True, True, False])


if __name__ == "__main__":
    unittest.main()
