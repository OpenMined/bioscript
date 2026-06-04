import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
UPSTREAM_ROOT = ROOT / "ports" / "vntyper" / "vntyper"
PORT_PATH = ROOT / "ports" / "vntyper" / "bioscript" / "vntyper_port.py"

spec = importlib.util.spec_from_file_location("vntyper_port", PORT_PATH)
vntyper_port = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vntyper_port)


def import_or_skip(module_name):
    try:
        return __import__(module_name, fromlist=["*"])
    except ModuleNotFoundError as exc:
        raise unittest.SkipTest(f"upstream VNtyper parity dependency missing: {exc.name}") from exc


class UpstreamScoringParityTests(unittest.TestCase):
    def test_scoring_confidence_and_alt_filter_subset_matches_upstream(self):
        pandas = import_or_skip("pandas")
        sys.path.insert(0, str(UPSTREAM_ROOT))
        scoring = import_or_skip("vntyper.scripts.scoring")
        confidence = import_or_skip("vntyper.scripts.confidence_assignment")
        variant_parsing = import_or_skip("vntyper.scripts.variant_parsing")

        rows = [
            {"REF": "C", "ALT": "CGGCA", "Sample": "Del:120:10000"},
            {"REF": "CGG", "ALT": "C", "Sample": "Del:21:4000"},
            {"REF": "C", "ALT": "CGG", "Sample": "Del:2:10000"},
        ]
        upstream = pandas.DataFrame(rows)
        upstream = scoring.split_depth_and_calculate_frame_score(upstream)
        upstream = scoring.split_frame_score(upstream)
        upstream = scoring.extract_frameshifts(upstream)
        upstream = confidence.calculate_depth_score_and_assign_confidence(
            upstream,
            vntyper_port.DEFAULT_KESTREL_CONFIG,
        )
        upstream = variant_parsing.filter_by_alt_values_and_finalize(
            upstream,
            vntyper_port.DEFAULT_KESTREL_CONFIG,
        )

        port = vntyper_port.split_depth_and_calculate_frame_score(rows)
        port = vntyper_port.split_frame_score(port)
        port = vntyper_port.extract_frameshifts(port)
        port = vntyper_port.calculate_depth_score_and_assign_confidence(
            port,
            vntyper_port.DEFAULT_KESTREL_CONFIG,
        )
        port = vntyper_port.filter_by_alt_values_and_finalize(
            port,
            vntyper_port.DEFAULT_KESTREL_CONFIG,
        )

        for index, port_row in enumerate(port):
            upstream_row = upstream.iloc[index]
            self.assertEqual(port_row["is_valid_frameshift"], bool(upstream_row["is_valid_frameshift"]))
            self.assertEqual(port_row["Confidence"], upstream_row["Confidence"])
            self.assertAlmostEqual(port_row["Depth_Score"], float(upstream_row["Depth_Score"]))
            self.assertEqual(port_row["alt_filter_pass"], bool(upstream_row["alt_filter_pass"]))


if __name__ == "__main__":
    unittest.main()
