from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assay_intermediate import load_assay_package_intermediate


class AssayIntermediateTest(unittest.TestCase):
    def test_apol1_intermediate_contains_runnable_variants(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "risk" / "APOL1"
        intermediate = load_assay_package_intermediate(assay_root)

        self.assertEqual(intermediate["schema"], "bioscript:assay-intermediate")
        self.assertEqual(intermediate["assay"]["id"], "apol1_risk")
        self.assertEqual(intermediate["implementation"]["kind"], "script")
        self.assertEqual(len(intermediate["runnable_variants"]), 3)
        self.assertEqual(len(intermediate["unsupported_variants"]), 0)

    def test_glp1_intermediate_contains_unsupported_reason(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "pgx" / "GLP1"
        intermediate = load_assay_package_intermediate(assay_root)

        reasons = {entry["reason"] for entry in intermediate["unsupported_variants"]}
        self.assertIn("indels not yet supported by bioscript runtime", reasons)


if __name__ == "__main__":
    unittest.main()
