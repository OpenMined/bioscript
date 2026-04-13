from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assay_loader import load_assay_package


class AssayLoaderTest(unittest.TestCase):
    def test_apol1_package_classifies_runnable_and_unsupported_variants(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "risk" / "APOL1"
        package = load_assay_package(assay_root)

        self.assertEqual(package.implementation_kind, "script")
        self.assertEqual(len(package.runnable_variants), 3)
        self.assertEqual(len(package.variants), 3)
        self.assertEqual(len(package.unsupported_variants), 0)

    def test_glp1_package_reports_unsupported_variants(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "pgx" / "GLP1"
        package = load_assay_package(assay_root)

        self.assertGreaterEqual(len(package.runnable_variants), 1)
        self.assertGreaterEqual(len(package.unsupported_variants), 1)
        reasons = {entry.reason for entry in package.unsupported_variants}
        self.assertIn("indels not yet supported by bioscript runtime", reasons)


if __name__ == "__main__":
    unittest.main()
