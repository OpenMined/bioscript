from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))

from assay_compiled import load_assay_package_compiled, write_assay_package_compiled


class AssayCompiledTest(unittest.TestCase):
    def test_apol1_compiled_contains_runnable_variants(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "risk" / "APOL1"
        compiled = load_assay_package_compiled(assay_root)

        self.assertEqual(compiled["schema"], "bioscript:assay-compiled")
        self.assertEqual(compiled["assay"]["id"], "apol1_risk")
        self.assertEqual(compiled["implementation"]["kind"], "script")
        self.assertEqual(len(compiled["runnable_variants"]), 3)
        self.assertEqual(len(compiled["unsupported_variants"]), 0)

    def test_glp1_compiled_contains_unsupported_reason(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "pgx" / "GLP1"
        compiled = load_assay_package_compiled(assay_root)

        reasons = {entry["reason"] for entry in compiled["unsupported_variants"]}
        self.assertIn("indels not yet supported by bioscript runtime", reasons)

    def test_write_assay_package_compiled_writes_generated_yaml(self) -> None:
        assay_root = Path(__file__).resolve().parent.parent / "assays" / "risk" / "APOL1"
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_path = Path(tmp_dir) / "assay.compiled.yaml"
            written_path = write_assay_package_compiled(assay_root, output_path)
            self.assertEqual(written_path, output_path)

            payload = yaml.safe_load(output_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "bioscript:assay-compiled")
            self.assertEqual(payload["assay"]["id"], "apol1_risk")


if __name__ == "__main__":
    unittest.main()
