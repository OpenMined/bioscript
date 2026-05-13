from __future__ import annotations

import os
import unittest
from unittest.mock import patch

from bioscript import bcftools, kestrel, pyfaidx, pysam, samtools
from bioscript.runtime import BackendMode, ModuleBackendPolicy, selected_backend


class BackendPolicyTests(unittest.TestCase):
    def test_backend_policy_is_explicit_for_each_module(self) -> None:
        modules = [bcftools, kestrel, pyfaidx, pysam, samtools]

        for module in modules:
            with self.subTest(module=module.__name__):
                policy = module.BACKEND_POLICY
                self.assertIsInstance(policy, ModuleBackendPolicy)
                self.assertTrue(policy.auto)
                self.assertTrue(policy.python)
                self.assertTrue(policy.rust)

        self.assertIn("bcftools-rs", bcftools.BACKEND_POLICY.rust)
        self.assertIn("kestrel-rs", kestrel.BACKEND_POLICY.rust)
        self.assertIn("samtools-rs", samtools.BACKEND_POLICY.rust)
        self.assertIn("real pysam", pysam.BACKEND_POLICY.python)
        self.assertIn("pure Python FASTA fallback", pyfaidx.BACKEND_POLICY.auto)

    def test_selected_backend_reports_invalid_values(self) -> None:
        with patch.dict(os.environ, {"BIOSCRIPT_BACKEND": "bad"}):
            with self.assertRaisesRegex(ValueError, "auto, python, rust"):
                selected_backend()

    def test_selected_backend_defaults_to_auto(self) -> None:
        env = {key: value for key, value in os.environ.items() if key != "BIOSCRIPT_BACKEND"}
        with patch.dict(os.environ, env, clear=True):
            self.assertEqual(selected_backend(), BackendMode.AUTO)


if __name__ == "__main__":
    unittest.main()
