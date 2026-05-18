from __future__ import annotations

import importlib.util
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bioscript import pyfaidx


class RealLibraryParityTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("pyfaidx"), "real pyfaidx is not installed")
    def test_pyfaidx_slice_matches_real_library_when_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ref.fa"
            path.write_text(">chr_test\nACGT\n")
            env = {
                key: value for key, value in os.environ.items() if key != "BIOSCRIPT_BACKEND"
            }
            with patch.dict(os.environ, {**env, "BIOSCRIPT_BACKEND": "python"}, clear=True):
                real_result = str(pyfaidx.Fasta(path)["chr_test"][:4])
            with patch.dict(os.environ, env, clear=True):
                shim_result = str(pyfaidx.Fasta(path)["chr_test"][:4])
            self.assertEqual(shim_result, real_result)

    @unittest.skipUnless(importlib.util.find_spec("pysam"), "real pysam is not installed")
    def test_pysam_real_library_available_for_future_alignment_parity(self) -> None:
        self.assertIsNotNone(importlib.util.find_spec("pysam"))
