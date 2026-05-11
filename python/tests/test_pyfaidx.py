from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from bioscript import pyfaidx


class PyfaidxTests(unittest.TestCase):
    def test_pyfaidx_auto_backend_has_pure_python_fallback(self) -> None:
        env = {key: value for key, value in os.environ.items() if key != "BIOSCRIPT_BACKEND"}
        with tempfile.TemporaryDirectory() as tmp, patch.dict(os.environ, env, clear=True):
            path = Path(tmp) / "ref.fa"
            path.write_text(">chr_test\nACGT\n")

            fasta = pyfaidx.Fasta(path)
            self.assertEqual(str(fasta["chr_test"][0:0]), "")
            self.assertEqual(str(fasta["chr_test"][:4]), "ACGT")
