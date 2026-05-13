from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace

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

    def test_pyfaidx_rust_backend_delegates_slice_to_native_extension(self) -> None:
        calls = []

        def fetch(path: str, contig: str, start: int, stop: int) -> str:
            calls.append((path, contig, start, stop))
            return "CG"

        fake_native = SimpleNamespace(pyfaidx_fetch_native=fetch)
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ref.fa"
            path.write_text(">chr_test\nACGT\n", encoding="utf-8")
            with patch.dict("sys.modules", {"bioscript._native": fake_native}), patch.dict(
                os.environ,
                {"BIOSCRIPT_BACKEND": "rust"},
            ):
                fasta = pyfaidx.Fasta(path)
                self.assertEqual(str(fasta["chr_test"][1:3]), "CG")

        self.assertEqual(calls, [(str(path), "chr_test", 1, 3)])

    def test_pyfaidx_rust_backend_requires_native_extension(self) -> None:
        with patch.dict(os.environ, {"BIOSCRIPT_BACKEND": "rust"}), patch.dict(
            "sys.modules",
            {"bioscript._native": None},
        ):
            with self.assertRaises(NotImplementedError):
                pyfaidx.Fasta("ref.fa")
