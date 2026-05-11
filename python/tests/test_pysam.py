from __future__ import annotations

import os
import unittest
from unittest.mock import patch

from bioscript import pysam


class PysamTests(unittest.TestCase):
    def test_pysam_rust_backend_reports_pending_native_extension(self) -> None:
        with patch.dict(os.environ, {"BIOSCRIPT_BACKEND": "rust"}):
            with self.assertRaises(NotImplementedError):
                pysam.AlignmentFile("sample.cram", "rc")
