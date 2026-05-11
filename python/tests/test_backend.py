from __future__ import annotations

import os
import unittest
from unittest.mock import patch

from bioscript.runtime import BackendMode, selected_backend


class BackendTests(unittest.TestCase):
    def test_selected_backend_defaults_to_auto(self) -> None:
        env = {key: value for key, value in os.environ.items() if key != "BIOSCRIPT_BACKEND"}
        with patch.dict(os.environ, env, clear=True):
            self.assertEqual(selected_backend(), BackendMode.AUTO)

    def test_selected_backend_rejects_unknown_value(self) -> None:
        with patch.dict(os.environ, {"BIOSCRIPT_BACKEND": "bad"}):
            with self.assertRaises(ValueError):
                selected_backend()
