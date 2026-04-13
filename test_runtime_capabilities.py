from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from runtime_capabilities import assess_variant_runtime_support


class RuntimeCapabilitiesTest(unittest.TestCase):
    def test_snp_is_supported(self) -> None:
        assessment = assess_variant_runtime_support({"kind": "snv", "ref": "A", "alt": "G"})
        self.assertTrue(assessment.supported)
        self.assertEqual(assessment.kind, "snp")
        self.assertIsNone(assessment.reason)

    def test_deletion_requires_deletion_length(self) -> None:
        unsupported = assess_variant_runtime_support({"kind": "deletion", "ref": "TACATG", "alt": "T"})
        self.assertFalse(unsupported.supported)
        self.assertEqual(unsupported.reason, "deletion missing deletion_length")

        supported = assess_variant_runtime_support(
            {"kind": "deletion", "ref": "TACATG", "alt": "T", "deletion_length": 5}
        )
        self.assertTrue(supported.supported)
        self.assertIsNone(supported.reason)

    def test_insertions_and_indels_are_not_yet_supported(self) -> None:
        insertion = assess_variant_runtime_support({"kind": "insertion"})
        self.assertFalse(insertion.supported)
        self.assertEqual(insertion.reason, "insertions not yet supported by bioscript runtime")

        indel = assess_variant_runtime_support({"kind": "indel"})
        self.assertFalse(indel.supported)
        self.assertEqual(indel.reason, "indels not yet supported by bioscript runtime")


if __name__ == "__main__":
    unittest.main()
