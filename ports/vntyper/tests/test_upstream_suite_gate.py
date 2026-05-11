import importlib.util
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
UPSTREAM_ROOT = ROOT / "ports" / "vntyper" / "vntyper"


class UpstreamVNtyperSuiteGateTests(unittest.TestCase):
    def test_upstream_unit_subset_runs_when_dependencies_are_installed(self):
        if importlib.util.find_spec("pytest") is None:
            self.skipTest("pytest is not installed for upstream VNtyper reference tests")
        if importlib.util.find_spec("pandas") is None:
            self.skipTest("pandas is not installed for upstream VNtyper reference tests")

        tests = [
            "tests/unit/test_scoring.py",
            "tests/unit/test_confidence_assignment.py",
            "tests/unit/test_variant_parsing.py",
            "tests/unit/test_region_utils.py",
            "tests/unit/test_reference_registry.py",
            "tests/unit/test_chromosome_utils.py",
        ]
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "-q", *tests],
            cwd=UPSTREAM_ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout)


if __name__ == "__main__":
    unittest.main()
