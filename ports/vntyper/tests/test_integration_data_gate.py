import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"


spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_manifest)


class VntyperIntegrationDataGateTests(unittest.TestCase):
    def setUp(self):
        try:
            self.manifest = data_manifest.require_test_data(check_md5=False)
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_large_data_manifest_is_available_for_integration_tests(self):
        self.assertGreater(self.manifest["present"], 0)
        self.assertEqual(self.manifest["missing"], [])


if __name__ == "__main__":
    unittest.main()
