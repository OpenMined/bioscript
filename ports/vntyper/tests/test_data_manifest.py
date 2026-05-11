import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"


spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_manifest)


class VntyperDataManifestTests(unittest.TestCase):
    def test_manifest_maps_upstream_test_data_into_port_tree(self):
        entries = data_manifest.load_manifest()
        self.assertGreater(len(entries), 0)
        first = entries[0]
        self.assertTrue(str(first["path"]).startswith(str(data_manifest.DATA_ROOT)))
        self.assertEqual(first["filename"], "example_6449_hg19_subset.bam")

    def test_validator_sees_copied_data_without_md5_scan(self):
        result = data_manifest.validate_manifest(check_md5=False)
        self.assertGreater(result["present"], 0)
        self.assertEqual(result["missing"], [])


if __name__ == "__main__":
    unittest.main()
