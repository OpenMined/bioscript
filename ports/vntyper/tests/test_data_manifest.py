import importlib.util
import os
import tempfile
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
        try:
            result = data_manifest.require_test_data(check_md5=False)
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))
        self.assertGreater(result["present"], 0)
        self.assertEqual(result["missing"], [])

    def test_validator_skip_message_names_data_drop_when_absent(self):
        missing = {
            "present": 0,
            "missing": [str(data_manifest.DATA_ROOT / "missing.bam")],
            "mismatched": [],
        }
        original = data_manifest.validate_manifest
        data_manifest.validate_manifest = lambda check_md5=False: missing
        try:
            with self.assertRaisesRegex(unittest.SkipTest, "ports/vntyper/test-data"):
                data_manifest.require_test_data(check_md5=False)
        finally:
            data_manifest.validate_manifest = original

    def test_kestrel_jar_can_be_overridden_by_environment(self):
        with tempfile.TemporaryDirectory() as tmp:
            jar = Path(tmp) / "kestrel.jar"
            jar.write_text("jar", encoding="utf-8")
            original = os.environ.get("BIOSCRIPT_KESTREL_JAR")
            os.environ["BIOSCRIPT_KESTREL_JAR"] = str(jar)
            try:
                self.assertEqual(data_manifest.resolve_kestrel_jar(), jar)
            finally:
                if original is None:
                    os.environ.pop("BIOSCRIPT_KESTREL_JAR", None)
                else:
                    os.environ["BIOSCRIPT_KESTREL_JAR"] = original

    def test_native_bam_skip_message_names_missing_opt_in_environment(self):
        original_env = os.environ.get("BIOSCRIPT_RUN_NATIVE_BAM_PARITY")
        os.environ.pop("BIOSCRIPT_RUN_NATIVE_BAM_PARITY", None)
        original_require = data_manifest.require_test_data
        original_import = data_manifest.import_native_module
        try:
            data_manifest.require_test_data = lambda check_md5=False: {
                "present": 1,
                "missing": [],
                "mismatched": [],
            }
            data_manifest.import_native_module = lambda: None
            with self.assertRaisesRegex(
                unittest.SkipTest,
                "BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1",
            ):
                data_manifest.require_all_native_bam_pipeline_prerequisites()
        finally:
            data_manifest.require_test_data = original_require
            data_manifest.import_native_module = original_import
            if original_env is None:
                os.environ.pop("BIOSCRIPT_RUN_NATIVE_BAM_PARITY", None)
            else:
                os.environ["BIOSCRIPT_RUN_NATIVE_BAM_PARITY"] = original_env

    def test_native_fastq_skip_message_names_missing_opt_in_environment(self):
        original_env = os.environ.get("BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY")
        os.environ.pop("BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY", None)
        original_require = data_manifest.require_test_data
        original_import = data_manifest.import_native_module
        try:
            data_manifest.require_test_data = lambda check_md5=False: {
                "present": 1,
                "missing": [],
                "mismatched": [],
            }
            data_manifest.import_native_module = lambda: None
            with self.assertRaisesRegex(
                unittest.SkipTest,
                "BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1",
            ):
                data_manifest.require_native_fastq_pipeline_prerequisites()
        finally:
            data_manifest.require_test_data = original_require
            data_manifest.import_native_module = original_import
            if original_env is None:
                os.environ.pop("BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY", None)
            else:
                os.environ["BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY"] = original_env

    def test_samtools_oracle_skip_message_names_missing_opt_in_environment(self):
        original_env = os.environ.get("BIOSCRIPT_RUN_SAMTOOLS_ORACLE")
        os.environ.pop("BIOSCRIPT_RUN_SAMTOOLS_ORACLE", None)
        original_require = data_manifest.require_test_data
        original_import = data_manifest.import_native_module
        try:
            data_manifest.require_test_data = lambda check_md5=False: {
                "present": 1,
                "missing": [],
                "mismatched": [],
            }
            data_manifest.import_native_module = lambda: None
            with self.assertRaisesRegex(
                unittest.SkipTest,
                "BIOSCRIPT_RUN_SAMTOOLS_ORACLE=1",
            ):
                data_manifest.require_samtools_fastq_oracle_prerequisites()
        finally:
            data_manifest.require_test_data = original_require
            data_manifest.import_native_module = original_import
            if original_env is None:
                os.environ.pop("BIOSCRIPT_RUN_SAMTOOLS_ORACLE", None)
            else:
                os.environ["BIOSCRIPT_RUN_SAMTOOLS_ORACLE"] = original_env


if __name__ == "__main__":
    unittest.main()
