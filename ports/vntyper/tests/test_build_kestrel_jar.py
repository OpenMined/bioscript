import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BUILDER_PATH = ROOT / "ports" / "vntyper" / "tests" / "build_kestrel_jar.py"


spec = importlib.util.spec_from_file_location("build_kestrel_jar", BUILDER_PATH)
build_kestrel_jar = importlib.util.module_from_spec(spec)
spec.loader.exec_module(build_kestrel_jar)


class BuildKestrelJarTests(unittest.TestCase):
    def test_discovers_vendored_sources_and_dependency_classpath(self):
        sources = build_kestrel_jar.source_files()

        self.assertGreater(len(sources), 50)
        self.assertTrue(any(source.endswith("edu/gatech/kestrel/clui/Main.java") for source in sources))
        self.assertFalse(any("/test/" in source for source in sources))
        classpath = build_kestrel_jar.classpath()
        self.assertIn("kanalyze.jar", classpath)
        self.assertIn("logback-classic-1.1.3.jar", classpath)

    def test_default_output_uses_ignored_test_data_tools_directory(self):
        self.assertIn("ports/vntyper/test-data/tools/kestrel", str(build_kestrel_jar.DEFAULT_OUTPUT))

    def test_manifest_uses_relative_lib_paths_for_kestrel_root_output(self):
        manifest = build_kestrel_jar.manifest_content(
            build_kestrel_jar.KESTREL_ROOT / "kestrel.jar"
        )

        self.assertIn("Main-Class: edu.gatech.kestrel.clui.Main", manifest)
        self.assertIn("Class-Path: lib/kanalyze.jar", manifest)

    def test_manifest_wraps_long_attribute_lines(self):
        attribute = build_kestrel_jar.manifest_attribute("Class-Path", "x" * 150)

        lines = attribute.splitlines()
        self.assertGreater(len(lines), 1)
        self.assertTrue(lines[1].startswith(" "))
        self.assertTrue(all(len(line) <= 70 for line in lines))


if __name__ == "__main__":
    unittest.main()
