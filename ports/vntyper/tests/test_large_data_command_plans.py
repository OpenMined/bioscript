import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
COMMANDS_PATH = BIOSCRIPT_PORT / "vntyper_commands.py"
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

manifest_spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(manifest_spec)
manifest_spec.loader.exec_module(data_manifest)

commands_spec = importlib.util.spec_from_file_location("vntyper_commands", COMMANDS_PATH)
vntyper_commands = importlib.util.module_from_spec(commands_spec)
sys.modules["vntyper_commands"] = vntyper_commands
commands_spec.loader.exec_module(vntyper_commands)


class LargeDataCommandPlanTests(unittest.TestCase):
    def setUp(self):
        try:
            data_manifest.require_test_data(check_md5=False)
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_representative_hg19_bams_plan_pre_kestrel_commands(self):
        samples = ["example_6449_hg19_subset", "example_66bf_hg19_subset"]
        for sample in samples:
            with self.subTest(sample=sample):
                bam = data_manifest.DATA_ROOT / f"{sample}.bam"
                bai = data_manifest.DATA_ROOT / f"{sample}.bam.bai"
                self.assertTrue(bam.exists())
                self.assertTrue(bai.exists())
                plan = vntyper_commands.plan_bam_pipeline(
                    str(bam),
                    sample,
                    assembly="hg19",
                    work_dir=f"work/{sample}",
                )
                self.assertEqual(plan.bam_region, "chr1:155158000-155163000")
                self.assertEqual(plan.vntr_region, "chr1:155160500-155162000")
                self.assertEqual(plan.samtools_view_command[0:3], ["samtools", "view", "-b"])
                self.assertEqual(plan.samtools_fastq_command[0], "samtools")
                self.assertIn(f"-s{sample}", plan.kestrel_command)

    def test_representative_fastq_pair_is_available_but_bwa_path_is_deferred(self):
        sample = "example_6449_hg19_subset"
        r1 = data_manifest.DATA_ROOT / f"{sample}_R1.fastq.gz"
        r2 = data_manifest.DATA_ROOT / f"{sample}_R2.fastq.gz"
        self.assertTrue(r1.exists())
        self.assertTrue(r2.exists())
        self.assertFalse(hasattr(vntyper_commands, "plan_fastq_alignment"))


if __name__ == "__main__":
    unittest.main()
