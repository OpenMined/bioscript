import importlib.util
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
COMMANDS_PATH = BIOSCRIPT_PORT / "vntyper_commands.py"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

spec = importlib.util.spec_from_file_location("vntyper_commands", COMMANDS_PATH)
vntyper_commands = importlib.util.module_from_spec(spec)
sys.modules["vntyper_commands"] = vntyper_commands
spec.loader.exec_module(vntyper_commands)


class VntyperCommandPlanTests(unittest.TestCase):
    def test_bam_pipeline_plan_uses_region_slice_before_fastq(self):
        plan = vntyper_commands.plan_bam_pipeline(
            "sample.bam",
            "sample1",
            assembly="hg19",
            work_dir="work",
        )
        self.assertEqual(plan.bam_region, "chr1:155158000-155163000")
        self.assertEqual(plan.vntr_region, "chr1:155160500-155162000")
        self.assertEqual(
            plan.samtools_view_command,
            ["samtools", "view", "-b", "sample.bam", "chr1:155158000-155163000", "-o", "work/sample1_sliced.bam"],
        )
        self.assertEqual(
            plan.samtools_fastq_command,
            [
                "samtools",
                "fastq",
                "-1",
                "work/sample1_R1.fastq.gz",
                "-2",
                "work/sample1_R2.fastq.gz",
                "work/sample1_sliced.bam",
            ],
        )
        self.assertIn("-ssample1", plan.kestrel_command)
        self.assertEqual(plan.bcftools_index_command, ["bcftools", "index", "-t", "work/kestrel/output.sorted.vcf.gz"])

    def test_bam_pipeline_can_plan_ncbi_regions(self):
        plan = vntyper_commands.plan_bam_pipeline(
            "sample.bam",
            "sample1",
            assembly="hg38_ncbi",
        )
        self.assertEqual(plan.bam_region, "NC_000001.11:155184000-155194000")
        self.assertEqual(plan.vntr_region, "NC_000001.11:155188000-155192500")

    def test_rejects_path_like_sample_names(self):
        with self.assertRaises(ValueError):
            vntyper_commands.plan_bam_pipeline("sample.bam", "../sample")


if __name__ == "__main__":
    unittest.main()
