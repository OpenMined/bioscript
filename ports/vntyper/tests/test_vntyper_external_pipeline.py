import importlib.util
import json
import shutil
import sys
import tempfile
import unittest
from types import SimpleNamespace
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
PIPELINE_PATH = BIOSCRIPT_PORT / "vntyper_external_pipeline.py"
FIXTURE_VCF = ROOT / "ports" / "vntyper" / "tests" / "fixtures" / "kestrel_minimal.vcf"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

spec = importlib.util.spec_from_file_location("vntyper_external_pipeline", PIPELINE_PATH)
vntyper_external_pipeline = importlib.util.module_from_spec(spec)
sys.modules["vntyper_external_pipeline"] = vntyper_external_pipeline
spec.loader.exec_module(vntyper_external_pipeline)


class VntyperExternalPipelineTests(unittest.TestCase):
    def test_dry_run_returns_ordered_external_commands(self):
        result = vntyper_external_pipeline.run_bam_pipeline(
            "sample.bam",
            "sample1",
            "work/sample1",
            dry_run=True,
        )

        self.assertEqual([command[0] for command in result.commands], ["samtools", "samtools", "samtools", "samtools", "java", "bcftools", "bcftools"])
        self.assertEqual(result.kestrel_vcf, "work/sample1/kestrel/output.vcf")
        self.assertEqual(result.kestrel_tsv, "work/sample1/kestrel/kestrel_result.tsv")
        self.assertEqual(result.report_json, "work/sample1/report.json")

    def test_dry_run_can_plan_native_samtools_bam_path(self):
        result = vntyper_external_pipeline.run_bam_pipeline(
            "sample.bam",
            "sample1",
            "work/sample1",
            dry_run=True,
            use_native_samtools=True,
        )

        self.assertEqual(
            [command[0] for command in result.commands],
            [
                "bioscript.samtools.view_region_native",
                "bioscript.samtools.fastq_native",
                "bioscript.samtools.depth_native",
                "java",
            ],
        )
        self.assertNotIn("bcftools", [command[0] for command in result.commands])
        self.assertEqual(result.commands[0][-1], "sample.bam.bai")

    def test_runner_materializes_kestrel_tsv_and_report_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls = []

            def fake_runner(command, check, **kwargs):
                calls.append(command)
                if command[0] == "samtools" and command[1] == "view":
                    Path(command[command.index("-o") + 1]).write_bytes(b"bam")
                if command[0] == "samtools" and command[1] == "fastq":
                    Path(command[command.index("-1") + 1]).write_bytes(b"r1")
                    Path(command[command.index("-2") + 1]).write_bytes(b"r2")
                if command[0] == "samtools" and command[1] == "depth":
                    self.assertTrue(kwargs["capture_output"])
                    return SimpleNamespace(stdout="chr1\t100\t10\nchr1\t101\t0\nchr1\t102\t20\n")
                if command[0] == "java":
                    shutil.copyfile(FIXTURE_VCF, command[command.index("-o") + 1])
                    Path(command[command.index("-p") + 1]).write_text("@HD\n", encoding="utf-8")
                return SimpleNamespace(stdout="")

            result = vntyper_external_pipeline.run_bam_pipeline(
                "sample.bam",
                "sample1",
                str(Path(tmp) / "sample1"),
                runner=fake_runner,
            )

            self.assertEqual([command[0] for command in calls], ["samtools", "samtools", "samtools", "samtools", "java", "bcftools", "bcftools"])
            self.assertTrue(Path(result.kestrel_tsv).exists())
            self.assertTrue(Path(result.report_json).exists())
            with open(result.kestrel_tsv, "r", encoding="utf-8") as handle:
                tsv = handle.read()
            self.assertIn("Depth_Score", tsv)
            self.assertIn("High_Precision", tsv)
            with open(result.report_json, "r", encoding="utf-8") as handle:
                report = json.load(handle)
            self.assertEqual(report["sample_name"], "sample1")
            self.assertEqual(report["metadata"]["alignment_pipeline"], "external samtools/kestrel")
            self.assertEqual(report["coverage"]["mean"], 10.0)
            self.assertEqual(report["coverage"]["median"], 10)
            self.assertEqual(report["coverage"]["min"], 0)
            self.assertEqual(report["coverage"]["max"], 20)
            self.assertEqual(report["coverage"]["region_length"], 3)
            self.assertEqual(report["coverage"]["uncovered_bases"], 1)
            self.assertEqual(len(report["pipeline_log"]), 7)

    def test_native_samtools_runner_materializes_bam_path_without_bcftools(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls = []

            class FakeNativeSamtools:
                def view_region_native(self, bam, region, output_bam, index=None):
                    calls.append(("view", bam, region, output_bam, index))
                    Path(output_bam).write_bytes(b"bam")
                    return 1

                def fastq_native(self, bam, region, fastq_1, fastq_2, index=None):
                    calls.append(("fastq", bam, region, fastq_1, fastq_2, index))
                    Path(fastq_1).write_bytes(b"r1")
                    Path(fastq_2).write_bytes(b"r2")
                    return {"read1_records": 1, "read2_records": 1, "skipped_records": 0}

                def depth_native(self, bam, region, index=None):
                    calls.append(("depth", bam, region, index))
                    return {
                        "mean": 10.0,
                        "median": 10.0,
                        "stdev": 8.16496580927726,
                        "min": 0,
                        "max": 20,
                        "region_length": 3,
                        "uncovered_bases": 1,
                        "percent_uncovered": 33.33333333333333,
                    }

            def fake_runner(command, check):
                calls.append(("kestrel", command))
                shutil.copyfile(FIXTURE_VCF, command[command.index("-o") + 1])
                Path(command[command.index("-p") + 1]).write_text("@HD\n", encoding="utf-8")

            result = vntyper_external_pipeline.run_bam_pipeline(
                "sample.bam",
                "sample1",
                str(Path(tmp) / "sample1"),
                runner=fake_runner,
                use_native_samtools=True,
                native_samtools=FakeNativeSamtools(),
            )

            self.assertEqual([call[0] for call in calls], ["view", "fastq", "depth", "kestrel"])
            self.assertTrue(Path(result.kestrel_tsv).exists())
            with open(result.report_json, "r", encoding="utf-8") as handle:
                report = json.load(handle)
            self.assertEqual(
                report["metadata"]["alignment_pipeline"],
                "native bioscript samtools/kestrel",
            )
            self.assertEqual(report["coverage"]["mean"], 10.0)
            self.assertEqual(len(report["pipeline_log"]), 4)
            self.assertEqual(
                report["pipeline_log"][0]["command"][0],
                "bioscript.samtools.view_region_native",
            )

    def test_coverage_from_depth_ignores_malformed_lines(self):
        coverage = vntyper_external_pipeline.coverage_from_depth(
            "chr1\t10\t5\nbad\nchr1\t11\tNA\nchr1\t12\t15\n"
        )

        self.assertEqual(coverage["mean"], 10.0)
        self.assertEqual(coverage["median"], 10.0)
        self.assertEqual(coverage["region_length"], 2)

    def test_fastq_kestrel_runner_materializes_outputs_without_samtools(self):
        with tempfile.TemporaryDirectory() as tmp:
            calls = []

            def fake_runner(command, check):
                calls.append(command)
                shutil.copyfile(FIXTURE_VCF, command[command.index("-o") + 1])
                Path(command[command.index("-p") + 1]).write_text("@HD\n", encoding="utf-8")

            result = vntyper_external_pipeline.run_fastq_kestrel(
                "sample_R1.fastq.gz",
                "sample_R2.fastq.gz",
                "sample1",
                str(Path(tmp) / "sample1"),
                runner=fake_runner,
            )

            self.assertEqual(len(calls), 1)
            self.assertEqual(calls[0][0], "java")
            self.assertTrue(Path(result.kestrel_tsv).exists())
            self.assertTrue(Path(result.report_json).exists())
            with open(result.report_json, "r", encoding="utf-8") as handle:
                report = json.load(handle)
            self.assertEqual(report["input_files"]["fastq_1"], "sample_R1.fastq.gz")
            self.assertEqual(report["input_files"]["fastq_2"], "sample_R2.fastq.gz")
            self.assertEqual(report["metadata"]["alignment_pipeline"], "external kestrel from FASTQ")


if __name__ == "__main__":
    unittest.main()
