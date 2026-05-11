from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

from bioscript import bcftools, kestrel, samtools


class ToolCommandTests(unittest.TestCase):
    def test_kestrel_build_command_matches_vntyper_defaults(self) -> None:
        self.assertEqual(
            kestrel.build_command(
                "kestrel.jar",
                "muc1.fa",
                "out.vcf",
                "out.sam",
                "tmp",
                "sample1",
                "r1.fastq.gz",
                "r2.fastq.gz",
            ),
            [
                "java",
                "-Xmx12g",
                "-jar",
                "kestrel.jar",
                "-k",
                "20",
                "--maxalignstates",
                "40",
                "--maxhapstates",
                "40",
                "-r",
                "muc1.fa",
                "-o",
                "out.vcf",
                "-ssample1",
                "r1.fastq.gz",
                "r2.fastq.gz",
                "--hapfmt",
                "sam",
                "-p",
                "out.sam",
                "--logstderr",
                "--logstdout",
                "--loglevel",
                "INFO",
                "--temploc",
                "tmp",
            ],
        )

    def test_kestrel_rejects_shell_program(self) -> None:
        with self.assertRaises(ValueError):
            kestrel.build_command(
                "kestrel.jar",
                "muc1.fa",
                "out.vcf",
                "out.sam",
                "tmp",
                "sample1",
                "r1.fastq.gz",
                "r2.fastq.gz",
                java_program="java;rm",
            )

    def test_kestrel_native_sequences_wrapper_delegates_to_extension(self) -> None:
        calls = []

        def call_sequences(*args):
            calls.append(args)
            return "##fileformat=VCF4.2\n"

        fake_native = SimpleNamespace(kestrel_call_sequences_native=call_sequences)
        with patch.dict("sys.modules", {"bioscript._native": fake_native}):
            self.assertEqual(
                kestrel.call_sequences_native(
                    "MUC1",
                    "ACGT",
                    ["ACGT"],
                    3,
                    sample_name="sample1",
                    minimum_difference=1,
                    difference_quantile=0.0,
                    locus_depth=10,
                ),
                "##fileformat=VCF4.2\n",
            )
        self.assertEqual(
            calls,
            [
                (
                    "MUC1",
                    "ACGT",
                    ["ACGT"],
                    3,
                    "sample1",
                    "native",
                    ".",
                    1,
                    0.0,
                    True,
                    1,
                    40,
                    500,
                    0,
                    40,
                    10,
                )
            ],
        )

    def test_kestrel_native_fastq_wrapper_delegates_to_extension(self) -> None:
        calls = []

        def call_fastq(*args):
            calls.append(args)
            return "##fileformat=VCF4.2\n"

        fake_native = SimpleNamespace(kestrel_call_fastq_native=call_fastq)
        with patch.dict("sys.modules", {"bioscript._native": fake_native}):
            self.assertEqual(
                kestrel.call_fastq_native(
                    "MUC1",
                    "ACGT",
                    ["reads.fastq"],
                    3,
                    sample_name="sample1",
                    minimum_difference=1,
                    difference_quantile=0.0,
                    locus_depth=10,
                ),
                "##fileformat=VCF4.2\n",
            )
        self.assertEqual(calls[0][0:5], ("MUC1", "ACGT", ["reads.fastq"], 3, "sample1"))

    def test_kestrel_native_sequences_wrapper_reports_missing_extension(self) -> None:
        with patch.dict("sys.modules", {"bioscript._native": None}):
            with self.assertRaises(NotImplementedError):
                kestrel.call_sequences_native("MUC1", "ACGT", ["ACGT"], 3)

    def test_samtools_fastq_and_view_region(self) -> None:
        self.assertEqual(
            samtools.fastq("slice.bam", "r1.fastq.gz", "r2.fastq.gz"),
            ["samtools", "fastq", "-1", "r1.fastq.gz", "-2", "r2.fastq.gz", "slice.bam"],
        )
        self.assertEqual(
            samtools.view_region("sample.bam", "chr1:1-10", "slice.bam"),
            ["samtools", "view", "-b", "sample.bam", "chr1:1-10", "-o", "slice.bam"],
        )
        self.assertEqual(
            samtools.depth("slice.bam", "chr1:1-10", include_zero=True),
            ["samtools", "depth", "-a", "-r", "chr1:1-10", "slice.bam"],
        )

    def test_samtools_native_wrappers_delegate_to_extension(self) -> None:
        calls = []

        def view_region_native(bam, index, region, output):
            calls.append((bam, index, region, output))
            return 7

        fake_native = SimpleNamespace(
            samtools_view_region_native=view_region_native,
            samtools_depth_native=lambda bam, index, region: {"mean": 2.5},
            samtools_fastq_native=lambda bam, index, region, fastq_1, fastq_2: {
                "read1_records": 3,
                "read2_records": 3,
                "skipped_records": 1,
            },
        )
        with patch.dict("sys.modules", {"bioscript._native": fake_native}):
            self.assertEqual(
                samtools.view_region_native(
                    "sample.bam",
                    "chr1:1-10",
                    "slice.bam",
                    index="sample.bam.bai",
                ),
                7,
            )
            self.assertEqual(
                calls,
                [("sample.bam", "sample.bam.bai", "chr1:1-10", "slice.bam")],
            )
            self.assertEqual(samtools.depth_native("slice.bam", "chr1:1-10"), {"mean": 2.5})
            self.assertEqual(
                samtools.fastq_native(
                    "slice.bam",
                    "chr1:1-10",
                    "r1.fastq.gz",
                    "r2.fastq.gz",
                ),
                {"read1_records": 3, "read2_records": 3, "skipped_records": 1},
            )

    def test_samtools_native_wrappers_report_missing_extension(self) -> None:
        with patch.dict("sys.modules", {"bioscript._native": None}):
            with self.assertRaises(NotImplementedError):
                samtools.fastq_native(
                    "slice.bam",
                    "chr1:1-10",
                    "r1.fastq.gz",
                    "r2.fastq.gz",
                )

    def test_bcftools_vcf_helpers(self) -> None:
        self.assertEqual(
            bcftools.sort("calls.vcf", "calls.vcf.gz"),
            ["bcftools", "sort", "-Oz", "-o", "calls.vcf.gz", "calls.vcf"],
        )
        self.assertEqual(
            bcftools.view_filter("calls.vcf", "pass.vcf.gz", 'FILTER="PASS"'),
            ["bcftools", "view", "-i", 'FILTER="PASS"', "-Oz", "-o", "pass.vcf.gz", "calls.vcf"],
        )


if __name__ == "__main__":
    unittest.main()
