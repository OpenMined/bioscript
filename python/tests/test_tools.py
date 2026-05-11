from __future__ import annotations

import unittest

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
