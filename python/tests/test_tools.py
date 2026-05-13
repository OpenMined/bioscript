from __future__ import annotations

import importlib
import sys
import tempfile
import unittest
from pathlib import Path
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
                    0.55,
                    0.8,
                    7,
                    7.0,
                    None,
                    True,
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

    def test_kestrel_native_multireference_fastq_wrapper_delegates_to_extension(self) -> None:
        calls = []

        def call_fastq_references(*args):
            calls.append(args)
            return "##fileformat=VCF4.2\n"

        fake_native = SimpleNamespace(kestrel_call_fastq_references_native=call_fastq_references)
        with patch.dict("sys.modules", {"bioscript._native": fake_native}):
            self.assertEqual(
                kestrel.call_fastq_references_native(
                    [("REF1", "ACGT", "md5-1"), ("REF2", "TGCA", "md5-2")],
                    ["reads.fastq"],
                    3,
                    sample_name="sample1",
                    minimum_difference=1,
                    difference_quantile=0.0,
                    locus_depth=10,
                ),
                "##fileformat=VCF4.2\n",
            )
        self.assertEqual(
            calls[0][0:5],
            (
                [("REF1", "ACGT", "md5-1"), ("REF2", "TGCA", "md5-2")],
                ["reads.fastq"],
                3,
                "sample1",
                "native",
            ),
        )
        self.assertEqual(calls[0][-1], 10)

    def test_kestrel_load_reference_regions_reads_fasta_with_md5(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "refs.fa"
            path.write_text(">REF1 description\nAAAA\nCCCC\n>REF2\nACAGTCCGTAAG\n", encoding="utf-8")

            self.assertEqual(
                kestrel.load_reference_regions(str(path)),
                [
                    ("REF1", "AAAACCCC", "7b0d393d76107409cd695d4a86386703"),
                    ("REF2", "ACAGTCCGTAAG", "f17cc056a4c30b8661b5585d2641a37a"),
                ],
            )

    def test_kestrel_load_reference_regions_rejects_empty_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "refs.fa"
            path.write_text("\n", encoding="utf-8")

            with self.assertRaises(ValueError):
                kestrel.load_reference_regions(str(path))

    def test_kestrel_native_sequences_wrapper_reports_missing_extension(self) -> None:
        with patch.dict("sys.modules", {"bioscript._native": None}):
            with self.assertRaises(NotImplementedError):
                kestrel.call_sequences_native("MUC1", "ACGT", ["ACGT"], 3)
            with self.assertRaises(NotImplementedError):
                kestrel.call_fastq_native("MUC1", "ACGT", ["reads.fastq"], 4)
            with self.assertRaises(NotImplementedError):
                kestrel.call_fastq_references_native([("MUC1", "ACGT", "md5")], ["reads.fastq"], 4)

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
                samtools.view_region_native("slice.bam", "chr1:1-10", "out.bam")
            with self.assertRaises(NotImplementedError):
                samtools.depth_native("slice.bam", "chr1:1-10")
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

    def test_bcftools_native_view_header_wrapper_delegates_to_extension(self) -> None:
        calls = []

        def view_header(input_vcf, output_vcf):
            calls.append((input_vcf, output_vcf))

        def view(input_vcf, output_vcf, output_type):
            calls.append((input_vcf, output_vcf, output_type))

        def index(input_vcf, output_index, tbi, force):
            calls.append((input_vcf, output_index, tbi, force))

        fake_native = SimpleNamespace(
            bcftools_view_header_native=view_header,
            bcftools_view_native=view,
            bcftools_index_native=index,
        )
        with patch.dict("sys.modules", {"bioscript._native": fake_native}):
            self.assertIsNone(bcftools.view_header_native("calls.vcf", "header.vcf"))
            self.assertIsNone(bcftools.view_native("calls.vcf", "calls.vcf.gz", output_type="z"))
            self.assertIsNone(
                bcftools.index_native(
                    "calls.vcf.gz",
                    output_index="calls.vcf.gz.tbi",
                    tbi=True,
                    force=False,
                )
            )

        self.assertEqual(
            calls,
            [
                ("calls.vcf", "header.vcf"),
                ("calls.vcf", "calls.vcf.gz", "z"),
                ("calls.vcf.gz", "calls.vcf.gz.tbi", True, False),
            ],
        )

    def test_bcftools_native_view_header_reports_missing_extension(self) -> None:
        with patch.dict("sys.modules", {"bioscript._native": None}):
            with self.assertRaises(NotImplementedError):
                bcftools.view_header_native("calls.vcf", "header.vcf")
            with self.assertRaises(NotImplementedError):
                bcftools.view_native("calls.vcf", "calls.vcf.gz", output_type="z")
            with self.assertRaises(NotImplementedError):
                bcftools.index_native("calls.vcf.gz", "calls.vcf.gz.tbi")

    def test_bcftools_native_view_header_real_extension_extracts_header(self) -> None:
        try:
            import bioscript as bioscript_package

            native = importlib.import_module("bioscript._native")
        except ImportError as exc:
            self.skipTest(f"BioScript native extension is not installed: {exc}")

        try:
            with tempfile.TemporaryDirectory() as tmp:
                input_vcf = Path(tmp) / "input.vcf"
                output_vcf = Path(tmp) / "header.vcf"
                input_vcf.write_text(
                    "##fileformat=VCFv4.2\n"
                    "##contig=<ID=chr1,length=16>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t5\t.\tC\tT\t.\tPASS\t.\n",
                    encoding="utf-8",
                )

                bcftools.view_header_native(str(input_vcf), str(output_vcf))

                header = output_vcf.read_text(encoding="utf-8")
        finally:
            if getattr(bioscript_package, "_native", None) is native:
                delattr(bioscript_package, "_native")
            sys.modules.pop("bioscript._native", None)

        self.assertIn("##fileformat=VCFv4.2\n", header)
        self.assertIn("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", header)
        self.assertNotIn("chr1\t5\t.\tC\tT", header)

    def test_bcftools_native_view_and_index_real_extension(self) -> None:
        try:
            import bioscript as bioscript_package

            native = importlib.import_module("bioscript._native")
        except ImportError as exc:
            self.skipTest(f"BioScript native extension is not installed: {exc}")

        try:
            with tempfile.TemporaryDirectory() as tmp:
                input_vcf = Path(tmp) / "input.vcf"
                output_vcf = Path(tmp) / "output.vcf"
                output_gz = Path(tmp) / "output.vcf.gz"
                output_tbi = Path(tmp) / "output.vcf.gz.tbi"
                input_vcf.write_text(
                    "##fileformat=VCFv4.2\n"
                    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
                    "##contig=<ID=chr1,length=1000>\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                    "chr1\t5\t.\tC\tT\t.\tPASS\t.\n",
                    encoding="utf-8",
                )

                bcftools.view_native(str(input_vcf), str(output_vcf))
                bcftools.view_native(str(input_vcf), str(output_gz), output_type="z")
                bcftools.index_native(str(output_gz), str(output_tbi))

                text = output_vcf.read_text(encoding="utf-8")
                index_size = output_tbi.stat().st_size
        finally:
            if getattr(bioscript_package, "_native", None) is native:
                delattr(bioscript_package, "_native")
            sys.modules.pop("bioscript._native", None)

        self.assertIn("chr1\t5\t.\tC\tT", text)
        self.assertGreater(index_size, 0)


if __name__ == "__main__":
    unittest.main()
