import csv
import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
PYTHON_ROOT = ROOT / "python"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))


spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_manifest)

import vntyper_port


class VntyperFastqExpectedOutputsTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_fastq_kestrel_expected_outputs()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_fastq_kestrel_outputs_are_parseable_for_representative_samples(self):
        for label in ["positive", "negative"]:
            with self.subTest(label=label):
                root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                vcf = root / "kestrel" / "output.vcf"
                tsv = root / "kestrel" / "kestrel_result.tsv"
                report_json = root / "report.json"

                self.assertGreater(vcf.stat().st_size, 0)
                with tsv.open("r", encoding="utf-8", newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                with report_json.open("r", encoding="utf-8") as handle:
                    report = json.load(handle)

                self.assertGreater(len(rows), 0)
                self.assertEqual(len(report["kestrel_variants"]), len(rows))
                self.assertIn(
                    report["algorithm_results"]["kestrel"],
                    ["negative", "Low_Precision", "High_Precision", "High_Precision_flagged"],
                )
                self.assertEqual(report["metadata"]["alignment_pipeline"], "external kestrel from FASTQ")

    def test_reprocessed_java_kestrel_vcf_matches_expected_classification(self):
        for label in ["positive", "negative"]:
            with self.subTest(label=label):
                root = data_manifest.EXPECTED_OUTPUT_ROOT / label
                rows = vntyper_port.process_kestrel_vcf(str(root / "kestrel" / "output.vcf"))
                with (root / "report.json").open("r", encoding="utf-8") as handle:
                    report = json.load(handle)
                rebuilt = vntyper_port.build_report_json(
                    sample_name=report["sample_name"],
                    input_files=report["input_files"],
                    kestrel_rows=rows,
                    metadata=report["metadata"],
                )

                self.assertEqual(
                    rebuilt["algorithm_results"]["kestrel"],
                    report["algorithm_results"]["kestrel"],
                )

    def test_native_kestrel_rs_adapter_emits_expected_tiny_variant(self):
        try:
            from bioscript import kestrel

            data_manifest.import_native_module()
        except Exception as exc:
            self.skipTest(f"bioscript native extension is unavailable: {exc}")

        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "reads.fastq.gz"
            with gzip.open(fastq, "wt", encoding="utf-8") as handle:
                for index in range(5):
                    handle.write(f"@r{index}\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n")

            vcf = kestrel.call_fastq_references_native(
                [("chr1", "AAAACCCCGGGGTTTT", "2a9fd43653a81f9ec44e34c7ec038636")],
                [str(fastq)],
                4,
                sample_name="sample1",
                minimum_difference=1,
                max_haplotypes=4,
                max_saved_states=4,
            )

        self.assertIn("##fileformat=VCFv4.2\n", vcf)
        self.assertIn("##contig=<ID=chr1,length=16", vcf)
        self.assertIn("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1", vcf)
        self.assertIn("chr1\t5\t.\tC\tT", vcf)


if __name__ == "__main__":
    unittest.main()
