import gzip
import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PYTHON_ROOT = ROOT / "python"
BIOSCRIPT_PORT = ROOT / "ports" / "vntyper" / "bioscript"
MANIFEST_PATH = ROOT / "ports" / "vntyper" / "tests" / "data_manifest.py"

sys.path.insert(0, str(PYTHON_ROOT))
sys.path.insert(0, str(BIOSCRIPT_PORT))

manifest_spec = importlib.util.spec_from_file_location("data_manifest", MANIFEST_PATH)
data_manifest = importlib.util.module_from_spec(manifest_spec)
manifest_spec.loader.exec_module(data_manifest)

from bioscript import samtools  # noqa: E402

try:
    import vntyper_regions  # noqa: E402
except ImportError:
    from ports.vntyper.bioscript import vntyper_regions  # noqa: E402


class SamtoolsFastqOracleTests(unittest.TestCase):
    def setUp(self):
        try:
            self.prereqs = data_manifest.require_samtools_fastq_oracle_prerequisites()
        except unittest.SkipTest as skip:
            self.skipTest(str(skip))

    def test_native_fastq_counts_match_samtools_name_sorted_pair_extraction(self):
        region = vntyper_regions.region_string("hg19", "bam_region_coords")
        for label, bam in self.prereqs["bam_cases"].items():
            with self.subTest(label=label):
                with tempfile.TemporaryDirectory() as tmp:
                    tmp = Path(tmp)
                    native_r1 = tmp / "native_R1.fastq.gz"
                    native_r2 = tmp / "native_R2.fastq.gz"
                    native_summary = samtools.fastq_native(
                        bam,
                        region,
                        str(native_r1),
                        str(native_r2),
                        index=f"{bam}.bai",
                    )

                    oracle_counts = run_samtools_oracle(bam, region, tmp)

                self.assertEqual(native_summary["read1_records"], oracle_counts["read1_records"])
                self.assertEqual(native_summary["read2_records"], oracle_counts["read2_records"])


def run_samtools_oracle(bam, region, tmp):
    sliced = tmp / "slice.bam"
    sorted_bam = tmp / "slice.name.bam"
    read1 = tmp / "samtools_R1.fastq.gz"
    read2 = tmp / "samtools_R2.fastq.gz"
    other = tmp / "samtools_other.fastq.gz"
    singleton = tmp / "samtools_single.fastq.gz"

    env = os.environ.copy()
    env["PATH"] = f"{data_manifest.LOCAL_TOOL_BIN}{os.pathsep}{env.get('PATH', '')}"
    subprocess.run(
        ["samtools", "view", "-P", "-b", bam, region, "-o", str(sliced)],
        check=True,
        env=env,
    )
    subprocess.run(
        ["samtools", "sort", "-n", "-o", str(sorted_bam), str(sliced)],
        check=True,
        env=env,
    )
    subprocess.run(
        [
            "samtools",
            "fastq",
            str(sorted_bam),
            "-1",
            str(read1),
            "-2",
            str(read2),
            "-0",
            str(other),
            "-s",
            str(singleton),
        ],
        check=True,
        env=env,
    )
    return {
        "read1_records": count_fastq_records(read1),
        "read2_records": count_fastq_records(read2),
        "other_records": count_fastq_records(other),
        "singleton_records": count_fastq_records(singleton),
    }


def count_fastq_records(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return sum(1 for index, _ in enumerate(handle, start=1) if index % 4 == 1)


if __name__ == "__main__":
    unittest.main()
