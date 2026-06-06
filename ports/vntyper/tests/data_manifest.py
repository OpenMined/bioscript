"""VNtyper test-data manifest and validator.

The copied large data lives in `ports/vntyper/test-data`. Upstream VNtyper's
manifest expects paths under `tests/data`, so this helper remaps those entries
into the BioScript port tree and can optionally verify MD5 checksums.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
UPSTREAM_CONFIG = ROOT / "ports" / "vntyper" / "vntyper" / "tests" / "test_data_config.json"
DATA_ROOT = ROOT / "ports" / "vntyper" / "test-data"
DEFAULT_KESTREL_JAR = ROOT / "ports" / "vntyper" / "kestrel" / "kestrel.jar"
TEST_DATA_KESTREL_JAR = DATA_ROOT / "tools" / "kestrel" / "kestrel.jar"
MUC1_REFERENCE = (
    ROOT
    / "ports"
    / "vntyper"
    / "vntyper"
    / "reference"
    / "All_Pairwise_and_Self_Merged_MUC1_motifs_filtered.fa"
)
EXPECTED_OUTPUT_ROOT = DATA_ROOT / "expected"
LOCAL_TOOL_BIN = DATA_ROOT / "tools" / "local" / "bin"
EXPECTED_OUTPUTS = [
    EXPECTED_OUTPUT_ROOT / "positive" / "kestrel" / "output.vcf",
    EXPECTED_OUTPUT_ROOT / "positive" / "kestrel" / "kestrel_result.tsv",
    EXPECTED_OUTPUT_ROOT / "positive" / "report.json",
    EXPECTED_OUTPUT_ROOT / "negative" / "kestrel" / "output.vcf",
    EXPECTED_OUTPUT_ROOT / "negative" / "kestrel" / "kestrel_result.tsv",
    EXPECTED_OUTPUT_ROOT / "negative" / "report.json",
]
REPRESENTATIVE_BAM_CASES = {
    "positive": DATA_ROOT / "example_6449_hg19_subset.bam",
    "negative": DATA_ROOT / "example_7a61_hg19_subset.bam",
}
REPRESENTATIVE_FASTQ_CASES = {
    "positive": (
        DATA_ROOT / "example_6449_hg19_subset_R1.fastq.gz",
        DATA_ROOT / "example_6449_hg19_subset_R2.fastq.gz",
    ),
    "negative": (
        DATA_ROOT / "example_7a61_hg19_subset_R1.fastq.gz",
        DATA_ROOT / "example_7a61_hg19_subset_R2.fastq.gz",
    ),
}

def resolve_kestrel_jar():
    env_path = os.environ.get("BIOSCRIPT_KESTREL_JAR")
    candidates = [
        Path(env_path) if env_path else None,
        TEST_DATA_KESTREL_JAR,
        DEFAULT_KESTREL_JAR,
    ]
    return next(
        (path for path in candidates if path is not None and path.exists()),
        TEST_DATA_KESTREL_JAR,
    )


KESTREL_JAR = resolve_kestrel_jar()


def require_test_data(check_md5=False):
    """Skip an integration test unless the ignored VNtyper data drop is present."""
    result = validate_manifest(check_md5=check_md5)
    if result["missing"]:
        preview = ", ".join(result["missing"][:3])
        remaining = len(result["missing"]) - min(len(result["missing"]), 3)
        suffix = f", plus {remaining} more" if remaining else ""
        raise unittest.SkipTest(
            "VNtyper integration data is absent from ports/vntyper/test-data: "
            f"{preview}{suffix}"
        )
    if result["mismatched"]:
        first = result["mismatched"][0]
        raise unittest.SkipTest(
            "VNtyper integration data checksum mismatch: "
            f"{first['path']} expected {first['expected']} got {first['actual']}"
        )
    return result


def require_full_pipeline_prerequisites():
    """Skip full external pipeline tests unless tools, data, and expected outputs exist."""
    manifest = require_test_data(check_md5=False)
    missing = []
    samtools_path = which_tool("samtools")
    bcftools_path = which_tool("bcftools")
    if samtools_path is None:
        missing.append("samtools on PATH or in ports/vntyper/test-data/tools/local/bin")
    if bcftools_path is None:
        missing.append("bcftools on PATH or in ports/vntyper/test-data/tools/local/bin")
    if shutil.which("java") is None:
        missing.append("java on PATH")
    if not KESTREL_JAR.exists():
        missing.append(str(KESTREL_JAR))
    if not MUC1_REFERENCE.exists():
        missing.append(str(MUC1_REFERENCE))
    missing_outputs = [str(path) for path in EXPECTED_OUTPUTS if not path.exists()]
    if missing_outputs:
        preview = ", ".join(missing_outputs[:3])
        remaining = len(missing_outputs) - min(len(missing_outputs), 3)
        suffix = f", plus {remaining} more" if remaining else ""
        missing.append(f"expected pipeline outputs: {preview}{suffix}")
    if missing:
        raise unittest.SkipTest(
            "VNtyper full pipeline prerequisites are missing: " + "; ".join(missing)
        )
    return {
        "manifest": manifest,
        "samtools": samtools_path,
        "bcftools": bcftools_path,
        "java": shutil.which("java"),
        "tool_path": str(LOCAL_TOOL_BIN),
        "kestrel_jar": str(KESTREL_JAR),
        "muc1_reference": str(MUC1_REFERENCE),
        "expected_outputs": [str(path) for path in EXPECTED_OUTPUTS],
    }


def require_external_bam_pipeline_prerequisites():
    """Skip unless the external samtools/bcftools BAM path is explicitly enabled."""
    prereqs = require_full_pipeline_prerequisites()
    missing = []
    if os.environ.get("BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY") != "1":
        missing.append("BIOSCRIPT_RUN_EXTERNAL_BAM_PARITY=1")
    missing_cases = [
        str(path)
        for bam in REPRESENTATIVE_BAM_CASES.values()
        for path in [bam, Path(f"{bam}.bai")]
        if not path.exists()
    ]
    missing.extend(missing_cases)
    if missing:
        raise unittest.SkipTest(
            "VNtyper external BAM pipeline prerequisites are missing: " + "; ".join(missing)
        )
    return {
        **prereqs,
        "bam_cases": {label: str(path) for label, path in REPRESENTATIVE_BAM_CASES.items()},
    }


def require_fastq_kestrel_expected_outputs():
    """Skip unless FASTQ-generated Kestrel expected outputs are present."""
    manifest = require_test_data(check_md5=False)
    missing = []
    if shutil.which("java") is None:
        missing.append("java on PATH")
    if not KESTREL_JAR.exists():
        missing.append(str(KESTREL_JAR))
    if not MUC1_REFERENCE.exists():
        missing.append(str(MUC1_REFERENCE))
    missing_outputs = [str(path) for path in EXPECTED_OUTPUTS if not path.exists()]
    if missing_outputs:
        preview = ", ".join(missing_outputs[:3])
        remaining = len(missing_outputs) - min(len(missing_outputs), 3)
        suffix = f", plus {remaining} more" if remaining else ""
        missing.append(f"FASTQ Kestrel expected outputs: {preview}{suffix}")
    if missing:
        raise unittest.SkipTest(
            "VNtyper FASTQ Kestrel expected outputs are missing: " + "; ".join(missing)
        )
    return {
        "manifest": manifest,
        "java": shutil.which("java"),
        "kestrel_jar": str(KESTREL_JAR),
        "muc1_reference": str(MUC1_REFERENCE),
        "expected_outputs": [str(path) for path in EXPECTED_OUTPUTS],
    }


def require_native_bam_pipeline_prerequisites():
    """Skip unless the native-samtools BAM path can run against copied data."""
    missing = []
    if shutil.which("java") is None:
        missing.append("java on PATH")
    if not KESTREL_JAR.exists():
        missing.append(str(KESTREL_JAR))
    try:
        prereqs = require_all_native_bam_pipeline_prerequisites()
    except unittest.SkipTest as skip:
        missing.append(str(skip))
        prereqs = {}
    if missing:
        raise unittest.SkipTest(
            "VNtyper native BAM pipeline prerequisites are missing: " + "; ".join(missing)
        )
    return {
        **prereqs,
        "java": shutil.which("java"),
        "kestrel_jar": str(KESTREL_JAR),
    }


def require_all_native_bam_pipeline_prerequisites():
    """Skip unless the all-native BAM path can run against copied data."""
    manifest = require_test_data(check_md5=False)
    missing = []
    if os.environ.get("BIOSCRIPT_RUN_NATIVE_BAM_PARITY") != "1":
        missing.append("BIOSCRIPT_RUN_NATIVE_BAM_PARITY=1")
    if not MUC1_REFERENCE.exists():
        missing.append(str(MUC1_REFERENCE))
    missing_cases = [
        str(path)
        for bam in REPRESENTATIVE_BAM_CASES.values()
        for path in [bam, Path(f"{bam}.bai")]
        if not path.exists()
    ]
    missing.extend(missing_cases)
    missing_outputs = [str(path) for path in EXPECTED_OUTPUTS if not path.exists()]
    if missing_outputs:
        preview = ", ".join(missing_outputs[:3])
        remaining = len(missing_outputs) - min(len(missing_outputs), 3)
        suffix = f", plus {remaining} more" if remaining else ""
        missing.append(f"native BAM expected outputs: {preview}{suffix}")
    try:
        import_native_module()
    except Exception as exc:
        missing.append(f"bioscript._native importable ({exc})")
    if missing:
        raise unittest.SkipTest(
            "VNtyper native BAM pipeline prerequisites are missing: " + "; ".join(missing)
        )
    return {
        "manifest": manifest,
        "muc1_reference": str(MUC1_REFERENCE),
        "expected_outputs": [str(path) for path in EXPECTED_OUTPUTS],
        "bam_cases": {label: str(path) for label, path in REPRESENTATIVE_BAM_CASES.items()},
    }


def require_native_fastq_pipeline_prerequisites():
    """Skip unless the native-Kestrel FASTQ path can run against copied data."""
    manifest = require_test_data(check_md5=False)
    missing = []
    if os.environ.get("BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY") != "1":
        missing.append("BIOSCRIPT_RUN_NATIVE_FASTQ_PARITY=1")
    if not MUC1_REFERENCE.exists():
        missing.append(str(MUC1_REFERENCE))
    missing_cases = [
        str(path)
        for pair in REPRESENTATIVE_FASTQ_CASES.values()
        for path in pair
        if not path.exists()
    ]
    missing.extend(missing_cases)
    missing_outputs = [str(path) for path in EXPECTED_OUTPUTS if not path.exists()]
    if missing_outputs:
        preview = ", ".join(missing_outputs[:3])
        remaining = len(missing_outputs) - min(len(missing_outputs), 3)
        suffix = f", plus {remaining} more" if remaining else ""
        missing.append(f"native FASTQ expected outputs: {preview}{suffix}")
    try:
        import_native_module()
    except Exception as exc:
        missing.append(f"bioscript._native importable ({exc})")
    if missing:
        raise unittest.SkipTest(
            "VNtyper native FASTQ pipeline prerequisites are missing: " + "; ".join(missing)
        )
    return {
        "manifest": manifest,
        "muc1_reference": str(MUC1_REFERENCE),
        "expected_outputs": [str(path) for path in EXPECTED_OUTPUTS],
        "fastq_cases": {
            label: (str(pair[0]), str(pair[1]))
            for label, pair in REPRESENTATIVE_FASTQ_CASES.items()
        },
    }


def require_samtools_fastq_oracle_prerequisites():
    """Skip unless native FASTQ extraction can be compared against samtools."""
    manifest = require_test_data(check_md5=False)
    missing = []
    if os.environ.get("BIOSCRIPT_RUN_SAMTOOLS_ORACLE") != "1":
        missing.append("BIOSCRIPT_RUN_SAMTOOLS_ORACLE=1")
    samtools_path = which_tool("samtools")
    if samtools_path is None:
        missing.append("samtools on PATH or in ports/vntyper/test-data/tools/local/bin")
    missing_cases = [
        str(path)
        for bam in REPRESENTATIVE_BAM_CASES.values()
        for path in [bam, Path(f"{bam}.bai")]
        if not path.exists()
    ]
    missing.extend(missing_cases)
    try:
        import_native_module()
    except Exception as exc:
        missing.append(f"bioscript._native importable ({exc})")
    if missing:
        raise unittest.SkipTest(
            "VNtyper samtools FASTQ oracle prerequisites are missing: " + "; ".join(missing)
        )
    return {
        "manifest": manifest,
        "samtools": samtools_path,
        "tool_path": str(LOCAL_TOOL_BIN),
        "bam_cases": {label: str(path) for label, path in REPRESENTATIVE_BAM_CASES.items()},
    }


def which_tool(name):
    path = shutil.which(name)
    if path is not None:
        return path
    local = LOCAL_TOOL_BIN / name
    if local.exists() and os.access(local, os.X_OK):
        return str(local)
    return None


def import_native_module():
    python_root = ROOT / "python"
    if str(python_root) not in sys.path:
        sys.path.insert(0, str(python_root))
    import bioscript._native as native

    return native


def load_manifest():
    with UPSTREAM_CONFIG.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    entries = []
    for resource in config.get("file_resources", []):
        local_path = resource["local_path"]
        prefix = "tests/data"
        if local_path == prefix:
            relative_dir = Path()
        elif local_path.startswith(prefix + "/"):
            relative_dir = Path(local_path[len(prefix) + 1 :])
        else:
            relative_dir = Path(local_path)
        entries.append(
            {
                "path": DATA_ROOT / relative_dir / resource["filename"],
                "filename": resource["filename"],
                "md5sum": resource["md5sum"],
            }
        )
    return entries


def validate_manifest(check_md5=False):
    missing = []
    mismatched = []
    present = 0
    for entry in load_manifest():
        path = entry["path"]
        if not path.exists():
            missing.append(str(path))
            continue
        present += 1
        if check_md5:
            actual = md5(path)
            if actual.lower() != entry["md5sum"].lower():
                mismatched.append(
                    {
                        "path": str(path),
                        "expected": entry["md5sum"],
                        "actual": actual,
                    }
                )
    return {
        "present": present,
        "missing": missing,
        "mismatched": mismatched,
    }


def md5(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    result = validate_manifest(check_md5=False)
    print(json.dumps(result, indent=2, sort_keys=True))
