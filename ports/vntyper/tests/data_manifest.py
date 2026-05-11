"""VNtyper test-data manifest and validator.

The copied large data lives in `ports/vntyper/test-data`. Upstream VNtyper's
manifest expects paths under `tests/data`, so this helper remaps those entries
into the BioScript port tree and can optionally verify MD5 checksums.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
UPSTREAM_CONFIG = ROOT / "ports" / "vntyper" / "vntyper" / "tests" / "test_data_config.json"
DATA_ROOT = ROOT / "ports" / "vntyper" / "test-data"


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
