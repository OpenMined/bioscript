"""Authoritative VNtyper fixture expectations, sourced from upstream.

Upstream's `tests/test_data_config.json` (`integration_tests.bam_tests`)
defines, per fixture, the `kestrel_assertions`: expected Confidence plus
Alt/ActiveRegion depth and Depth_Score with a tolerance percentage. This
module reads that file directly so the harness can never drift from
upstream's own assertions, and applies upstream's exact comparison rules
(see `tests/integration/test_pipeline_integration.py`):

- Confidence "Negative": no positive call expected.
- Confidence ending in "*": actual must start with the prefix.
- otherwise: exact Confidence match.
- depth fields: |actual - expected| <= |expected| * tol%/100.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
UPSTREAM_CONFIG = (
    ROOT / "ports" / "vntyper" / "vntyper" / "tests" / "test_data_config.json"
)
TEST_DATA = ROOT / "ports" / "vntyper" / "test-data"


def _bam_basename(bam_path):
    return Path(str(bam_path)).name


def load_expectations():
    """Return {fixture_stem: expectation dict} for every asserted fixture.

    fixture_stem is e.g. "example_66bf_hg19_subset" (no .bam). Only fixtures
    whose BAM is present under test-data are included.
    """
    config = json.loads(UPSTREAM_CONFIG.read_text(encoding="utf-8"))
    out = {}
    for test in config.get("integration_tests", {}).get("bam_tests", []):
        bam_name = _bam_basename(test.get("bam", ""))
        if not bam_name.endswith(".bam"):
            continue
        stem = bam_name[: -len(".bam")]
        bam_path = TEST_DATA / bam_name
        if not bam_path.exists():
            continue
        ka = test.get("kestrel_assertions", {})
        out[stem] = {
            "test_name": test.get("test_name"),
            "reference_assembly": test.get("reference_assembly", "hg19"),
            "bam": str(bam_path),
            "confidence": ka.get("Confidence"),
            "alt_depth": _assertion(ka.get("Estimated_Depth_AlternateVariant")),
            "region_depth": _assertion(
                ka.get("Estimated_Depth_Variant_ActiveRegion")
            ),
            "depth_score": _assertion(ka.get("Depth_Score")),
            "is_negative": ka.get("Confidence") == "Negative",
        }
    return out


def _assertion(value):
    """Normalize an upstream assertion into {value, tol} or {value: None}."""
    if isinstance(value, dict):
        v = value.get("value")
        if isinstance(v, str) and v == "None":
            return {"value": None, "tol": value.get("tolerance_percentage", 5)}
        return {"value": v, "tol": value.get("tolerance_percentage", 5)}
    if isinstance(value, str) and value == "None":
        return {"value": None, "tol": 5}
    if value is None:
        return None
    return {"value": value, "tol": 5}


def _within(actual, expected_value, tol_pct):
    if expected_value is None:
        return actual is None
    if actual is None:
        return False
    allowed = abs(float(expected_value)) * (float(tol_pct) / 100.0)
    return abs(float(actual) - float(expected_value)) <= allowed


def confidence_matches(expected, actual):
    """Upstream rule: '*' suffix → prefix match; else exact. Negative → n/a."""
    if expected is None:
        return True
    if expected == "Negative":
        return True  # handled separately via "no positive call"
    if actual is None:
        return False
    if expected.endswith("*"):
        return str(actual).startswith(expected[:-1])
    return str(actual) == expected


def evaluate(expectation, called):
    """Compare a called variant (or None) against one fixture expectation.

    `called` is the top passing row dict (or None when no positive call).
    Returns (ok: bool, reasons: list[str]).
    """
    reasons = []
    if expectation["is_negative"]:
        if called is not None:
            reasons.append(
                f"expected Negative (no call) but got "
                f"{called.get('Confidence')} "
                f"alt={called.get('Estimated_Depth_AlternateVariant')}"
            )
        return (not reasons, reasons)

    # Positive fixture: a variant must be called.
    if called is None:
        reasons.append(
            f"expected {expectation['confidence']} call but no variant passed"
        )
        return (False, reasons)

    if not confidence_matches(expectation["confidence"], called.get("Confidence")):
        reasons.append(
            f"Confidence: expected {expectation['confidence']!r}, "
            f"got {called.get('Confidence')!r}"
        )

    checks = (
        ("alt_depth", "Estimated_Depth_AlternateVariant"),
        ("region_depth", "Estimated_Depth_Variant_ActiveRegion"),
        ("depth_score", "Depth_Score"),
    )
    for key, field in checks:
        spec = expectation.get(key)
        if not spec:
            continue
        actual = _to_float(called.get(field))
        if not _within(actual, spec["value"], spec["tol"]):
            reasons.append(
                f"{field}: expected ~{spec['value']} "
                f"(±{spec['tol']}%), got {actual}"
            )
    return (not reasons, reasons)


def _to_float(value):
    try:
        if value is None or value == "" or value == "None":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None
