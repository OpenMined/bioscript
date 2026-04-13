from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

CAPABILITIES_SPEC_PATH = Path(__file__).resolve().with_name("runtime_capabilities.json")
CAPABILITIES_SPEC = json.loads(CAPABILITIES_SPEC_PATH.read_text(encoding="utf-8"))


@dataclass(frozen=True)
class RuntimeCapabilityAssessment:
    kind: str
    reason: str | None
    supported: bool


def normalize_variant_kind(value: Any) -> str:
    kind = str(value or "").strip().lower()
    aliases = CAPABILITIES_SPEC.get("kindAliases", {})
    if not isinstance(aliases, dict):
        return kind
    return str(aliases.get(kind, kind))


def assess_variant_runtime_support(fields: dict[str, Any]) -> RuntimeCapabilityAssessment:
    kind = normalize_variant_kind(fields.get("kind"))
    kinds = CAPABILITIES_SPEC.get("kinds", {})
    kind_spec = kinds.get(kind) if isinstance(kinds, dict) else None

    if not isinstance(kind_spec, dict):
        template = str(CAPABILITIES_SPEC.get("defaultUnsupportedReason", "unsupported variant kind: {kind}"))
        return RuntimeCapabilityAssessment(
            kind=kind,
            supported=False,
            reason=template.format(kind=kind),
        )

    required_fields = kind_spec.get("requiredFields", [])
    if isinstance(required_fields, list):
        for field in required_fields:
            if not fields.get(field):
                template = str(kind_spec.get("missingFieldReason", "{kind} missing {field}"))
                return RuntimeCapabilityAssessment(
                    kind=kind,
                    supported=False,
                    reason=template.format(kind=kind, field=field),
                )

    supported = bool(kind_spec.get("supported"))
    if supported:
        return RuntimeCapabilityAssessment(kind=kind, supported=True, reason=None)

    return RuntimeCapabilityAssessment(
        kind=kind,
        supported=False,
        reason=str(kind_spec.get("reason", CAPABILITIES_SPEC.get("defaultUnsupportedReason", "unsupported variant kind: {kind}"))).format(kind=kind),
    )
