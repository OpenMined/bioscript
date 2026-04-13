from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import yaml

from assay_loader import AssayPackage, VariantDefinition, load_assay_package


def as_string(value: Any) -> str | None:
    return value if isinstance(value, str) and value.strip() else None


def as_string_list(value: Any) -> list[str]:
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, str) and item.strip()]


def variant_definition_to_record(variant: VariantDefinition) -> dict[str, Any]:
    data = yaml.safe_load(variant.source)
    if not isinstance(data, dict):
        raise RuntimeError(f"{variant.name} did not deserialize to a YAML mapping")

    identifiers = data.get("identifiers", {}) or {}
    coordinates = data.get("coordinates", {}) or {}
    alleles = data.get("alleles", {}) or {}
    findings = data.get("findings", []) or []
    first_finding = findings[0] if isinstance(findings, list) and findings and isinstance(findings[0], dict) else {}

    return {
        "name": variant.name,
        "gene": as_string(data.get("gene")),
        "summary": as_string(data.get("summary")),
        "rsids": as_string_list(identifiers.get("rsids")),
        "kind": as_string(alleles.get("kind")) or "snv",
        "ref": as_string(alleles.get("ref")),
        "alts": as_string_list(alleles.get("alts")),
        "deletion_length": alleles.get("deletion_length"),
        "grch37": coordinates.get("grch37"),
        "grch38": coordinates.get("grch38"),
        "note": as_string(first_finding.get("notes")),
        "fields": variant.fields,
    }


def assay_package_to_intermediate(package: AssayPackage) -> dict[str, Any]:
    manifest = package.manifest
    metadata = manifest.get("metadata", {}) if isinstance(manifest.get("metadata"), dict) else {}
    package_block = manifest.get("package", {}) if isinstance(manifest.get("package"), dict) else {}
    ui = manifest.get("ui", {}) if isinstance(manifest.get("ui"), dict) else {}
    compatibility = manifest.get("compatibility", {}) if isinstance(manifest.get("compatibility"), dict) else {}
    privacy = manifest.get("privacy", {}) if isinstance(manifest.get("privacy"), dict) else {}

    return {
        "schema": "bioscript:assay-intermediate",
        "version": "1.0",
        "assay": {
            "id": as_string(manifest.get("assay_id")),
            "label": as_string(manifest.get("label")),
            "summary": as_string(manifest.get("summary")),
            "category": as_string(metadata.get("category")),
            "tags": as_string_list(metadata.get("tags")),
            "disclaimer": as_string(metadata.get("disclaimer")),
            "package_version": as_string(package_block.get("assay_version")) or as_string(manifest.get("version")),
            "source_of_truth": as_string(package_block.get("source_of_truth")),
        },
        "implementation": {
            "kind": package.implementation_kind,
            "script_path": str(package.script_path.relative_to(package.root)) if package.script_path else None,
        },
        "ui": {
            "template": as_string(ui.get("template")),
            "version": as_string(ui.get("version")),
        },
        "compatibility": {
            "works_with": as_string_list(compatibility.get("works_with")),
            "assemblies": as_string_list(compatibility.get("assemblies")),
            "notes": as_string_list(compatibility.get("notes"))
            if isinstance(compatibility.get("notes"), list)
            else ([compatibility["notes"]] if as_string(compatibility.get("notes")) else []),
        },
        "privacy": {
            "mode": as_string(privacy.get("mode")),
            "uploads_data": bool(privacy.get("uploads_data")),
            "stores_results_locally": bool(privacy.get("stores_results_locally")),
            "external_urls": as_string_list(privacy.get("external_urls")),
        },
        "runnable_variants": [variant_definition_to_record(entry.variant) for entry in package.runnable_variants],
        "unsupported_variants": [
            {
                **variant_definition_to_record(entry.variant),
                "reason": entry.reason,
            }
            for entry in package.unsupported_variants
        ],
    }


def load_assay_package_intermediate(assay_path: Path) -> dict[str, Any]:
    return assay_package_to_intermediate(load_assay_package(assay_path))


def load_assay_package_intermediate_json(assay_path: Path) -> str:
    return json.dumps(load_assay_package_intermediate(assay_path), indent=2, sort_keys=True)


def write_assay_package_intermediate(assay_path: Path, output_path: Path | None = None) -> Path:
    target_path = output_path if output_path is not None else assay_path / "assay.intermediate.json"
    target_path.write_text(load_assay_package_intermediate_json(assay_path) + "\n", encoding="utf-8")
    return target_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate a bioscript assay intermediate artifact.")
    parser.add_argument("assay_path", type=Path, help="Path to the assay package directory.")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Optional output path. Defaults to <assay_path>/assay.intermediate.json.",
    )
    args = parser.parse_args()

    output_path = write_assay_package_intermediate(args.assay_path.resolve(), args.output.resolve() if args.output else None)
    print(output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
