from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass
class VariantDefinition:
    name: str
    fields: dict[str, Any]
    source: str


@dataclass
class AssayPackage:
    root: Path
    manifest: dict[str, Any]
    implementation_kind: str
    script_path: Path | None
    variants: list[VariantDefinition]
    unsupported_variants: list[tuple[VariantDefinition, str]]


def read_yaml_dict(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise RuntimeError(f"{path} did not contain a YAML mapping")
    return data


def fmt_coord(coords: dict[str, Any]) -> str | None:
    if not coords:
        return None
    chrom = coords.get("chrom")
    if not chrom:
        return None
    pos = coords.get("pos")
    start = coords.get("start")
    end = coords.get("end")
    if pos is not None:
        return f"{chrom}:{pos}-{pos}"
    if start is not None and end is not None:
        return f"{chrom}:{start}-{end}"
    if start is not None:
        return f"{chrom}:{start}-{start}"
    return None


def yaml_to_variant_definition(yaml_path: Path) -> VariantDefinition | None:
    data = read_yaml_dict(yaml_path)
    schema = data.get("schema", "")
    if not isinstance(schema, str) or "variant" not in schema:
        return None

    coords = data.get("coordinates", {}) or {}
    alleles_block = data.get("alleles", {}) or {}
    identifiers = data.get("identifiers", {}) or {}
    rsids = identifiers.get("rsids", []) or []
    rsid = rsids[0] if rsids else None

    fields: dict[str, Any] = {}
    if rsid:
        fields["rsid"] = rsid
    grch37 = fmt_coord(coords.get("grch37", {}))
    grch38 = fmt_coord(coords.get("grch38", {}))
    if grch37:
        fields["grch37"] = grch37
    if grch38:
        fields["grch38"] = grch38

    ref = alleles_block.get("ref")
    if ref:
        fields["ref"] = ref

    alts = alleles_block.get("alts", [])
    if alts:
        canonical_alt = alleles_block.get("canonical_alt")
        if canonical_alt in alts:
            fields["alt"] = canonical_alt
        else:
            fields["alt"] = alts[0]

    kind = alleles_block.get("kind")
    if kind:
        kind_map = {
            "snv": "snp",
            "deletion": "deletion",
            "insertion": "insertion",
            "indel": "indel",
        }
        fields["kind"] = kind_map.get(kind, kind)

    deletion_length = alleles_block.get("deletion_length")
    if deletion_length is not None:
        fields["deletion_length"] = deletion_length

    motifs = alleles_block.get("motifs")
    if motifs:
        fields["motifs"] = motifs

    name_str = data.get("name", "") or yaml_path.stem
    safe_name = name_str.replace("-", "_").replace(".", "_").replace(" ", "_")
    return VariantDefinition(
        name=safe_name,
        fields=fields,
        source=yaml_path.read_text(encoding="utf-8"),
    )


def extract_yaml_variant_definitions(directory: Path) -> list[VariantDefinition]:
    variants: list[VariantDefinition] = []
    for yaml_file in sorted(directory.glob("*.yaml")):
        variant = yaml_to_variant_definition(yaml_file)
        if variant:
            variants.append(variant)
    return variants


def resolve_catalogue_variants(assay_root: Path, manifest: dict[str, Any]) -> list[VariantDefinition]:
    inputs = manifest.get("inputs", {}) or {}
    catalogue_ref = inputs.get("catalogue")
    if not isinstance(catalogue_ref, str) or not catalogue_ref:
        return extract_yaml_variant_definitions(assay_root)

    catalogue_path = (assay_root / catalogue_ref).resolve()
    catalogue = read_yaml_dict(catalogue_path)
    variants: list[VariantDefinition] = []
    for entry in catalogue.get("variants", []) or []:
        if not isinstance(entry, dict):
            continue
        variant_ref = entry.get("path")
        if not isinstance(variant_ref, str) or not variant_ref:
            continue
        variant_path = (catalogue_path.parent / variant_ref).resolve()
        variant = yaml_to_variant_definition(variant_path)
        if variant:
            variants.append(variant)
    return variants


def runtime_supports_variant(variant: VariantDefinition) -> tuple[bool, str]:
    kind = str(variant.fields.get("kind", "") or "").lower()
    if kind in {"", "snp"}:
        return True, ""
    if kind == "deletion":
        if variant.fields.get("deletion_length"):
            return True, ""
        return False, "deletion missing deletion_length"
    if kind == "insertion":
        return False, "insertions not yet supported by bioscript runtime"
    if kind == "indel":
        return False, "indels not yet supported by bioscript runtime"
    return False, f"unsupported variant kind: {kind}"


def load_assay_package(assay_path: Path) -> AssayPackage:
    manifest_path = assay_path / "assay.yaml"
    manifest = read_yaml_dict(manifest_path)
    implementation = manifest.get("implementation", {}) or {}
    implementation_kind = str(implementation.get("kind", "panel") or "panel")
    script_path: Path | None = None
    if implementation_kind == "script":
        script_ref = implementation.get("path")
        if not isinstance(script_ref, str) or not script_ref:
            raise RuntimeError(f"{manifest_path} missing implementation.path for script assay")
        script_path = (assay_path / script_ref).resolve()

    all_variants = resolve_catalogue_variants(assay_path, manifest)
    supported_variants: list[VariantDefinition] = []
    unsupported_variants: list[tuple[VariantDefinition, str]] = []
    for variant in all_variants:
        supported, reason = runtime_supports_variant(variant)
        if supported:
            supported_variants.append(variant)
        else:
            unsupported_variants.append((variant, reason))

    return AssayPackage(
        root=assay_path,
        manifest=manifest,
        implementation_kind=implementation_kind,
        script_path=script_path,
        variants=supported_variants,
        unsupported_variants=unsupported_variants,
    )
