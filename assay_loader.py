from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from runtime_capabilities import assess_variant_runtime_support


@dataclass
class VariantDefinition:
    name: str
    fields: dict[str, Any]
    source: str


@dataclass
class RunnableVariantEntry:
    variant: VariantDefinition


@dataclass
class UnsupportedVariantEntry:
    reason: str
    variant: VariantDefinition


@dataclass
class AssayPackage:
    implementation_kind: str
    manifest: dict[str, Any]
    root: Path
    runnable_variants: list[RunnableVariantEntry]
    script_path: Path | None
    unsupported_variants: list[UnsupportedVariantEntry]

    @property
    def variants(self) -> list[VariantDefinition]:
        return [entry.variant for entry in self.runnable_variants]


VALID_IMPLEMENTATION_KINDS = {"panel", "script"}


def read_yaml_dict(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise RuntimeError(f"{path} does not exist")
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


def validate_manifest(manifest_path: Path, manifest: dict[str, Any]) -> None:
    schema = manifest.get("schema")
    if schema != "bioscript:assay":
        raise RuntimeError(f"{manifest_path} must declare schema: bioscript:assay")

    assay_id = manifest.get("assay_id")
    if not isinstance(assay_id, str) or not assay_id.strip():
        raise RuntimeError(f"{manifest_path} missing assay_id")

    implementation = manifest.get("implementation")
    if not isinstance(implementation, dict):
        raise RuntimeError(f"{manifest_path} missing implementation block")

    kind = implementation.get("kind")
    if not isinstance(kind, str) or kind not in VALID_IMPLEMENTATION_KINDS:
        allowed = ", ".join(sorted(VALID_IMPLEMENTATION_KINDS))
        raise RuntimeError(f"{manifest_path} implementation.kind must be one of: {allowed}")


def validate_catalogue(catalogue_path: Path, catalogue: dict[str, Any]) -> list[dict[str, Any]]:
    schema = catalogue.get("schema")
    if schema != "bioscript:catalogue":
        raise RuntimeError(f"{catalogue_path} must declare schema: bioscript:catalogue")

    entries = catalogue.get("variants")
    if not isinstance(entries, list):
        raise RuntimeError(f"{catalogue_path} missing variants list")
    if not entries:
        raise RuntimeError(f"{catalogue_path} variants list is empty")

    validated: list[dict[str, Any]] = []
    for idx, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise RuntimeError(f"{catalogue_path} variants[{idx}] must be a mapping")
        variant_ref = entry.get("path")
        if not isinstance(variant_ref, str) or not variant_ref:
            raise RuntimeError(f"{catalogue_path} variants[{idx}] missing path")
        validated.append(entry)
    return validated


def validate_variant_fields(variant_path: Path, variant: VariantDefinition) -> None:
    fields = variant.fields
    has_rsid = bool(fields.get("rsid"))
    has_coords = bool(fields.get("grch37") or fields.get("grch38"))
    if not has_rsid and not has_coords:
        raise RuntimeError(f"{variant_path} must declare at least one rsid or genomic coordinate")

    kind = str(fields.get("kind", "") or "").lower()
    if kind in {"", "snp"}:
        if not fields.get("ref"):
            raise RuntimeError(f"{variant_path} SNP variant missing ref allele")
        if not fields.get("alt"):
            raise RuntimeError(f"{variant_path} SNP variant missing alt allele")
        return

    if kind == "deletion":
        if not fields.get("ref"):
            raise RuntimeError(f"{variant_path} deletion variant missing ref allele")
        if not fields.get("alt"):
            raise RuntimeError(f"{variant_path} deletion variant missing alt allele")
        if not fields.get("deletion_length"):
            raise RuntimeError(f"{variant_path} deletion variant missing deletion_length")
        return

    if kind in {"insertion", "indel", "other"}:
        return

    raise RuntimeError(f"{variant_path} declares unsupported variant kind: {kind}")


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
        variants = extract_yaml_variant_definitions(assay_root)
        if not variants:
            raise RuntimeError(f"{assay_root} does not declare a catalogue and contains no root-level variant YAML files")
        return variants

    catalogue_path = (assay_root / catalogue_ref).resolve()
    catalogue = read_yaml_dict(catalogue_path)
    entries = validate_catalogue(catalogue_path, catalogue)
    variants: list[VariantDefinition] = []
    for entry in entries:
        variant_ref = entry.get("path")
        variant_path = (catalogue_path.parent / variant_ref).resolve()
        if not variant_path.exists():
            raise RuntimeError(f"{catalogue_path} references missing variant file: {variant_ref}")
        variant = yaml_to_variant_definition(variant_path)
        if variant is None:
            raise RuntimeError(f"{variant_path} is not a bioscript:variant record")
        validate_variant_fields(variant_path, variant)
        variants.append(variant)
    return variants


def load_assay_package(assay_path: Path) -> AssayPackage:
    manifest_path = assay_path / "assay.yaml"
    manifest = read_yaml_dict(manifest_path)
    validate_manifest(manifest_path, manifest)
    implementation = manifest.get("implementation", {}) or {}
    implementation_kind = str(implementation.get("kind", "panel") or "panel")
    script_path: Path | None = None
    if implementation_kind == "script":
        script_ref = implementation.get("path")
        if not isinstance(script_ref, str) or not script_ref:
            raise RuntimeError(f"{manifest_path} missing implementation.path for script assay")
        script_path = (assay_path / script_ref).resolve()
        if not script_path.exists():
            raise RuntimeError(f"{manifest_path} references missing script: {script_ref}")

    all_variants = resolve_catalogue_variants(assay_path, manifest)
    supported_variants: list[RunnableVariantEntry] = []
    unsupported_variants: list[UnsupportedVariantEntry] = []
    for variant in all_variants:
        assessment = assess_variant_runtime_support(variant.fields)
        if assessment.supported:
            supported_variants.append(RunnableVariantEntry(variant=variant))
        else:
            unsupported_variants.append(
                UnsupportedVariantEntry(
                    variant=variant,
                    reason=assessment.reason or "unsupported by bioscript runtime",
                )
            )

    return AssayPackage(
        implementation_kind=implementation_kind,
        manifest=manifest,
        root=assay_path,
        runnable_variants=supported_variants,
        script_path=script_path,
        unsupported_variants=unsupported_variants,
    )
