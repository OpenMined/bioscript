use std::{fs, path::PathBuf};

use bioscript_core::{Assembly, VariantKind, VariantSpec};
use bioscript_formats::{GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore};
use bioscript_schema::load_variant_manifest_text_for_lookup;

use crate::types::{
    NamedVariantSpec, RunVariantYamlRequest, RunVariantYamlResult, observation_result,
};

/// Runs a BioScript variant YAML assay against a supported genome file.
///
/// The native desktop/mobile path uses this instead of the web WASM exports.
/// It intentionally mirrors the web variant YAML flow: compile YAML through
/// `bioscript-schema`, choose the preferred assembly-specific variant, and
/// execute lookup through `bioscript-formats`.
pub fn run_variant_yaml_request(
    request: RunVariantYamlRequest,
) -> Result<RunVariantYamlResult, String> {
    let yaml_path = PathBuf::from(&request.yaml_path);
    let yaml_text = fs::read_to_string(&yaml_path)
        .map_err(|err| format!("failed to read YAML assay {}: {err}", yaml_path.display()))?;
    let variants = compile_variant_yaml_named(
        yaml_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("variant.yaml"),
        &yaml_text,
    )?;
    let selected = select_preferred_assembly_variants(&request.genome_path, variants);
    let loader = variant_loader(&request)?;

    let genome_path = PathBuf::from(&request.genome_path);
    let store = GenotypeStore::from_file_with_options(&genome_path, &loader)
        .map_err(|err| format!("failed to load genome {}: {err}", genome_path.display()))?;
    let specs = selected
        .iter()
        .map(|variant| variant.spec.clone())
        .collect::<Vec<_>>();
    let observations = store
        .lookup_variants(&specs)
        .map_err(|err| format!("variant lookup failed: {err}"))?;

    Ok(RunVariantYamlResult {
        observations: selected
            .into_iter()
            .zip(observations)
            .map(|(variant, observation)| observation_result(variant, observation))
            .collect(),
    })
}

fn variant_loader(request: &RunVariantYamlRequest) -> Result<GenotypeLoadOptions, String> {
    let mut loader = GenotypeLoadOptions::default();
    if let Some(value) = request.input_format.as_deref() {
        if value.eq_ignore_ascii_case("auto") {
            loader.format = None;
        } else {
            loader.format = Some(
                value
                    .parse::<GenotypeSourceFormat>()
                    .map_err(|err| format!("invalid inputFormat value {value}: {err}"))?,
            );
        }
    }
    loader.input_index = request.input_index.clone().map(PathBuf::from);
    loader.reference_file = request.reference_file.clone().map(PathBuf::from);
    loader.reference_index = request.reference_index.clone().map(PathBuf::from);
    loader.allow_reference_md5_mismatch = request.allow_md5_mismatch.unwrap_or(false);
    Ok(loader)
}

fn compile_variant_yaml_named(name: &str, text: &str) -> Result<Vec<NamedVariantSpec>, String> {
    let manifest = load_variant_manifest_text_for_lookup(name, text)
        .map_err(|err| format!("compile variant YAML failed: {err}"))?;
    let mut out = Vec::new();
    if let Some(locus) = manifest.spec.grch38.clone() {
        let mut spec = manifest.spec.clone();
        spec.grch37 = None;
        spec.grch38 = Some(locus);
        out.push(NamedVariantSpec {
            name: manifest.name.clone(),
            spec,
        });
    }
    if let Some(locus) = manifest.spec.grch37.clone() {
        let mut spec = manifest.spec;
        spec.grch37 = Some(locus);
        spec.grch38 = None;
        out.push(NamedVariantSpec {
            name: if out.is_empty() {
                manifest.name.clone()
            } else {
                format!("{}_grch37", manifest.name)
            },
            spec,
        });
    }
    if out.is_empty() {
        return Err(format!("variant {} has no coordinates", manifest.name));
    }
    Ok(out)
}

fn select_preferred_assembly_variants(
    genome_name: &str,
    variants: Vec<NamedVariantSpec>,
) -> Vec<NamedVariantSpec> {
    let target = infer_genome_assembly(genome_name).unwrap_or(Assembly::Grch38);
    let mut groups: Vec<(String, Vec<NamedVariantSpec>)> = Vec::new();
    for variant in variants {
        let key = variant_group_key(&variant);
        if let Some((_, values)) = groups.iter_mut().find(|(candidate, _)| candidate == &key) {
            values.push(variant);
        } else {
            groups.push((key, vec![variant]));
        }
    }

    groups
        .into_iter()
        .flat_map(|(_, group)| select_preferred_group(group, target))
        .collect()
}

fn select_preferred_group(group: Vec<NamedVariantSpec>, target: Assembly) -> Vec<NamedVariantSpec> {
    let has_multiple_assemblies = group
        .iter()
        .filter_map(|variant| variant_assembly(&variant.spec))
        .fold(Vec::new(), |mut assemblies, assembly| {
            if !assemblies.contains(&assembly) {
                assemblies.push(assembly);
            }
            assemblies
        })
        .len()
        > 1;
    if group.len() <= 1 || !has_multiple_assemblies {
        return group;
    }
    let fallback = group[0].clone();
    vec![
        group
            .iter()
            .find(|variant| variant_assembly(&variant.spec) == Some(target))
            .or_else(|| {
                group
                    .iter()
                    .find(|variant| variant_assembly(&variant.spec) == Some(Assembly::Grch38))
            })
            .cloned()
            .unwrap_or(fallback),
    ]
}

fn infer_genome_assembly(name: &str) -> Option<Assembly> {
    let lower = name.to_ascii_lowercase();
    if lower.contains("grch38") || lower.contains("hg38") {
        return Some(Assembly::Grch38);
    }
    if lower.contains("grch37") || lower.contains("hg19") {
        return Some(Assembly::Grch37);
    }
    None
}

fn variant_group_key(variant: &NamedVariantSpec) -> String {
    if let Some(rsid) = variant.spec.rsids.first() {
        return format!(
            "{}|{}|{}|{}",
            rsid.to_ascii_lowercase(),
            variant
                .spec
                .reference
                .as_deref()
                .unwrap_or("")
                .to_ascii_uppercase(),
            variant
                .spec
                .alternate
                .as_deref()
                .unwrap_or("")
                .to_ascii_uppercase(),
            variant_kind_label(variant.spec.kind)
        );
    }
    variant
        .name
        .trim_end_matches("_grch37")
        .trim_end_matches("_grch38")
        .to_ascii_lowercase()
}

fn variant_assembly(spec: &VariantSpec) -> Option<Assembly> {
    if spec.grch37.is_some() {
        return Some(Assembly::Grch37);
    }
    if spec.grch38.is_some() {
        return Some(Assembly::Grch38);
    }
    None
}

pub(crate) fn assembly_label(assembly: Assembly) -> String {
    match assembly {
        Assembly::Grch37 => "grch37",
        Assembly::Grch38 => "grch38",
    }
    .to_owned()
}

fn variant_kind_label(kind: Option<VariantKind>) -> &'static str {
    match kind {
        Some(VariantKind::Snp) => "snv",
        Some(VariantKind::Insertion) => "insertion",
        Some(VariantKind::Deletion) => "deletion",
        Some(VariantKind::Indel) => "indel",
        Some(VariantKind::Other) => "other",
        None => "",
    }
}
