use bioscript_schema::load_variant_manifest_text_for_lookup;
use serde::Serialize;
use wasm_bindgen::prelude::*;

#[derive(Serialize)]
struct CompiledVariantSpecJs {
    name: String,
    chrom: String,
    start: i64,
    end: i64,
    #[serde(rename = "ref")]
    ref_base: String,
    #[serde(rename = "alt")]
    alt_base: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    rsid: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    kind: Option<String>,
}

#[wasm_bindgen(js_name = compileVariantYamlText)]
pub fn compile_variant_yaml_text(name: &str, text: &str) -> Result<String, JsError> {
    let manifest = load_variant_manifest_text_for_lookup(name, text)
        .map_err(|err| JsError::new(&format!("compile variant YAML failed: {err}")))?;
    let spec = manifest.spec;
    let ref_base = spec
        .reference
        .clone()
        .ok_or_else(|| JsError::new(&format!("variant {}: alleles.ref missing", manifest.name)))?;
    let alt_base = spec
        .alternate
        .clone()
        .ok_or_else(|| JsError::new(&format!("variant {}: alleles.alts missing", manifest.name)))?;
    let rsid = spec.rsids.first().cloned();
    let kind = spec.kind.map(|kind| {
        match kind {
            bioscript_core::VariantKind::Snp => "snv",
            bioscript_core::VariantKind::Insertion => "insertion",
            bioscript_core::VariantKind::Deletion => "deletion",
            bioscript_core::VariantKind::Indel => "indel",
            bioscript_core::VariantKind::Other => "other",
        }
        .to_owned()
    });
    let mut out = Vec::new();
    if let Some(locus) = spec.grch38 {
        out.push(CompiledVariantSpecJs {
            name: manifest.name.clone(),
            chrom: locus.chrom,
            start: locus.start,
            end: locus.end,
            ref_base: ref_base.clone(),
            alt_base: alt_base.clone(),
            rsid: rsid.clone(),
            assembly: Some("grch38".to_owned()),
            kind: kind.clone(),
        });
    }
    if let Some(locus) = spec.grch37 {
        out.push(CompiledVariantSpecJs {
            name: if out.is_empty() {
                manifest.name.clone()
            } else {
                format!("{}_grch37", manifest.name)
            },
            chrom: locus.chrom,
            start: locus.start,
            end: locus.end,
            ref_base,
            alt_base,
            rsid,
            assembly: Some("grch37".to_owned()),
            kind,
        });
    }
    if out.is_empty() {
        return Err(JsError::new(&format!(
            "variant {} has no coordinates",
            manifest.name
        )));
    }
    serde_json::to_string(&out)
        .map_err(|err| JsError::new(&format!("failed to encode compiled variant: {err}")))
}
