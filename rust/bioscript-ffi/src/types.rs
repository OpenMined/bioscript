use bioscript_core::{VariantObservation, VariantSpec};
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunFileRequest {
    pub script_path: String,
    pub root: Option<String>,
    pub input_file: Option<String>,
    pub output_file: Option<String>,
    pub participant_id: Option<String>,
    pub trace_report_path: Option<String>,
    pub timing_report_path: Option<String>,
    pub input_format: Option<String>,
    pub input_index: Option<String>,
    pub reference_file: Option<String>,
    pub reference_index: Option<String>,
    pub allow_md5_mismatch: Option<bool>,
    pub auto_index: Option<bool>,
    pub cache_dir: Option<String>,
    pub max_duration_ms: Option<u64>,
    pub max_memory_bytes: Option<usize>,
    pub max_allocations: Option<usize>,
    pub max_recursion_depth: Option<usize>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunFileResult {
    pub ok: bool,
}

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunVariantYamlRequest {
    pub yaml_path: String,
    pub genome_path: String,
    pub input_format: Option<String>,
    pub input_index: Option<String>,
    pub reference_file: Option<String>,
    pub reference_index: Option<String>,
    pub allow_md5_mismatch: Option<bool>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunVariantYamlResult {
    pub observations: Vec<VariantObservationResult>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VariantObservationResult {
    pub name: String,
    pub backend: String,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    pub reference: Option<String>,
    #[serde(rename = "alt", skip_serializing_if = "Option::is_none")]
    pub alternate: Option<String>,
    #[serde(rename = "matchedRsid", skip_serializing_if = "Option::is_none")]
    pub matched_rsid: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype: Option<String>,
    #[serde(rename = "refCount", skip_serializing_if = "Option::is_none")]
    pub ref_count: Option<u32>,
    #[serde(rename = "altCount", skip_serializing_if = "Option::is_none")]
    pub alt_count: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub depth: Option<u32>,
    #[serde(rename = "rawCounts")]
    pub raw_counts: std::collections::BTreeMap<String, u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub decision: Option<String>,
    pub evidence: Vec<String>,
}

#[derive(Debug, Clone)]
pub(crate) struct NamedVariantSpec {
    pub(crate) name: String,
    pub(crate) spec: VariantSpec,
}

pub(crate) fn observation_result(
    variant: NamedVariantSpec,
    observation: VariantObservation,
) -> VariantObservationResult {
    VariantObservationResult {
        name: variant.name,
        backend: observation.backend,
        reference: variant.spec.reference,
        alternate: variant.spec.alternate,
        matched_rsid: observation.matched_rsid,
        assembly: observation
            .assembly
            .map(super::variant_yaml::assembly_label),
        genotype: observation.genotype,
        ref_count: observation.ref_count,
        alt_count: observation.alt_count,
        depth: observation.depth,
        raw_counts: observation.raw_counts,
        decision: observation.decision,
        evidence: observation.evidence,
    }
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub(crate) struct FfiResult<T> {
    pub(crate) ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) value: Option<T>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) error: Option<String>,
}
