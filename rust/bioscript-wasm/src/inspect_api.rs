use std::path::PathBuf;

use bioscript_formats::{
    inspect_bytes as inspect_bytes_rs, DetectedKind, DetectionConfidence, FileContainer,
    FileInspection, InspectOptions, SourceMetadata,
};
use bioscript_schema::resolve_remote_resource_text as resolve_remote_resource_text_rs;
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[derive(Default, Deserialize)]
struct InspectOptionsJs {
    input_index: Option<String>,
    reference_file: Option<String>,
    reference_index: Option<String>,
}

/// Classify bytes as a known genomic file. Mirrors `bioscript-formats::inspect::inspect_bytes`.
/// Returns JSON matching the `Inspection` shape the app already uses.
#[wasm_bindgen(js_name = inspectBytes)]
pub fn inspect_bytes(
    name: &str,
    bytes: &[u8],
    options_json: Option<String>,
) -> Result<String, JsError> {
    let options_js: InspectOptionsJs = match options_json {
        Some(text) if !text.is_empty() => serde_json::from_str(&text)
            .map_err(|err| JsError::new(&format!("invalid InspectOptions JSON: {err}")))?,
        _ => InspectOptionsJs::default(),
    };
    let options = InspectOptions {
        input_index: options_js.input_index.map(PathBuf::from),
        reference_file: options_js.reference_file.map(PathBuf::from),
        reference_index: options_js.reference_index.map(PathBuf::from),
    };

    let inspection = inspect_bytes_rs(name, bytes, &options)
        .map_err(|err| JsError::new(&format!("inspect_bytes failed: {err:?}")))?;

    let resp = InspectionJs::from(inspection);
    serde_json::to_string(&resp)
        .map_err(|err| JsError::new(&format!("failed to encode response: {err}")))
}

/// Classify a fetched remote resource and return dependency requirements.
///
/// Network access stays in the host app so each platform can prompt before
/// fetching. The schema/type/dependency logic lives here so web, mobile,
/// desktop, and CLI share one implementation.
#[wasm_bindgen(js_name = resolveRemoteResourceText)]
pub fn resolve_remote_resource_text(
    source_url: &str,
    name: &str,
    text: &str,
) -> Result<String, JsError> {
    let resolved = resolve_remote_resource_text_rs(source_url, name, text)
        .map_err(|err| JsError::new(&format!("resolve remote resource failed: {err}")))?;
    serde_json::to_string(&resolved)
        .map_err(|err| JsError::new(&format!("failed to encode response: {err}")))
}

// Wire types — we flatten the Rust FileInspection into the shape the app's
// TS Inspection type already expects (matches widgets/FilePicker/types.ts).
#[derive(Serialize)]
struct InspectionJs {
    #[serde(rename = "fileName")]
    file_name: String,
    container: &'static str,
    #[serde(rename = "detectedKind")]
    detected_kind: &'static str,
    confidence: &'static str,
    assembly: Option<&'static str>,
    phased: Option<bool>,
    source: Option<SourceJs>,
    #[serde(rename = "selectedEntry", skip_serializing_if = "Option::is_none")]
    selected_entry: Option<String>,
    #[serde(rename = "hasIndex", skip_serializing_if = "Option::is_none")]
    has_index: Option<bool>,
    #[serde(rename = "referenceMatches", skip_serializing_if = "Option::is_none")]
    reference_matches: Option<bool>,
    evidence: Vec<String>,
    warnings: Vec<String>,
    #[serde(rename = "durationMs")]
    duration_ms: u128,
}

#[derive(Serialize)]
struct SourceJs {
    vendor: String,
    #[serde(rename = "platformVersion", skip_serializing_if = "Option::is_none")]
    platform_version: Option<String>,
    confidence: &'static str,
    evidence: Vec<String>,
}

impl From<FileInspection> for InspectionJs {
    fn from(i: FileInspection) -> Self {
        InspectionJs {
            file_name: i.path.display().to_string(),
            container: render_container(i.container),
            detected_kind: render_kind(i.detected_kind),
            confidence: render_confidence(i.confidence),
            assembly: i.assembly.map(render_assembly),
            phased: i.phased,
            source: i.source.map(SourceJs::from),
            selected_entry: i.selected_entry,
            has_index: i.has_index,
            reference_matches: i.reference_matches,
            evidence: i.evidence,
            warnings: i.warnings,
            duration_ms: i.duration_ms,
        }
    }
}

impl From<SourceMetadata> for SourceJs {
    fn from(s: SourceMetadata) -> Self {
        SourceJs {
            vendor: s.vendor.unwrap_or_default(),
            platform_version: s.platform_version,
            confidence: render_confidence(s.confidence),
            evidence: s.evidence,
        }
    }
}

fn render_container(c: FileContainer) -> &'static str {
    match c {
        FileContainer::Plain => "plain",
        FileContainer::Zip => "zip",
    }
}

fn render_kind(k: DetectedKind) -> &'static str {
    match k {
        DetectedKind::GenotypeText => "genotype_text",
        DetectedKind::Vcf => "vcf",
        DetectedKind::AlignmentCram => "alignment_cram",
        DetectedKind::AlignmentBam => "alignment_bam",
        DetectedKind::ReferenceFasta => "reference_fasta",
        DetectedKind::Unknown => "unknown",
    }
}

fn render_confidence(c: DetectionConfidence) -> &'static str {
    match c {
        DetectionConfidence::Authoritative => "authoritative",
        DetectionConfidence::StrongHeuristic => "strong_heuristic",
        DetectionConfidence::WeakHeuristic => "weak_heuristic",
        DetectionConfidence::Unknown => "unknown",
    }
}

pub(crate) fn render_assembly(a: bioscript_core::Assembly) -> &'static str {
    match a {
        bioscript_core::Assembly::Grch37 => "grch37",
        bioscript_core::Assembly::Grch38 => "grch38",
    }
}
