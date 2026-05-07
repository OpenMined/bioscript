use std::{
    collections::BTreeMap,
    fmt::Write as _,
    path::{Path, PathBuf},
    time::Duration,
};

use bioscript_core::{Assembly, OBSERVATION_TSV_HEADERS, VariantObservation};
use bioscript_formats::{
    GenotypeLoadOptions, GenotypeStore, InferredSex, InspectOptions, SexDetectionConfidence,
    SexInference, inspect_bytes as inspect_bytes_rs,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use bioscript_schema::{
    AssayManifest, PanelInterpretation, PanelManifest, VariantManifest, load_assay_manifest_text,
    load_panel_manifest_text, load_variant_manifest_text,
};
use monty::{MontyObject, ResourceLimits};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[path = "report_render.rs"]
mod report_render;
#[path = "report_helpers.rs"]
mod report_helpers;
#[path = "report_workspace.rs"]
mod report_workspace;

use report_helpers::*;
use report_render::{app_report_json, match_app_findings, render_app_html_document, AppReportJsonInput};
use report_workspace::PackageWorkspace;

include!("../../bioscript-cli/src/report_matching.rs");
include!("../../bioscript-cli/src/report_html_sections.rs");
include!("../../bioscript-cli/src/report_html_analysis.rs");
include!("../../bioscript-cli/src/report_html_provenance.rs");
include!("../../bioscript-cli/src/report_html_observations.rs");
include!("../../bioscript-cli/src/report_html_pgx.rs");
include!("../../bioscript-cli/src/report_html_helpers.rs");

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
struct PackageFileInput {
    path: String,
    contents: String,
    #[serde(default)]
    source_url: Option<String>,
}

#[derive(Default, Deserialize)]
#[serde(rename_all = "camelCase")]
struct ReportOptionsInput {
    #[serde(default = "default_analysis_max_duration_ms")]
    analysis_max_duration_ms: u64,
    #[serde(default)]
    detect_sex: bool,
    #[serde(default)]
    filters: Vec<String>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ReportArtifactOutput {
    name: String,
    path: String,
    mime_type: String,
    text: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct ReportRunOutput {
    artifacts: Vec<ReportArtifactOutput>,
    duration_ms: u128,
    text_output: String,
}

#[wasm_bindgen(js_name = runPackageReportBytes)]
pub fn run_package_report_bytes(
    manifest_path: &str,
    package_files_json: &str,
    input_name: &str,
    input_bytes: &[u8],
    options_json: Option<String>,
) -> Result<String, JsError> {
    let started_ms = js_sys::Date::now();
    let package_files: Vec<PackageFileInput> = serde_json::from_str(package_files_json)
        .map_err(|err| JsError::new(&format!("invalid package files JSON: {err}")))?;
    let options = match options_json {
        Some(text) if !text.is_empty() => serde_json::from_str(&text)
            .map_err(|err| JsError::new(&format!("invalid report options JSON: {err}")))?,
        _ => ReportOptionsInput::default(),
    };
    let workspace = PackageWorkspace::new(package_files)?;
    let participant_id = participant_id_from_name(input_name);
    let assay_id = app_assay_id(Path::new(manifest_path))?;
    let manifest_metadata = workspace.report_manifest_metadata(manifest_path)?;
    let findings = workspace.load_manifest_findings(manifest_path)?;
    let provenance = workspace.load_manifest_provenance_links(manifest_path)?;
    let inspect_options = InspectOptions {
        input_index: None,
        reference_file: None,
        reference_index: None,
        detect_sex: options.detect_sex,
    };
    let input_inspection = inspect_bytes_rs(input_name, input_bytes, &inspect_options)
        .map_err(|err| JsError::new(&format!("inspect input failed: {err:?}")))?;
    let mut loader = GenotypeLoadOptions::default();
    loader.assembly = input_inspection.assembly;
    loader.inferred_sex = input_inspection
        .inferred_sex
        .as_ref()
        .map(|inference| inference.sex);
    let store = GenotypeStore::from_bytes(input_name, input_bytes)
        .map_err(|err| JsError::new(&format!("load genotypes failed: {err:?}")))?;
    let rows = workspace.run_manifest_rows(manifest_path, &store, &participant_id, &options.filters)?;
    let observations = rows
        .iter()
        .map(|row| {
            workspace.app_observation_from_manifest_row(
                row,
                &assay_id,
                input_inspection.inferred_sex.as_ref(),
                input_inspection.assembly,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let analyses = workspace.run_manifest_analyses(
        manifest_path,
        input_name,
        input_bytes,
        &participant_id,
        &loader,
        &options,
    )?;
    let matched_findings = match_app_findings(&findings, &observations, &analyses);
    let reports = vec![app_report_json(AppReportJsonInput {
        assay_id: &assay_id,
        participant_id: &participant_id,
        input_file_name: input_name,
        observations: &observations,
        analyses: &analyses,
        findings: &matched_findings,
        provenance: &provenance,
        input_inspection: Some(&input_inspection),
        manifest_metadata: &manifest_metadata,
    })];
    let observations_tsv = render_app_observations_tsv(&observations)?;
    let analysis_jsonl = render_jsonl(&analyses)?;
    let reports_jsonl = render_jsonl(&reports)?;
    let html = render_app_html_document(&observations, &reports)?;
    let text_output = format!(
        "observations: observations.tsv\nanalysis: analysis.jsonl\nreports: reports.jsonl\nhtml: index.html\n"
    );
    serde_json::to_string(&ReportRunOutput {
        artifacts: vec![
            artifact("observations.tsv", "text/tab-separated-values", observations_tsv),
            artifact("analysis.jsonl", "application/jsonl", analysis_jsonl),
            artifact("reports.jsonl", "application/jsonl", reports_jsonl),
            artifact("index.html", "text/html", html),
        ],
        duration_ms: (js_sys::Date::now() - started_ms).max(0.0) as u128,
        text_output,
    })
    .map_err(|err| JsError::new(&format!("failed to encode report output: {err}")))
}

