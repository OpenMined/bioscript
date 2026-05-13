use std::{
    collections::BTreeMap,
    fmt::Write as _,
    path::{Path, PathBuf},
    time::Duration,
};

use bioscript_core::{
    Assembly, GenomicLocus, OBSERVATION_TSV_HEADERS, RuntimeError, VariantKind, VariantObservation,
    VariantSpec,
};
use bioscript_formats::{
    DetectedKind, DetectionConfidence, FileContainer, FileInspection, GenotypeLoadOptions,
    GenotypeStore, InferredSex, InspectOptions, SexDetectionConfidence, SexInference,
    inspect_bytes as inspect_bytes_rs,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use bioscript_schema::{
    AssayManifest, PanelInterpretation, PanelManifest, VariantManifest, load_assay_manifest_text,
    load_panel_manifest_text, load_variant_manifest_text,
};
use monty::{MontyObject, ResourceLimits};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[path = "report_helpers.rs"]
mod report_helpers;
#[path = "report_input_inspection.rs"]
mod report_input_inspection;
#[path = "report_lookup.rs"]
mod report_lookup;
#[path = "report_render.rs"]
mod report_render;
#[path = "report_workspace.rs"]
mod report_workspace;

use report_helpers::*;
use report_input_inspection::{
    decompress_vcf_head_lines, explicit_sex_from_options, inspect_head_via_js_reader,
    vcf_sex_via_tabix,
};
use report_lookup::{BamReportLookup, CramReportLookup, VcfReportLookup};
use report_render::{
    AppReportJsonInput, app_report_json, match_app_findings, render_app_html_document,
};
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
pub(super) struct PackageFileInput {
    path: String,
    contents: String,
    #[serde(default)]
    source_url: Option<String>,
}

#[derive(Default, Deserialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct ReportOptionsInput {
    #[serde(default = "default_analysis_max_duration_ms")]
    analysis_max_duration_ms: u64,
    #[serde(default)]
    detect_sex: bool,
    #[serde(default)]
    filters: Vec<String>,
    #[serde(default)]
    output_dir: Option<String>,
    /// Optional explicit sample sex (mirrors the CLI's `--sample-sex` flag).
    /// When set, takes precedence over inference: the report carries
    /// `method=explicit_sample_sex` like the CLI.
    #[serde(default)]
    sample_sex: Option<String>,
}

fn analysis_cache_observations(
    manifest_observations: &[VariantObservation],
    app_observations: &[serde_json::Value],
) -> Vec<VariantObservation> {
    manifest_observations
        .iter()
        .map(|observation| {
            let mut observation = observation.clone();
            if let Some(app_observation) = matching_app_observation(&observation, app_observations)
                && let Some(genotype_display) = app_observation
                    .get("genotype_display")
                    .and_then(serde_json::Value::as_str)
                    .filter(|value| !value.is_empty() && *value != "??")
            {
                observation.genotype = Some(genotype_display.to_owned());
            }
            observation
        })
        .collect()
}

fn matching_app_observation<'a>(
    observation: &VariantObservation,
    app_observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let matched_rsid = observation.matched_rsid.as_deref()?;
    app_observations.iter().find(|app_observation| {
        app_observation
            .get("rsid")
            .and_then(serde_json::Value::as_str)
            == Some(matched_rsid)
    })
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
    let assay_id = app_assay_id_from_workspace(&workspace, manifest_path)?;
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
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &store, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
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
        &[],
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
            artifact(
                "observations.tsv",
                "text/tab-separated-values",
                observations_tsv,
            ),
            artifact("analysis.jsonl", "application/jsonl", analysis_jsonl),
            artifact("reports.jsonl", "application/jsonl", reports_jsonl),
            artifact("index.html", "text/html", html),
        ],
        duration_ms: (js_sys::Date::now() - started_ms).max(0.0) as u128,
        text_output,
    })
    .map_err(|err| JsError::new(&format!("failed to encode report output: {err}")))
}

/// Mirrors `runPackageReportBytes` but for CRAM input. The CRAM body and
/// FASTA reference are streamed via JS-supplied `readAt` callbacks so the
/// browser doesn't have to load multi-GB genomes into wasm memory. The CRAI
/// and FAI indexes are passed inline.
///
/// Analyses run against the observations produced from the CRAM lookup. The
/// per-script Python interpreter still receives `input_bytes` as a virtual
/// file; for CRAM that's an empty buffer because typical PGx analysis scripts
/// (apoe, mthfr, apol1, …) read observation rows rather than raw genome bytes.
#[wasm_bindgen(js_name = runPackageReportFromCram)]
#[allow(clippy::too_many_arguments)]
pub fn run_package_report_from_cram(
    manifest_path: &str,
    package_files_json: &str,
    input_name: &str,
    cram_read_at: js_sys::Function,
    cram_len: f64,
    crai_bytes: &[u8],
    fasta_read_at: js_sys::Function,
    fasta_len: f64,
    fai_bytes: &[u8],
    options_json: Option<String>,
) -> Result<String, JsError> {
    use crate::js_reader::JsReader;
    use std::io::BufReader;
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
    let assay_id = app_assay_id_from_workspace(&workspace, manifest_path)?;
    let manifest_metadata = workspace.report_manifest_metadata(manifest_path)?;
    let findings = workspace.load_manifest_findings(manifest_path)?;
    let provenance = workspace.load_manifest_provenance_links(manifest_path)?;
    let mut head_inspection = inspect_head_via_js_reader(
        &cram_read_at,
        cram_len as u64,
        input_name,
        false, // sex detection runs separately below via the indexed reader
    );

    let crai_index = bioscript_formats::alignment::parse_crai_bytes(crai_bytes)
        .map_err(|err| JsError::new(&format!("parse crai: {err:?}")))?;
    let fai_index = bioscript_formats::alignment::parse_fai_bytes(fai_bytes)
        .map_err(|err| JsError::new(&format!("parse fai: {err:?}")))?;
    let fasta_reader = BufReader::new(JsReader::new(fasta_read_at, fasta_len as u64, "fasta"));
    let repository = bioscript_formats::alignment::build_reference_repository_from_readers(
        fasta_reader,
        fai_index,
    );
    let cram_reader = JsReader::new(cram_read_at, cram_len as u64, "cram");
    let indexed = bioscript_formats::alignment::build_cram_indexed_reader_from_reader(
        cram_reader,
        crai_index,
        repository,
    )
    .map_err(|err| JsError::new(&format!("build cram reader: {err:?}")))?;

    let lookup = CramReportLookup {
        reader: std::cell::RefCell::new(indexed),
        label: input_name.to_owned(),
    };

    // CRAM sex detection: explicit override wins, otherwise alignment Y/X
    // coverage analysis through the same reader the variant lookup will use.
    if let Some(explicit) = explicit_sex_from_options(&options) {
        head_inspection.inferred_sex = Some(explicit);
    } else if options.detect_sex {
        let mut reader_borrow = lookup.reader.borrow_mut();
        match bioscript_formats::infer_sex_from_alignment_reader(
            &mut reader_borrow,
            &lookup.label,
            true,
        ) {
            Ok(inference) => head_inspection.inferred_sex = Some(inference),
            Err(err) => {
                head_inspection
                    .evidence
                    .push(format!("alignment sex detection failed: {err:?}"));
            }
        }
    }

    let mut loader = GenotypeLoadOptions::default();
    loader.format = Some(bioscript_formats::GenotypeSourceFormat::Cram);
    loader.allow_reference_md5_mismatch = true;
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &lookup, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
        .iter()
        .map(|row| workspace.app_observation_from_manifest_row(row, &assay_id, None, None))
        .collect::<Result<Vec<_>, _>>()?;
    let analysis_observations =
        analysis_cache_observations(&manifest_output.observations, &observations);
    // Analysis scripts call `bioscript.load_genotypes(input_file)` then rsid
    // lookups via `genotypes.lookup_variants(plan)`. The runtime now layers a
    // pre-resolved-observation cache over whatever the input file resolves
    // to (Plan B in genotype/types.rs:QueryBackend::Cached), so for CRAM the
    // cache hits and we skip re-walking the genome. The input bytes can be
    // empty since every spec the panel/assays declared is in the cache.
    let analyses = workspace.run_manifest_analyses(
        manifest_path,
        input_name,
        &[],
        &analysis_observations,
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
        input_inspection: Some(&head_inspection),
        manifest_metadata: &manifest_metadata,
    })];
    let observations_tsv = render_app_observations_tsv(&observations)?;
    let analysis_jsonl = render_jsonl(&analyses)?;
    let reports_jsonl = render_jsonl(&reports)?;
    let html = render_app_html_document(&observations, &reports)?;
    let text_output = "observations: observations.tsv\nanalysis: analysis.jsonl\nreports: reports.jsonl\nhtml: index.html\n".to_owned();
    serde_json::to_string(&ReportRunOutput {
        artifacts: vec![
            artifact(
                "observations.tsv",
                "text/tab-separated-values",
                observations_tsv,
            ),
            artifact("analysis.jsonl", "application/jsonl", analysis_jsonl),
            artifact("reports.jsonl", "application/jsonl", reports_jsonl),
            artifact("index.html", "text/html", html),
        ],
        duration_ms: (js_sys::Date::now() - started_ms).max(0.0) as u128,
        text_output,
    })
    .map_err(|err| JsError::new(&format!("failed to encode report output: {err}")))
}

/// Mirrors `runPackageReportBytes` but for BAM input. The BAM body is streamed
/// via a JS-supplied `readAt` callback and the BAI index is passed inline.
#[wasm_bindgen(js_name = runPackageReportFromBam)]
#[allow(clippy::too_many_arguments)]
pub fn run_package_report_from_bam(
    manifest_path: &str,
    package_files_json: &str,
    input_name: &str,
    bam_read_at: js_sys::Function,
    bam_len: f64,
    bai_bytes: &[u8],
    options_json: Option<String>,
) -> Result<String, JsError> {
    use crate::js_reader::JsReader;
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
    let assay_id = app_assay_id_from_workspace(&workspace, manifest_path)?;
    let manifest_metadata = workspace.report_manifest_metadata(manifest_path)?;
    let findings = workspace.load_manifest_findings(manifest_path)?;
    let provenance = workspace.load_manifest_provenance_links(manifest_path)?;
    let mut head_inspection =
        inspect_head_via_js_reader(&bam_read_at, bam_len as u64, input_name, false);

    let bai_index = bioscript_formats::alignment::parse_bai_bytes(bai_bytes)
        .map_err(|err| JsError::new(&format!("parse bai: {err:?}")))?;
    let bam_reader = JsReader::new(bam_read_at, bam_len as u64, "bam");
    let indexed = bioscript_formats::alignment::build_bam_indexed_reader_from_reader(
        bam_reader,
        bai_index,
    )
    .map_err(|err| JsError::new(&format!("build bam reader: {err:?}")))?;

    let lookup = BamReportLookup {
        reader: std::cell::RefCell::new(indexed),
        label: input_name.to_owned(),
    };

    if let Some(explicit) = explicit_sex_from_options(&options) {
        head_inspection.inferred_sex = Some(explicit);
    }

    let mut loader = GenotypeLoadOptions::default();
    loader.format = Some(bioscript_formats::GenotypeSourceFormat::Bam);
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &lookup, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
        .iter()
        .map(|row| workspace.app_observation_from_manifest_row(row, &assay_id, None, None))
        .collect::<Result<Vec<_>, _>>()?;
    let analysis_observations =
        analysis_cache_observations(&manifest_output.observations, &observations);
    let analyses = workspace.run_manifest_analyses(
        manifest_path,
        input_name,
        &[],
        &analysis_observations,
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
        input_inspection: Some(&head_inspection),
        manifest_metadata: &manifest_metadata,
    })];
    let observations_tsv = render_app_observations_tsv(&observations)?;
    let analysis_jsonl = render_jsonl(&analyses)?;
    let reports_jsonl = render_jsonl(&reports)?;
    let html = render_app_html_document(&observations, &reports)?;
    let text_output = "observations: observations.tsv\nanalysis: analysis.jsonl\nreports: reports.jsonl\nhtml: index.html\n".to_owned();
    serde_json::to_string(&ReportRunOutput {
        artifacts: vec![
            artifact(
                "observations.tsv",
                "text/tab-separated-values",
                observations_tsv,
            ),
            artifact("analysis.jsonl", "application/jsonl", analysis_jsonl),
            artifact("reports.jsonl", "application/jsonl", reports_jsonl),
            artifact("index.html", "text/html", html),
        ],
        duration_ms: (js_sys::Date::now() - started_ms).max(0.0) as u128,
        text_output,
    })
    .map_err(|err| JsError::new(&format!("failed to encode report output: {err}")))
}

/// Mirrors `runPackageReportBytes` but for a bgzipped, tabix-indexed VCF
/// streamed via JS-supplied `readAt` callbacks. The TBI is passed inline.
#[wasm_bindgen(js_name = runPackageReportFromVcf)]
#[allow(clippy::too_many_arguments)]
pub fn run_package_report_from_vcf(
    manifest_path: &str,
    package_files_json: &str,
    input_name: &str,
    vcf_read_at: js_sys::Function,
    vcf_len: f64,
    tbi_bytes: &[u8],
    options_json: Option<String>,
) -> Result<String, JsError> {
    use crate::js_reader::JsReader;
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
    let assay_id = app_assay_id_from_workspace(&workspace, manifest_path)?;
    let manifest_metadata = workspace.report_manifest_metadata(manifest_path)?;
    let findings = workspace.load_manifest_findings(manifest_path)?;
    let provenance = workspace.load_manifest_provenance_links(manifest_path)?;
    // Inspect format/source/assembly from the head, but skip the byte-stream
    // sex detection — we'll do that via tabix-targeted X non-PAR queries
    // below, which works on indexed VCFs of any size.
    let mut head_inspection =
        inspect_head_via_js_reader(&vcf_read_at, vcf_len as u64, input_name, false);
    // Decompress the head once to grab the VCF header lines (## meta + #CHROM
    // column header) — these are needed by `infer_sex_from_text_lines` to
    // figure out delimiter / column indexes for the data lines we'll pull
    // via tabix below.
    let head_lines = decompress_vcf_head_lines(&vcf_read_at, vcf_len as u64);

    let tabix_index = bioscript_formats::alignment::parse_tbi_bytes(tbi_bytes)
        .map_err(|err| JsError::new(&format!("parse tbi: {err:?}")))?;
    let vcf_reader = JsReader::new(vcf_read_at, vcf_len as u64, "vcf");
    let indexed = noodles::csi::io::IndexedReader::new(vcf_reader, tabix_index);

    let lookup = VcfReportLookup {
        reader: std::cell::RefCell::new(indexed),
        label: input_name.to_owned(),
        detected_assembly: head_inspection.assembly,
    };

    if let Some(explicit) = explicit_sex_from_options(&options) {
        head_inspection.inferred_sex = Some(explicit);
    } else if options.detect_sex {
        let mut reader_borrow = lookup.reader.borrow_mut();
        if let Some(inference) = vcf_sex_via_tabix(&mut reader_borrow, &head_lines) {
            head_inspection.inferred_sex = Some(inference);
        }
    }

    let mut loader = GenotypeLoadOptions::default();
    loader.format = Some(bioscript_formats::GenotypeSourceFormat::Vcf);
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &lookup, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
        .iter()
        .map(|row| workspace.app_observation_from_manifest_row(row, &assay_id, None, None))
        .collect::<Result<Vec<_>, _>>()?;
    let analysis_observations =
        analysis_cache_observations(&manifest_output.observations, &observations);
    // Pre-resolved observation cache replaces the synth approach: analysis
    // scripts hit the cache via QueryBackend::Cached and skip re-opening the
    // VCF. See report_api.rs:run_package_report_from_cram for the same
    // pattern and bioscript-formats::genotype::types::QueryBackend::Cached
    // for the dispatch.
    let analyses = workspace.run_manifest_analyses(
        manifest_path,
        input_name,
        &[],
        &analysis_observations,
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
        input_inspection: Some(&head_inspection),
        manifest_metadata: &manifest_metadata,
    })];
    let observations_tsv = render_app_observations_tsv(&observations)?;
    let analysis_jsonl = render_jsonl(&analyses)?;
    let reports_jsonl = render_jsonl(&reports)?;
    let html = render_app_html_document(&observations, &reports)?;
    let text_output = "observations: observations.tsv\nanalysis: analysis.jsonl\nreports: reports.jsonl\nhtml: index.html\n".to_owned();
    serde_json::to_string(&ReportRunOutput {
        artifacts: vec![
            artifact(
                "observations.tsv",
                "text/tab-separated-values",
                observations_tsv,
            ),
            artifact("analysis.jsonl", "application/jsonl", analysis_jsonl),
            artifact("reports.jsonl", "application/jsonl", reports_jsonl),
            artifact("index.html", "text/html", html),
        ],
        duration_ms: (js_sys::Date::now() - started_ms).max(0.0) as u128,
        text_output,
    })
    .map_err(|err| JsError::new(&format!("failed to encode report output: {err}")))
}
