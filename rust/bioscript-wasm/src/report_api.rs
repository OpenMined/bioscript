use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    time::Duration,
};

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};
use bioscript_formats::{
    GenotypeLoadOptions, GenotypeStore, InspectOptions, SexInference,
    inspect_bytes as inspect_bytes_rs,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use bioscript_schema::PanelInterpretation;
use monty::{MontyObject, ResourceLimits};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[path = "report_api/analysis_cache.rs"]
mod analysis_cache;
#[path = "report_helpers.rs"]
mod report_helpers;
#[path = "report_input_inspection.rs"]
mod report_input_inspection;
#[path = "report_lookup.rs"]
mod report_lookup;
#[path = "report_workspace.rs"]
mod report_workspace;

use analysis_cache::analysis_cache_observations;
use report_helpers::*;
use report_input_inspection::{
    decompress_vcf_head_lines, explicit_sex_from_options, inspect_head_via_js_reader,
    vcf_sex_via_tabix,
};
use report_lookup::{BamReportLookup, CramReportLookup, VcfReportLookup};
use report_workspace::PackageWorkspace;

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
    let manifest_context = workspace.report_manifest_context(manifest_path)?;
    let inspect_options = InspectOptions {
        input_index: None,
        reference_file: None,
        reference_index: None,
        detect_sex: options.detect_sex,
    };
    let input_inspection = inspect_bytes_rs(input_name, input_bytes, &inspect_options)
        .map_err(|err| JsError::new(&format!("inspect input failed: {err:?}")))?;
    let loader = GenotypeLoadOptions {
        assembly: input_inspection.assembly,
        inferred_sex: input_inspection
            .inferred_sex
            .as_ref()
            .map(|inference| inference.sex),
        ..Default::default()
    };
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
                &manifest_context.assay_id,
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
    let artifacts = bioscript_reporting::render_input_report_artifact_texts(
        bioscript_reporting::AppInputReportInput {
            assay_id: &manifest_context.assay_id,
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path: input_name,
            observations: &observations,
            analyses: &analyses,
            findings: &manifest_context.findings,
            provenance: &manifest_context.provenance,
            input_inspection: Some(&input_inspection),
            manifest_metadata: &manifest_context.manifest_metadata,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, artifacts)
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
    let manifest_context = workspace.report_manifest_context(manifest_path)?;
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

    let loader = GenotypeLoadOptions {
        format: Some(bioscript_formats::GenotypeSourceFormat::Cram),
        allow_reference_md5_mismatch: true,
        ..Default::default()
    };
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &lookup, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
        .iter()
        .map(|row| {
            workspace.app_observation_from_manifest_row(row, &manifest_context.assay_id, None, None)
        })
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
    let artifacts = bioscript_reporting::render_input_report_artifact_texts(
        bioscript_reporting::AppInputReportInput {
            assay_id: &manifest_context.assay_id,
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path: input_name,
            observations: &observations,
            analyses: &analyses,
            findings: &manifest_context.findings,
            provenance: &manifest_context.provenance,
            input_inspection: Some(&head_inspection),
            manifest_metadata: &manifest_context.manifest_metadata,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, artifacts)
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
    let manifest_context = workspace.report_manifest_context(manifest_path)?;
    let mut head_inspection =
        inspect_head_via_js_reader(&bam_read_at, bam_len as u64, input_name, false);

    let bai_index = bioscript_formats::alignment::parse_bai_bytes(bai_bytes)
        .map_err(|err| JsError::new(&format!("parse bai: {err:?}")))?;
    let bam_reader = JsReader::new(bam_read_at, bam_len as u64, "bam");
    let indexed =
        bioscript_formats::alignment::build_bam_indexed_reader_from_reader(bam_reader, bai_index)
            .map_err(|err| JsError::new(&format!("build bam reader: {err:?}")))?;

    let lookup = BamReportLookup {
        reader: std::cell::RefCell::new(indexed),
        label: input_name.to_owned(),
    };

    if let Some(explicit) = explicit_sex_from_options(&options) {
        head_inspection.inferred_sex = Some(explicit);
    }

    let loader = GenotypeLoadOptions {
        format: Some(bioscript_formats::GenotypeSourceFormat::Bam),
        ..Default::default()
    };
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &lookup, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
        .iter()
        .map(|row| {
            workspace.app_observation_from_manifest_row(row, &manifest_context.assay_id, None, None)
        })
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
    let artifacts = bioscript_reporting::render_input_report_artifact_texts(
        bioscript_reporting::AppInputReportInput {
            assay_id: &manifest_context.assay_id,
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path: input_name,
            observations: &observations,
            analyses: &analyses,
            findings: &manifest_context.findings,
            provenance: &manifest_context.provenance,
            input_inspection: Some(&head_inspection),
            manifest_metadata: &manifest_context.manifest_metadata,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, artifacts)
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
    let manifest_context = workspace.report_manifest_context(manifest_path)?;
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

    let loader = GenotypeLoadOptions {
        format: Some(bioscript_formats::GenotypeSourceFormat::Vcf),
        ..Default::default()
    };
    let manifest_output =
        workspace.run_manifest_rows(manifest_path, &lookup, &participant_id, &options.filters)?;
    let observations = manifest_output
        .rows
        .iter()
        .map(|row| {
            workspace.app_observation_from_manifest_row(row, &manifest_context.assay_id, None, None)
        })
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
    let artifacts = bioscript_reporting::render_input_report_artifact_texts(
        bioscript_reporting::AppInputReportInput {
            assay_id: &manifest_context.assay_id,
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path: input_name,
            observations: &observations,
            analyses: &analyses,
            findings: &manifest_context.findings,
            provenance: &manifest_context.provenance,
            input_inspection: Some(&head_inspection),
            manifest_metadata: &manifest_context.manifest_metadata,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, artifacts)
}
