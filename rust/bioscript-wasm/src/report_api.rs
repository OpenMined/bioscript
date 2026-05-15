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
use bioscript_schema::{PanelInterpretation, VariantManifest, load_variant_manifest_text};
use monty::{MontyObject, ResourceLimits};
use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

#[path = "report_helpers.rs"]
mod report_helpers;
#[path = "report_input_inspection.rs"]
mod report_input_inspection;
#[path = "report_lookup.rs"]
mod report_lookup;
#[path = "report_workspace.rs"]
mod report_workspace;

use report_helpers::*;
use report_input_inspection::{
    detect_vcf_assembly_via_js_reader, explicit_sex_from_options, inspect_head_via_js_reader,
    vcf_sex_via_js_reader,
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
    #[serde(default)]
    input_file_path: Option<String>,
    #[serde(default)]
    input_index_path: Option<String>,
    #[serde(default)]
    reference_file_path: Option<String>,
    #[serde(default)]
    reference_index_path: Option<String>,
    /// Optional explicit sample sex (mirrors the CLI's `--sample-sex` flag).
    /// When set, takes precedence over inference: the report carries
    /// `method=explicit_sample_sex` like the CLI.
    #[serde(default)]
    sample_sex: Option<String>,
}

impl ReportOptionsInput {
    fn inspect_options(&self, detect_sex: bool) -> InspectOptions {
        InspectOptions {
            input_index: self.input_index_path.as_ref().map(PathBuf::from),
            reference_file: self.reference_file_path.as_ref().map(PathBuf::from),
            reference_index: self.reference_index_path.as_ref().map(PathBuf::from),
            detect_sex,
        }
    }
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
    let input_file_path = options.input_file_path.as_deref().unwrap_or(input_name);
    let inspect_options = options.inspect_options(options.detect_sex);
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
    let analysis_runner = report_workspace::WasmReportAnalysisRunner {
        workspace: &workspace,
        input_name,
        input_bytes,
        participant_id: &participant_id,
        loader: &loader,
        options: &options,
    };
    let run = bioscript_reporting::run_report(
        &workspace,
        manifest_path,
        &store,
        &analysis_runner,
        bioscript_reporting::ReportInputContext {
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path,
            input_inspection: Some(&input_inspection),
        },
        bioscript_reporting::ReportRunOptions {
            filters: &options.filters,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, run.artifacts)
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
    let input_file_path = options.input_file_path.as_deref().unwrap_or(input_name);
    let mut head_inspection = inspect_head_via_js_reader(
        &cram_read_at,
        cram_len as u64,
        input_name,
        &options.inspect_options(false),
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
        label: input_file_path.to_owned(),
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
    let analysis_runner = report_workspace::WasmReportAnalysisRunner {
        workspace: &workspace,
        input_name,
        input_bytes: &[],
        participant_id: &participant_id,
        loader: &loader,
        options: &options,
    };
    let run = bioscript_reporting::run_report(
        &workspace,
        manifest_path,
        &lookup,
        &analysis_runner,
        bioscript_reporting::ReportInputContext {
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path,
            input_inspection: Some(&head_inspection),
        },
        bioscript_reporting::ReportRunOptions {
            filters: &options.filters,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, run.artifacts)
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
    let input_file_path = options.input_file_path.as_deref().unwrap_or(input_name);
    let mut head_inspection = inspect_head_via_js_reader(
        &bam_read_at,
        bam_len as u64,
        input_name,
        &options.inspect_options(false),
        false,
    );

    let bai_index = bioscript_formats::alignment::parse_bai_bytes(bai_bytes)
        .map_err(|err| JsError::new(&format!("parse bai: {err:?}")))?;
    let bam_reader = JsReader::new(bam_read_at, bam_len as u64, "bam");
    let indexed =
        bioscript_formats::alignment::build_bam_indexed_reader_from_reader(bam_reader, bai_index)
            .map_err(|err| JsError::new(&format!("build bam reader: {err:?}")))?;

    let lookup = BamReportLookup {
        reader: std::cell::RefCell::new(indexed),
        label: input_file_path.to_owned(),
    };

    if let Some(explicit) = explicit_sex_from_options(&options) {
        head_inspection.inferred_sex = Some(explicit);
    }

    let loader = GenotypeLoadOptions {
        format: Some(bioscript_formats::GenotypeSourceFormat::Bam),
        ..Default::default()
    };
    let analysis_runner = report_workspace::WasmReportAnalysisRunner {
        workspace: &workspace,
        input_name,
        input_bytes: &[],
        participant_id: &participant_id,
        loader: &loader,
        options: &options,
    };
    let run = bioscript_reporting::run_report(
        &workspace,
        manifest_path,
        &lookup,
        &analysis_runner,
        bioscript_reporting::ReportInputContext {
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path,
            input_inspection: Some(&head_inspection),
        },
        bioscript_reporting::ReportRunOptions {
            filters: &options.filters,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, run.artifacts)
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
    let input_file_path = options.input_file_path.as_deref().unwrap_or(input_name);
    // Inspect format/source/assembly from the head, but skip the byte-stream
    // sex detection — we'll do that via tabix-targeted X non-PAR queries
    // below, which works on indexed VCFs of any size.
    let mut head_inspection = inspect_head_via_js_reader(
        &vcf_read_at,
        vcf_len as u64,
        input_name,
        &options.inspect_options(false),
        false,
    );
    if head_inspection.assembly.is_none() {
        head_inspection.assembly =
            detect_vcf_assembly_via_js_reader(vcf_read_at.clone(), vcf_len as u64, input_name);
    }
    let tabix_index = bioscript_formats::alignment::parse_tbi_bytes(tbi_bytes)
        .map_err(|err| JsError::new(&format!("parse tbi: {err:?}")))?;
    let vcf_reader = JsReader::new(vcf_read_at.clone(), vcf_len as u64, "vcf");
    let indexed = noodles::csi::io::IndexedReader::new(vcf_reader, tabix_index);

    let lookup = VcfReportLookup {
        reader: std::cell::RefCell::new(indexed),
        label: input_file_path.to_owned(),
        detected_assembly: head_inspection.assembly,
    };

    if let Some(explicit) = explicit_sex_from_options(&options) {
        head_inspection.inferred_sex = Some(explicit);
    } else if options.detect_sex
        && let Some(inference) = vcf_sex_via_js_reader(vcf_read_at, vcf_len as u64, input_name)
    {
        head_inspection.inferred_sex = Some(inference);
    }

    let loader = GenotypeLoadOptions {
        format: Some(bioscript_formats::GenotypeSourceFormat::Vcf),
        ..Default::default()
    };
    let analysis_runner = report_workspace::WasmReportAnalysisRunner {
        workspace: &workspace,
        input_name,
        input_bytes: &[],
        participant_id: &participant_id,
        loader: &loader,
        options: &options,
    };
    let run = bioscript_reporting::run_report(
        &workspace,
        manifest_path,
        &lookup,
        &analysis_runner,
        bioscript_reporting::ReportInputContext {
            participant_id: &participant_id,
            input_file_name: input_name,
            input_file_path,
            input_inspection: Some(&head_inspection),
        },
        bioscript_reporting::ReportRunOptions {
            filters: &options.filters,
        },
    )
    .map_err(|err| JsError::new(&err))?;
    encode_report_run_output(started_ms, run.artifacts)
}
