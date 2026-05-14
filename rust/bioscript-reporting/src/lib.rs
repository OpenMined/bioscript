#![allow(clippy::missing_errors_doc, clippy::must_use_candidate)]

mod analysis;
mod artifacts;
mod html;
mod manifest;
mod matching;
mod observation;
mod report_json;
mod rows;

pub use analysis::{
    AnalysisOutputFormat, AnalysisOutputJsonInput, analysis_observations_relative_file,
    analysis_output_format, analysis_output_json, analysis_output_relative_file,
    parse_analysis_output_text, participant_id_from_name, participant_id_from_path,
    render_analysis_observations_tsv, validate_bioscript_interpretation,
};
pub use artifacts::{
    ReportArtifactTexts, json_field_as_tsv, render_input_report_artifact_texts, render_jsonl,
    render_observations_tsv, render_report_artifact_texts, standard_text_output,
};
pub use html::render_app_html_document;
pub use manifest::{
    AnalysisManifestTask, ExecutableAssayMember, ExecutablePanelMember,
    FilesystemManifestWorkspace, ManifestWorkspace, ReportManifestContext, ReportManifestKind,
    VariantManifestTask, assay_executable_member, assay_executable_member_path,
    collect_analysis_manifest_tasks, collect_manifest_provenance_entries,
    collect_variant_manifest_tasks, load_manifest_findings, load_manifest_provenance_links,
    load_report_manifest_context, matches_analysis_path_filters, matches_variant_manifest_filters,
    panel_executable_member, panel_executable_member_path, report_assay_id, report_manifest_kind,
    report_manifest_metadata, report_manifest_schema, resolve_filesystem_manifest_path,
};
pub use matching::match_app_findings;
pub use observation::{AppObservationInput, app_observation_from_manifest_row};
pub use report_json::{
    AppInputReportInput, AppReportJsonInput, app_input_report_json, app_report_json,
};
pub use rows::{
    MANIFEST_ROW_TSV_HEADERS, render_manifest_rows_tsv, render_manifest_trace_tsv, variant_row,
};
