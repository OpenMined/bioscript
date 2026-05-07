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

mod report_render;
mod report_workspace;

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

#[derive(Clone, Copy)]
struct AppReportJsonInput<'a> {
    assay_id: &'a str,
    participant_id: &'a str,
    input_file_name: &'a str,
    observations: &'a [serde_json::Value],
    analyses: &'a [serde_json::Value],
    findings: &'a [serde_json::Value],
    provenance: &'a [serde_json::Value],
    input_inspection: Option<&'a bioscript_formats::FileInspection>,
    manifest_metadata: &'a serde_json::Value,
}

fn app_report_json(input: AppReportJsonInput<'_>) -> serde_json::Value {
    let called = input
        .observations
        .iter()
        .filter(|item| item.get("call_status").and_then(serde_json::Value::as_str) == Some("called"))
        .count();
    serde_json::json!({
        "schema": "bioscript:report:1.0",
        "version": "1.0",
        "participant_id": input.participant_id,
        "assay_id": input.assay_id,
        "assay_version": "1.0",
        "manifest": input.manifest_metadata,
        "input": {
            "file_name": input.input_file_name,
            "file_path": input.input_file_name,
            "debug": input.input_inspection.map(input_inspection_json),
        },
        "report_status": if called == input.observations.len() { "complete" } else { "partial" },
        "derived_from": input.observations.iter().filter_map(|item| item.get("variant_key").cloned()).collect::<Vec<_>>(),
        "analyses": input.analyses,
        "findings": input.findings,
        "provenance": input.provenance,
        "metrics": {
            "n_sites_tested": input.observations.len(),
            "n_sites_called": called,
            "n_sites_missing": input.observations.len().saturating_sub(called),
            "n_analyses": input.analyses.len(),
            "n_findings_matched": input.findings.len(),
        }
    })
}

fn match_app_findings(
    findings: &[serde_json::Value],
    observations: &[serde_json::Value],
    analyses: &[serde_json::Value],
) -> Vec<serde_json::Value> {
    let mut matched = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for finding in findings {
        if let Some(effects) = finding.get("effects").and_then(serde_json::Value::as_array) {
            for effect in effects {
                if let Some(observation) = app_finding_match_observation(effect, observations) {
                    let mut item = finding.clone();
                    if let Some(object) = item.as_object_mut() {
                        object.remove("effects");
                        object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                        object.insert("matched_effect".to_owned(), effect.clone());
                        object.insert(
                            "matched_observation".to_owned(),
                            app_finding_observation_context(observation),
                        );
                    }
                    if seen.insert(app_finding_dedupe_key(&item)) {
                        matched.push(item);
                    }
                } else if let Some(analysis) = app_finding_match_analysis(effect, analyses) {
                    let mut item = finding.clone();
                    if let Some(object) = item.as_object_mut() {
                        object.remove("effects");
                        object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                        object.insert("matched_effect".to_owned(), effect.clone());
                        object.insert("matched_analysis".to_owned(), analysis);
                    }
                    if seen.insert(app_finding_dedupe_key(&item)) {
                        matched.push(item);
                    }
                }
            }
        } else if let Some(observation) = app_finding_match_observation(finding, observations) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert(
                    "matched_observation".to_owned(),
                    app_finding_observation_context(observation),
                );
            }
            if seen.insert(app_finding_dedupe_key(&item)) {
                matched.push(item);
            }
        } else if let Some(analysis) = app_finding_match_analysis(finding, analyses) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert("matched_analysis".to_owned(), analysis);
            }
            if seen.insert(app_finding_dedupe_key(&item)) {
                matched.push(item);
            }
        }
    }
    matched
}

fn render_app_html_document(
    observations: &[serde_json::Value],
    reports: &[serde_json::Value],
) -> Result<String, JsError> {
    let mut out = String::from(
        r##"<!doctype html><meta charset="utf-8"><title>BioScript report</title><style>body{font-family:system-ui,sans-serif;margin:0;background:#f7f8fa;color:#1f2933}.wrap{max-width:1440px;margin:0 auto;padding:24px}h1{margin:0 0 10px}h2{margin:32px 0 10px;scroll-margin-top:82px}.nav{position:sticky;top:0;z-index:20;display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:16px -24px 22px;padding:10px 24px;background:rgba(247,248,250,.96);border-block:1px solid #d8dee6;backdrop-filter:saturate(160%) blur(8px)}.nav a{border:1px solid #cbd5df;background:#fff;color:#1f2933;text-decoration:none;padding:7px 10px;border-radius:6px}.nav a:hover{background:#eef2f6}.table-tools,.level-filter{display:flex;justify-content:space-between;gap:12px;align-items:flex-start;margin:6px 0}.table-tools input{width:min(420px,100%);border:1px solid #cbd5df;border-radius:6px;padding:7px 9px;font:inherit;background:#fff}.table-tools label,.level-filter label{display:flex;gap:5px;align-items:center}.level-filter{justify-content:flex-start;flex-wrap:wrap}.filter-scale{display:flex;flex-direction:column;gap:6px}.filter-scale-title{font-weight:700;color:#344054;display:flex;gap:6px;align-items:center}.filter-options{display:flex;flex-wrap:wrap;gap:6px 10px}.filter-actions{display:flex;gap:8px;align-self:flex-end}.level-filter a{display:inline-grid;place-items:center;width:20px;height:20px;border:1px solid #cbd5df;border-radius:50%;text-decoration:none;color:#1f2933;background:#fff;font-weight:700}.filter-action,.pgx-tabs button{border:1px solid #cbd5df;background:#fff;color:#1f2933;border-radius:6px;padding:4px 7px;font:inherit;cursor:pointer}.filter-action:hover,.pgx-tabs button:hover{background:#eef2f6}.pgx-tabs{display:flex;gap:8px;margin:10px 0}.pgx-tabs button.active{background:#1f2933;border-color:#1f2933;color:#fff}.table-wrap{overflow:auto;border:1px solid #d8dee6;background:white;border-radius:8px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #e5e9ef;padding:6px 8px;text-align:left;vertical-align:top}th{position:sticky;top:0;background:#eef2f6;z-index:1;white-space:nowrap;cursor:pointer;user-select:none}.sort-mark{font-size:10px;color:#667085;margin-left:3px}.row-variant td{background:#fff7cc}.row-reference td{background:#eaf7ee}.debug-hidden .debug-col{display:none}.participant-filter{display:flex;gap:8px;align-items:center;margin:12px 0}.participant-filter select{border:1px solid #cbd5df;border-radius:6px;padding:6px 8px;background:#fff;font:inherit}.genotype-hit{font-weight:700}.allele-hit{color:#075985;background:#dff4ff;border-radius:3px;padding:0 2px}.level-badge,.pgx-badge{display:inline-block;min-width:2.2em;text-align:center;border-radius:999px;padding:2px 8px;color:#fff;font-weight:700}.level-1,.pgx-informative{background:#0abc72}.level-2,.pgx-actionable{background:#2a74df}.level-3,.pgx-recommended{background:#ffc107;color:#1f2933}.level-4,.pgx-required{background:#c53b3b}.pgx-no-clinical,.pgx-criteria,.pgx-unknown,.level-unknown{background:#667085}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.muted{color:#667085}.effect{max-width:760px;min-width:360px}.analysis-kv{display:grid;grid-template-columns:max-content minmax(0,1fr);gap:6px 14px;background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 12px;margin:8px 0 10px}.analysis-kv dt{font-weight:700;color:#344054}.analysis-kv dd{margin:0;min-width:0}.analysis-badge{display:inline-block;border:1px solid #cbd5df;border-radius:999px;background:#eef2f6;padding:1px 8px;font-weight:700}.analysis-badge-normal{background:#eaf7ee;border-color:#9fd6ad;color:#14532d}.analysis-badge-variant{background:#fff7cc;border-color:#f0d66a;color:#713f12}.analysis-badge-unknown{background:#eef2f6;border-color:#cbd5df;color:#475467}.analysis-card{background:#fff;border:1px solid #cbd5df;border-radius:8px;margin:14px 0;overflow:hidden}.analysis-card summary{display:flex;align-items:center;gap:10px;padding:11px 13px;background:#eef2f6;cursor:pointer;font-weight:700}.analysis-card-body{padding:12px 14px}.analysis-card-title{font-size:16px}.logic-note,.analysis-notes{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 12px;margin:8px 0 10px}.logic-note h4,.analysis-notes h4,h4{margin:0 0 6px}.logic-note p,.analysis-notes p{margin:0 0 6px}.analysis-notes p:last-child{margin-bottom:0}.provenance-list{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 18px}.provenance-list li{margin:8px 0}pre{white-space:pre-wrap;background:#fff;padding:12px;border:1px solid #d8dee6;border-radius:8px;max-height:520px;overflow:auto}</style><script>const sortState={};function cellText(row,i){const cell=row.cells[i];return(cell?.dataset.sort||cell?.innerText||"").trim()}function cmp(a,b){const an=Number(a),bn=Number(b);if(a!==""&&b!==""&&!Number.isNaN(an)&&!Number.isNaN(bn))return an-bn;return a.localeCompare(b,undefined,{numeric:true,sensitivity:"base"})}function sortTable(id,col){const table=document.getElementById(id);const tbody=table.tBodies[0];const key=id+":"+col;const dir=sortState[key]==="asc"?"desc":"asc";sortState[key]=dir;table.querySelectorAll(".sort-mark").forEach(s=>s.textContent="");table.tHead.rows[0].cells[col].querySelector(".sort-mark").textContent=dir==="asc"?"^":"v";Array.from(tbody.rows).sort((a,b)=>{const v=cmp(cellText(a,col),cellText(b,col));return dir==="asc"?v:-v}).forEach(r=>tbody.appendChild(r))}let selectedParticipant="";function participantOk(row){return !selectedParticipant||!row.dataset.participant||row.dataset.participant===selectedParticipant}function pgxLevelOk(row){if(!row.dataset.pgxLevel&&!row.dataset.level)return true;const key=row.dataset.pgxLevel||row.dataset.level||"unknown";const input=document.querySelector('[data-pgx-any-level-filter="'+key+'"]')||document.querySelector('[data-pgx-any-level-filter="unknown"]');return !!input&&input.checked}function pgxOutcomeOk(row){const key=row.dataset.pgxOutcome;if(!key)return true;const input=document.querySelector('[data-pgx-outcome-filter="'+key+'"]');return !input||input.checked}function observationOk(row){const key=row.dataset.observation;if(!key)return true;const input=document.querySelector('[data-observation-filter="'+key+'"]');return !input||input.checked}function applyTableFilters(id){const q=(document.querySelector('[data-filter-for="'+id+'"]')?.value||"").toLowerCase();document.querySelectorAll("#"+id+" tbody tr").forEach(row=>{let ok=row.innerText.toLowerCase().includes(q)&&participantOk(row)&&pgxLevelOk(row)&&pgxOutcomeOk(row)&&observationOk(row);row.style.display=ok?"":"none"})}function applyPgxFilters(){document.querySelectorAll('#pgx table[id]').forEach(table=>applyTableFilters(table.id));document.querySelectorAll('#pgx .pgx-drug-group').forEach(group=>{const visible=Array.from(group.querySelectorAll('tbody tr')).some(row=>row.style.display!=='none');group.style.display=visible?'':'none'})}function setPgxFilterGroup(checked){document.querySelectorAll("[data-pgx-any-level-filter]").forEach(input=>input.checked=checked);applyPgxFilters()}function setObservationFilterGroup(checked){document.querySelectorAll("[data-observation-filter]").forEach(input=>input.checked=checked);applyTableFilters("observations-table")}function setPgxView(view){document.getElementById("pgx-view-variant").hidden=view!=="variant";document.getElementById("pgx-view-drug").hidden=view!=="drug";document.getElementById("pgx-tab-variant")?.classList.toggle("active",view==="variant");document.getElementById("pgx-tab-drug")?.classList.toggle("active",view==="drug")}function setParticipant(value){selectedParticipant=value;document.querySelectorAll("table[id]").forEach(table=>applyTableFilters(table.id))}function toggleDebug(show){document.getElementById("report-wrap").classList.toggle("debug-hidden",!show)}</script><div class="wrap debug-hidden" id="report-wrap">"##,
    );
    let label_findings = collect_report_findings(reports, "bioscript:pgx-label:1.0");
    let summary_findings = collect_report_findings(reports, "bioscript:pgx-summary:1.0");
    let analysis_outputs = collect_report_analyses(reports);
    let participants = collect_report_participants(reports);
    render_report_manifest_header(&mut out, reports);
    let _ = write!(out, "<div class=\"muted\">{} observation(s), {} analysis output(s), {} PGx label finding(s), {} PGx summary finding(s)</div>", observations.len(), analysis_outputs.len(), label_findings.len(), summary_findings.len());
    render_participant_filter(&mut out, &participants);
    out.push_str("<nav class=\"nav\"><a href=\"#input-info\">Input</a><a href=\"#observations\">Observations</a><a href=\"#analysis\">Analysis</a><a href=\"#pgx\">PGx</a><a href=\"#provenance\">Provenance</a><a href=\"#source\">Source</a><a href=\"#json\">Raw JSON</a></nav>");
    out.push_str("<section id=\"input-info\"><h2>Input</h2>");
    render_input_debug(&mut out, reports, participants.len() > 1);
    out.push_str("</section><section id=\"observations\"><h2>Observations</h2>");
    render_observation_table(&mut out, observations, participants.len() > 1);
    out.push_str("</section><section id=\"analysis\"><h2>Analysis</h2>");
    render_analysis_tables(&mut out, &analysis_outputs, participants.len() > 1);
    out.push_str("</section><section id=\"pgx\"><h2>PGx</h2>");
    render_pgx_table(&mut out, &label_findings, &summary_findings);
    out.push_str("</section><section id=\"provenance\"><h2>Provenance</h2>");
    render_provenance_links(&mut out, reports);
    out.push_str("</section><section id=\"source\"><h2>Source</h2>");
    render_report_source_section(&mut out, reports);
    out.push_str("</section><section id=\"json\"><h2>Raw Reports JSON</h2><details><summary>Show raw report JSON</summary>");
    for report in reports {
        let text = serde_json::to_string_pretty(report).map_err(|err| JsError::new(&err.to_string()))?;
        let _ = write!(out, "<pre>{}</pre>", html_escape(&text));
    }
    out.push_str("</details></section></div>");
    Ok(out)
}

fn artifact(name: &str, mime_type: &str, text: String) -> ReportArtifactOutput {
    ReportArtifactOutput {
        name: name.to_owned(),
        path: name.to_owned(),
        mime_type: mime_type.to_owned(),
        text,
    }
}

fn variant_row(
    path: &str,
    name: &str,
    tags: &[String],
    observation: &VariantObservation,
    participant_id: &str,
) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("kind".to_owned(), "variant".to_owned());
    row.insert("name".to_owned(), name.to_owned());
    row.insert("path".to_owned(), path.to_owned());
    row.insert("tags".to_owned(), tags.join(","));
    row.insert("backend".to_owned(), observation.backend.clone());
    row.insert("participant_id".to_owned(), participant_id.to_owned());
    row.insert("matched_rsid".to_owned(), observation.matched_rsid.clone().unwrap_or_default());
    row.insert("assembly".to_owned(), observation.assembly.map(assembly_row_value).unwrap_or_default());
    row.insert("genotype".to_owned(), observation.genotype.clone().unwrap_or_default());
    row.insert("ref_count".to_owned(), observation.ref_count.map_or_else(String::new, |value| value.to_string()));
    row.insert("alt_count".to_owned(), observation.alt_count.map_or_else(String::new, |value| value.to_string()));
    row.insert("depth".to_owned(), observation.depth.map_or_else(String::new, |value| value.to_string()));
    row.insert("raw_counts".to_owned(), serde_json::to_string(&observation.raw_counts).unwrap_or_default());
    row.insert("evidence".to_owned(), observation.evidence.join(" | "));
    row
}

fn render_app_observations_tsv(observations: &[serde_json::Value]) -> Result<String, JsError> {
    let mut out = OBSERVATION_TSV_HEADERS.join("\t");
    out.push('\n');
    for observation in observations {
        let line = OBSERVATION_TSV_HEADERS
            .iter()
            .map(|header| json_field_as_tsv(observation.get(*header)))
            .collect::<Vec<_>>()
            .join("\t");
        out.push_str(&line);
        out.push('\n');
    }
    Ok(out)
}

fn render_jsonl(rows: &[serde_json::Value]) -> Result<String, JsError> {
    let mut out = String::new();
    for row in rows {
        out.push_str(&serde_json::to_string(row).map_err(|err| JsError::new(&err.to_string()))?);
        out.push('\n');
    }
    Ok(out)
}

fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}

fn normalize_package_path(path: &str) -> Result<String, JsError> {
    let mut out = PathBuf::new();
    for component in Path::new(path).components() {
        match component {
            std::path::Component::Normal(value) => out.push(value),
            std::path::Component::CurDir => {}
            _ => return Err(JsError::new(&format!("unsafe package path: {path}"))),
        }
    }
    Ok(out.display().to_string().replace('\\', "/"))
}

fn default_analysis_max_duration_ms() -> u64 {
    30_000
}

fn participant_id_from_name(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or(path)
        .replace([' ', '\t', '\n'], "_")
}

fn app_assay_id(path: &Path) -> Result<String, JsError> {
    path.file_stem()
        .and_then(|value| value.to_str())
        .map(ToOwned::to_owned)
        .ok_or_else(|| JsError::new(&format!("failed to derive assay id from {}", path.display())))
}

fn matches_filters(manifest: &VariantManifest, path: &str, filters: &[String]) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("kind", value)) => value == "variant",
        Some(("name", value)) => manifest.name.contains(value),
        Some(("path", value)) => path.contains(value),
        Some(("tag", value)) => manifest.tags.iter().any(|tag| tag == value),
        Some(_) | None => false,
    })
}

fn parse_analysis_output_text(
    text: &str,
    format: &str,
) -> Result<(Vec<serde_json::Value>, Vec<String>), JsError> {
    match format {
        "tsv" => Ok(parse_analysis_tsv(text)),
        "json" => {
            let value: serde_json::Value = serde_json::from_str(text)
                .map_err(|err| JsError::new(&format!("failed to parse analysis JSON: {err}")))?;
            let rows = match value {
                serde_json::Value::Array(rows) => rows,
                serde_json::Value::Object(mut object) => object
                    .remove("rows")
                    .and_then(|rows| rows.as_array().cloned())
                    .unwrap_or_else(|| vec![serde_json::Value::Object(object)]),
                other => vec![other],
            };
            let row_headers = rows
                .iter()
                .find_map(|row| row.as_object())
                .map(|object| object.keys().cloned().collect())
                .unwrap_or_default();
            Ok((rows, row_headers))
        }
        "jsonl" => {
            let mut rows: Vec<serde_json::Value> = Vec::new();
            for line in text.lines().filter(|line| !line.trim().is_empty()) {
                rows.push(serde_json::from_str(line)
                    .map_err(|err| JsError::new(&format!("failed to parse analysis JSONL: {err}")))?);
            }
            let row_headers = rows
                .iter()
                .find_map(|row| row.as_object())
                .map(|object| object.keys().cloned().collect())
                .unwrap_or_default();
            Ok((rows, row_headers))
        }
        other => Err(JsError::new(&format!("unsupported analysis output_format '{other}'"))),
    }
}

fn parse_analysis_tsv(text: &str) -> (Vec<serde_json::Value>, Vec<String>) {
    let mut lines = text.lines();
    let headers = lines
        .next()
        .map(|line| line.split('\t').map(ToOwned::to_owned).collect::<Vec<_>>())
        .unwrap_or_default();
    let rows = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let object = headers
                .iter()
                .enumerate()
                .map(|(index, header)| {
                    (
                        header.clone(),
                        serde_json::Value::String(fields.get(index).copied().unwrap_or_default().to_owned()),
                    )
                })
                .collect();
            serde_json::Value::Object(object)
        })
        .collect();
    (rows, headers)
}

fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, JsError> {
    serde_json::to_value(value).map_err(|err| JsError::new(&format!("failed to convert YAML to JSON: {err}")))
}

fn collect_manifest_provenance_entries(
    value: &serde_yaml::Value,
    links: &mut BTreeMap<String, serde_json::Value>,
) -> Result<(), JsError> {
    if let Some(sources) = value
        .get("provenance")
        .and_then(|provenance| provenance.get("sources"))
        .and_then(serde_yaml::Value::as_sequence)
    {
        for source in sources {
            let json = yaml_to_json(source.clone())?;
            if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
                links.entry(url.to_owned()).or_insert(json);
            }
        }
    }
    if let Some(source) = value.get("source") {
        let json = yaml_to_json(source.clone())?;
        if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
            links.entry(url.to_owned()).or_insert(json);
        }
    }
    Ok(())
}

fn input_inspection_json(inspection: &bioscript_formats::FileInspection) -> serde_json::Value {
    serde_json::json!({
        "container": match inspection.container {
            bioscript_formats::FileContainer::Plain => "plain",
            bioscript_formats::FileContainer::Zip => "zip",
        },
        "format": match inspection.detected_kind {
            bioscript_formats::DetectedKind::GenotypeText => "genotype_text",
            bioscript_formats::DetectedKind::Vcf => "vcf",
            bioscript_formats::DetectedKind::AlignmentCram => "alignment_cram",
            bioscript_formats::DetectedKind::AlignmentBam => "alignment_bam",
            bioscript_formats::DetectedKind::ReferenceFasta => "reference_fasta",
            bioscript_formats::DetectedKind::Unknown => "unknown",
        },
        "format_confidence": match inspection.confidence {
            bioscript_formats::DetectionConfidence::Authoritative => "authoritative",
            bioscript_formats::DetectionConfidence::StrongHeuristic => "strong_heuristic",
            bioscript_formats::DetectionConfidence::WeakHeuristic => "weak_heuristic",
            bioscript_formats::DetectionConfidence::Unknown => "unknown",
        },
        "assembly": inspection.assembly.map(|assembly| match assembly {
            Assembly::Grch37 => "grch37",
            Assembly::Grch38 => "grch38",
        }),
        "selected_entry": inspection.selected_entry,
        "source": inspection.source.as_ref().map(|source| serde_json::json!({
            "vendor": source.vendor,
            "platform_version": source.platform_version,
            "evidence": source.evidence,
        })),
        "inferred_sex": inspection.inferred_sex.as_ref().map(|sex| serde_json::json!({
            "sex": inferred_sex_name(sex.sex),
            "confidence": sex_detection_confidence_name(sex.confidence),
            "method": sex.method,
            "evidence": sex.evidence,
        })),
        "evidence": inspection.evidence,
        "warnings": inspection.warnings,
        "duration_ms": inspection.duration_ms,
    })
}

fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value.get(key).and_then(serde_yaml::Value::as_str).map(ToOwned::to_owned)
}

fn yaml_string_sequence(value: &serde_yaml::Value, key: &str) -> Vec<serde_json::Value> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_sequence)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_yaml::Value::as_str)
                .map(serde_json::Value::from)
                .collect()
        })
        .unwrap_or_default()
}

fn yaml_mapping_string(mapping: &serde_yaml::Mapping, key: &str) -> Option<String> {
    mapping
        .get(serde_yaml::Value::String(key.to_owned()))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

fn variant_primary_source_from_yaml(value: &serde_yaml::Value) -> serde_json::Value {
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    let _ = collect_manifest_provenance_entries(value, &mut links);
    if let Some(rsid) = value
        .get("identifiers")
        .and_then(|identifiers| identifiers.get("rsids"))
        .and_then(serde_yaml::Value::as_sequence)
        .and_then(|items| items.iter().find_map(serde_yaml::Value::as_str))
    {
        return serde_json::json!({
            "kind": "database",
            "label": "dbSNP / NCBI SNP",
            "url": format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"),
            "fields": ["identifiers.rsids"],
        });
    }
    links.into_values().next().unwrap_or(serde_json::Value::Null)
}

fn normalize_app_genotype(
    display: &str,
    ref_allele: &str,
    alt_allele: &str,
    chrom: &str,
    inferred_sex: Option<&SexInference>,
) -> (String, String) {
    if display.is_empty() {
        return ("./.".to_owned(), "unknown".to_owned());
    }
    let alleles: Vec<char> = display.chars().filter(char::is_ascii_alphabetic).collect();
    if ref_allele.len() != 1 || alt_allele.len() != 1 {
        return (display.to_owned(), "unknown".to_owned());
    }
    let ref_ch = ref_allele.chars().next().unwrap_or_default();
    let alt_ch = alt_allele.chars().next().unwrap_or_default();
    if is_confident_male_sex_chromosome(chrom, inferred_sex) && alleles.len() == 2 && alleles[0] == alleles[1] {
        let allele = alleles[0];
        if allele == ref_ch {
            return ("0".to_owned(), "hem_ref".to_owned());
        }
        if allele == alt_ch {
            return ("1".to_owned(), "hem_alt".to_owned());
        }
    }
    if alleles.len() != 2 {
        return (display.to_owned(), "unknown".to_owned());
    }
    let alt_count = alleles.iter().filter(|allele| **allele == alt_ch).count();
    let ref_count = alleles.iter().filter(|allele| **allele == ref_ch).count();
    match (ref_count, alt_count) {
        (2, 0) => ("0/0".to_owned(), "hom_ref".to_owned()),
        (1, 1) => ("0/1".to_owned(), "het".to_owned()),
        (0, 2) => ("1/1".to_owned(), "hom_alt".to_owned()),
        _ => (display.to_owned(), "unknown".to_owned()),
    }
}

fn is_confident_male_sex_chromosome(chrom: &str, inferred_sex: Option<&SexInference>) -> bool {
    matches!(
        chrom.trim().trim_start_matches("chr").to_ascii_uppercase().as_str(),
        "X" | "Y" | "23" | "24"
    ) && inferred_sex.is_some_and(|sex| {
        sex.sex == InferredSex::Male
            && matches!(sex.confidence, SexDetectionConfidence::High | SexDetectionConfidence::Medium)
    })
}

fn assembly_row_value(assembly: Assembly) -> String {
    match assembly {
        Assembly::Grch37 => "grch37".to_owned(),
        Assembly::Grch38 => "grch38".to_owned(),
    }
}

fn inferred_sex_name(value: InferredSex) -> &'static str {
    match value {
        InferredSex::Male => "male",
        InferredSex::Female => "female",
        InferredSex::Unknown => "unknown",
    }
}

fn sex_detection_confidence_name(value: SexDetectionConfidence) -> &'static str {
    match value {
        SexDetectionConfidence::High => "high",
        SexDetectionConfidence::Medium => "medium",
        SexDetectionConfidence::Low => "low",
    }
}

fn parse_optional_u32(value: Option<&String>) -> Option<u32> {
    value.and_then(|value| value.parse::<u32>().ok())
}
