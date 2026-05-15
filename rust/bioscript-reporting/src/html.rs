use std::fmt::Write as _;

mod analysis;
mod helpers;
mod observations;
mod pgx;
mod provenance;
mod sections;

use analysis::render_analysis_tables;
use helpers::html_escape;
use observations::render_observation_table;
use pgx::render_pgx_table;
use provenance::render_provenance_links;
use sections::{
    collect_report_analyses, collect_report_findings, collect_report_participants,
    render_input_debug, render_participant_filter, render_report_manifest_header,
    render_report_source_section,
};

pub fn render_app_html_document(
    observations: &[serde_json::Value],
    reports: &[serde_json::Value],
) -> Result<String, String> {
    let mut out = String::from(
        r##"<!doctype html><meta charset="utf-8"><title>BioScript report</title><style>body{font-family:system-ui,sans-serif;margin:0;background:#f7f8fa;color:#1f2933}.wrap{max-width:1440px;margin:0 auto;padding:24px}h1{margin:0 0 10px}h2{margin:32px 0 10px;scroll-margin-top:82px}.nav{position:sticky;top:0;z-index:20;display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:16px -24px 22px;padding:10px 24px;background:rgba(247,248,250,.96);border-block:1px solid #d8dee6;backdrop-filter:saturate(160%) blur(8px)}.nav a{border:1px solid #cbd5df;background:#fff;color:#1f2933;text-decoration:none;padding:7px 10px;border-radius:6px}.nav a:hover{background:#eef2f6}.table-tools,.level-filter{display:flex;justify-content:space-between;gap:12px;align-items:flex-start;margin:6px 0}.table-tools input{width:min(420px,100%);border:1px solid #cbd5df;border-radius:6px;padding:7px 9px;font:inherit;background:#fff}.table-tools label,.level-filter label{display:flex;gap:5px;align-items:center}.level-filter{justify-content:flex-start;flex-wrap:wrap}.filter-scale{display:flex;flex-direction:column;gap:6px}.filter-scale-title{font-weight:700;color:#344054;display:flex;gap:6px;align-items:center}.filter-options{display:flex;flex-wrap:wrap;gap:6px 10px}.filter-actions{display:flex;gap:8px;align-self:flex-end}.level-filter a{display:inline-grid;place-items:center;width:20px;height:20px;border:1px solid #cbd5df;border-radius:50%;text-decoration:none;color:#1f2933;background:#fff;font-weight:700}.filter-action,.pgx-tabs button{border:1px solid #cbd5df;background:#fff;color:#1f2933;border-radius:6px;padding:4px 7px;font:inherit;cursor:pointer}.filter-action:hover,.pgx-tabs button:hover{background:#eef2f6}.pgx-tabs{display:flex;gap:8px;margin:10px 0}.pgx-tabs button.active{background:#1f2933;border-color:#1f2933;color:#fff}.table-wrap{overflow:auto;border:1px solid #d8dee6;background:white;border-radius:8px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #e5e9ef;padding:6px 8px;text-align:left;vertical-align:top}th{position:sticky;top:0;background:#eef2f6;z-index:1;white-space:nowrap;cursor:pointer;user-select:none}.sort-mark{font-size:10px;color:#667085;margin-left:3px}.row-variant td{background:#fff7cc}.row-reference td{background:#eaf7ee}.debug-hidden .debug-col{display:none}.participant-filter{display:flex;gap:8px;align-items:center;margin:12px 0}.participant-filter select{border:1px solid #cbd5df;border-radius:6px;padding:6px 8px;background:#fff;font:inherit}.genotype-hit{font-weight:700}.allele-hit{color:#075985;background:#dff4ff;border-radius:3px;padding:0 2px}.level-badge,.pgx-badge{display:inline-block;min-width:2.2em;text-align:center;border-radius:999px;padding:2px 8px;color:#fff;font-weight:700}.level-1,.pgx-informative{background:#0abc72}.level-2,.pgx-actionable{background:#2a74df}.level-3,.pgx-recommended{background:#ffc107;color:#1f2933}.level-4,.pgx-required{background:#c53b3b}.pgx-no-clinical,.pgx-criteria,.pgx-unknown,.level-unknown{background:#667085}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.muted{color:#667085}.effect{max-width:760px;min-width:360px}.analysis-kv{display:grid;grid-template-columns:max-content minmax(0,1fr);gap:6px 14px;background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 12px;margin:8px 0 10px}.analysis-kv dt{font-weight:700;color:#344054}.analysis-kv dd{margin:0;min-width:0}.analysis-badge{display:inline-block;border:1px solid #cbd5df;border-radius:999px;background:#eef2f6;padding:1px 8px;font-weight:700}.analysis-badge-normal{background:#eaf7ee;border-color:#9fd6ad;color:#14532d}.analysis-badge-variant{background:#fff7cc;border-color:#f0d66a;color:#713f12}.analysis-badge-unknown{background:#eef2f6;border-color:#cbd5df;color:#475467}.analysis-card{background:#fff;border:1px solid #cbd5df;border-radius:8px;margin:14px 0;overflow:hidden}.analysis-card summary{display:flex;align-items:center;gap:10px;padding:11px 13px;background:#eef2f6;cursor:pointer;font-weight:700}.analysis-card-body{padding:12px 14px}.analysis-card-title{font-size:16px}.logic-note,.analysis-notes{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 12px;margin:8px 0 10px}.logic-note h4,.analysis-notes h4,h4{margin:0 0 6px}.logic-note p,.analysis-notes p{margin:0 0 6px}.analysis-notes p:last-child{margin-bottom:0}.provenance-list{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 18px}.provenance-list li{margin:8px 0}pre{white-space:pre-wrap;background:#fff;padding:12px;border:1px solid #d8dee6;border-radius:8px;max-height:520px;overflow:auto}</style><script>const sortState={};function cellText(row,i){const cell=row.cells[i];return(cell?.dataset.sort||cell?.innerText||"").trim()}function cmp(a,b){const an=Number(a),bn=Number(b);if(a!==""&&b!==""&&!Number.isNaN(an)&&!Number.isNaN(bn))return an-bn;return a.localeCompare(b,undefined,{numeric:true,sensitivity:"base"})}function sortTable(id,col){const table=document.getElementById(id);const tbody=table.tBodies[0];const key=id+":"+col;const dir=sortState[key]==="asc"?"desc":"asc";sortState[key]=dir;table.querySelectorAll(".sort-mark").forEach(s=>s.textContent="");table.tHead.rows[0].cells[col].querySelector(".sort-mark").textContent=dir==="asc"?"^":"v";Array.from(tbody.rows).sort((a,b)=>{const v=cmp(cellText(a,col),cellText(b,col));return dir==="asc"?v:-v}).forEach(r=>tbody.appendChild(r))}let selectedParticipant="";function participantOk(row){return !selectedParticipant||!row.dataset.participant||row.dataset.participant===selectedParticipant}function pgxLevelOk(row){if(!row.dataset.pgxLevel&&!row.dataset.level)return true;const key=row.dataset.pgxLevel||row.dataset.level||"unknown";const input=document.querySelector('[data-pgx-any-level-filter="'+key+'"]')||document.querySelector('[data-pgx-any-level-filter="unknown"]');return !!input&&input.checked}function pgxOutcomeOk(row){const key=row.dataset.pgxOutcome;if(!key)return true;const input=document.querySelector('[data-pgx-outcome-filter="'+key+'"]');return !input||input.checked}function observationOk(row){const key=row.dataset.observation;if(!key)return true;const input=document.querySelector('[data-observation-filter="'+key+'"]');return !input||input.checked}function applyTableFilters(id){const q=(document.querySelector('[data-filter-for="'+id+'"]')?.value||"").toLowerCase();document.querySelectorAll("#"+id+" tbody tr").forEach(row=>{let ok=row.innerText.toLowerCase().includes(q)&&participantOk(row)&&pgxLevelOk(row)&&pgxOutcomeOk(row)&&observationOk(row);row.style.display=ok?"":"none"})}function applyPgxFilters(){document.querySelectorAll('#pgx table[id]').forEach(table=>applyTableFilters(table.id));document.querySelectorAll('#pgx .pgx-drug-group').forEach(group=>{const visible=Array.from(group.querySelectorAll('tbody tr')).some(row=>row.style.display!=='none');group.style.display=visible?'':'none'})}function setPgxFilterGroup(checked){document.querySelectorAll("[data-pgx-any-level-filter]").forEach(input=>input.checked=checked);applyPgxFilters()}function setObservationFilterGroup(checked){document.querySelectorAll("[data-observation-filter]").forEach(input=>input.checked=checked);applyTableFilters("observations-table")}function setPgxView(view){document.getElementById("pgx-view-variant").hidden=view!=="variant";document.getElementById("pgx-view-drug").hidden=view!=="drug";document.getElementById("pgx-tab-variant")?.classList.toggle("active",view==="variant");document.getElementById("pgx-tab-drug")?.classList.toggle("active",view==="drug")}function setParticipant(value){selectedParticipant=value;document.querySelectorAll("table[id]").forEach(table=>applyTableFilters(table.id))}function toggleDebug(show){document.getElementById("report-wrap").classList.toggle("debug-hidden",!show)}</script><div class="wrap debug-hidden" id="report-wrap">"##,
    );
    out.push_str(
        r#"<style>.ref-alt-cell{max-width:38rem}.ref-alt-cell .truncate-cell{display:inline-block;max-width:min(38rem,62vw);overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:bottom}</style>"#,
    );
    let label_findings = collect_report_findings(reports, "bioscript:pgx-label:1.0");
    let summary_findings = collect_report_findings(reports, "bioscript:pgx-summary:1.0");
    let has_pgx_findings = !label_findings.is_empty() || !summary_findings.is_empty();
    let analysis_outputs = collect_report_analyses(reports);
    let participants = collect_report_participants(reports);
    render_report_manifest_header(&mut out, reports);
    let _ = write!(
        out,
        "<div class=\"muted\">{} observation(s), {} analysis output(s), {} PGx label finding(s), {} PGx summary finding(s)</div>",
        observations.len(),
        analysis_outputs.len(),
        label_findings.len(),
        summary_findings.len()
    );
    render_participant_filter(&mut out, &participants);
    out.push_str("<nav class=\"nav\"><a href=\"#input-info\">Input</a><a href=\"#observations\">Observations</a><a href=\"#analysis\">Analysis</a>");
    if has_pgx_findings {
        out.push_str("<a href=\"#pgx\">PGx</a>");
    }
    out.push_str("<a href=\"#provenance\">Provenance</a><a href=\"#source\">Source</a><a href=\"#json\">Raw JSON</a></nav>");
    out.push_str("<section id=\"input-info\"><h2>Input</h2>");
    render_input_debug(&mut out, reports, participants.len() > 1);
    out.push_str("</section>");
    out.push_str("<section id=\"observations\"><h2>Observations</h2>");
    render_observation_table(&mut out, observations, participants.len() > 1);
    out.push_str("</section>");
    out.push_str("<section id=\"analysis\"><h2>Analysis</h2>");
    render_analysis_tables(
        &mut out,
        &analysis_outputs,
        observations,
        participants.len() > 1,
    );
    out.push_str("</section>");
    if has_pgx_findings {
        out.push_str("<section id=\"pgx\"><h2>PGx</h2>");
        render_pgx_table(&mut out, &label_findings, &summary_findings);
        out.push_str("</section>");
    }
    out.push_str("<section id=\"provenance\"><h2>Provenance</h2>");
    render_provenance_links(&mut out, reports);
    out.push_str("</section>");
    out.push_str("<section id=\"source\"><h2>Source</h2>");
    render_report_source_section(&mut out, reports);
    out.push_str("</section>");
    out.push_str("<section id=\"json\"><h2>Raw Reports JSON</h2><details><summary>Show raw report JSON</summary>");
    for report in reports {
        let text = serde_json::to_string_pretty(report).map_err(|err| err.to_string())?;
        let _ = write!(out, "<pre>{}</pre>", html_escape(&text));
    }
    out.push_str("</details></section></div>");
    Ok(out)
}
