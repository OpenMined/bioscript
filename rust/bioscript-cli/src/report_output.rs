fn app_report_json(
    assay_id: &str,
    participant_id: &str,
    input_file: &Path,
    observations: &[serde_json::Value],
    analyses: &[serde_json::Value],
    findings: &[serde_json::Value],
    provenance: &[serde_json::Value],
) -> serde_json::Value {
    let called = observations
        .iter()
        .filter(|item| {
            item.get("call_status").and_then(serde_json::Value::as_str) == Some("called")
        })
        .count();
    serde_json::json!({
        "schema": "bioscript:report:1.0",
        "version": "1.0",
        "participant_id": participant_id,
        "assay_id": assay_id,
        "assay_version": "1.0",
        "input": {
            "file_name": input_file.file_name().and_then(|value| value.to_str()).unwrap_or_default(),
            "file_path": input_file.display().to_string(),
        },
        "report_status": if called == observations.len() { "complete" } else { "partial" },
        "derived_from": observations.iter().filter_map(|item| item.get("variant_key").cloned()).collect::<Vec<_>>(),
        "analyses": analyses,
        "findings": findings,
        "provenance": provenance,
        "metrics": {
            "n_sites_tested": observations.len(),
            "n_sites_called": called,
            "n_sites_missing": observations.len().saturating_sub(called),
            "n_analyses": analyses.len(),
            "n_findings_matched": findings.len(),
        }
    })
}

fn write_app_observations(
    output_dir: &Path,
    observations: &[serde_json::Value],
    format: AppOutputFormat,
) -> Result<(), String> {
    if matches!(format, AppOutputFormat::Tsv | AppOutputFormat::Both) {
        let mut out = bioscript_core::OBSERVATION_TSV_HEADERS.join("\t");
        out.push('\n');
        for observation in observations {
            let line = bioscript_core::OBSERVATION_TSV_HEADERS
                .iter()
                .map(|header| json_field_as_tsv(observation.get(*header)))
                .collect::<Vec<_>>()
                .join("\t");
            out.push_str(&line);
            out.push('\n');
        }
        fs::write(output_dir.join("observations.tsv"), out)
            .map_err(|err| format!("failed to write observations.tsv: {err}"))?;
    }
    if matches!(format, AppOutputFormat::Jsonl | AppOutputFormat::Both) {
        write_jsonl(&output_dir.join("observations.jsonl"), observations)?;
    }
    if matches!(format, AppOutputFormat::Json) {
        write_json_pretty(
            &output_dir.join("observations.json"),
            &serde_json::json!({"observations": observations}),
        )?;
    }
    Ok(())
}

fn write_app_analyses(output_dir: &Path, analyses: &[serde_json::Value]) -> Result<(), String> {
    write_jsonl(&output_dir.join("analysis.jsonl"), analyses)
}

fn write_app_reports(
    output_dir: &Path,
    reports: &[serde_json::Value],
    format: AppOutputFormat,
) -> Result<(), String> {
    if matches!(format, AppOutputFormat::Jsonl | AppOutputFormat::Both) {
        write_jsonl(&output_dir.join("reports.jsonl"), reports)?;
    }
    if matches!(format, AppOutputFormat::Json | AppOutputFormat::Both) {
        write_json_pretty(
            &output_dir.join("reports.json"),
            &serde_json::json!({
                "schema": "bioscript:report-set:1.0",
                "version": "1.0",
                "reports": reports,
            }),
        )?;
    }
    Ok(())
}

fn write_jsonl(path: &Path, rows: &[serde_json::Value]) -> Result<(), String> {
    let mut out = String::new();
    for row in rows {
        let line = serde_json::to_string(row).map_err(|err| err.to_string())?;
        out.push_str(&line);
        out.push('\n');
    }
    fs::write(path, out).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn write_json_pretty(path: &Path, value: &serde_json::Value) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|err| err.to_string())?;
    fs::write(path, text).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}

fn write_app_html(
    output_dir: &Path,
    observations: &[serde_json::Value],
    reports: &[serde_json::Value],
) -> Result<(), String> {
    let mut out = String::from(
        r##"<!doctype html><meta charset="utf-8"><title>BioScript report</title><style>body{font-family:system-ui,sans-serif;margin:0;background:#f7f8fa;color:#1f2933}.wrap{max-width:1440px;margin:0 auto;padding:24px}h1{margin:0 0 10px}h2{margin:32px 0 10px;scroll-margin-top:82px}.nav{position:sticky;top:0;z-index:20;display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:16px -24px 22px;padding:10px 24px;background:rgba(247,248,250,.96);border-block:1px solid #d8dee6;backdrop-filter:saturate(160%) blur(8px)}.nav a{border:1px solid #cbd5df;background:#fff;color:#1f2933;text-decoration:none;padding:7px 10px;border-radius:6px}.nav a:hover{background:#eef2f6}.table-tools,.level-filter{display:flex;justify-content:space-between;gap:12px;align-items:center;margin:6px 0}.table-tools input{width:min(420px,100%);border:1px solid #cbd5df;border-radius:6px;padding:7px 9px;font:inherit;background:#fff}.table-tools label,.level-filter label{display:flex;gap:5px;align-items:center}.level-filter{justify-content:flex-start;flex-wrap:wrap}.level-filter a{display:inline-grid;place-items:center;width:20px;height:20px;border:1px solid #cbd5df;border-radius:50%;text-decoration:none;color:#1f2933;background:#fff;font-weight:700}.filter-action{border:1px solid #cbd5df;background:#fff;color:#1f2933;border-radius:6px;padding:4px 7px;font:inherit;cursor:pointer}.filter-action:hover{background:#eef2f6}.table-wrap{overflow:auto;border:1px solid #d8dee6;background:white;border-radius:8px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #e5e9ef;padding:6px 8px;text-align:left;vertical-align:top}th{position:sticky;top:0;background:#eef2f6;z-index:1;white-space:nowrap;cursor:pointer;user-select:none}.sort-mark{font-size:10px;color:#667085;margin-left:3px}.row-variant td{background:#fff7cc}.row-reference td{background:#eaf7ee}.refs-hidden .row-reference{display:none}.genotype-hit{font-weight:700}.allele-hit{color:#075985;background:#dff4ff;border-radius:3px;padding:0 2px}.level-badge,.pgx-badge{display:inline-block;min-width:2.2em;text-align:center;border-radius:999px;padding:2px 8px;color:#fff;font-weight:700}.level-1,.pgx-informative{background:#0abc72}.level-2,.pgx-actionable{background:#2a74df}.level-3,.pgx-recommended{background:#ffc107;color:#1f2933}.level-4,.pgx-required{background:#c53b3b}.pgx-no-clinical,.pgx-criteria{background:#667085}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.muted{color:#667085}.effect{max-width:760px;min-width:360px}.logic-note{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 12px;margin:8px 0 10px}.logic-note p{margin:0 0 6px}.provenance-list{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 18px}.provenance-list li{margin:8px 0}pre{white-space:pre-wrap;background:#fff;padding:12px;border:1px solid #d8dee6;border-radius:8px;max-height:520px;overflow:auto}</style><script>const sortState={};function cellText(row,i){return(row.cells[i]?.innerText||"").trim()}function cmp(a,b){const an=Number(a),bn=Number(b);if(a!==""&&b!==""&&!Number.isNaN(an)&&!Number.isNaN(bn))return an-bn;return a.localeCompare(b,undefined,{numeric:true,sensitivity:"base"})}function sortTable(id,col){const table=document.getElementById(id);const tbody=table.tBodies[0];const key=id+":"+col;const dir=sortState[key]==="asc"?"desc":"asc";sortState[key]=dir;table.querySelectorAll(".sort-mark").forEach(s=>s.textContent="");table.tHead.rows[0].cells[col].querySelector(".sort-mark").textContent=dir==="asc"?"^":"v";Array.from(tbody.rows).sort((a,b)=>{const v=cmp(cellText(a,col),cellText(b,col));return dir==="asc"?v:-v}).forEach(r=>tbody.appendChild(r))}function applyTableFilters(id){const q=(document.querySelector('[data-filter-for="'+id+'"]')?.value||"").toLowerCase();document.querySelectorAll("#"+id+" tbody tr").forEach(row=>{let ok=row.innerText.toLowerCase().includes(q);if(id==="summaries-table"&&row.dataset.level){ok=ok&&!!document.querySelector('[data-level-filter="'+row.dataset.level+'"]:checked')}if(id==="labels-table"&&row.dataset.pgxLevel){ok=ok&&!!document.querySelector('[data-pgx-level-filter="'+row.dataset.pgxLevel+'"]:checked')}row.style.display=ok?"":"none"})}function setFilterGroup(group,id,checked){document.querySelectorAll(group==="level"?"[data-level-filter]":"[data-pgx-level-filter]").forEach(input=>input.checked=checked);applyTableFilters(id)}function toggleRefs(show){document.getElementById("report-wrap").classList.toggle("refs-hidden",!show)}</script><div class="wrap refs-hidden" id="report-wrap"><h1>BioScript Report</h1>"##,
    );
    let label_findings = collect_report_findings(reports, "bioscript:pgx-label:1.0");
    let summary_findings = collect_report_findings(reports, "bioscript:pgx-summary:1.0");
    let analysis_outputs = collect_report_analyses(reports);
    let _ = write!(
        out,
        "<div class=\"muted\">{} observation(s), {} analysis output(s), {} PGx label finding(s), {} PGx summary finding(s)</div>",
        observations.len(),
        analysis_outputs.len(),
        label_findings.len(),
        summary_findings.len()
    );
    out.push_str("<nav class=\"nav\"><a href=\"#observations\">Observations</a><a href=\"#analysis\">Analysis</a><a href=\"#labels\">PGx Labels</a><a href=\"#summaries\">PGx Summaries</a><a href=\"#provenance\">Provenance</a><a href=\"#json\">Raw JSON</a></nav>");
    out.push_str("<section id=\"observations\"><h2>Observations</h2>");
    render_observation_table(&mut out, observations);
    out.push_str("</section>");
    out.push_str("<section id=\"analysis\"><h2>Analysis</h2>");
    render_analysis_tables(&mut out, &analysis_outputs);
    out.push_str("</section>");
    out.push_str("<section id=\"labels\"><h2>PGx Label Annotations</h2>");
    render_pgx_label_table(&mut out, &label_findings);
    out.push_str("</section>");
    out.push_str("<section id=\"summaries\"><h2>PGx Summary Annotations</h2>");
    render_pgx_summary_table(&mut out, &summary_findings);
    out.push_str("</section>");
    out.push_str("<section id=\"provenance\"><h2>Provenance</h2>");
    render_provenance_links(&mut out, reports);
    out.push_str("</section>");
    out.push_str("<section id=\"json\"><h2>Raw Reports JSON</h2>");
    for report in reports {
        let text = serde_json::to_string_pretty(report).map_err(|err| err.to_string())?;
        let _ = write!(out, "<pre>{}</pre>", html_escape(&text));
    }
    out.push_str("</section></div>");
    fs::write(output_dir.join("index.html"), out)
        .map_err(|err| format!("failed to write index.html: {err}"))
}

