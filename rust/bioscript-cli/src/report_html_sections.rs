fn collect_report_analyses(reports: &[serde_json::Value]) -> Vec<serde_json::Value> {
    reports
        .iter()
        .filter_map(|report| report.get("analyses").and_then(serde_json::Value::as_array))
        .flat_map(|analyses| analyses.iter())
        .cloned()
        .collect()
}

fn collect_report_findings(reports: &[serde_json::Value], schema: &str) -> Vec<serde_json::Value> {
    reports
        .iter()
        .filter_map(|report| report.get("findings").and_then(serde_json::Value::as_array))
        .flat_map(|findings| findings.iter())
        .filter(|finding| finding.get("schema").and_then(serde_json::Value::as_str) == Some(schema))
        .cloned()
        .collect()
}

fn render_analysis_tables(out: &mut String, analyses: &[serde_json::Value]) {
    if analyses.is_empty() {
        out.push_str("<p class=\"muted\">No analysis outputs.</p>");
        return;
    }
    for (index, analysis) in analyses.iter().enumerate() {
        let table_id = format!("analysis-table-{index}");
        let title = format!(
            "{} / {}",
            value_str(analysis, "participant_id"),
            value_str(analysis, "analysis_id")
        );
        let _ = write!(out, "<h3>{}</h3>", html_escape(&title));
        render_analysis_logic(out, analysis);
        let rows = analysis
            .get("rows")
            .and_then(serde_json::Value::as_array)
            .cloned()
            .unwrap_or_default();
        if rows.is_empty() {
            out.push_str("<p class=\"muted\">No rows emitted.</p>");
            continue;
        }
        let headers = analysis_row_headers(&rows);
        let header_refs = headers.iter().map(String::as_str).collect::<Vec<_>>();
        render_table_start(out, &table_id, &header_refs);
        for row in rows {
            out.push_str("<tr>");
            for header in &headers {
                table_cell(out, &json_field_as_tsv(row.get(header)));
            }
            out.push_str("</tr>");
        }
        render_table_end(out);
    }
}

fn analysis_row_headers(rows: &[serde_json::Value]) -> Vec<String> {
    let mut headers = Vec::new();
    for row in rows {
        let Some(object) = row.as_object() else {
            continue;
        };
        for key in object.keys() {
            if !headers.contains(key) {
                headers.push(key.clone());
            }
        }
    }
    headers
}

fn render_analysis_logic(out: &mut String, analysis: &serde_json::Value) {
    let Some(logic) = analysis.get("logic") else {
        return;
    };
    if logic.is_null() {
        return;
    }
    let description = logic
        .get("description")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let source = logic.get("source").unwrap_or(&serde_json::Value::Null);
    let source_name = source
        .get("name")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("source");
    let source_url = source
        .get("url")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    out.push_str("<div class=\"logic-note\">");
    if !description.is_empty() {
        let _ = write!(out, "<p>{}</p>", html_escape(description));
    }
    if !source_url.is_empty() {
        let _ = write!(
            out,
            "<p class=\"muted\">Logic source: <a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">{}</a></p>",
            html_escape(source_url),
            html_escape(source_name)
        );
    }
    out.push_str("</div>");
}

fn render_provenance_links(out: &mut String, reports: &[serde_json::Value]) {
    let mut links = BTreeMap::<String, String>::new();
    for report in reports {
        collect_provenance_links_from_value(report, &mut links);
    }
    if links.is_empty() {
        out.push_str("<p class=\"muted\">No provenance links.</p>");
        return;
    }
    out.push_str("<ul class=\"provenance-list\">");
    for (url, label) in links {
        let display = if label.is_empty() { &url } else { &label };
        let _ = write!(
            out,
            "<li><a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">{}</a><div class=\"muted mono\">{}</div></li>",
            html_escape(&url),
            html_escape(display),
            html_escape(&url)
        );
    }
    out.push_str("</ul>");
}

fn collect_provenance_links_from_value(
    value: &serde_json::Value,
    links: &mut BTreeMap<String, String>,
) {
    match value {
        serde_json::Value::Object(object) => {
            if let Some(url) = object.get("url").and_then(serde_json::Value::as_str)
                && url.starts_with("http")
            {
                let label = object
                    .get("name")
                    .or_else(|| object.get("label"))
                    .or_else(|| object.get("source"))
                    .and_then(value_as_string)
                    .unwrap_or_default();
                links.entry(url.to_owned()).or_insert(label);
            }
            for child in object.values() {
                collect_provenance_links_from_value(child, links);
            }
        }
        serde_json::Value::Array(items) => {
            for item in items {
                collect_provenance_links_from_value(item, links);
            }
        }
        _ => {}
    }
}

fn render_observation_table(out: &mut String, observations: &[serde_json::Value]) {
    let headers = [
        "participant_id",
        "rsid",
        "ref",
        "alt",
        "genotype_display",
        "genotype",
        "zygosity",
        "outcome",
        "match_status",
        "coverage_status",
        "call_status",
        "assembly",
        "chrom",
        "pos_start",
        "pos_end",
        "kind",
        "ref_count",
        "alt_count",
        "depth",
        "genotype_quality",
        "allele_balance",
        "evidence_type",
        "evidence_raw",
        "facets",
        "assay_id",
        "assay_version",
        "variant_key",
    ];
    render_table_start(out, "observations-table", &headers);
    for observation in observations {
        let _ = write!(out, "<tr class=\"{}\">", observation_row_class(observation));
        for header in headers {
            render_observation_cell(out, observation, header);
        }
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></div>");
}

fn observation_row_class(observation: &serde_json::Value) -> &'static str {
    match observation
        .get("outcome")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "variant" => "row-variant",
        "reference" => "row-reference",
        _ => "",
    }
}

fn render_observation_cell(out: &mut String, observation: &serde_json::Value, header: &str) {
    if header == "genotype_display" {
        let outcome = observation
            .get("outcome")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let value = json_field_as_tsv(observation.get(header));
        if outcome == "variant" {
            let alt = observation
                .get("alt")
                .and_then(serde_json::Value::as_str)
                .unwrap_or_default();
            let _ = write!(
                out,
                "<td class=\"genotype-hit\">{}</td>",
                highlight_allele(&value, alt)
            );
            return;
        }
    }
    let _ = write!(
        out,
        "<td>{}</td>",
        html_escape(&json_field_as_tsv(observation.get(header)))
    );
}

fn highlight_allele(value: &str, allele: &str) -> String {
    if value.is_empty() || allele.is_empty() {
        return html_escape(value);
    }
    if allele.chars().count() == 1 {
        let target = allele
            .chars()
            .next()
            .unwrap_or_default()
            .to_ascii_uppercase();
        let mut out = String::new();
        for ch in value.chars() {
            let escaped = html_escape(&ch.to_string());
            if ch.to_ascii_uppercase() == target {
                let _ = write!(out, "<span class=\"allele-hit\">{escaped}</span>");
            } else {
                out.push_str(&escaped);
            }
        }
        return out;
    }
    let escaped_value = html_escape(value);
    let escaped_allele = html_escape(allele);
    escaped_value.replace(
        &escaped_allele,
        &format!("<span class=\"allele-hit\">{escaped_allele}</span>"),
    )
}

