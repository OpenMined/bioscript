use super::helpers::{
    html_escape, json_field_as_tsv, render_table_end, render_table_start, table_cell,
    table_header_label, value_str,
};
use std::fmt::Write as _;
pub(super) fn render_analysis_tables(
    out: &mut String,
    analyses: &[serde_json::Value],
    observations: &[serde_json::Value],
    show_participant_id: bool,
) {
    if analyses.is_empty() {
        out.push_str("<p class=\"muted\">No analysis outputs.</p>");
        return;
    }
    for (index, analysis) in analyses.iter().enumerate() {
        let table_id = format!("analysis-table-{index}");
        let title = analysis_title(analysis);
        let weak_indel_dependency = analysis_depends_on_weak_observation(analysis, observations);
        let row_count = analysis
            .get("rows")
            .and_then(serde_json::Value::as_array)
            .map_or(0, Vec::len);
        let _ = write!(
            out,
            "<details class=\"analysis-card\" open><summary><span class=\"analysis-card-title\">{}</span><span class=\"muted\">{} result row(s)</span></summary><div class=\"analysis-card-body\">",
            html_escape(&title),
            row_count
        );
        render_analysis_logic(out, analysis);
        let rows = analysis
            .get("rows")
            .and_then(serde_json::Value::as_array)
            .cloned()
            .unwrap_or_default();
        if rows.is_empty() {
            out.push_str("<p class=\"muted\">No rows emitted.</p>");
            out.push_str("</div></details>");
            continue;
        }
        let headers = analysis_row_headers(analysis, &rows, show_participant_id);
        let notes = analysis_notes(&rows);
        if rows.len() == 1 {
            out.push_str("<h4>Results</h4>");
            render_analysis_key_values(out, analysis, &rows[0], &headers);
            render_analysis_notes(out, &notes);
            render_weak_indel_analysis_note(out, weak_indel_dependency);
            out.push_str("</div></details>");
            continue;
        }
        out.push_str("<h4>Results</h4>");
        let header_refs = headers.iter().map(String::as_str).collect::<Vec<_>>();
        render_table_start(out, &table_id, &header_refs);
        let participant = value_str(analysis, "participant_id");
        for row in rows {
            let _ = write!(
                out,
                "<tr data-participant=\"{}\">",
                html_escape(participant)
            );
            for header in &headers {
                table_cell(out, &json_field_as_tsv(row.get(header)));
            }
            out.push_str("</tr>");
        }
        render_table_end(out);
        render_analysis_notes(out, &notes);
        render_weak_indel_analysis_note(out, weak_indel_dependency);
        out.push_str("</div></details>");
    }
}

pub(super) fn analysis_title(analysis: &serde_json::Value) -> String {
    let label = value_str(analysis, "analysis_label");
    if label.is_empty() {
        value_str(analysis, "analysis_id").to_owned()
    } else {
        label.to_owned()
    }
}

pub(super) fn analysis_row_headers(
    analysis: &serde_json::Value,
    rows: &[serde_json::Value],
    show_participant_id: bool,
) -> Vec<String> {
    let mut headers = Vec::new();
    if let Some(emits) = analysis.get("emits").and_then(serde_json::Value::as_array) {
        for key in emits
            .iter()
            .filter_map(|emit| emit.get("key"))
            .filter_map(serde_json::Value::as_str)
        {
            if should_show_analysis_header(key, show_participant_id)
                && rows.iter().any(|row| row.get(key).is_some())
                && !headers.iter().any(|item| item == key)
            {
                headers.push(key.to_owned());
            }
        }
    }
    if let Some(row_headers) = analysis
        .get("row_headers")
        .and_then(serde_json::Value::as_array)
    {
        for header in row_headers.iter().filter_map(serde_json::Value::as_str) {
            if should_show_analysis_header(header, show_participant_id)
                && !headers.iter().any(|item| item == header)
            {
                headers.push(header.to_owned());
            }
        }
    }
    for row in rows {
        let Some(object) = row.as_object() else {
            continue;
        };
        for key in object.keys() {
            if !should_show_analysis_header(key, show_participant_id) {
                continue;
            }
            if !headers.contains(key) {
                headers.push(key.clone());
            }
        }
    }
    headers
}

pub(super) fn should_show_analysis_header(key: &str, show_participant_id: bool) -> bool {
    (show_participant_id || key != "participant_id") && key != "notes" && key != "report_notes"
}

pub(super) fn analysis_notes(rows: &[serde_json::Value]) -> Vec<String> {
    let mut notes = Vec::new();
    for row in rows {
        if let Some(note) = row
            .get("notes")
            .or_else(|| row.get("report_notes"))
            .and_then(serde_json::Value::as_str)
            && !note.trim().is_empty()
            && !notes.iter().any(|item: &String| item == note)
        {
            notes.push(note.to_owned());
        }
    }
    notes
}

pub(super) fn render_analysis_key_values(
    out: &mut String,
    analysis: &serde_json::Value,
    row: &serde_json::Value,
    headers: &[String],
) {
    out.push_str("<dl class=\"analysis-kv\">");
    for header in headers {
        let value = json_field_as_tsv(row.get(header));
        let _ = write!(
            out,
            "<dt>{}</dt><dd>{}</dd>",
            html_escape(&analysis_header_label(analysis, header)),
            render_analysis_value(header, &value)
        );
    }
    out.push_str("</dl>");
}

pub(super) fn analysis_header_label(analysis: &serde_json::Value, key: &str) -> String {
    analysis
        .get("emits")
        .and_then(serde_json::Value::as_array)
        .and_then(|emits| {
            emits.iter().find_map(|emit| {
                if emit.get("key").and_then(serde_json::Value::as_str) == Some(key) {
                    emit.get("label")
                        .and_then(serde_json::Value::as_str)
                        .map(ToOwned::to_owned)
                } else {
                    None
                }
            })
        })
        .unwrap_or_else(|| table_header_label(key))
}

pub(super) fn render_analysis_value(key: &str, value: &str) -> String {
    if value.starts_with("http://") || value.starts_with("https://") {
        return format!(
            "<a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">Source</a>",
            html_escape(value)
        );
    }
    if is_analysis_badge_key(key) {
        format!(
            "<span class=\"analysis-badge {}\">{}</span>",
            analysis_badge_class(value),
            html_escape(value)
        )
    } else {
        html_escape(value)
    }
}

pub(super) fn is_analysis_badge_key(key: &str) -> bool {
    key.ends_with("_status") || key.ends_with("_outcome")
}

pub(super) fn analysis_badge_class(value: &str) -> &'static str {
    match value {
        "normal" | "reference" => "analysis-badge-normal",
        "variant" => "analysis-badge-variant",
        "unknown" | "unresolved_missing_variant" => "analysis-badge-unknown",
        _ => "",
    }
}

pub(super) fn render_analysis_notes(out: &mut String, notes: &[String]) {
    if notes.is_empty() {
        return;
    }
    out.push_str("<div class=\"analysis-notes\">");
    out.push_str("<h4>Notes</h4>");
    for note in notes {
        let _ = write!(out, "<p>{}</p>", html_escape(note));
    }
    out.push_str("</div>");
}

pub(super) fn render_weak_indel_analysis_note(out: &mut String, weak_indel_dependency: bool) {
    if !weak_indel_dependency {
        return;
    }
    out.push_str("<div class=\"analysis-notes\"><h4>Notes</h4><p>* Result depends on a weak indel match from a consumer genotype file.</p></div>");
}

pub(super) fn analysis_depends_on_weak_observation(
    analysis: &serde_json::Value,
    observations: &[serde_json::Value],
) -> bool {
    let weak_paths = observations
        .iter()
        .filter(|observation| analysis_observation_is_weak_indel_match(observation))
        .filter_map(|observation| {
            observation
                .get("variant_path")
                .and_then(serde_json::Value::as_str)
        })
        .collect::<Vec<_>>();
    if weak_paths.is_empty() {
        return false;
    }
    analysis
        .get("derived_from")
        .and_then(serde_json::Value::as_array)
        .into_iter()
        .flatten()
        .filter_map(serde_json::Value::as_str)
        .any(|derived| {
            weak_paths
                .iter()
                .any(|path| path.ends_with(derived) || derived.ends_with(path))
        })
}

pub(super) fn analysis_observation_is_weak_indel_match(observation: &serde_json::Value) -> bool {
    observation
        .get("match_quality")
        .and_then(serde_json::Value::as_str)
        == Some("weak")
}

pub(super) fn render_analysis_logic(out: &mut String, analysis: &serde_json::Value) {
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
    out.push_str("<h4>Description</h4>");
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
