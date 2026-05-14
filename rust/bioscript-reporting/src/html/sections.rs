use super::helpers::*;
use std::fmt::Write as _;
pub(super) fn collect_report_analyses(reports: &[serde_json::Value]) -> Vec<serde_json::Value> {
    reports
        .iter()
        .filter_map(|report| report.get("analyses").and_then(serde_json::Value::as_array))
        .flat_map(|analyses| analyses.iter())
        .cloned()
        .collect()
}

pub(super) fn collect_report_findings(
    reports: &[serde_json::Value],
    schema: &str,
) -> Vec<serde_json::Value> {
    reports
        .iter()
        .filter_map(|report| report.get("findings").and_then(serde_json::Value::as_array))
        .flat_map(|findings| findings.iter())
        .filter(|finding| finding.get("schema").and_then(serde_json::Value::as_str) == Some(schema))
        .cloned()
        .collect()
}

pub(super) fn collect_report_participants(reports: &[serde_json::Value]) -> Vec<String> {
    let mut participants = Vec::new();
    for report in reports {
        let participant = value_str(report, "participant_id");
        if !participant.is_empty() && !participants.iter().any(|item| item == participant) {
            participants.push(participant.to_owned());
        }
    }
    participants
}

pub(super) fn render_report_manifest_header(out: &mut String, reports: &[serde_json::Value]) {
    let manifest = reports
        .first()
        .and_then(|report| report.get("manifest"))
        .unwrap_or(&serde_json::Value::Null);
    let title = value_str(manifest, "label");
    let fallback = value_str(manifest, "name");
    let title = if title.is_empty() { fallback } else { title };
    let title = if title.is_empty() {
        "BioScript Report"
    } else {
        title
    };
    let _ = write!(out, "<h1>{}</h1>", html_escape(title));
    out.push_str("<div class=\"logic-note\"><p><strong>Disclaimer:</strong> This is not medical or clinical advice, only for research purposes. Always consult a licensed professional to interpret medical information.</p><p>This report was generated offline on your system.</p><p>For more information see <a href=\"https://app.biovault.net\" target=\"_blank\" rel=\"noopener noreferrer\">https://app.biovault.net</a></p></div>");
}

pub(super) fn render_report_source_section(out: &mut String, reports: &[serde_json::Value]) {
    let manifest = reports
        .first()
        .and_then(|report| report.get("manifest"))
        .unwrap_or(&serde_json::Value::Null);
    out.push_str("<dl class=\"analysis-kv\">");
    report_manifest_kv(out, "Schema", value_str(manifest, "schema"));
    report_manifest_kv(out, "Version", value_str(manifest, "version"));
    report_manifest_kv(out, "Name", value_str(manifest, "name"));
    report_manifest_kv(out, "Label", value_str(manifest, "label"));
    report_manifest_kv(out, "Tags", &manifest_tags(manifest));
    report_manifest_kv(out, "Members", &manifest_member_summary(manifest));
    out.push_str("</dl>");
    render_manifest_members(out, manifest);
}

pub(super) fn report_manifest_kv(out: &mut String, key: &str, value: &str) {
    if value.is_empty() {
        return;
    }
    let _ = write!(
        out,
        "<dt>{}</dt><dd>{}</dd>",
        html_escape(key),
        html_escape(value)
    );
}

pub(super) fn manifest_tags(manifest: &serde_json::Value) -> String {
    manifest
        .get("tags")
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_json::Value::as_str)
                .collect::<Vec<_>>()
                .join(", ")
        })
        .unwrap_or_default()
}

pub(super) fn manifest_member_summary(manifest: &serde_json::Value) -> String {
    let Some(members) = manifest
        .get("members")
        .and_then(serde_json::Value::as_array)
    else {
        return String::new();
    };
    let preview = members
        .iter()
        .take(5)
        .filter_map(|member| {
            let kind = value_str(member, "kind");
            let path = value_str(member, "path");
            if kind.is_empty() && path.is_empty() {
                None
            } else if path.is_empty() {
                Some(kind.to_owned())
            } else {
                Some(format!("{kind}: {path}"))
            }
        })
        .collect::<Vec<_>>();
    let remaining = members.len().saturating_sub(preview.len());
    if remaining == 0 {
        preview.join("; ")
    } else {
        format!("{}; +{} more", preview.join("; "), remaining)
    }
}

pub(super) fn render_manifest_members(out: &mut String, manifest: &serde_json::Value) {
    let Some(members) = manifest
        .get("members")
        .and_then(serde_json::Value::as_array)
    else {
        return;
    };
    if members.is_empty() {
        return;
    }
    out.push_str("<details><summary>Show panel members</summary><div class=\"table-wrap\"><table id=\"manifest-members-table\"><thead><tr><th>Kind</th><th>Path</th><th>Version</th></tr></thead><tbody>");
    for member in members {
        out.push_str("<tr>");
        table_cell(out, value_str(member, "kind"));
        table_cell(out, value_str(member, "path"));
        table_cell(out, value_str(member, "version"));
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></div></details>");
}

pub(super) fn render_participant_filter(out: &mut String, participants: &[String]) {
    if participants.len() <= 1 {
        return;
    }
    out.push_str("<div class=\"participant-filter\"><label class=\"muted\" for=\"participant-filter\">Participant</label><select id=\"participant-filter\" onchange=\"setParticipant(this.value)\"><option value=\"\">All participants</option>");
    for participant in participants {
        let _ = write!(
            out,
            "<option value=\"{}\">{}</option>",
            html_escape(participant),
            html_escape(participant)
        );
    }
    out.push_str("</select></div>");
}

pub(super) fn render_input_debug(
    out: &mut String,
    reports: &[serde_json::Value],
    show_participant_id: bool,
) {
    if reports.is_empty() {
        out.push_str("<p class=\"muted\">No input metadata.</p>");
        return;
    }
    if reports.len() == 1 {
        render_input_debug_key_values(out, &reports[0]);
        return;
    }
    out.push_str("<div class=\"table-wrap\"><table id=\"input-debug-table\"><thead><tr>");
    let mut headers = Vec::new();
    if show_participant_id {
        headers.push("Participant");
    }
    headers.extend([
        "File",
        "Format",
        "Source",
        "Assembly",
        "Inferred Sex",
        "VCF Ref Imputation",
        "Evidence",
    ]);
    for (idx, header) in headers.iter().enumerate() {
        let _ = write!(
            out,
            "<th onclick=\"sortTable('input-debug-table',{})\">{}<span class=\"sort-mark\"></span></th>",
            idx,
            html_escape(header)
        );
    }
    out.push_str("</tr></thead><tbody>");
    for report in reports {
        let participant = value_str(report, "participant_id");
        let input = report.get("input").unwrap_or(&serde_json::Value::Null);
        let debug = input.get("debug").unwrap_or(&serde_json::Value::Null);
        let source = debug.get("source").unwrap_or(&serde_json::Value::Null);
        let sex = debug
            .get("inferred_sex")
            .unwrap_or(&serde_json::Value::Null);
        let _ = write!(
            out,
            "<tr data-participant=\"{}\">",
            html_escape(participant)
        );
        if show_participant_id {
            table_cell(out, participant);
        }
        table_cell(out, value_str(input, "file_name"));
        table_cell(
            out,
            &compact_join(&[
                value_str(debug, "format"),
                value_str(debug, "format_confidence"),
            ]),
        );
        table_cell(
            out,
            &compact_join(&[
                value_str(source, "vendor"),
                value_str(source, "platform_version"),
            ]),
        );
        table_cell(out, value_str(debug, "assembly"));
        table_cell(
            out,
            &compact_join(&[
                value_str(sex, "sex"),
                value_str(sex, "confidence"),
                value_str(sex, "method"),
            ]),
        );
        table_cell(out, input_debug_vcf_imputation(debug));
        table_cell(out, &input_debug_evidence(debug));
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></div>");
}

pub(super) fn render_input_debug_key_values(out: &mut String, report: &serde_json::Value) {
    let input = report.get("input").unwrap_or(&serde_json::Value::Null);
    let debug = input.get("debug").unwrap_or(&serde_json::Value::Null);
    let source = debug.get("source").unwrap_or(&serde_json::Value::Null);
    let sex = debug
        .get("inferred_sex")
        .unwrap_or(&serde_json::Value::Null);
    out.push_str("<dl class=\"analysis-kv\">");
    input_debug_kv(out, "File", value_str(input, "file_name"));
    input_debug_kv(
        out,
        "Format",
        &compact_join(&[
            value_str(debug, "format"),
            value_str(debug, "format_confidence"),
        ]),
    );
    input_debug_kv(
        out,
        "Source",
        &compact_join(&[
            value_str(source, "vendor"),
            value_str(source, "platform_version"),
        ]),
    );
    input_debug_kv(out, "Assembly", value_str(debug, "assembly"));
    input_debug_kv(
        out,
        "Inferred sex",
        &compact_join(&[
            value_str(sex, "sex"),
            value_str(sex, "confidence"),
            value_str(sex, "method"),
        ]),
    );
    input_debug_kv(out, "VCF ref imputation", input_debug_vcf_imputation(debug));
    input_debug_kv(out, "Evidence", &input_debug_evidence(debug));
    out.push_str("</dl>");
}

pub(super) fn input_debug_kv(out: &mut String, key: &str, value: &str) {
    let _ = write!(
        out,
        "<dt>{}</dt><dd>{}</dd>",
        html_escape(key),
        html_escape(value)
    );
}

pub(super) fn compact_join(values: &[&str]) -> String {
    values
        .iter()
        .filter(|value| !value.is_empty())
        .copied()
        .collect::<Vec<_>>()
        .join(" / ")
}

pub(super) fn input_debug_vcf_imputation(debug: &serde_json::Value) -> &'static str {
    if debug
        .get("vcf_missing_reference_imputation")
        .and_then(serde_json::Value::as_bool)
        == Some(true)
    {
        "used"
    } else {
        ""
    }
}

pub(super) fn input_debug_evidence(debug: &serde_json::Value) -> String {
    let mut evidence = Vec::new();
    collect_string_array(debug.get("evidence"), &mut evidence);
    collect_string_array(debug.get("warnings"), &mut evidence);
    if let Some(source) = debug.get("source") {
        collect_string_array(source.get("evidence"), &mut evidence);
    }
    if let Some(sex) = debug.get("inferred_sex") {
        collect_string_array(sex.get("evidence"), &mut evidence);
    }
    evidence.join(" | ")
}

pub(super) fn collect_string_array(value: Option<&serde_json::Value>, out: &mut Vec<String>) {
    if let Some(items) = value.and_then(serde_json::Value::as_array) {
        for item in items.iter().filter_map(serde_json::Value::as_str) {
            if !item.is_empty() && !out.iter().any(|existing| existing == item) {
                out.push(item.to_owned());
            }
        }
    }
}
