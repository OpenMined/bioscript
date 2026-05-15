use std::fmt::Write as _;
pub(super) fn render_table_start(out: &mut String, table_id: &str, headers: &[&str]) {
    let escaped_id = html_escape(table_id);
    let _ = write!(
        out,
        "<div class=\"table-tools\"><input type=\"search\" placeholder=\"Filter table\" data-filter-for=\"{escaped_id}\" oninput=\"scheduleTableFilter('{escaped_id}')\"></div><div class=\"table-wrap\"><table id=\"{escaped_id}\"><thead><tr>"
    );
    for (index, header) in headers.iter().enumerate() {
        let _ = write!(
            out,
            "<th class=\"{}\" onclick=\"sortTable('{}',{})\">{}<span class=\"sort-mark\"></span></th>",
            table_column_class(header),
            escaped_id,
            index,
            html_escape(&table_header_label(header))
        );
    }
    out.push_str("</tr></thead><tbody>");
}

pub(super) fn table_column_class(header: &str) -> &'static str {
    if is_debug_column(header) {
        "debug-col"
    } else {
        ""
    }
}

pub(super) fn is_debug_column(header: &str) -> bool {
    matches!(
        header,
        "allele_balance"
            | "genotype_quality"
            | "evidence_type"
            | "evidence_raw"
            | "match_quality"
            | "match_notes"
            | "assay_id"
            | "assay_version"
            | "variant_key"
            | "match_status"
            | "coverage_status"
            | "call_status"
    )
}

pub(super) fn table_header_label(header: &str) -> String {
    match header {
        "participant_id" => "Participant ID".to_owned(),
        "rsid" => "RSID".to_owned(),
        "gene" => "Gene".to_owned(),
        "ref" => "Ref".to_owned(),
        "alt" => "Alt".to_owned(),
        "ref_alt" | "Ref/Alt" => "Ref / Alt".to_owned(),
        "genotype_display" => "Genotype".to_owned(),
        "genotype" => "GT".to_owned(),
        "zygosity" => "Zygosity".to_owned(),
        "outcome" => "Outcome".to_owned(),
        "match_status" => "Match Status".to_owned(),
        "coverage_status" => "Coverage Status".to_owned(),
        "call_status" => "Call Status".to_owned(),
        "assembly" => "Assembly".to_owned(),
        "chrom" => "Chrom".to_owned(),
        "pos_start" => "Start".to_owned(),
        "pos_end" => "End".to_owned(),
        "kind" => "Kind".to_owned(),
        "ref_count" => "Ref Count".to_owned(),
        "alt_count" => "Alt Count".to_owned(),
        "depth" => "Depth".to_owned(),
        "genotype_quality" => "Genotype Quality".to_owned(),
        "allele_balance" => "Allele Balance".to_owned(),
        "evidence_type" => "Evidence Type".to_owned(),
        "source" => "Source".to_owned(),
        "evidence_raw" => "Evidence Raw".to_owned(),
        "match_quality" => "Match Quality".to_owned(),
        "match_notes" => "Match Notes".to_owned(),
        "facets" => "Facets".to_owned(),
        "assay_id" => "Assay ID".to_owned(),
        "assay_version" => "Assay Version".to_owned(),
        "variant_key" => "Variant Key".to_owned(),
        other => other.to_owned(),
    }
}

pub(super) fn render_table_end(out: &mut String) {
    out.push_str("</tbody></table></div>");
}

pub(super) fn table_cell(out: &mut String, value: &str) {
    class_cell(out, value, "");
}

pub(super) fn class_cell(out: &mut String, value: &str, class_name: &str) {
    if class_name.is_empty() {
        let _ = write!(out, "<td>{}</td>", html_escape(value));
    } else {
        let _ = write!(
            out,
            "<td class=\"{}\">{}</td>",
            class_name,
            html_escape(value)
        );
    }
}

pub(super) fn link_cell(out: &mut String, url: &str) {
    if url.is_empty() {
        out.push_str("<td></td>");
    } else {
        let escaped = html_escape(url);
        let _ = write!(
            out,
            "<td><a href=\"{escaped}\" target=\"_blank\" rel=\"noopener noreferrer\">source</a></td>"
        );
    }
}

pub(super) fn value_str<'a>(value: &'a serde_json::Value, key: &str) -> &'a str {
    value
        .get(key)
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
}

pub(super) fn join_string_array(value: Option<&serde_json::Value>) -> String {
    value
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

pub(super) fn join_drugs(finding: &serde_json::Value) -> String {
    finding
        .get("drugs")
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(|drug| drug.get("name").and_then(serde_json::Value::as_str))
                .collect::<Vec<_>>()
                .join(", ")
        })
        .unwrap_or_default()
}

pub(super) fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}

pub(super) fn html_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
