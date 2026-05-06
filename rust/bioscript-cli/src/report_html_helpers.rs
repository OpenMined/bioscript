fn render_table_start(out: &mut String, table_id: &str, headers: &[&str]) {
    let escaped_id = html_escape(table_id);
    let refs_control = if table_id == "observations-table" {
        "<label><input type=\"checkbox\" onchange=\"toggleRefs(this.checked)\"> Show refs</label>"
    } else {
        ""
    };
    let _ = write!(
        out,
        "<div class=\"table-tools\"><input type=\"search\" placeholder=\"Filter table\" data-filter-for=\"{escaped_id}\" oninput=\"applyTableFilters('{escaped_id}')\">{refs_control}</div><div class=\"table-wrap\"><table id=\"{escaped_id}\"><thead><tr>"
    );
    for (index, header) in headers.iter().enumerate() {
        let _ = write!(
            out,
            "<th onclick=\"sortTable('{}',{})\">{}<span class=\"sort-mark\"></span></th>",
            escaped_id,
            index,
            html_escape(header)
        );
    }
    out.push_str("</tr></thead><tbody>");
}

fn render_table_end(out: &mut String) {
    out.push_str("</tbody></table></div>");
}

fn table_cell(out: &mut String, value: &str) {
    class_cell(out, value, "");
}

fn class_cell(out: &mut String, value: &str, class_name: &str) {
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

fn link_cell(out: &mut String, url: &str) {
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

fn value_str<'a>(value: &'a serde_json::Value, key: &str) -> &'a str {
    value
        .get(key)
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
}

fn join_string_array(value: Option<&serde_json::Value>) -> String {
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

fn join_drugs(finding: &serde_json::Value) -> String {
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

fn html_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

