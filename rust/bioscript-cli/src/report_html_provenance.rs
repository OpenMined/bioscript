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
    for (domain, domain_links) in group_provenance_links_by_domain(links) {
        let _ = write!(
            out,
            "<li><details><summary><strong>{}</strong> <span class=\"muted\">({} links)</span></summary><ul>",
            html_escape(&domain),
            domain_links.len()
        );
        for (url, label) in domain_links {
            let display = if label.is_empty() { &url } else { &label };
            let _ = write!(
                out,
                "<li><a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">{}</a><div class=\"muted mono\">{}</div></li>",
                html_escape(&url),
                html_escape(display),
                html_escape(&url)
            );
        }
        out.push_str("</ul></details></li>");
    }
    out.push_str("</ul>");
}

fn group_provenance_links_by_domain(
    links: BTreeMap<String, String>,
) -> BTreeMap<String, BTreeMap<String, String>> {
    let mut grouped = BTreeMap::<String, BTreeMap<String, String>>::new();
    for (url, label) in links {
        grouped
            .entry(domain_from_url(&url).unwrap_or_else(|| "other".to_owned()))
            .or_default()
            .insert(url, label);
    }
    grouped
}

fn domain_from_url(url: &str) -> Option<String> {
    let without_scheme = url.split_once("://")?.1;
    let host = without_scheme.split(['/', '?', '#']).next()?.trim();
    if host.is_empty() {
        None
    } else {
        Some(host.to_ascii_lowercase())
    }
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

