fn render_pgx_label_table(out: &mut String, findings: &[serde_json::Value]) {
    let headers = [
        "Variant",
        "Ref/Alt",
        "Genes",
        "Drug(s)",
        "Regulator",
        "Action",
        "Label",
        "Evidence",
    ];
    render_pgx_label_filters(out);
    render_table_start(out, "labels-table", &headers);
    for finding in findings {
        let evidence = finding.get("evidence");
        let url = evidence
            .and_then(|value| value.get("url"))
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let pgx_level = value_str(finding, "pgx_action_level");
        let _ = write!(
            out,
            "<tr data-pgx-level=\"{}\">",
            html_escape(&pgx_level_slug(pgx_level))
        );
        table_cell(out, value_str(finding, "variant"));
        class_cell(out, &matched_ref_alt(finding), "mono");
        table_cell(out, &join_string_array(finding.get("genes")));
        table_cell(out, &join_drugs(finding));
        table_cell(out, &join_string_array(finding.get("regulatory_sources")));
        pgx_level_cell(out, pgx_level);
        table_cell(out, value_str(finding, "label"));
        link_cell(out, url);
        out.push_str("</tr>");
    }
    render_table_end(out);
}

fn render_pgx_summary_table(out: &mut String, findings: &[serde_json::Value]) {
    let headers = [
        "Variant",
        "Ref/Alt",
        "Genotype",
        "Drug(s)",
        "Category",
        "Level",
        "Phenotype",
        "Effect",
        "Evidence",
    ];
    render_evidence_level_filters(out);
    render_table_start(out, "summaries-table", &headers);
    for finding in findings {
        let effect = finding
            .get("matched_effect")
            .unwrap_or(&serde_json::Value::Null);
        let evidence = finding.get("evidence");
        let url = evidence
            .and_then(|value| value.get("url"))
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let evidence_level = value_str(finding, "evidence_level");
        let _ = write!(
            out,
            "<tr data-level=\"{}\">",
            html_escape(&evidence_level_group(evidence_level))
        );
        table_cell(out, value_str(finding, "variant"));
        class_cell(out, &matched_ref_alt(finding), "mono");
        table_cell(out, value_str(effect, "label"));
        table_cell(out, &join_drugs(finding));
        table_cell(out, &join_string_array(finding.get("phenotype_categories")));
        evidence_level_cell(out, evidence_level);
        table_cell(out, &join_string_array(finding.get("phenotypes")));
        class_cell(out, value_str(effect, "text"), "effect");
        link_cell(out, url);
        out.push_str("</tr>");
    }
    render_table_end(out);
}

fn render_evidence_level_filters(out: &mut String) {
    out.push_str("<div class=\"level-filter\"><span class=\"muted\">Evidence:</span>");
    for (level, label) in [
        ("1", "Level 1"),
        ("1a", "Level 1A"),
        ("1b", "Level 1B"),
        ("2", "Level 2"),
        ("2a", "Level 2A"),
        ("2b", "Level 2B"),
        ("3", "Level 3"),
        ("4", "Level 4"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-level-filter=\"{level}\" checked onchange=\"applyTableFilters('summaries-table')\"> {label}</label>"
        );
    }
    out.push_str("<button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('level','summaries-table',true)\">Show all</button><button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('level','summaries-table',false)\">Hide all</button>");
    out.push_str("<a href=\"https://www.clinpgx.org/page/clinAnnLevels\" title=\"ClinPGx levels of evidence\">i</a></div>");
}

fn render_pgx_label_filters(out: &mut String) {
    out.push_str("<div class=\"level-filter\"><span class=\"muted\">PGx level:</span>");
    for (level, label) in [
        ("required", "Testing Required"),
        ("recommended", "Testing Recommended"),
        ("actionable", "Actionable PGx"),
        ("informative", "Informative PGx"),
        ("no-clinical", "No Clinical PGx"),
        ("criteria", "Criteria Not Met"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-pgx-level-filter=\"{level}\" checked onchange=\"applyTableFilters('labels-table')\"> {label}</label>"
        );
    }
    out.push_str("<button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('pgx-level','labels-table',true)\">Show all</button><button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('pgx-level','labels-table',false)\">Hide all</button>");
    out.push_str("<a href=\"https://www.clinpgx.org/page/drugLabelLegend#pgx-level\" title=\"ClinPGx drug label PGx levels\">i</a></div>");
}

fn matched_ref_alt(finding: &serde_json::Value) -> String {
    let Some(observation) = finding.get("matched_observation") else {
        return String::new();
    };
    let ref_allele = observation
        .get("ref")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let alt_allele = observation
        .get("alt")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if ref_allele.is_empty() && alt_allele.is_empty() {
        String::new()
    } else {
        let alt_display = alt_allele.replace(',', "/");
        format!("{ref_allele}->{alt_display}")
    }
}

fn evidence_level_group(level: &str) -> String {
    let normalized = level.trim().to_ascii_lowercase();
    if normalized.starts_with("1a") {
        "1a".to_owned()
    } else if normalized.starts_with("1b") {
        "1b".to_owned()
    } else if normalized.starts_with('1') {
        "1".to_owned()
    } else if normalized.starts_with("2a") {
        "2a".to_owned()
    } else if normalized.starts_with("2b") {
        "2b".to_owned()
    } else if normalized.starts_with('2') {
        "2".to_owned()
    } else if normalized.starts_with('3') {
        "3".to_owned()
    } else if normalized.starts_with('4') {
        "4".to_owned()
    } else {
        "unknown".to_owned()
    }
}

fn evidence_level_color_group(level: &str) -> String {
    level
        .chars()
        .find(char::is_ascii_digit)
        .map_or_else(|| "unknown".to_owned(), |ch| ch.to_string())
}

fn evidence_level_cell(out: &mut String, level: &str) {
    if level.is_empty() {
        out.push_str("<td></td>");
        return;
    }
    let group = evidence_level_color_group(level);
    let _ = write!(
        out,
        "<td><span class=\"level-badge level-{}\">{}</span></td>",
        html_escape(&group),
        html_escape(level)
    );
}

fn pgx_level_slug(level: &str) -> String {
    let normalized = level.to_ascii_lowercase();
    if normalized.contains("required") {
        "required".to_owned()
    } else if normalized.contains("recommended") {
        "recommended".to_owned()
    } else if normalized.contains("actionable") {
        "actionable".to_owned()
    } else if normalized.contains("informative") {
        "informative".to_owned()
    } else if normalized.contains("no clinical") {
        "no-clinical".to_owned()
    } else if normalized.contains("criteria") {
        "criteria".to_owned()
    } else {
        "unknown".to_owned()
    }
}

fn pgx_level_cell(out: &mut String, level: &str) {
    if level.is_empty() {
        out.push_str("<td></td>");
        return;
    }
    let slug = pgx_level_slug(level);
    let _ = write!(
        out,
        "<td><span class=\"pgx-badge pgx-{}\">{}</span></td>",
        html_escape(&slug),
        html_escape(level)
    );
}

