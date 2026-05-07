fn render_pgx_table(
    out: &mut String,
    label_findings: &[serde_json::Value],
    summary_findings: &[serde_json::Value],
) {
    let mut findings = Vec::new();
    findings.extend(label_findings.iter());
    findings.extend(summary_findings.iter());
    if findings.is_empty() {
        out.push_str("<p class=\"muted\">No PGx findings.</p>");
        return;
    }
    render_pgx_filters(out);
    out.push_str("<div class=\"pgx-tabs\"><button id=\"pgx-tab-variant\" class=\"active\" type=\"button\" onclick=\"setPgxView('variant')\">Variant</button><button id=\"pgx-tab-drug\" type=\"button\" onclick=\"setPgxView('drug')\">Drug</button></div>");
    out.push_str("<div id=\"pgx-view-variant\">");
    let headers = [
        "Type",
        "RSID",
        "Gene",
        "Ref / Alt",
        "Genotype",
        "Drug(s)",
        "Level",
        "Category",
        "Phenotype",
        "Finding",
        "Source",
    ];
    render_table_start(out, "pgx-variant-table", &headers);
    for finding in &findings {
        render_pgx_row(out, finding, true);
    }
    render_table_end(out);
    out.push_str("</div><div id=\"pgx-view-drug\" hidden>");
    let mut drugs = pgx_drug_names(&findings);
    if drugs.is_empty() {
        drugs.push("unknown".to_owned());
    }
    let drug_headers = [
        "Type",
        "RSID",
        "Gene",
        "Ref / Alt",
        "Genotype",
        "Level",
        "Category",
        "Phenotype",
        "Finding",
        "Source",
    ];
    for (index, drug) in drugs.iter().enumerate() {
        let _ = write!(
            out,
            "<div class=\"pgx-drug-group\"><h3>{}</h3>",
            html_escape(drug)
        );
        let table_id = format!("pgx-drug-table-{index}");
        render_table_start(out, &table_id, &drug_headers);
        for finding in &findings {
            let names = finding_drug_names(finding);
            if (names.is_empty() && drug == "unknown") || names.iter().any(|name| name == drug) {
                render_pgx_row(out, finding, false);
            }
        }
        render_table_end(out);
        out.push_str("</div>");
    }
    out.push_str("</div>");
}

fn render_pgx_row(out: &mut String, finding: &serde_json::Value, show_drug: bool) {
    let source_type = pgx_source_type(finding);
    let level = pgx_level_value(finding);
    let level_slug = pgx_level_filter_slug(finding);
    let outcome = pgx_outcome_filter_slug(finding);
    let _ = write!(
        out,
        "<tr data-level=\"{}\" data-pgx-level=\"{}\" data-pgx-outcome=\"{}\" data-participant=\"{}\">",
        html_escape(&evidence_level_group(level)),
        html_escape(&level_slug),
        html_escape(outcome),
        html_escape(&finding_participant(finding))
    );
    table_cell(out, source_type);
    table_cell(out, &finding_rsid(finding));
    table_cell(out, &finding_gene(finding));
    class_cell(out, &matched_ref_alt(finding), "mono");
    pgx_genotype_cell(out, finding);
    if show_drug {
        table_cell(out, &join_drugs(finding));
    }
    pgx_any_level_cell(out, finding);
    table_cell(out, &pgx_category(finding));
    table_cell(out, &join_string_array(finding.get("phenotypes")));
    class_cell(out, &pgx_finding_text(finding), "effect");
    link_cell(out, pgx_evidence_url(finding));
    out.push_str("</tr>");
}

fn pgx_source_type(finding: &serde_json::Value) -> &str {
    match value_str(finding, "schema") {
        "bioscript:pgx-label:1.0" => "Drug Label",
        _ => "Summary",
    }
}

fn pgx_level_value(finding: &serde_json::Value) -> &str {
    if pgx_source_type(finding) == "Drug Label" {
        value_str(finding, "pgx_action_level")
    } else {
        value_str(finding, "evidence_level")
    }
}

fn pgx_level_filter_slug(finding: &serde_json::Value) -> String {
    if pgx_source_type(finding) == "Drug Label" {
        format!("drug-{}", pgx_level_slug(pgx_level_value(finding)))
    } else {
        format!("summary-{}", evidence_level_group(pgx_level_value(finding)))
    }
}

fn pgx_category(finding: &serde_json::Value) -> String {
    if pgx_source_type(finding) == "Drug Label" {
        let actions = join_string_array(finding.get("prescribing_actions"));
        let sources = join_string_array(finding.get("regulatory_sources"));
        if actions.is_empty() {
            sources
        } else if sources.is_empty() {
            actions
        } else {
            format!("{sources}; {actions}")
        }
    } else {
        join_string_array(finding.get("phenotype_categories"))
    }
}

fn pgx_finding_text(finding: &serde_json::Value) -> String {
    if pgx_source_type(finding) == "Drug Label" {
        for key in ["prescribing_information", "summary", "notes", "label"] {
            let value = value_str(finding, key);
            if !value.is_empty() {
                return value.to_owned();
            }
        }
        return String::new();
    }
    let effect = finding
        .get("matched_effect")
        .unwrap_or(&serde_json::Value::Null);
    let text = value_str(effect, "text");
    if !text.is_empty() {
        return text.to_owned();
    }
    value_str(finding, "notes").to_owned()
}

fn pgx_evidence_url(finding: &serde_json::Value) -> &str {
    finding
        .get("evidence")
        .and_then(|value| value.get("url"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
}

fn pgx_drug_names(findings: &[&serde_json::Value]) -> Vec<String> {
    let mut drugs = Vec::new();
    for finding in findings {
        for drug in finding_drug_names(finding) {
            if !drugs.contains(&drug) {
                drugs.push(drug);
            }
        }
    }
    drugs.sort();
    drugs
}

fn finding_drug_names(finding: &serde_json::Value) -> Vec<String> {
    finding
        .get("drugs")
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(|drug| drug.get("name").and_then(serde_json::Value::as_str))
                .map(ToOwned::to_owned)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

fn render_pgx_filters(out: &mut String) {
    out.push_str("<div class=\"level-filter\"><div class=\"filter-scale\"><div class=\"filter-scale-title\">Drug Label PGx Level <a href=\"https://www.clinpgx.org/page/drugLabelLegend#pgx-level\" target=\"_blank\" rel=\"noopener noreferrer\" title=\"ClinPGx drug label PGx levels\">i</a></div><div class=\"filter-options\">");
    for (level, label) in [
        ("required", "Testing Required"),
        ("recommended", "Testing Recommended"),
        ("actionable", "Actionable PGx"),
        ("informative", "Informative PGx"),
        ("no-clinical", "No Clinical PGx"),
        ("criteria", "Criteria Not Met"),
        ("unknown", "Unknown"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-pgx-any-level-filter=\"drug-{level}\" checked onchange=\"applyPgxFilters()\"> <span class=\"pgx-badge pgx-{level}\">{label}</span></label>"
        );
    }
    out.push_str("</div></div><div class=\"filter-scale\"><div class=\"filter-scale-title\">Summary Evidence Level <a href=\"https://www.clinpgx.org/page/clinAnnLevels\" target=\"_blank\" rel=\"noopener noreferrer\" title=\"ClinPGx summary annotation evidence levels\">i</a></div><div class=\"filter-options\">");
    for (level, label, class_level) in [
        ("1", "Level 1"),
        ("1a", "Level 1A"),
        ("1b", "Level 1B"),
        ("2", "Level 2"),
        ("2a", "Level 2A"),
        ("2b", "Level 2B"),
        ("3", "Level 3"),
        ("4", "Level 4"),
        ("unknown", "Unknown"),
    ]
    .map(|(level, label)| (level, label, evidence_level_color_group(level)))
    {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-pgx-any-level-filter=\"summary-{level}\" checked onchange=\"applyPgxFilters()\"> <span class=\"level-badge level-{}\">{label}</span></label>",
            html_escape(&class_level)
        );
    }
    out.push_str("</div></div><div class=\"filter-scale\"><div class=\"filter-scale-title\">Result</div><div class=\"filter-options\">");
    for (outcome, label, class_name) in [
        ("variant", "Variant", "analysis-badge-variant"),
        ("reference", "Normal", "analysis-badge-normal"),
        ("missing", "Missing", "analysis-badge-unknown"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-pgx-outcome-filter=\"{outcome}\" checked onchange=\"applyPgxFilters()\"> <span class=\"analysis-badge {class_name}\">{label}</span></label>"
        );
    }
    out.push_str("</div></div><div class=\"filter-actions\"><button class=\"filter-action\" type=\"button\" onclick=\"setPgxFilterGroup(true)\">Show all</button><button class=\"filter-action\" type=\"button\" onclick=\"setPgxFilterGroup(false)\">Hide all</button></div></div>");
}

fn pgx_any_level_cell(out: &mut String, finding: &serde_json::Value) {
    let level = pgx_level_value(finding);
    if pgx_source_type(finding) == "Drug Label" {
        pgx_level_cell(out, level);
    } else {
        evidence_level_cell(out, level);
    }
}

fn pgx_genotype_cell(out: &mut String, finding: &serde_json::Value) {
    let value = finding
        .get("matched_observation")
        .and_then(|observation| observation.get("genotype_display"))
        .and_then(serde_json::Value::as_str)
        .or_else(|| {
            finding
                .get("matched_effect")
                .and_then(|effect| effect.get("label"))
                .and_then(serde_json::Value::as_str)
        })
        .unwrap_or_default();
    let alt = finding
        .get("matched_observation")
        .and_then(|observation| observation.get("alt"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if alt.is_empty() {
        class_cell(out, value, "mono");
    } else {
        let _ = write!(
            out,
            "<td class=\"mono genotype-hit\">{}</td>",
            highlight_allele(value, alt)
        );
    }
}

fn finding_participant(finding: &serde_json::Value) -> String {
    finding
        .get("matched_observation")
        .or_else(|| finding.get("matched_analysis"))
        .and_then(|value| value.get("participant_id"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
        .to_owned()
}

fn pgx_outcome_filter_slug(finding: &serde_json::Value) -> &'static str {
    match finding
        .get("matched_observation")
        .and_then(|observation| observation.get("outcome"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "variant" | "observed_alt" | "unknown_alt" => "variant",
        "reference" => "reference",
        _ => "missing",
    }
}

fn finding_rsid(finding: &serde_json::Value) -> String {
    finding
        .get("matched_observation")
        .and_then(|observation| observation.get("rsid"))
        .and_then(serde_json::Value::as_str)
        .or_else(|| finding.get("rsid").and_then(serde_json::Value::as_str))
        .unwrap_or_default()
        .to_owned()
}

fn finding_gene(finding: &serde_json::Value) -> String {
    finding
        .get("matched_observation")
        .and_then(|observation| observation.get("gene"))
        .and_then(serde_json::Value::as_str)
        .map(ToOwned::to_owned)
        .or_else(|| {
            let genes = join_string_array(finding.get("genes"));
            if genes.is_empty() {
                None
            } else {
                Some(genes)
            }
        })
        .unwrap_or_default()
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
    let display = if level.is_empty() { "Unknown" } else { level };
    let group = evidence_level_color_group(display);
    let _ = write!(
        out,
        "<td data-sort=\"{}\"><span class=\"level-badge level-{}\">{}</span></td>",
        evidence_level_sort_rank(display),
        html_escape(&group),
        html_escape(display)
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
    let display = if level.is_empty() { "Unknown" } else { level };
    let slug = pgx_level_slug(display);
    let _ = write!(
        out,
        "<td data-sort=\"{}\"><span class=\"pgx-badge pgx-{}\">{}</span></td>",
        pgx_level_sort_rank(display),
        html_escape(&slug),
        html_escape(display)
    );
}

fn pgx_level_sort_rank(level: &str) -> u8 {
    match pgx_level_slug(level).as_str() {
        "required" => 1,
        "recommended" => 2,
        "actionable" => 3,
        "informative" => 4,
        "no-clinical" => 5,
        "criteria" => 6,
        _ => 7,
    }
}

fn evidence_level_sort_rank(level: &str) -> u8 {
    match evidence_level_group(level).as_str() {
        "1a" => 11,
        "1b" => 12,
        "1" => 13,
        "2a" => 21,
        "2b" => 22,
        "2" => 23,
        "3" => 30,
        "4" => 40,
        _ => 99,
    }
}
