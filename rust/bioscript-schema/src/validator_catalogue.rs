fn validate_variant_catalogue_root(root: &Value, issues: &mut Vec<Issue>) {
    require_const(root, &["schema"], "bioscript:variant-catalogue:1.0", issues);
    require_const(root, &["version"], "1.0", issues);
    validate_optional_strings(root, &["name", "label", "summary"], issues);
    require_path(root, &["name"], issues);
    validate_tags(root, issues);
    validate_catalogue_table(root, "variants", true, issues);
    validate_catalogue_table(root, "findings", false, issues);
    validate_catalogue_provenance(root, issues);
}

fn validate_catalogue_table(
    root: &Value,
    field: &str,
    required: bool,
    issues: &mut Vec<Issue>,
) {
    let Some(value) = value_at(root, &[field]) else {
        if required {
            issues.push(Issue {
                severity: Severity::Error,
                path: field.to_owned(),
                message: "missing required table declaration".to_owned(),
            });
        }
        return;
    };
    let Some(mapping) = value.as_mapping() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: field.to_owned(),
            message: "expected mapping".to_owned(),
        });
        return;
    };
    validate_required_non_empty_mapping_string(mapping, field, "source", issues);
    validate_optional_catalogue_mapping_string(mapping, field, "format", issues);
    validate_optional_catalogue_mapping_string(mapping, field, "key", issues);

    if let Some(source) = mapping
        .get(Value::String("source".to_owned()))
        .and_then(Value::as_str)
        && !std::path::Path::new(source)
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("tsv"))
    {
        issues.push(Issue {
            severity: Severity::Warning,
            path: format!("{field}.source"),
            message: "expected a .tsv source for auditable tabular catalogue data".to_owned(),
        });
    }
    if let Some(format) = mapping
        .get(Value::String("format".to_owned()))
        .and_then(Value::as_str)
        && format != "tsv"
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{field}.format"),
            message: "expected 'tsv'".to_owned(),
        });
    }
}

fn validate_catalogue_provenance(root: &Value, issues: &mut Vec<Issue>) {
    let Some(sources) = value_at(root, &["provenance", "sources"]).and_then(Value::as_sequence)
    else {
        return;
    };
    let mut ids = BTreeSet::new();
    for (idx, source) in sources.iter().enumerate() {
        let Some(mapping) = source.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("provenance.sources[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        for field in ["id", "kind", "label"] {
            validate_required_non_empty_mapping_string(
                mapping,
                &format!("provenance.sources[{idx}]"),
                field,
                issues,
            );
        }
        if let Some(id) = mapping
            .get(Value::String("id".to_owned()))
            .and_then(Value::as_str)
            && !ids.insert(id.to_owned())
        {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("provenance.sources[{idx}].id"),
                message: format!("duplicate provenance source id '{id}'"),
            });
        }

        let url = mapping
            .get(Value::String("url".to_owned()))
            .and_then(Value::as_str);
        let url_template = mapping
            .get(Value::String("url_template".to_owned()))
            .and_then(Value::as_str);
        match (url, url_template) {
            (Some(value), _) => validate_url_string(
                value,
                &format!("provenance.sources[{idx}].url"),
                false,
                issues,
            ),
            (None, Some(value)) => validate_url_template(
                value,
                &format!("provenance.sources[{idx}].url_template"),
                issues,
            ),
            (None, None) => issues.push(Issue {
                severity: Severity::Error,
                path: format!("provenance.sources[{idx}]"),
                message: "expected url or url_template".to_owned(),
            }),
        }
    }
}

fn validate_url_template(value: &str, path: &str, issues: &mut Vec<Issue>) {
    let sample = value
        .replace("{rsid}", "rs1")
        .replace("{genomic_hgvs}", "NG_000001.1%3Ag.1A%3ET")
        .replace("{variant_id}", "rs1");
    validate_url_string(&sample, path, false, issues);
}

fn validate_required_non_empty_mapping_string(
    mapping: &Mapping,
    parent: &str,
    field: &str,
    issues: &mut Vec<Issue>,
) {
    match mapping
        .get(Value::String(field.to_owned()))
        .and_then(Value::as_str)
    {
        Some(text) if !text.trim().is_empty() => {}
        Some(_) => issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.{field}"),
            message: "empty string".to_owned(),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.{field}"),
            message: "missing required field".to_owned(),
        }),
    }
}

fn validate_optional_catalogue_mapping_string(
    mapping: &Mapping,
    parent: &str,
    field: &str,
    issues: &mut Vec<Issue>,
) {
    let Some(value) = mapping.get(Value::String(field.to_owned())) else {
        return;
    };
    if !matches!(value, Value::String(text) if !text.trim().is_empty()) {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.{field}"),
            message: "expected non-empty string".to_owned(),
        });
    }
}
