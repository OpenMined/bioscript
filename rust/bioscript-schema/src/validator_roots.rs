fn validate_variant_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(
        root,
        "bioscript:variant:1.0",
        Some("bioscript:variant"),
        issues,
    );
    validate_optional_strings(root, &["name", "label", "gene", "summary"], issues);
    validate_tags(root, issues);
    validate_identifiers(root, issues);
    validate_coordinates(root, issues);
    validate_alleles(root, issues);
    validate_findings(root, issues);
    validate_provenance(root, issues);

    let has_identifiers = value_at(root, &["identifiers"])
        .and_then(Value::as_mapping)
        .is_some_and(|mapping| !mapping.is_empty());
    let has_coordinates = ["grch37", "grch38"]
        .iter()
        .any(|assembly| value_at(root, &["coordinates", assembly]).is_some());
    if !has_identifiers && !has_coordinates {
        issues.push(Issue {
            severity: Severity::Error,
            path: "identifiers/coordinates".to_owned(),
            message: "expected at least one identifier block or one coordinate block".to_owned(),
        });
    }
}

fn validate_panel_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(root, "bioscript:panel:1.0", None, issues);
    validate_optional_strings(root, &["name", "label", "summary"], issues);
    validate_tags(root, issues);
    validate_permissions(root, issues);
    validate_downloads(root, issues);
    validate_panel_members(root, &["variant", "variant-catalogue", "assay"], issues);
    validate_panel_interpretations(root, issues);
    validate_findings(root, issues);
}

fn validate_assay_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(root, "bioscript:assay:1.0", None, issues);
    validate_optional_strings(root, &["name", "label", "summary"], issues);
    validate_tags(root, issues);
    validate_panel_members(root, &["variant", "variant-catalogue"], issues);
    validate_panel_interpretations(root, issues);
    validate_findings(root, issues);
}

fn validate_pgx_findings_root(root: &Value, issues: &mut Vec<Issue>) {
    require_const(root, &["schema"], "bioscript:pgx-findings:1.0", issues);
    require_const(root, &["version"], "1.0", issues);
    validate_optional_strings(root, &["variant", "gene", "rsid", "variant_pa_id"], issues);
    if value_at(root, &["variant"]).is_none() && value_at(root, &["rsid"]).is_none() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "variant/rsid".to_owned(),
            message: "expected at least one variant identifier".to_owned(),
        });
    }
    match value_at(root, &["findings"]) {
        Some(Value::Sequence(_)) => {}
        Some(_) => issues.push(Issue {
            severity: Severity::Error,
            path: "findings".to_owned(),
            message: "expected a sequence of findings".to_owned(),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: "findings".to_owned(),
            message: "missing required field".to_owned(),
        }),
    }
    validate_findings(root, issues);
}

fn validate_schema_and_identity(
    root: &Value,
    canonical_schema: &str,
    legacy_schema: Option<&str>,
    issues: &mut Vec<Issue>,
) {
    let schema = scalar_at(root, &["schema"]);
    let valid_schema = schema
        .as_deref()
        .is_some_and(|value| value == canonical_schema || legacy_schema == Some(value));
    if !valid_schema {
        issues.push(Issue {
            severity: Severity::Error,
            path: "schema".to_owned(),
            message: format!("expected schema to be '{canonical_schema}'"),
        });
    }
    if let Some(legacy_schema) = legacy_schema
        && matches!(schema.as_deref(), Some(value) if value == legacy_schema)
    {
        issues.push(Issue {
            severity: Severity::Warning,
            path: "schema".to_owned(),
            message: format!("legacy schema value '{legacy_schema}'; prefer '{canonical_schema}'"),
        });
    }
    require_const(root, &["version"], "1.0", issues);
    match scalar_at(root, &["name"]) {
        Some(name) if !name.trim().is_empty() => {}
        Some(_) => issues.push(Issue {
            severity: Severity::Error,
            path: "name".to_owned(),
            message: "empty string".to_owned(),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: "name".to_owned(),
            message: "missing required field".to_owned(),
        }),
    }
    if value_at(root, &["variant_id"]).is_some() {
        issues.push(Issue {
            severity: Severity::Warning,
            path: "variant_id".to_owned(),
            message: "variant_id is legacy; prefer name".to_owned(),
        });
    }
}

fn validate_optional_strings(root: &Value, fields: &[&str], issues: &mut Vec<Issue>) {
    for field in fields {
        if let Some(value) = value_at(root, &[*field]) {
            match value.as_str() {
                Some(text) if !text.trim().is_empty() => {}
                Some(_) => issues.push(Issue {
                    severity: Severity::Warning,
                    path: (*field).to_owned(),
                    message: "empty string".to_owned(),
                }),
                None => issues.push(Issue {
                    severity: Severity::Error,
                    path: (*field).to_owned(),
                    message: "expected string".to_owned(),
                }),
            }
        }
    }
}

fn validate_tags(root: &Value, issues: &mut Vec<Issue>) {
    let Some(value) = value_at(root, &["tags"]) else {
        return;
    };
    let Some(items) = value.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: "tags".to_owned(),
            message: "expected a sequence of strings".to_owned(),
        });
        return;
    };

    for (idx, item) in items.iter().enumerate() {
        let Some(tag) = item.as_str() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("tags[{idx}]"),
                message: "expected string".to_owned(),
            });
            continue;
        };
        if tag.trim().is_empty() {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("tags[{idx}]"),
                message: "empty tag string".to_owned(),
            });
        }
    }
}

fn validate_identifiers(root: &Value, issues: &mut Vec<Issue>) {
    for field in ["rsids", "aliases"] {
        let Some(values) = value_at(root, &["identifiers", field]) else {
            continue;
        };
        let Some(items) = values.as_sequence() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("identifiers.{field}"),
                message: "expected a sequence of strings".to_owned(),
            });
            continue;
        };
        let mut seen = BTreeSet::new();
        for (idx, item) in items.iter().enumerate() {
            let Some(value) = item.as_str() else {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: "expected string".to_owned(),
                });
                continue;
            };
            if !is_rsid(value) {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: format!("expected rsid like rs123, found '{value}'"),
                });
            }
            if !seen.insert(value.to_owned()) {
                issues.push(Issue {
                    severity: Severity::Warning,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: format!("duplicate identifier '{value}'"),
                });
            }
        }
    }
}

fn validate_coordinates(root: &Value, issues: &mut Vec<Issue>) {
    for assembly in ["grch37", "grch38"] {
        let Some(coord) = mapping_at(root, &["coordinates", assembly]) else {
            continue;
        };

        let Some(chrom) = coord
            .get(Value::String("chrom".to_owned()))
            .and_then(Value::as_str)
        else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.chrom"),
                message: "missing chrom".to_owned(),
            });
            continue;
        };
        if !is_allowed_chromosome(chrom) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.chrom"),
                message: format!("invalid chromosome '{chrom}'; expected 1-22, X, Y, or MT"),
            });
        }

        let has_pos = coord.contains_key(Value::String("pos".to_owned()));
        let has_start = coord.contains_key(Value::String("start".to_owned()));
        let has_end = coord.contains_key(Value::String("end".to_owned()));
        if has_pos && (has_start || has_end) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}"),
                message: "use either pos or start/end, not both".to_owned(),
            });
            continue;
        }
        if !(has_pos || has_start && has_end) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}"),
                message: "expected either pos or start/end".to_owned(),
            });
            continue;
        }

        if has_pos {
            validate_coordinate_pos(coord, assembly, issues);
        } else {
            validate_coordinate_range(coord, assembly, issues);
        }
    }
}

fn validate_coordinate_pos(coord: &Mapping, assembly: &str, issues: &mut Vec<Issue>) {
    if let Some(pos) = i64_at_mapping(coord, "pos") {
        if pos < 1 {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.pos"),
                message: "expected integer >= 1".to_owned(),
            });
        }
    } else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.pos"),
            message: "expected integer".to_owned(),
        });
    }
}

fn validate_coordinate_range(coord: &Mapping, assembly: &str, issues: &mut Vec<Issue>) {
    let start = i64_at_mapping(coord, "start");
    let end = i64_at_mapping(coord, "end");
    match (start, end) {
        (Some(start), Some(end)) => validate_coordinate_range_values(start, end, assembly, issues),
        _ => issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}"),
            message: "expected integer start/end".to_owned(),
        }),
    }
}

fn validate_coordinate_range_values(start: i64, end: i64, assembly: &str, issues: &mut Vec<Issue>) {
    if start < 1 {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.start"),
            message: "expected integer >= 1".to_owned(),
        });
    }
    if end < 1 {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.end"),
            message: "expected integer >= 1".to_owned(),
        });
    }
    if end < start {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.end"),
            message: "expected end >= start".to_owned(),
        });
    }
    if start == end {
        issues.push(Issue {
            severity: Severity::Warning,
            path: format!("coordinates.{assembly}"),
            message: "single-position coordinate uses start/end; prefer pos".to_owned(),
        });
    }
}
