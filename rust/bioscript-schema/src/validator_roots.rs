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
    validate_panel_members(root, &["variant", "assay"], issues);
    validate_panel_interpretations(root, issues);
    validate_findings(root, issues);
}

fn validate_assay_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(root, "bioscript:assay:1.0", None, issues);
    validate_optional_strings(root, &["name", "label", "summary"], issues);
    validate_tags(root, issues);
    validate_panel_members(root, &["variant"], issues);
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

#[cfg(test)]
mod root_validator_tests {
    use super::*;

    fn yaml(text: &str) -> Value {
        serde_yaml::from_str(text).unwrap()
    }

    fn messages(issues: &[Issue]) -> Vec<String> {
        issues
            .iter()
            .map(|issue| format!("{}:{}", issue.path, issue.message))
            .collect()
    }

    #[test]
    fn variant_root_reports_identity_optional_tag_identifier_and_coordinate_edges() {
        let root = yaml(
            r#"
schema: bioscript:variant
version: "1.0"
name: ""
label: ""
gene: [not, string]
tags: ["", 7]
variant_id: legacy
identifiers:
  rsids: [bad, rs1, rs1, 9]
  aliases: bad-shape
coordinates:
  grch37:
    chrom: 99
    pos: 0
  grch38:
    chrom: "1"
    pos: 1
    start: 1
    end: 1
"#,
        );
        let mut issues = Vec::new();
        validate_variant_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");

        assert!(rendered.contains("schema:legacy schema value"));
        assert!(rendered.contains("name:empty string"));
        assert!(rendered.contains("label:empty string"));
        assert!(rendered.contains("gene:expected string"));
        assert!(rendered.contains("tags[0]:empty tag string"));
        assert!(rendered.contains("tags[1]:expected string"));
        assert!(rendered.contains("variant_id:variant_id is legacy"));
        assert!(rendered.contains("identifiers.rsids[0]:expected rsid"));
        assert!(rendered.contains("identifiers.rsids[2]:duplicate identifier"));
        assert!(rendered.contains("identifiers.rsids[3]:expected string"));
        assert!(rendered.contains("identifiers.aliases:expected a sequence"));
        assert!(rendered.contains("coordinates.grch37.chrom:missing chrom"));
        assert!(rendered.contains("coordinates.grch38:use either pos or start/end"));
    }

    #[test]
    fn coordinate_range_validation_reports_missing_non_integer_bounds_and_single_position() {
        let root = yaml(
            r#"
schema: bioscript:variant:1.0
version: "1.0"
name: coordinate-errors
coordinates:
  grch37:
    chrom: Z
    start: zero
    end: 1
  grch38:
    chrom: X
    start: 5
    end: 5
alleles:
  kind: snv
  ref: A
  alts: [G]
"#,
        );
        let mut issues = Vec::new();
        validate_variant_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");

        assert!(rendered.contains("coordinates.grch37.chrom:invalid chromosome"));
        assert!(rendered.contains("coordinates.grch37:expected integer start/end"));
        assert!(rendered.contains("coordinates.grch38:single-position coordinate"));

        let root = yaml(
            r#"
schema: bioscript:variant:1.0
version: "1.0"
name: bad-range
coordinates:
  grch38:
    chrom: MT
    start: 0
    end: -1
alleles:
  kind: snv
  ref: A
  alts: [G]
"#,
        );
        let mut issues = Vec::new();
        validate_variant_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");
        assert!(rendered.contains("coordinates.grch38.start:expected integer >= 1"));
        assert!(rendered.contains("coordinates.grch38.end:expected integer >= 1"));
        assert!(rendered.contains("coordinates.grch38.end:expected end >= start"));
    }

    #[test]
    fn panel_assay_and_pgx_roots_report_root_level_shape_errors() {
        let root = yaml(
            r#"
schema: bioscript:panel:1.0
version: "1.0"
name: panel
tags: not-a-list
"#,
        );
        let mut issues = Vec::new();
        validate_panel_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");
        assert!(rendered.contains("tags:expected a sequence of strings"));

        let root = yaml(
            r#"
schema: bioscript:assay:1.0
version: "1.0"
name: assay
tags: [ok, ""]
"#,
        );
        let mut issues = Vec::new();
        validate_assay_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");
        assert!(rendered.contains("tags[1]:empty tag string"));

        let root = yaml(
            r#"
schema: bioscript:pgx-findings:1.0
version: "1.0"
gene: ABC
findings: not-a-list
"#,
        );
        let mut issues = Vec::new();
        validate_pgx_findings_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");
        assert!(rendered.contains("variant/rsid:expected at least one variant identifier"));
        assert!(rendered.contains("findings:expected a sequence of findings"));

        let root = yaml(
            r#"
schema: bioscript:pgx-findings:1.0
version: "1.0"
rsid: rs1
"#,
        );
        let mut issues = Vec::new();
        validate_pgx_findings_root(&root, &mut issues);
        let rendered = messages(&issues).join("\n");
        assert!(rendered.contains("findings:missing required field"));
    }
}
