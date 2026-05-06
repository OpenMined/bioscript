fn validate_findings(root: &Value, issues: &mut Vec<Issue>) {
    let alts = seq_of_strings(root, &["alleles", "alts"]).unwrap_or_default();
    let Some(findings) = value_at(root, &["findings"]).and_then(Value::as_sequence) else {
        return;
    };

    for (idx, finding) in findings.iter().enumerate() {
        let Some(mapping) = finding.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };

        let Some(schema) = mapping
            .get(Value::String("schema".to_owned()))
            .and_then(Value::as_str)
        else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].schema"),
                message: "missing schema".to_owned(),
            });
            continue;
        };
        if schema.trim().is_empty() {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].schema"),
                message: "empty string".to_owned(),
            });
        }
        if schema == "bioscript:pgx:1.0" {
            issues.push(Issue {
                severity: Severity::Warning,
                path: format!("findings[{idx}].schema"),
                message: "legacy PGx finding schema; prefer bioscript:pgx-summary:1.0 or bioscript:pgx-label:1.0".to_owned(),
            });
        }
        if let Some(alt) = mapping
            .get(Value::String("alt".to_owned()))
            .and_then(Value::as_str)
            && !alts.is_empty()
            && alt != "*"
            && !alts.iter().any(|item| item == alt)
        {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].alt"),
                message: format!("finding alt '{alt}' is not present in alleles.alts {alts:?}"),
            });
        }
        let has_summary = mapping
            .get(Value::String("summary".to_owned()))
            .and_then(Value::as_str)
            .is_some_and(|value| !value.trim().is_empty());
        let has_notes = mapping
            .get(Value::String("notes".to_owned()))
            .and_then(Value::as_str)
            .is_some_and(|value| !value.trim().is_empty());
        if !has_summary && !has_notes {
            issues.push(Issue {
                severity: Severity::Warning,
                path: format!("findings[{idx}]"),
                message: "finding has neither summary nor notes".to_owned(),
            });
        }
        validate_finding_binding(&format!("findings[{idx}]"), mapping, issues);
        validate_finding_effects(idx, mapping, issues);
    }
}

fn validate_finding_effects(idx: usize, mapping: &Mapping, issues: &mut Vec<Issue>) {
    let Some(effects) = mapping.get(Value::String("effects".to_owned())) else {
        return;
    };
    let Some(effects) = effects.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("findings[{idx}].effects"),
            message: "expected a sequence of mappings".to_owned(),
        });
        return;
    };
    for (effect_idx, effect) in effects.iter().enumerate() {
        let Some(effect) = effect.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].effects[{effect_idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        validate_finding_binding(
            &format!("findings[{idx}].effects[{effect_idx}]"),
            effect,
            issues,
        );
    }
}

fn validate_finding_binding(parent: &str, mapping: &Mapping, issues: &mut Vec<Issue>) {
    let Some(binding) = mapping.get(Value::String("binding".to_owned())) else {
        return;
    };
    let Some(binding) = binding.as_mapping() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.binding"),
            message: "expected mapping".to_owned(),
        });
        return;
    };
    validate_required_mapping_string(binding, "source", &format!("{parent}.binding"), issues);
    validate_finding_binding_source(parent, binding, issues);
    validate_finding_binding_operator(parent, binding, issues);
}

fn validate_finding_binding_source(parent: &str, binding: &Mapping, issues: &mut Vec<Issue>) {
    let source = binding
        .get(Value::String("source".to_owned()))
        .and_then(Value::as_str);
    match source {
        Some("variant") => {
            if !binding.contains_key(Value::String("variant".to_owned()))
                && !binding.contains_key(Value::String("path".to_owned()))
            {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.variant"),
                    message: "variant findings require variant or path".to_owned(),
                });
            }
        }
        Some("analysis") => {
            validate_required_mapping_string(binding, "key", &format!("{parent}.binding"), issues);
            validate_required_mapping_string(
                binding,
                "analysis_id",
                &format!("{parent}.binding"),
                issues,
            );
        }
        Some(other) => issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.binding.source"),
            message: format!("unsupported source '{other}'"),
        }),
        None => {}
    }
}

fn validate_finding_binding_operator(parent: &str, binding: &Mapping, issues: &mut Vec<Issue>) {
    let operator = binding
        .get(Value::String("operator".to_owned()))
        .and_then(Value::as_str)
        .unwrap_or("equals");
    match operator {
        "equals" => {
            validate_required_mapping_string(binding, "key", &format!("{parent}.binding"), issues);
            if !binding.contains_key(Value::String("value".to_owned())) {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.value"),
                    message: "equals requires value".to_owned(),
                });
            }
        }
        "in" => {
            validate_required_mapping_string(binding, "key", &format!("{parent}.binding"), issues);
            let values = binding
                .get(Value::String("values".to_owned()))
                .and_then(Value::as_sequence);
            if values.is_none_or(Vec::is_empty) {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.values"),
                    message: "in requires non-empty values".to_owned(),
                });
            }
        }
        "dosage_equals" => {
            if binding
                .get(Value::String("allele".to_owned()))
                .and_then(Value::as_str)
                .is_none_or(|value| value.trim().is_empty())
            {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.allele"),
                    message: "dosage_equals requires allele".to_owned(),
                });
            }
            if binding
                .get(Value::String("value".to_owned()))
                .and_then(Value::as_i64)
                .is_none_or(|value| !(0..=2).contains(&value))
            {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.value"),
                    message: "dosage_equals requires integer value 0, 1, or 2".to_owned(),
                });
            }
        }
        "dosage_in" => {
            if binding
                .get(Value::String("allele".to_owned()))
                .and_then(Value::as_str)
                .is_none_or(|value| value.trim().is_empty())
            {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.allele"),
                    message: "dosage_in requires allele".to_owned(),
                });
            }
            let values = binding
                .get(Value::String("values".to_owned()))
                .and_then(Value::as_sequence);
            let invalid_values = match values {
                Some(items) if !items.is_empty() => items
                    .iter()
                    .any(|value| value.as_i64().is_none_or(|n| !(0..=2).contains(&n))),
                _ => true,
            };
            if invalid_values {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("{parent}.binding.values"),
                    message: "dosage_in requires integer values from 0 to 2".to_owned(),
                });
            }
        }
        other => issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.binding.operator"),
            message: format!(
                "unsupported operator '{other}'; expected 'equals', 'in', 'dosage_equals', or 'dosage_in'"
            ),
        }),
    }
}

