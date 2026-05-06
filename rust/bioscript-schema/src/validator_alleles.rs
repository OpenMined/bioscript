fn validate_alleles(root: &Value, issues: &mut Vec<Issue>) {
    require_path(root, &["alleles"], issues);
    require_path(root, &["alleles", "kind"], issues);
    require_path(root, &["alleles", "ref"], issues);
    require_path(root, &["alleles", "alts"], issues);

    let Some(kind) = scalar_at(root, &["alleles", "kind"]) else {
        return;
    };
    if !matches!(kind.as_str(), "snv" | "deletion" | "insertion" | "indel") {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.kind".to_owned(),
            message: "expected one of snv, deletion, insertion, indel".to_owned(),
        });
    }

    if value_at(root, &["alleles", "canonical_alt"]).is_some() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.canonical_alt".to_owned(),
            message: "canonical_alt is not part of the current schema".to_owned(),
        });
    }

    let Some(reference) = scalar_at(root, &["alleles", "ref"]) else {
        return;
    };
    if reference.trim().is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.ref".to_owned(),
            message: "empty string".to_owned(),
        });
    }

    let Some(alts_value) = value_at(root, &["alleles", "alts"]) else {
        return;
    };
    let Some(alts_seq) = alts_value.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.alts".to_owned(),
            message: "expected a non-empty sequence of strings".to_owned(),
        });
        return;
    };
    if alts_seq.is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.alts".to_owned(),
            message: "expected at least one alternate allele".to_owned(),
        });
        return;
    }

    let mut alts = Vec::new();
    for (idx, item) in alts_seq.iter().enumerate() {
        let Some(alt) = item.as_str() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message: "expected string".to_owned(),
            });
            continue;
        };
        if alt.trim().is_empty() {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message: "empty string".to_owned(),
            });
            continue;
        }
        alts.push(alt.to_owned());
    }
    let observed_alts = match seq_of_strings(root, &["alleles", "observed_alts"]) {
        Some(items) => {
            if items.is_empty() {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: "alleles.observed_alts".to_owned(),
                    message: "expected a non-empty sequence of strings when present".to_owned(),
                });
            }
            for alt in &alts {
                if !items.iter().any(|item| item == alt) {
                    issues.push(Issue {
                        severity: Severity::Error,
                        path: "alleles.observed_alts".to_owned(),
                        message: format!("significant alt '{alt}' is not present in observed_alts"),
                    });
                }
            }
            items
        }
        None => alts.clone(),
    };
    validate_symbolic_alleles(&reference, &observed_alts, issues);
    validate_snv_alleles(&kind, &reference, &observed_alts, issues);
}

fn validate_symbolic_alleles(reference: &str, alts: &[String], issues: &mut Vec<Issue>) {
    if reference == "I" || reference == "D" {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.ref".to_owned(),
            message: "symbolic I/D alleles are not allowed in stored YAML; use biological alleles"
                .to_owned(),
        });
    }
    for (idx, alt) in alts.iter().enumerate() {
        if alt == "I" || alt == "D" {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message:
                    "symbolic I/D alleles are not allowed in stored YAML; use biological alleles"
                        .to_owned(),
            });
        }
    }
}

fn validate_snv_alleles(kind: &str, reference: &str, alts: &[String], issues: &mut Vec<Issue>) {
    if kind != "snv" {
        return;
    }
    if !is_base_allele(reference) {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.ref".to_owned(),
            message: "snv ref must be one of A/C/G/T".to_owned(),
        });
    }
    for (idx, alt) in alts.iter().enumerate() {
        if !is_base_allele(alt) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message: "snv alt must be one of A/C/G/T".to_owned(),
            });
        }
    }
}

