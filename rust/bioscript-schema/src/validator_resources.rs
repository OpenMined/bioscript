fn validate_provenance(root: &Value, issues: &mut Vec<Issue>) {
    let Some(sources) = value_at(root, &["provenance", "sources"]).and_then(Value::as_sequence)
    else {
        return;
    };
    for (idx, source) in sources.iter().enumerate() {
        let Some(mapping) = source.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("provenance.sources[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        for field in ["kind", "label", "url"] {
            match mapping
                .get(Value::String(field.to_owned()))
                .and_then(Value::as_str)
            {
                Some(text) if !text.trim().is_empty() => {}
                Some(_) => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("provenance.sources[{idx}].{field}"),
                    message: "empty string".to_owned(),
                }),
                None => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("provenance.sources[{idx}].{field}"),
                    message: "missing required field".to_owned(),
                }),
            }
        }
        if let Some(url) = mapping
            .get(Value::String("url".to_owned()))
            .and_then(Value::as_str)
        {
            validate_url_string(
                url,
                &format!("provenance.sources[{idx}].url"),
                false,
                issues,
            );
        }
    }
}

fn validate_permissions(root: &Value, issues: &mut Vec<Issue>) {
    let Some(domains) = value_at(root, &["permissions", "domains"]) else {
        return;
    };
    let Some(items) = domains.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: "permissions.domains".to_owned(),
            message: "expected a sequence of origins".to_owned(),
        });
        return;
    };
    let mut seen = BTreeSet::new();
    for (idx, item) in items.iter().enumerate() {
        let Some(value) = item.as_str() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("permissions.domains[{idx}]"),
                message: "expected string".to_owned(),
            });
            continue;
        };
        match normalize_origin(value) {
            Ok(origin) => {
                if !seen.insert(origin.clone()) {
                    issues.push(Issue {
                        severity: Severity::Warning,
                        path: format!("permissions.domains[{idx}]"),
                        message: format!("duplicate origin '{origin}'"),
                    });
                }
            }
            Err(message) => issues.push(Issue {
                severity: Severity::Error,
                path: format!("permissions.domains[{idx}]"),
                message,
            }),
        }
    }
}

fn validate_downloads(root: &Value, issues: &mut Vec<Issue>) {
    let allowed_origins: BTreeSet<String> = seq_of_strings(root, &["permissions", "domains"])
        .unwrap_or_default()
        .into_iter()
        .filter_map(|domain| normalize_origin(&domain).ok())
        .collect();
    let Some(downloads) = value_at(root, &["downloads"]).and_then(Value::as_sequence) else {
        return;
    };
    let mut ids = BTreeSet::new();
    for (idx, item) in downloads.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("downloads[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        for field in ["id", "url", "sha256", "version"] {
            match mapping
                .get(Value::String(field.to_owned()))
                .and_then(Value::as_str)
            {
                Some(text) if !text.trim().is_empty() => {}
                Some(_) => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("downloads[{idx}].{field}"),
                    message: "empty string".to_owned(),
                }),
                None => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("downloads[{idx}].{field}"),
                    message: "missing required field".to_owned(),
                }),
            }
        }

        if let Some(id) = mapping
            .get(Value::String("id".to_owned()))
            .and_then(Value::as_str)
            && !ids.insert(id.to_owned())
        {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("downloads[{idx}].id"),
                message: format!("duplicate download id '{id}'"),
            });
        }
        if let Some(sha) = mapping
            .get(Value::String("sha256".to_owned()))
            .and_then(Value::as_str)
            && !is_sha256(sha)
        {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("downloads[{idx}].sha256"),
                message: "expected 64 lowercase hex characters".to_owned(),
            });
        }
        if let Some(url) = mapping
            .get(Value::String("url".to_owned()))
            .and_then(Value::as_str)
        {
            match normalize_download_url(url) {
                Ok(origin) => {
                    if !allowed_origins.is_empty() && !allowed_origins.contains(&origin) {
                        issues.push(Issue {
                            severity: Severity::Error,
                            path: format!("downloads[{idx}].url"),
                            message: format!(
                                "download origin '{origin}' is not listed in permissions.domains"
                            ),
                        });
                    }
                }
                Err(message) => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("downloads[{idx}].url"),
                    message,
                }),
            }
        }
    }
}

