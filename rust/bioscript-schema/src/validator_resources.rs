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

#[cfg(test)]
mod resource_validator_tests {
    use super::*;

    fn yaml(text: &str) -> Value {
        serde_yaml::from_str(text).unwrap()
    }

    fn issue_paths(issues: &[Issue]) -> Vec<String> {
        issues.iter().map(|issue| issue.path.clone()).collect()
    }

    #[test]
    fn provenance_validation_reports_source_shape_and_url_errors() {
        let root = yaml(
            r#"
provenance:
  sources:
    - not-a-map
    - kind: ""
      label: 3
      url: not-a-url
"#,
        );
        let mut issues = Vec::new();
        validate_provenance(&root, &mut issues);
        let paths = issue_paths(&issues);
        assert!(paths.contains(&"provenance.sources[0]".to_owned()));
        assert!(paths.contains(&"provenance.sources[1].kind".to_owned()));
        assert!(paths.contains(&"provenance.sources[1].label".to_owned()));
        assert!(paths.contains(&"provenance.sources[1].url".to_owned()));
    }

    #[test]
    fn permissions_validation_reports_origin_shape_duplicates_and_bad_urls() {
        let root = yaml(
            r#"
permissions:
  domains:
    - https://example.com
    - https://example.com
    - 3
    - ftp://example.com
"#,
        );
        let mut issues = Vec::new();
        validate_permissions(&root, &mut issues);
        let paths = issue_paths(&issues);
        assert!(paths.contains(&"permissions.domains[1]".to_owned()));
        assert!(paths.contains(&"permissions.domains[2]".to_owned()));
        assert!(paths.contains(&"permissions.domains[3]".to_owned()));

        let mut issues = Vec::new();
        validate_permissions(&yaml("permissions: {domains: nope}"), &mut issues);
        assert_eq!(issues[0].path, "permissions.domains");
    }

    #[test]
    fn downloads_validation_reports_required_hash_duplicate_and_origin_errors() {
        let root = yaml(
            r#"
permissions:
  domains: [https://allowed.example]
downloads:
  - not-a-map
  - id: dl
    url: https://blocked.example/file.txt
    sha256: BAD
    version: ""
  - id: dl
    url: not-a-url
    sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    version: "1"
  - id: ""
    sha256: ""
"#,
        );
        let mut issues = Vec::new();
        validate_downloads(&root, &mut issues);
        let paths = issue_paths(&issues);
        assert!(paths.contains(&"downloads[0]".to_owned()));
        assert!(paths.contains(&"downloads[1].url".to_owned()));
        assert!(paths.contains(&"downloads[1].sha256".to_owned()));
        assert!(paths.contains(&"downloads[1].version".to_owned()));
        assert!(paths.contains(&"downloads[2].id".to_owned()));
        assert!(paths.contains(&"downloads[2].url".to_owned()));
        assert!(paths.contains(&"downloads[3].id".to_owned()));
        assert!(paths.contains(&"downloads[3].url".to_owned()));
        assert!(paths.contains(&"downloads[3].version".to_owned()));
    }
}
