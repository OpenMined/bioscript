use std::collections::BTreeSet;

use serde_yaml::{Mapping, Value};

use super::{
    Download, Issue, PanelMember, Severity,
    common::{
        is_sha256, mapping_required_string, normalize_download_url, normalize_origin,
        seq_of_strings, validate_optional_strings, validate_schema_and_identity, validate_tags,
        value_at,
    },
};

pub(crate) fn validate_panel_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(root, "bioscript:panel:1.0", None, issues);
    validate_optional_strings(root, &["name", "label", "summary"], issues);
    validate_tags(root, issues);
    validate_permissions(root, issues);
    validate_downloads(root, issues);
    validate_panel_members(root, issues);
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

fn validate_panel_members(root: &Value, issues: &mut Vec<Issue>) {
    let Some(members) = value_at(root, &["members"]).and_then(Value::as_sequence) else {
        issues.push(Issue {
            severity: Severity::Error,
            path: "members".to_owned(),
            message: "missing required field".to_owned(),
        });
        return;
    };
    if members.is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "members".to_owned(),
            message: "expected at least one member".to_owned(),
        });
        return;
    }

    let download_ids = panel_download_ids(root);

    for (idx, item) in members.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("members[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        validate_panel_member(idx, mapping, &download_ids, issues);
    }
}

fn panel_download_ids(root: &Value) -> BTreeSet<String> {
    value_at(root, &["downloads"])
        .and_then(Value::as_sequence)
        .into_iter()
        .flatten()
        .filter_map(|item| {
            item.as_mapping()?
                .get(Value::String("id".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned)
        })
        .collect()
}

fn validate_panel_member(
    idx: usize,
    mapping: &Mapping,
    download_ids: &BTreeSet<String>,
    issues: &mut Vec<Issue>,
) {
    let kind = mapping
        .get(Value::String("kind".to_owned()))
        .and_then(Value::as_str);
    match kind {
        Some("variant") => {}
        Some(other) => issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].kind"),
            message: format!(
                "unsupported member kind '{other}'; panel support is currently variant-only"
            ),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].kind"),
            message: "missing required field".to_owned(),
        }),
    }

    let path_value = mapping
        .get(Value::String("path".to_owned()))
        .and_then(Value::as_str);
    let download_value = mapping
        .get(Value::String("download".to_owned()))
        .and_then(Value::as_str);
    if path_value.is_some() == download_value.is_some() {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}]"),
            message: "expected exactly one of path or download".to_owned(),
        });
    }
    validate_panel_member_path(idx, path_value, issues);
    validate_panel_member_download(idx, download_value, download_ids, issues);
    validate_panel_member_metadata(idx, mapping, issues);
}

fn validate_panel_member_path(idx: usize, path_value: Option<&str>, issues: &mut Vec<Issue>) {
    if let Some(path) = path_value
        && path.trim().is_empty()
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].path"),
            message: "empty string".to_owned(),
        });
    }
}

fn validate_panel_member_download(
    idx: usize,
    download_value: Option<&str>,
    download_ids: &BTreeSet<String>,
    issues: &mut Vec<Issue>,
) {
    let Some(download) = download_value else {
        return;
    };
    if download.trim().is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].download"),
            message: "empty string".to_owned(),
        });
    } else if !download_ids.contains(download) {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].download"),
            message: format!("unknown download id '{download}'"),
        });
    }
}

fn validate_panel_member_metadata(idx: usize, mapping: &Mapping, issues: &mut Vec<Issue>) {
    if let Some(version) = mapping
        .get(Value::String("version".to_owned()))
        .and_then(Value::as_str)
        && version.trim().is_empty()
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].version"),
            message: "empty string".to_owned(),
        });
    }
    if let Some(sha) = mapping
        .get(Value::String("sha256".to_owned()))
        .and_then(Value::as_str)
        && !is_sha256(sha)
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].sha256"),
            message: "expected 64 lowercase hex characters".to_owned(),
        });
    }
}

pub(crate) fn parse_downloads(root: &Value) -> Result<Vec<Download>, String> {
    let mut downloads = Vec::new();
    let Some(items) = value_at(root, &["downloads"]).and_then(Value::as_sequence) else {
        return Ok(downloads);
    };
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            return Err(format!("downloads[{idx}] must be a mapping"));
        };
        let id = mapping_required_string(mapping, "id", idx, "downloads")?;
        let url = mapping_required_string(mapping, "url", idx, "downloads")?;
        let sha256 = mapping_required_string(mapping, "sha256", idx, "downloads")?;
        let version = mapping_required_string(mapping, "version", idx, "downloads")?;
        let origin = normalize_download_url(&url)?;
        downloads.push(Download {
            id,
            url,
            origin,
            sha256,
            version,
        });
    }
    Ok(downloads)
}

pub(crate) fn parse_panel_members(root: &Value) -> Result<Vec<PanelMember>, String> {
    let mut members = Vec::new();
    let Some(items) = value_at(root, &["members"]).and_then(Value::as_sequence) else {
        return Ok(members);
    };
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            return Err(format!("members[{idx}] must be a mapping"));
        };
        members.push(PanelMember {
            kind: mapping_required_string(mapping, "kind", idx, "members")?,
            path: mapping
                .get(Value::String("path".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            download: mapping
                .get(Value::String("download".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            sha256: mapping
                .get(Value::String("sha256".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            version: mapping
                .get(Value::String("version".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
        });
    }
    Ok(members)
}
