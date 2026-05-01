use std::{
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
};

use serde_yaml::{Mapping, Value};
use url::Url;

use super::{Issue, Severity};

pub(crate) fn collect_yaml_files(path: &Path) -> Result<Vec<PathBuf>, String> {
    if path.is_file() {
        return Ok(vec![path.to_path_buf()]);
    }

    let mut files = Vec::new();
    collect_yaml_files_recursive(path, &mut files)?;
    files.sort();
    Ok(files)
}

fn collect_yaml_files_recursive(path: &Path, files: &mut Vec<PathBuf>) -> Result<(), String> {
    let entries = fs::read_dir(path)
        .map_err(|err| format!("failed to read directory {}: {err}", path.display()))?;
    for entry in entries {
        let entry = entry.map_err(|err| format!("failed to read directory entry: {err}"))?;
        let entry_path = entry.path();
        if entry_path.is_dir() {
            collect_yaml_files_recursive(&entry_path, files)?;
            continue;
        }
        if entry_path.extension().is_some_and(|extension| {
            ["yaml", "yml"]
                .iter()
                .any(|item| extension.eq_ignore_ascii_case(item))
        }) {
            files.push(entry_path);
        }
    }
    Ok(())
}

pub(crate) fn validate_schema_and_identity(
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

pub(crate) fn validate_optional_strings(root: &Value, fields: &[&str], issues: &mut Vec<Issue>) {
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

pub(crate) fn validate_tags(root: &Value, issues: &mut Vec<Issue>) {
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

pub(crate) fn mapping_required_string(
    mapping: &Mapping,
    field: &str,
    idx: usize,
    parent: &str,
) -> Result<String, String> {
    mapping
        .get(Value::String(field.to_owned()))
        .and_then(Value::as_str)
        .filter(|value| !value.trim().is_empty())
        .map(ToOwned::to_owned)
        .ok_or_else(|| format!("{parent}[{idx}].{field} missing or empty"))
}

pub(crate) fn validate_url_string(
    value: &str,
    path: &str,
    require_origin_only: bool,
    issues: &mut Vec<Issue>,
) {
    let normalized = if require_origin_only {
        normalize_origin(value)
    } else {
        normalize_download_url(value)
    };
    if let Err(message) = normalized {
        issues.push(Issue {
            severity: Severity::Error,
            path: path.to_owned(),
            message,
        });
    }
}

pub(crate) fn normalize_origin(value: &str) -> Result<String, String> {
    let url = Url::parse(value).map_err(|err| format!("invalid URL: {err}"))?;
    if !matches!(url.scheme(), "http" | "https") {
        return Err("expected http or https origin".to_owned());
    }
    if url.host_str().is_none() {
        return Err("origin is missing host".to_owned());
    }
    if url.path() != "/" || url.query().is_some() || url.fragment().is_some() {
        return Err("expected origin only, without path, query, or fragment".to_owned());
    }
    let mut origin = format!("{}://{}", url.scheme(), url.host_str().unwrap_or_default());
    if let Some(port) = url.port() {
        let _ = write!(origin, ":{port}");
    }
    Ok(origin)
}

pub(crate) fn normalize_download_url(value: &str) -> Result<String, String> {
    let url = Url::parse(value).map_err(|err| format!("invalid URL: {err}"))?;
    if !matches!(url.scheme(), "http" | "https") {
        return Err("expected http or https URL".to_owned());
    }
    if url.host_str().is_none() {
        return Err("URL is missing host".to_owned());
    }
    let mut origin = format!("{}://{}", url.scheme(), url.host_str().unwrap_or_default());
    if let Some(port) = url.port() {
        let _ = write!(origin, ":{port}");
    }
    Ok(origin)
}

pub(crate) fn is_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .chars()
            .all(|ch| ch.is_ascii_hexdigit() && !ch.is_ascii_uppercase())
}

pub(crate) fn i64_at_mapping(mapping: &Mapping, key: &str) -> Option<i64> {
    mapping
        .get(Value::String(key.to_owned()))
        .and_then(Value::as_i64)
}

pub(crate) fn required_non_empty_string(root: &Value, path: &[&str]) -> Result<String, String> {
    scalar_at(root, path)
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| format!("{} missing or empty", path.join(".")))
}

pub(crate) fn render_single_manifest_errors(path: &Path, issues: &[Issue]) -> String {
    let mut out = format!("invalid manifest {}:\n", path.display());
    for issue in issues {
        let _ = writeln!(
            out,
            "  - [{}] {}: {}",
            issue.severity, issue.path, issue.message
        );
    }
    out
}

pub(crate) fn load_yaml(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))
}

pub(crate) fn require_const(root: &Value, path: &[&str], expected: &str, issues: &mut Vec<Issue>) {
    match scalar_at(root, path) {
        Some(actual) if actual == expected => {}
        Some(actual) => issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: format!("expected '{expected}', found '{actual}'"),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: "missing required field".to_owned(),
        }),
    }
}

pub(crate) fn require_path(root: &Value, path: &[&str], issues: &mut Vec<Issue>) {
    if value_at(root, path).is_none() {
        issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: "missing required field".to_owned(),
        });
    }
}

pub(crate) fn value_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut current = root;
    for key in path {
        let mapping = current.as_mapping()?;
        current = mapping.get(Value::String((*key).to_owned()))?;
    }
    Some(current)
}

pub(crate) fn mapping_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Mapping> {
    value_at(root, path)?.as_mapping()
}

pub(crate) fn scalar_at(root: &Value, path: &[&str]) -> Option<String> {
    value_at(root, path).and_then(|value| match value {
        Value::String(text) => Some(text.clone()),
        Value::Number(number) => Some(number.to_string()),
        _ => None,
    })
}

pub(crate) fn seq_of_strings(root: &Value, path: &[&str]) -> Option<Vec<String>> {
    value_at(root, path)?.as_sequence().map(|items| {
        items
            .iter()
            .filter_map(|item| item.as_str().map(ToOwned::to_owned))
            .collect()
    })
}
