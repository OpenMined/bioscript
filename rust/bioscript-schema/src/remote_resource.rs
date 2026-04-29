use std::{fmt::Write as _, path::Path};

use serde::Serialize;
use serde_yaml::Value;
use sha2::{Digest, Sha256};
use url::Url;

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum RemoteResourceKind {
    Assay,
    Catalogue,
    Panel,
    Python,
    Unknown,
    Variant,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RemoteDependency {
    pub kind: String,
    pub label: String,
    pub optional: bool,
    pub url: String,
    pub version: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct RemoteResourceResolution {
    pub dependencies: Vec<RemoteDependency>,
    pub kind: RemoteResourceKind,
    pub name: String,
    pub schema: Option<String>,
    pub sha256: String,
    pub source_url: String,
    pub title: String,
    pub version: Option<String>,
}

/// Classify a fetched remote resource and discover dependency requirements.
///
/// This function intentionally does not perform network access. Every platform
/// asks for user consent, fetches the target bytes with its own transport, then
/// calls this shared parser before asking whether dependencies should be
/// fetched.
///
/// # Errors
///
/// Returns an error when YAML/JSON-like content cannot be parsed for recognized
/// file extensions, or when relative dependency URLs cannot be resolved.
pub fn resolve_remote_resource_text(
    source_url: &str,
    name: &str,
    text: &str,
) -> Result<RemoteResourceResolution, String> {
    let parsed = parse_structured_text(name, text)?;
    let schema = parsed
        .as_ref()
        .and_then(|value| scalar_at(value, &["schema"]));
    let kind = infer_kind(name, schema.as_deref(), parsed.as_ref());
    let version = parsed.as_ref().and_then(extract_version);
    let title = parsed
        .as_ref()
        .and_then(|value| {
            scalar_at(value, &["label"])
                .or_else(|| scalar_at(value, &["title"]))
                .or_else(|| scalar_at(value, &["name"]))
                .or_else(|| scalar_at(value, &["variant_id"]))
        })
        .unwrap_or_else(|| name.to_owned());
    let dependencies = parsed
        .as_ref()
        .map_or_else(Vec::new, |value| collect_dependencies(source_url, value));

    Ok(RemoteResourceResolution {
        dependencies,
        kind,
        name: name.to_owned(),
        schema,
        sha256: sha256_hex(text.as_bytes()),
        source_url: source_url.to_owned(),
        title,
        version,
    })
}

fn parse_structured_text(name: &str, text: &str) -> Result<Option<Value>, String> {
    if has_extension(name, &["yaml", "yml"]) {
        let value: Value = serde_yaml::from_str(text)
            .map_err(|err| format!("failed to parse YAML resource {name}: {err}"))?;
        return Ok(Some(value));
    }
    if has_extension(name, &["json"]) {
        let value: Value = serde_yaml::from_str(text)
            .map_err(|err| format!("failed to parse JSON resource {name}: {err}"))?;
        return Ok(Some(value));
    }
    Ok(None)
}

fn infer_kind(name: &str, schema: Option<&str>, parsed: Option<&Value>) -> RemoteResourceKind {
    if has_extension(name, &["py"]) {
        return RemoteResourceKind::Python;
    }
    if let Some(schema) = schema {
        let lower = schema.to_ascii_lowercase();
        if lower.contains("variant") {
            return RemoteResourceKind::Variant;
        }
        if lower.contains("panel") {
            return RemoteResourceKind::Panel;
        }
        if lower.contains("catalogue") || lower.contains("catalog") || lower.contains("index") {
            return RemoteResourceKind::Catalogue;
        }
        if lower.contains("assay") {
            return RemoteResourceKind::Assay;
        }
    }
    let Some(value) = parsed else {
        return RemoteResourceKind::Unknown;
    };
    if value_at(value, &["members"])
        .and_then(Value::as_sequence)
        .is_some()
        || value_at(value, &["variants"])
            .and_then(Value::as_sequence)
            .is_some()
    {
        return RemoteResourceKind::Panel;
    }
    if value_at(value, &["assays"])
        .and_then(Value::as_sequence)
        .is_some()
    {
        return RemoteResourceKind::Catalogue;
    }
    if value_at(value, &["assay"]).is_some() {
        return RemoteResourceKind::Assay;
    }
    if value_at(value, &["rsid"]).is_some()
        || value_at(value, &["variant_id"]).is_some()
        || value_at(value, &["coordinates"]).is_some()
        || value_at(value, &["alleles"]).is_some()
    {
        return RemoteResourceKind::Variant;
    }
    RemoteResourceKind::Unknown
}

fn collect_dependencies(source_url: &str, root: &Value) -> Vec<RemoteDependency> {
    let mut dependencies = Vec::new();
    collect_dependency_values(source_url, root, &mut Vec::new(), &mut dependencies);
    dependencies.sort_by(|left, right| left.url.cmp(&right.url));
    dependencies.dedup_by(|left, right| left.url == right.url);
    dependencies
}

fn collect_dependency_values(
    source_url: &str,
    value: &Value,
    path: &mut Vec<String>,
    dependencies: &mut Vec<RemoteDependency>,
) {
    match value {
        Value::Sequence(items) => {
            for (idx, item) in items.iter().enumerate() {
                path.push(idx.to_string());
                collect_dependency_values(source_url, item, path, dependencies);
                path.pop();
            }
        }
        Value::Mapping(mapping) => {
            for (key, entry) in mapping {
                let Some(key) = key.as_str() else {
                    continue;
                };
                path.push(key.to_owned());
                if let Some(text) = entry.as_str()
                    && looks_like_dependency_key(key)
                    && looks_like_resource_value(text)
                    && let Some(url) = resolve_resource_url(source_url, text)
                {
                    dependencies.push(RemoteDependency {
                        kind: dependency_kind_from_path(path),
                        label: path.join("."),
                        optional: false,
                        url,
                        version: sibling_scalar(mapping, "version"),
                    });
                }
                collect_dependency_values(source_url, entry, path, dependencies);
                path.pop();
            }
        }
        _ => {}
    }
}

fn dependency_kind_from_path(path: &[String]) -> String {
    if path.iter().any(|part| part == "members") {
        "member".to_owned()
    } else if path.iter().any(|part| part == "downloads") {
        "download".to_owned()
    } else {
        "file".to_owned()
    }
}

fn looks_like_dependency_key(key: &str) -> bool {
    matches!(
        key,
        "artifact_url"
            | "catalog"
            | "catalogue"
            | "compiled_path"
            | "dependency"
            | "download"
            | "file"
            | "index"
            | "panel"
            | "path"
            | "url"
            | "variant"
            | "variants"
    )
}

fn looks_like_resource_value(value: &str) -> bool {
    value.starts_with("http://")
        || value.starts_with("https://")
        || has_extension(value, &["yaml", "yml", "json", "py"])
}

fn resolve_resource_url(source_url: &str, value: &str) -> Option<String> {
    if let Ok(url) = Url::parse(value) {
        return Some(url.to_string());
    }
    let source = Url::parse(source_url).ok()?;
    if source.host_str() == Some("github.com") {
        let mut parts = source
            .path_segments()?
            .map(ToOwned::to_owned)
            .collect::<Vec<_>>();
        let blob_idx = parts.iter().position(|part| part == "blob")?;
        if value.starts_with('/') {
            let owner = parts.first()?.clone();
            let repo = parts.get(1)?.clone();
            let reference = parts.get(blob_idx + 1)?.clone();
            return Url::parse(&format!(
                "https://github.com/{owner}/{repo}/blob/{reference}{value}"
            ))
            .ok()
            .map(|url| url.to_string());
        }
        parts.pop();
        parts.push(value.to_owned());
        return Url::parse(&format!("https://github.com/{}", parts.join("/")))
            .ok()
            .map(|url| url.to_string());
    }
    source.join(value).ok().map(|url| url.to_string())
}

fn extract_version(root: &Value) -> Option<String> {
    scalar_at(root, &["version"])
        .or_else(|| scalar_at(root, &["package_version"]))
        .or_else(|| scalar_at(root, &["schema_version"]))
        .or_else(|| scalar_at(root, &["assay", "package_version"]))
        .or_else(|| scalar_at(root, &["assay", "version"]))
}

fn sibling_scalar(mapping: &serde_yaml::Mapping, field: &str) -> Option<String> {
    mapping
        .get(Value::String(field.to_owned()))
        .and_then(Value::as_str)
        .map(ToOwned::to_owned)
}

fn value_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut current = root;
    for segment in path {
        current = current
            .as_mapping()?
            .get(Value::String((*segment).to_owned()))?;
    }
    Some(current)
}

fn scalar_at(root: &Value, path: &[&str]) -> Option<String> {
    value_at(root, path)
        .and_then(Value::as_str)
        .map(ToOwned::to_owned)
}

fn sha256_hex(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    digest.iter().fold(String::new(), |mut output, byte| {
        let _ = write!(output, "{byte:02x}");
        output
    })
}

fn has_extension(value: &str, extensions: &[&str]) -> bool {
    Path::new(value).extension().is_some_and(|extension| {
        extensions
            .iter()
            .any(|item| extension.eq_ignore_ascii_case(item))
    })
}
