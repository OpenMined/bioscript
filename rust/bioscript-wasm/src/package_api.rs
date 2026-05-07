use std::{
    collections::BTreeSet,
    io::{Cursor, Read},
    path::{Path, PathBuf},
};

use bioscript_schema::{
    resolve_remote_resource_text, RemoteResourceKind, RemoteResourceResolution,
};
use serde::Serialize;
use sha2::{Digest, Sha256};
use wasm_bindgen::prelude::*;

const PACKAGE_DESCRIPTOR: &str = "manifest.yaml";
const LEGACY_PACKAGE_DESCRIPTOR: &str = "bioscript-package.yaml";
const MAX_PACKAGE_FILES: usize = 1000;
const MAX_PACKAGE_FILE_BYTES: u64 = 16 * 1024 * 1024;
const MAX_PACKAGE_TOTAL_BYTES: u64 = 64 * 1024 * 1024;

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PackageFileJs {
    path: String,
    contents: String,
    source_url: String,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PackageResourceJs {
    path: String,
    contents: String,
    resolution: RemoteResourceResolution,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PackageResolutionJs {
    entrypoint: String,
    files: Vec<PackageFileJs>,
    name: Option<String>,
    resources: Vec<PackageResourceJs>,
}

#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct PackageReleaseJs {
    artifact_sha256: Option<String>,
    artifact_size_bytes: Option<u64>,
    artifact_url: String,
    entrypoint: Option<String>,
    name: Option<String>,
    title: String,
    version: Option<String>,
}

struct PackageDescriptor {
    entrypoint: PathBuf,
    name: Option<String>,
}

struct ExtractedPackageFile {
    path: PathBuf,
    contents: String,
}

/// Resolve a BioScript package zip from bytes.
///
/// This mirrors the CLI package importer enough for browser/mobile callers:
/// path safety, package size limits, descriptor/entrypoint discovery, and
/// resource classification all stay in Rust.
#[wasm_bindgen(js_name = resolvePackageZipBytes)]
pub fn resolve_package_zip_bytes(
    source_url: &str,
    name: &str,
    bytes: &[u8],
) -> Result<String, JsError> {
    let files = extract_package_zip(name, bytes)
        .map_err(|err| JsError::new(&format!("resolve package zip failed: {err}")))?;
    let descriptor = load_package_descriptor(&files)
        .map_err(|err| JsError::new(&format!("resolve package descriptor failed: {err}")))?;
    let entrypoint = descriptor.entrypoint.display().to_string();
    let entry_file = files
        .iter()
        .find(|file| file.path == descriptor.entrypoint)
        .ok_or_else(|| JsError::new(&format!("package entrypoint not found: {entrypoint}")))?;
    let entry_resolution = resolve_remote_resource_text(
        &package_member_url(source_url, &descriptor.entrypoint),
        &entrypoint,
        &entry_file.contents,
    )
    .map_err(|err| JsError::new(&format!("resolve package entrypoint failed: {err}")))?;
    match entry_resolution.kind {
        RemoteResourceKind::Assay
        | RemoteResourceKind::Panel
        | RemoteResourceKind::Python
        | RemoteResourceKind::Variant => {}
        _ => {
            return Err(JsError::new(&format!(
                "package entrypoint has unsupported resource kind: {:?}",
                entry_resolution.kind
            )));
        }
    }

    let mut resources = Vec::new();
    for file in &files {
        if !is_resource_file(&file.path) {
            continue;
        }
        let path = file.path.display().to_string();
        let resolution = resolve_remote_resource_text(
            &package_member_url(source_url, &file.path),
            &path,
            &file.contents,
        )
        .map_err(|err| JsError::new(&format!("resolve package member {path} failed: {err}")))?;
        if matches!(
            resolution.kind,
            RemoteResourceKind::Assay
                | RemoteResourceKind::Panel
                | RemoteResourceKind::Python
                | RemoteResourceKind::Variant
        ) {
            resources.push(PackageResourceJs {
                path,
                contents: file.contents.clone(),
                resolution,
            });
        }
    }

    let files_js = files
        .iter()
        .map(|file| PackageFileJs {
            path: file.path.display().to_string(),
            contents: file.contents.clone(),
            source_url: package_member_url(source_url, &file.path),
        })
        .collect();

    serde_json::to_string(&PackageResolutionJs {
        entrypoint,
        files: files_js,
        name: descriptor.name,
        resources,
    })
    .map_err(|err| JsError::new(&format!("failed to encode package response: {err}")))
}

/// Resolve a BioScript package release YAML into the package zip artifact URL.
#[wasm_bindgen(js_name = resolvePackageReleaseText)]
pub fn resolve_package_release_text(
    source_url: &str,
    name: &str,
    text: &str,
) -> Result<String, JsError> {
    let value: serde_yaml::Value = serde_yaml::from_str(text)
        .map_err(|err| JsError::new(&format!("failed to parse package release {name}: {err}")))?;
    let schema = yaml_string(&value, "schema");
    if schema.as_deref() != Some("bioscript:package-release:1.0") {
        return Err(JsError::new(&format!(
            "{name} is not a bioscript:package-release:1.0 manifest"
        )));
    }
    let artifact = value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("artifact".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .ok_or_else(|| JsError::new(&format!("package release {name} is missing artifact")))?;
    let artifact_path = artifact
        .get(serde_yaml::Value::String("path".to_owned()))
        .and_then(serde_yaml::Value::as_str);
    let artifact_url = artifact
        .get(serde_yaml::Value::String("url".to_owned()))
        .and_then(serde_yaml::Value::as_str);
    let artifact_url = if let Some(url) = artifact_url {
        url.to_owned()
    } else if let Some(relative) = artifact_path {
        join_url(source_url, relative)
    } else {
        return Err(JsError::new(&format!(
            "package release {name} artifact needs path or url"
        )));
    };
    let artifact_sha256 = artifact
        .get(serde_yaml::Value::String("sha256".to_owned()))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned);
    let artifact_size_bytes = artifact
        .get(serde_yaml::Value::String("size_bytes".to_owned()))
        .and_then(serde_yaml::Value::as_u64);
    let title = scalar_at(&value, "label")
        .or_else(|| scalar_at(&value, "title"))
        .or_else(|| scalar_at(&value, "name"))
        .unwrap_or_else(|| name.to_owned());
    let release = PackageReleaseJs {
        artifact_sha256,
        artifact_size_bytes,
        artifact_url,
        entrypoint: scalar_at(&value, "entrypoint"),
        name: scalar_at(&value, "name"),
        title,
        version: scalar_at(&value, "package_version").or_else(|| scalar_at(&value, "version")),
    };
    serde_json::to_string(&release)
        .map_err(|err| JsError::new(&format!("failed to encode package release: {err}")))
}

/// Verify package artifact bytes against a package-release sha256 value.
#[wasm_bindgen(js_name = verifyPackageArtifactSha256)]
pub fn verify_package_artifact_sha256(
    name: &str,
    bytes: &[u8],
    expected: &str,
) -> Result<(), JsError> {
    let actual = sha256_hex(bytes);
    if actual != expected {
        return Err(JsError::new(&format!(
            "package artifact sha256 mismatch for {name}: expected {expected}, got {actual}"
        )));
    }
    Ok(())
}

fn extract_package_zip(name: &str, bytes: &[u8]) -> Result<Vec<ExtractedPackageFile>, String> {
    let mut archive = zip::ZipArchive::new(Cursor::new(bytes))
        .map_err(|err| format!("failed to read package zip {name}: {err}"))?;
    if archive.len() > MAX_PACKAGE_FILES {
        return Err(format!(
            "package has too many entries: {} > {MAX_PACKAGE_FILES}",
            archive.len()
        ));
    }

    let mut seen = BTreeSet::new();
    let mut total_size = 0_u64;
    let mut files = Vec::new();
    for idx in 0..archive.len() {
        let mut entry = archive
            .by_index(idx)
            .map_err(|err| format!("failed to read package zip entry {idx}: {err}"))?;
        let Some(enclosed) = entry.enclosed_name() else {
            return Err(format!(
                "package zip entry has unsafe path: {}",
                entry.name()
            ));
        };
        let relative = checked_relative_package_path(&enclosed.to_string_lossy())?;
        if entry.is_dir() {
            continue;
        }
        if entry
            .unix_mode()
            .is_some_and(|mode| mode & 0o170_000 == 0o120_000)
        {
            return Err(format!("package zip entry is a symlink: {}", entry.name()));
        }
        if !is_allowed_package_file(&relative) {
            return Err(format!(
                "package zip entry has unsupported extension: {}",
                relative.display()
            ));
        }
        if !seen.insert(relative.clone()) {
            return Err(format!(
                "package zip contains duplicate path: {}",
                relative.display()
            ));
        }
        let size = entry.size();
        if size > MAX_PACKAGE_FILE_BYTES {
            return Err(format!(
                "package file {} exceeds {} bytes",
                relative.display(),
                MAX_PACKAGE_FILE_BYTES
            ));
        }
        total_size = total_size.saturating_add(size);
        if total_size > MAX_PACKAGE_TOTAL_BYTES {
            return Err(format!(
                "package contents exceed {MAX_PACKAGE_TOTAL_BYTES} bytes"
            ));
        }
        let mut contents = String::new();
        entry
            .read_to_string(&mut contents)
            .map_err(|err| format!("failed to read package file {}: {err}", relative.display()))?;
        files.push(ExtractedPackageFile {
            path: relative,
            contents,
        });
    }
    Ok(files)
}

fn load_package_descriptor(files: &[ExtractedPackageFile]) -> Result<PackageDescriptor, String> {
    for name in [PACKAGE_DESCRIPTOR, LEGACY_PACKAGE_DESCRIPTOR] {
        if let Some(file) = files.iter().find(|file| file.path == Path::new(name)) {
            let value: serde_yaml::Value = serde_yaml::from_str(&file.contents)
                .map_err(|err| format!("failed to parse package descriptor {name}: {err}"))?;
            let schema = value
                .as_mapping()
                .and_then(|mapping| mapping.get(serde_yaml::Value::String("schema".to_owned())))
                .and_then(serde_yaml::Value::as_str)
                .ok_or_else(|| format!("package descriptor {name} is missing schema"))?;
            if matches!(
                schema,
                "bioscript:panel:1.0"
                    | "bioscript:assay:1.0"
                    | "bioscript:variant:1.0"
                    | "bioscript:variant"
            ) {
                let package_name = value
                    .as_mapping()
                    .and_then(|mapping| mapping.get(serde_yaml::Value::String("name".to_owned())))
                    .and_then(serde_yaml::Value::as_str)
                    .map(ToOwned::to_owned);
                return Ok(PackageDescriptor {
                    entrypoint: PathBuf::from(name),
                    name: package_name,
                });
            }
            if schema != "bioscript:package:1.0" {
                return Err(format!(
                    "package descriptor {name} has unsupported schema '{schema}'"
                ));
            }
            let entrypoint = value
                .as_mapping()
                .and_then(|mapping| mapping.get(serde_yaml::Value::String("entrypoint".to_owned())))
                .and_then(serde_yaml::Value::as_str)
                .ok_or_else(|| format!("package descriptor {name} is missing entrypoint"))?;
            let name = value
                .as_mapping()
                .and_then(|mapping| mapping.get(serde_yaml::Value::String("name".to_owned())))
                .and_then(serde_yaml::Value::as_str)
                .map(ToOwned::to_owned);
            return Ok(PackageDescriptor {
                entrypoint: checked_relative_package_path(entrypoint)?,
                name,
            });
        }
    }
    for candidate in ["panel.yaml", "assay.yaml", "variant.yaml"] {
        if files.iter().any(|file| file.path == Path::new(candidate)) {
            return Ok(PackageDescriptor {
                entrypoint: PathBuf::from(candidate),
                name: None,
            });
        }
    }
    Err(format!(
        "package does not contain {PACKAGE_DESCRIPTOR}, {LEGACY_PACKAGE_DESCRIPTOR}, panel.yaml, assay.yaml, or variant.yaml"
    ))
}

fn checked_relative_package_path(raw: &str) -> Result<PathBuf, String> {
    let path = Path::new(raw);
    if path.is_absolute() {
        return Err(format!("package path must be relative: {raw}"));
    }
    let mut out = PathBuf::new();
    for component in path.components() {
        match component {
            std::path::Component::Normal(part) => out.push(part),
            std::path::Component::CurDir => {}
            std::path::Component::ParentDir
            | std::path::Component::RootDir
            | std::path::Component::Prefix(_) => {
                return Err(format!("package path escapes package root: {raw}"));
            }
        }
    }
    if out.as_os_str().is_empty() {
        return Err("package path is empty".to_owned());
    }
    Ok(out)
}

fn is_allowed_package_file(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| name == PACKAGE_DESCRIPTOR || name == LEGACY_PACKAGE_DESCRIPTOR)
        || path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| {
                matches!(
                    ext.to_ascii_lowercase().as_str(),
                    "yaml" | "yml" | "py" | "md" | "txt" | "tsv" | "json" | "jsonl"
                )
            })
}

fn is_resource_file(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| matches!(ext.to_ascii_lowercase().as_str(), "yaml" | "yml" | "py"))
}

fn package_member_url(source_url: &str, path: &Path) -> String {
    format!(
        "{}/{}",
        source_url.trim_end_matches('/'),
        path.display().to_string().replace('\\', "/")
    )
}

fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String(key.to_owned())))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

fn scalar_at(value: &serde_yaml::Value, key: &str) -> Option<String> {
    yaml_string(value, key)
}

fn join_url(base_url: &str, relative: &str) -> String {
    if relative.starts_with("https://") || relative.starts_with("http://") {
        return relative.to_owned();
    }
    let base = base_url.split('?').next().unwrap_or(base_url);
    match base.rsplit_once('/') {
        Some((prefix, _)) => format!("{prefix}/{relative}"),
        None => relative.to_owned(),
    }
}

fn sha256_hex(bytes: &[u8]) -> String {
    let mut digest = Sha256::new();
    digest.update(bytes);
    format!("{:x}", digest.finalize())
}
