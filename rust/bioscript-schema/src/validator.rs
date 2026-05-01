use std::{
    fmt::{self, Write as _},
    path::{Path, PathBuf},
};

use bioscript_core::VariantSpec;
use serde_yaml::Value;

mod common;
mod panel;
mod spec;
mod variant;

use common::{
    collect_yaml_files, load_yaml, render_single_manifest_errors, required_non_empty_string,
    scalar_at, seq_of_strings, validate_schema_and_identity,
};
use panel::{parse_downloads, parse_panel_members, validate_panel_root};
use spec::variant_spec_from_root;
use variant::{
    validate_alleles, validate_coordinates, validate_identifiers, validate_variant_root,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
}

impl fmt::Display for Severity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Error => f.write_str("error"),
            Self::Warning => f.write_str("warning"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Issue {
    pub severity: Severity,
    pub path: String,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileReport {
    pub file: PathBuf,
    pub issues: Vec<Issue>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationReport {
    pub files_scanned: usize,
    pub reports: Vec<FileReport>,
}

impl ValidationReport {
    #[must_use]
    pub fn total_issues(&self) -> usize {
        self.reports.iter().map(|report| report.issues.len()).sum()
    }

    #[must_use]
    pub fn total_errors(&self) -> usize {
        self.reports
            .iter()
            .flat_map(|report| &report.issues)
            .filter(|issue| issue.severity == Severity::Error)
            .count()
    }

    #[must_use]
    pub fn total_warnings(&self) -> usize {
        self.reports
            .iter()
            .flat_map(|report| &report.issues)
            .filter(|issue| issue.severity == Severity::Warning)
            .count()
    }

    #[must_use]
    pub fn has_errors(&self) -> bool {
        self.total_errors() > 0
    }

    #[must_use]
    pub fn render_text(&self) -> String {
        let mut out = String::new();
        let _ = write!(
            out,
            "files_scanned: {}\nerrors: {}\nwarnings: {}\n",
            self.files_scanned,
            self.total_errors(),
            self.total_warnings()
        );
        for report in &self.reports {
            out.push('\n');
            let _ = writeln!(out, "file: {}", report.file.display());
            for issue in &report.issues {
                let _ = writeln!(
                    out,
                    "  - [{}] {}: {}",
                    issue.severity, issue.path, issue.message
                );
            }
        }
        out
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantManifest {
    pub path: PathBuf,
    pub name: String,
    pub tags: Vec<String>,
    pub spec: VariantSpec,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelManifest {
    pub path: PathBuf,
    pub name: String,
    pub tags: Vec<String>,
    pub permissions: Permissions,
    pub downloads: Vec<Download>,
    pub members: Vec<PanelMember>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Permissions {
    pub domains: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Download {
    pub id: String,
    pub url: String,
    pub origin: String,
    pub sha256: String,
    pub version: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelMember {
    pub kind: String,
    pub path: Option<String>,
    pub download: Option<String>,
    pub sha256: Option<String>,
    pub version: Option<String>,
}

/// Validate a variant file or directory of variant files.
///
/// # Errors
///
/// Returns an error when the input path cannot be read, traversed, or parsed
/// as YAML.
pub fn validate_variants_path(path: &Path) -> Result<ValidationReport, String> {
    validate_manifest_path(path, ManifestSelector::Variant)
}

/// Validate a panel file or directory of panel files.
///
/// # Errors
///
/// Returns an error when the input path cannot be read, traversed, or parsed
/// as YAML.
pub fn validate_panels_path(path: &Path) -> Result<ValidationReport, String> {
    validate_manifest_path(path, ManifestSelector::Panel)
}

/// Load a single variant manifest from YAML.
///
/// # Errors
///
/// Returns an error when the file does not parse or is not a valid variant
/// manifest.
pub fn load_variant_manifest(path: &Path) -> Result<VariantManifest, String> {
    let value = load_yaml(path)?;
    variant_manifest_from_root(path, &value)
}

/// Load a single variant manifest from YAML text.
///
/// # Errors
///
/// Returns an error when the text does not parse or is not a valid variant
/// manifest.
pub fn load_variant_manifest_text(name: &str, text: &str) -> Result<VariantManifest, String> {
    let value: Value =
        serde_yaml::from_str(text).map_err(|err| format!("failed to parse YAML {name}: {err}"))?;
    variant_manifest_from_root(Path::new(name), &value)
}

/// Compile a variant manifest from YAML text for lookup execution.
///
/// This validates the execution-critical fields only: identity, identifiers,
/// coordinates, and alleles. Full manifest validation still reports metadata
/// issues such as missing finding schemas, but those do not block local lookup.
///
/// # Errors
///
/// Returns an error when the text does not parse or the execution-critical
/// fields are invalid.
pub fn load_variant_manifest_text_for_lookup(
    name: &str,
    text: &str,
) -> Result<VariantManifest, String> {
    let value: Value =
        serde_yaml::from_str(text).map_err(|err| format!("failed to parse YAML {name}: {err}"))?;
    let path = Path::new(name);
    let mut issues = Vec::new();
    validate_schema_and_identity(
        &value,
        "bioscript:variant:1.0",
        Some("bioscript:variant"),
        &mut issues,
    );
    validate_identifiers(&value, &mut issues);
    validate_coordinates(&value, &mut issues);
    validate_alleles(&value, &mut issues);
    if issues.iter().any(|issue| issue.severity == Severity::Error) {
        return Err(render_single_manifest_errors(path, &issues));
    }

    Ok(VariantManifest {
        path: path.to_path_buf(),
        name: required_non_empty_string(&value, &["name"])?,
        tags: seq_of_strings(&value, &["tags"]).unwrap_or_default(),
        spec: variant_spec_from_root(&value)?,
    })
}

fn variant_manifest_from_root(path: &Path, value: &Value) -> Result<VariantManifest, String> {
    let mut issues = Vec::new();
    validate_variant_root(value, &mut issues);
    if issues.iter().any(|issue| issue.severity == Severity::Error) {
        return Err(render_single_manifest_errors(path, &issues));
    }

    Ok(VariantManifest {
        path: path.to_path_buf(),
        name: required_non_empty_string(value, &["name"])?,
        tags: seq_of_strings(value, &["tags"]).unwrap_or_default(),
        spec: variant_spec_from_root(value)?,
    })
}

/// Load a single panel manifest from YAML.
///
/// # Errors
///
/// Returns an error when the file does not parse or is not a valid panel
/// manifest.
pub fn load_panel_manifest(path: &Path) -> Result<PanelManifest, String> {
    let value = load_yaml(path)?;
    let mut issues = Vec::new();
    validate_panel_root(&value, &mut issues);
    if issues.iter().any(|issue| issue.severity == Severity::Error) {
        return Err(render_single_manifest_errors(path, &issues));
    }

    let permissions = Permissions {
        domains: seq_of_strings(&value, &["permissions", "domains"]).unwrap_or_default(),
    };
    let downloads = parse_downloads(&value)?;
    let members = parse_panel_members(&value)?;

    Ok(PanelManifest {
        path: path.to_path_buf(),
        name: required_non_empty_string(&value, &["name"])?,
        tags: seq_of_strings(&value, &["tags"]).unwrap_or_default(),
        permissions,
        downloads,
        members,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ManifestSelector {
    Variant,
    Panel,
}

fn validate_manifest_path(
    path: &Path,
    selector: ManifestSelector,
) -> Result<ValidationReport, String> {
    let files = collect_yaml_files(path)?;
    let mut reports = Vec::new();
    for file in &files {
        let report = match selector {
            ManifestSelector::Variant => validate_variant_file(file)?,
            ManifestSelector::Panel => validate_panel_file(file)?,
        };
        if !report.issues.is_empty() {
            reports.push(report);
        }
    }
    Ok(ValidationReport {
        files_scanned: files.len(),
        reports,
    })
}

fn validate_variant_file(path: &Path) -> Result<FileReport, String> {
    let value = load_yaml(path)?;
    let Some(schema) = scalar_at(&value, &["schema"]) else {
        return Ok(FileReport {
            file: path.to_path_buf(),
            issues: vec![Issue {
                severity: Severity::Error,
                path: "schema".to_owned(),
                message: "missing schema".to_owned(),
            }],
        });
    };
    if !schema.contains("variant") {
        return Ok(FileReport {
            file: path.to_path_buf(),
            issues: Vec::new(),
        });
    }

    let mut issues = Vec::new();
    validate_variant_root(&value, &mut issues);
    Ok(FileReport {
        file: path.to_path_buf(),
        issues,
    })
}

fn validate_panel_file(path: &Path) -> Result<FileReport, String> {
    let value = load_yaml(path)?;
    let Some(schema) = scalar_at(&value, &["schema"]) else {
        return Ok(FileReport {
            file: path.to_path_buf(),
            issues: vec![Issue {
                severity: Severity::Error,
                path: "schema".to_owned(),
                message: "missing schema".to_owned(),
            }],
        });
    };
    if !schema.contains("panel") {
        return Ok(FileReport {
            file: path.to_path_buf(),
            issues: Vec::new(),
        });
    }

    let mut issues = Vec::new();
    validate_panel_root(&value, &mut issues);
    Ok(FileReport {
        file: path.to_path_buf(),
        issues,
    })
}
