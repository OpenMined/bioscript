use std::{
    fmt::{self, Write as _},
    fs,
    path::{Path, PathBuf},
};

use serde_yaml::Value;

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

/// Validate a variant file or directory of variant files.
///
/// # Errors
///
/// Returns an error when the input path cannot be read, traversed, or parsed
/// as YAML.
pub fn validate_variants_path(path: &Path) -> Result<ValidationReport, String> {
    let files = collect_variant_files(path)?;
    let mut reports = Vec::new();
    for file in &files {
        let report = validate_variant_file(file)?;
        if !report.issues.is_empty() {
            reports.push(report);
        }
    }
    Ok(ValidationReport {
        files_scanned: files.len(),
        reports,
    })
}

fn collect_variant_files(path: &Path) -> Result<Vec<PathBuf>, String> {
    if path.is_file() {
        return Ok(vec![path.to_path_buf()]);
    }

    let mut files = Vec::new();
    collect_variant_files_recursive(path, &mut files)?;
    files.sort();
    Ok(files)
}

fn collect_variant_files_recursive(path: &Path, files: &mut Vec<PathBuf>) -> Result<(), String> {
    let entries = fs::read_dir(path)
        .map_err(|err| format!("failed to read directory {}: {err}", path.display()))?;
    for entry in entries {
        let entry = entry.map_err(|err| format!("failed to read directory entry: {err}"))?;
        let entry_path = entry.path();
        if entry_path.is_dir() {
            collect_variant_files_recursive(&entry_path, files)?;
            continue;
        }
        let Some(file_name) = entry_path.file_name().and_then(|name| name.to_str()) else {
            continue;
        };
        if matches!(file_name, "variant.yaml" | "variant.yml") {
            files.push(entry_path);
        }
    }
    Ok(())
}

fn validate_variant_file(path: &Path) -> Result<FileReport, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    let value: Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))?;

    let mut issues = Vec::new();
    validate_required_shape(&value, &mut issues);
    validate_kind_vs_tags(&value, &mut issues);
    validate_pgx_shape(&value, &mut issues);

    Ok(FileReport {
        file: path.to_path_buf(),
        issues,
    })
}

fn validate_required_shape(root: &Value, issues: &mut Vec<Issue>) {
    require_const(root, &["schema"], "bioscript:variant", issues);
    require_const(root, &["version"], "1.0", issues);
    require_path(root, &["variant_id"], issues);
    require_path(root, &["alleles"], issues);
    require_path(root, &["alleles", "kind"], issues);
    require_path(root, &["alleles", "ref"], issues);
    if value_at(root, &["alleles", "alts"]).is_none() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles".to_owned(),
            message: "missing allele definition; expected alleles.alts".to_owned(),
        });
    }

    if value_at(root, &["alleles", "alt"]).is_some() {
        issues.push(Issue {
            severity: Severity::Warning,
            path: "alleles.alt".to_owned(),
            message: "alleles.alt is legacy shape; prefer alleles.alts and optional alleles.canonical_alt".to_owned(),
        });
    }

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
    if let Some(canonical_alt) = scalar_at(root, &["alleles", "canonical_alt"]) {
        let alts = seq_at(root, &["alleles", "alts"]).unwrap_or_default();
        if !alts.iter().any(|alt| alt == &canonical_alt) {
            issues.push(Issue {
                severity: Severity::Error,
                path: "alleles.canonical_alt".to_owned(),
                message: format!(
                    "canonical_alt '{canonical_alt}' is not present in alleles.alts {alts:?}"
                ),
            });
        }
    }
}

fn validate_kind_vs_tags(root: &Value, issues: &mut Vec<Issue>) {
    let Some(kind) = scalar_at(root, &["alleles", "kind"]) else {
        return;
    };
    let Some(tags) = seq_at(root, &["research", "tags"]) else {
        return;
    };

    let has_legacy_snp_tag = tags.iter().any(|tag| tag == "snp");
    let has_preferred_snv_tag = tags.iter().any(|tag| tag == "snv");
    if kind == "snv" && has_legacy_snp_tag && !has_preferred_snv_tag {
        issues.push(Issue {
            severity: Severity::Warning,
            path: "research.tags".to_owned(),
            message: "alleles.kind is 'snv' but research.tags uses 'snp'; pick one vocabulary and use it consistently".to_owned(),
        });
    }
}

fn validate_pgx_shape(root: &Value, issues: &mut Vec<Issue>) {
    let Some(pgx) = mapping_at(root, &["clinical", "pgx"]) else {
        return;
    };

    for key in ["drug_labels", "annotations", "clinical_annotations"] {
        let Some(items) = pgx
            .get(Value::String(key.to_owned()))
            .and_then(Value::as_sequence)
        else {
            continue;
        };
        for (idx, item) in items.iter().enumerate() {
            let Some(mapping) = item.as_mapping() else {
                issues.push(Issue {
                    severity: Severity::Warning,
                    path: format!("clinical.pgx.{key}[{idx}]"),
                    message: "expected mapping".to_owned(),
                });
                continue;
            };

            if let Some(level) = mapping
                .get(Value::String("pgx_level".to_owned()))
                .and_then(Value::as_str)
                && level.trim().is_empty()
            {
                issues.push(Issue {
                    severity: Severity::Warning,
                    path: format!("clinical.pgx.{key}[{idx}].pgx_level"),
                    message: "empty pgx_level string; prefer null/omitted or a normalized controlled value".to_owned(),
                });
            }
        }
    }
}

fn require_const(root: &Value, path: &[&str], expected: &str, issues: &mut Vec<Issue>) {
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

fn require_path(root: &Value, path: &[&str], issues: &mut Vec<Issue>) {
    if value_at(root, path).is_none() {
        issues.push(Issue {
            severity: Severity::Error,
            path: path.join("."),
            message: "missing required field".to_owned(),
        });
    }
}

fn value_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut current = root;
    for key in path {
        let mapping = current.as_mapping()?;
        current = mapping.get(Value::String((*key).to_owned()))?;
    }
    Some(current)
}

fn mapping_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a serde_yaml::Mapping> {
    value_at(root, path)?.as_mapping()
}

fn scalar_at(root: &Value, path: &[&str]) -> Option<String> {
    value_at(root, path).and_then(|value| match value {
        Value::String(text) => Some(text.clone()),
        Value::Number(number) => Some(number.to_string()),
        _ => None,
    })
}

fn seq_at(root: &Value, path: &[&str]) -> Option<Vec<String>> {
    value_at(root, path)?.as_sequence().map(|items| {
        items
            .iter()
            .filter_map(|item| item.as_str().map(ToOwned::to_owned))
            .collect()
    })
}
