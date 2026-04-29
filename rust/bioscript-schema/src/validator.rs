use std::{
    collections::BTreeSet,
    fmt::{self, Write as _},
    fs,
    path::{Path, PathBuf},
};

use bioscript_core::{GenomicLocus, VariantKind, VariantSpec};
use serde_yaml::{Mapping, Value};
use url::Url;

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

fn collect_yaml_files(path: &Path) -> Result<Vec<PathBuf>, String> {
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

fn validate_variant_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(
        root,
        "bioscript:variant:1.0",
        Some("bioscript:variant"),
        issues,
    );
    validate_optional_strings(root, &["name", "label", "gene", "summary"], issues);
    validate_tags(root, issues);
    validate_identifiers(root, issues);
    validate_coordinates(root, issues);
    validate_alleles(root, issues);
    validate_findings(root, issues);
    validate_provenance(root, issues);

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
}

fn validate_panel_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(root, "bioscript:panel:1.0", None, issues);
    validate_optional_strings(root, &["name", "label", "summary"], issues);
    validate_tags(root, issues);
    validate_permissions(root, issues);
    validate_downloads(root, issues);
    validate_panel_members(root, issues);
}

fn validate_schema_and_identity(
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

fn validate_optional_strings(root: &Value, fields: &[&str], issues: &mut Vec<Issue>) {
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

fn validate_tags(root: &Value, issues: &mut Vec<Issue>) {
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

fn validate_identifiers(root: &Value, issues: &mut Vec<Issue>) {
    for field in ["rsids", "aliases"] {
        let Some(values) = value_at(root, &["identifiers", field]) else {
            continue;
        };
        let Some(items) = values.as_sequence() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("identifiers.{field}"),
                message: "expected a sequence of strings".to_owned(),
            });
            continue;
        };
        let mut seen = BTreeSet::new();
        for (idx, item) in items.iter().enumerate() {
            let Some(value) = item.as_str() else {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: "expected string".to_owned(),
                });
                continue;
            };
            if !is_rsid(value) {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: format!("expected rsid like rs123, found '{value}'"),
                });
            }
            if !seen.insert(value.to_owned()) {
                issues.push(Issue {
                    severity: Severity::Warning,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: format!("duplicate identifier '{value}'"),
                });
            }
        }
    }
}

fn validate_coordinates(root: &Value, issues: &mut Vec<Issue>) {
    for assembly in ["grch37", "grch38"] {
        let Some(coord) = mapping_at(root, &["coordinates", assembly]) else {
            continue;
        };

        let Some(chrom) = coord
            .get(Value::String("chrom".to_owned()))
            .and_then(Value::as_str)
        else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.chrom"),
                message: "missing chrom".to_owned(),
            });
            continue;
        };
        if !is_allowed_chromosome(chrom) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.chrom"),
                message: format!("invalid chromosome '{chrom}'; expected 1-22, X, Y, or MT"),
            });
        }

        let has_pos = coord.contains_key(Value::String("pos".to_owned()));
        let has_start = coord.contains_key(Value::String("start".to_owned()));
        let has_end = coord.contains_key(Value::String("end".to_owned()));
        if has_pos && (has_start || has_end) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}"),
                message: "use either pos or start/end, not both".to_owned(),
            });
            continue;
        }
        if !(has_pos || has_start && has_end) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}"),
                message: "expected either pos or start/end".to_owned(),
            });
            continue;
        }

        if has_pos {
            validate_coordinate_pos(coord, assembly, issues);
        } else {
            validate_coordinate_range(coord, assembly, issues);
        }
    }
}

fn validate_coordinate_pos(coord: &Mapping, assembly: &str, issues: &mut Vec<Issue>) {
    if let Some(pos) = i64_at_mapping(coord, "pos") {
        if pos < 1 {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.pos"),
                message: "expected integer >= 1".to_owned(),
            });
        }
    } else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.pos"),
            message: "expected integer".to_owned(),
        });
    }
}

fn validate_coordinate_range(coord: &Mapping, assembly: &str, issues: &mut Vec<Issue>) {
    let start = i64_at_mapping(coord, "start");
    let end = i64_at_mapping(coord, "end");
    match (start, end) {
        (Some(start), Some(end)) => validate_coordinate_range_values(start, end, assembly, issues),
        _ => issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}"),
            message: "expected integer start/end".to_owned(),
        }),
    }
}

fn validate_coordinate_range_values(start: i64, end: i64, assembly: &str, issues: &mut Vec<Issue>) {
    if start < 1 {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.start"),
            message: "expected integer >= 1".to_owned(),
        });
    }
    if end < 1 {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.end"),
            message: "expected integer >= 1".to_owned(),
        });
    }
    if end < start {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.end"),
            message: "expected end >= start".to_owned(),
        });
    }
    if start == end {
        issues.push(Issue {
            severity: Severity::Warning,
            path: format!("coordinates.{assembly}"),
            message: "single-position coordinate uses start/end; prefer pos".to_owned(),
        });
    }
}

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
    validate_symbolic_alleles(&reference, &alts, issues);
    validate_snv_alleles(&kind, &reference, &alts, issues);
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
        if let Some(alt) = mapping
            .get(Value::String("alt".to_owned()))
            .and_then(Value::as_str)
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
    }
}

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

fn variant_spec_from_root(root: &Value) -> Result<VariantSpec, String> {
    let rsids = seq_of_strings(root, &["identifiers", "rsids"]).unwrap_or_default();
    let grch37 = locus_from_root(root, "grch37")?;
    let grch38 = locus_from_root(root, "grch38")?;
    let reference = scalar_at(root, &["alleles", "ref"]);
    let alternate =
        seq_of_strings(root, &["alleles", "alts"]).and_then(|alts| alts.first().cloned());
    let deletion_length = value_at(root, &["alleles", "deletion_length"])
        .and_then(Value::as_u64)
        .and_then(|value| usize::try_from(value).ok());
    let motifs = seq_of_strings(root, &["alleles", "motifs"]).unwrap_or_default();
    let kind = scalar_at(root, &["alleles", "kind"]).map(|kind| match kind.as_str() {
        "snv" => VariantKind::Snp,
        "deletion" => VariantKind::Deletion,
        "insertion" => VariantKind::Insertion,
        "indel" => VariantKind::Indel,
        _ => VariantKind::Other,
    });

    Ok(VariantSpec {
        rsids,
        grch37,
        grch38,
        reference,
        alternate,
        kind,
        deletion_length,
        motifs,
    })
}

fn locus_from_root(root: &Value, assembly: &str) -> Result<Option<GenomicLocus>, String> {
    let Some(mapping) = mapping_at(root, &["coordinates", assembly]) else {
        return Ok(None);
    };
    let chrom = mapping
        .get(Value::String("chrom".to_owned()))
        .and_then(Value::as_str)
        .ok_or_else(|| format!("coordinates.{assembly}.chrom missing"))?;
    let (start, end) = if let Some(pos) = i64_at_mapping(mapping, "pos") {
        (pos, pos)
    } else {
        let start = i64_at_mapping(mapping, "start")
            .ok_or_else(|| format!("coordinates.{assembly}.start missing"))?;
        let end = i64_at_mapping(mapping, "end")
            .ok_or_else(|| format!("coordinates.{assembly}.end missing"))?;
        (start, end)
    };
    Ok(Some(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    }))
}

fn parse_downloads(root: &Value) -> Result<Vec<Download>, String> {
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

fn parse_panel_members(root: &Value) -> Result<Vec<PanelMember>, String> {
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

fn mapping_required_string(
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

fn validate_url_string(
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

fn normalize_origin(value: &str) -> Result<String, String> {
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

fn normalize_download_url(value: &str) -> Result<String, String> {
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

fn is_allowed_chromosome(value: &str) -> bool {
    matches!(value, "X" | "Y" | "MT")
        || value
            .parse::<u8>()
            .is_ok_and(|chrom| (1..=22).contains(&chrom))
}

fn is_base_allele(value: &str) -> bool {
    matches!(value, "A" | "C" | "G" | "T")
}

fn is_rsid(value: &str) -> bool {
    value.starts_with("rs") && value[2..].chars().all(|ch| ch.is_ascii_digit())
}

fn is_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .chars()
            .all(|ch| ch.is_ascii_hexdigit() && !ch.is_ascii_uppercase())
}

fn i64_at_mapping(mapping: &Mapping, key: &str) -> Option<i64> {
    mapping
        .get(Value::String(key.to_owned()))
        .and_then(Value::as_i64)
}

fn required_non_empty_string(root: &Value, path: &[&str]) -> Result<String, String> {
    scalar_at(root, path)
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| format!("{} missing or empty", path.join(".")))
}

fn render_single_manifest_errors(path: &Path, issues: &[Issue]) -> String {
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

fn load_yaml(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))
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

fn mapping_at<'a>(root: &'a Value, path: &[&str]) -> Option<&'a Mapping> {
    value_at(root, path)?.as_mapping()
}

fn scalar_at(root: &Value, path: &[&str]) -> Option<String> {
    value_at(root, path).and_then(|value| match value {
        Value::String(text) => Some(text.clone()),
        Value::Number(number) => Some(number.to_string()),
        _ => None,
    })
}

fn seq_of_strings(root: &Value, path: &[&str]) -> Option<Vec<String>> {
    value_at(root, path)?.as_sequence().map(|items| {
        items
            .iter()
            .filter_map(|item| item.as_str().map(ToOwned::to_owned))
            .collect()
    })
}
