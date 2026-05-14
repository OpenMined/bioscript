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
/// Validate a panel file or directory of panel files.
///
/// # Errors
///
/// Returns an error when the input path cannot be read, traversed, or parsed
/// as YAML.
pub fn validate_panels_path(path: &Path) -> Result<ValidationReport, String> {
    validate_manifest_path(path, ManifestSelector::Panel)
}

/// Validate an assay file or directory of assay files.
///
/// # Errors
///
/// Returns an error when the input path cannot be read, traversed, or parsed
/// as YAML.
/// Validate an assay file or directory of assay files.
///
/// # Errors
///
/// Returns an error when the input path cannot be read, traversed, or parsed
/// as YAML.
pub fn validate_assays_path(path: &Path) -> Result<ValidationReport, String> {
    validate_manifest_path(path, ManifestSelector::Assay)
}

/// Load a single variant manifest from YAML.
///
/// # Errors
///
/// Returns an error when the file does not parse or is not a valid variant
/// manifest.
/// Load a variant manifest from a YAML file.
///
/// # Errors
///
/// Returns an error when the file cannot be read, parsed, or converted into a
/// valid variant manifest shape.
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
/// Load a variant manifest from YAML text.
///
/// # Errors
///
/// Returns an error when the text cannot be parsed or converted into a valid
/// variant manifest shape.
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
/// Load a panel manifest from a YAML file.
///
/// # Errors
///
/// Returns an error when the file cannot be read, parsed, or converted into a
/// valid panel manifest shape.
pub fn load_panel_manifest(path: &Path) -> Result<PanelManifest, String> {
    let value = load_yaml(path)?;
    panel_manifest_from_root(path, &value)
}

/// Load a panel manifest from YAML text.
///
/// # Errors
///
/// Returns an error when the text cannot be parsed or converted into a valid
/// panel manifest shape.
pub fn load_panel_manifest_text(name: &str, text: &str) -> Result<PanelManifest, String> {
    let value: Value =
        serde_yaml::from_str(text).map_err(|err| format!("failed to parse YAML {name}: {err}"))?;
    panel_manifest_from_root(Path::new(name), &value)
}

fn panel_manifest_from_root(path: &Path, value: &Value) -> Result<PanelManifest, String> {
    let mut issues = Vec::new();
    validate_panel_root(value, &mut issues);
    if issues.iter().any(|issue| issue.severity == Severity::Error) {
        return Err(render_single_manifest_errors(path, &issues));
    }

    let permissions = Permissions {
        domains: seq_of_strings(value, &["permissions", "domains"]).unwrap_or_default(),
    };
    let downloads = parse_downloads(value)?;
    let members = parse_panel_members(value)?;
    let interpretations = parse_panel_interpretations(value)?;

    Ok(PanelManifest {
        path: path.to_path_buf(),
        name: required_non_empty_string(value, &["name"])?,
        label: scalar_at(value, &["label"]),
        tags: seq_of_strings(value, &["tags"]).unwrap_or_default(),
        permissions,
        downloads,
        members,
        interpretations,
    })
}

/// Load a single assay manifest from YAML.
///
/// # Errors
///
/// Returns an error when the file does not parse or is not a valid assay
/// manifest.
/// Load an assay manifest from a YAML file.
///
/// # Errors
///
/// Returns an error when the file cannot be read, parsed, or converted into a
/// valid assay manifest shape.
pub fn load_assay_manifest(path: &Path) -> Result<AssayManifest, String> {
    let value = load_yaml(path)?;
    assay_manifest_from_root(path, &value)
}

/// Load an assay manifest from YAML text.
///
/// # Errors
///
/// Returns an error when the text cannot be parsed or converted into a valid
/// assay manifest shape.
pub fn load_assay_manifest_text(name: &str, text: &str) -> Result<AssayManifest, String> {
    let value: Value =
        serde_yaml::from_str(text).map_err(|err| format!("failed to parse YAML {name}: {err}"))?;
    assay_manifest_from_root(Path::new(name), &value)
}

fn assay_manifest_from_root(path: &Path, value: &Value) -> Result<AssayManifest, String> {
    let mut issues = Vec::new();
    validate_assay_root(value, &mut issues);
    if issues.iter().any(|issue| issue.severity == Severity::Error) {
        return Err(render_single_manifest_errors(path, &issues));
    }

    Ok(AssayManifest {
        path: path.to_path_buf(),
        name: required_non_empty_string(value, &["name"])?,
        tags: seq_of_strings(value, &["tags"]).unwrap_or_default(),
        members: parse_panel_members(value)?,
        interpretations: parse_panel_interpretations(value)?,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ManifestSelector {
    Assay,
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
            ManifestSelector::Assay => validate_assay_file(file)?,
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

fn validate_assay_file(path: &Path) -> Result<FileReport, String> {
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
    if !schema.contains("assay") {
        return Ok(FileReport {
            file: path.to_path_buf(),
            issues: Vec::new(),
        });
    }

    let mut issues = Vec::new();
    validate_assay_root(&value, &mut issues);
    Ok(FileReport {
        file: path.to_path_buf(),
        issues,
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
        if schema == "bioscript:pgx-findings:1.0" {
            let mut issues = Vec::new();
            validate_pgx_findings_root(&value, &mut issues);
            return Ok(FileReport {
                file: path.to_path_buf(),
                issues,
            });
        }
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

#[cfg(test)]
mod load_validator_tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "bioscript-schema-load-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn variant_yaml(name: &str) -> String {
        format!(
            r#"
schema: bioscript:variant:1.0
version: "1.0"
name: {name}
tags: [tag:test]
identifiers:
  rsids: [rs1]
coordinates:
  grch38:
    chrom: "1"
    pos: 100
alleles:
  kind: snv
  ref: A
  alts: [G]
"#
        )
    }

    #[test]
    fn manifest_loaders_parse_variant_panel_and_assay_text() {
        let variant = load_variant_manifest_text("variant.yaml", &variant_yaml("rs1")).unwrap();
        assert_eq!(variant.name, "rs1");
        assert_eq!(variant.tags, vec!["tag:test"]);

        let lookup = load_variant_manifest_text_for_lookup(
            "legacy.yaml",
            &variant_yaml("rs1").replace("bioscript:variant:1.0", "bioscript:variant"),
        )
        .unwrap();
        assert_eq!(lookup.spec.rsids, vec!["rs1"]);

        let panel = load_panel_manifest_text(
            "panel.yaml",
            r#"
schema: bioscript:panel:1.0
version: "1.0"
name: panel
label: Panel
permissions:
  domains: [https://example.test]
downloads:
  - id: dl
    url: https://example.test/file.yaml
    sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    version: "1"
members:
  - kind: variant
    path: variant.yaml
analyses:
  - id: a
    kind: bioscript
    path: analysis.bs
    output_format: json
    derived_from: [variant.yaml]
"#,
        )
        .unwrap();
        assert_eq!(panel.downloads.len(), 1);
        assert_eq!(panel.members.len(), 1);
        assert_eq!(panel.interpretations.len(), 1);

        let assay = load_assay_manifest_text(
            "assay.yaml",
            r#"
schema: bioscript:assay:1.0
version: "1.0"
name: assay
members:
  - kind: variant
    path: variant.yaml
"#,
        )
        .unwrap();
        assert_eq!(assay.members.len(), 1);
    }

    #[test]
    fn manifest_loaders_report_parse_and_validation_errors() {
        assert!(load_variant_manifest_text("bad.yaml", "{")
            .unwrap_err()
            .contains("failed to parse YAML"));
        assert!(load_variant_manifest_text("bad.yaml", "schema: bioscript:variant:1.0\n")
            .unwrap_err()
            .contains("missing required field"));
        assert!(load_variant_manifest_text_for_lookup("bad.yaml", "schema: bad\n")
            .unwrap_err()
            .contains("expected schema"));
        assert!(load_panel_manifest_text("bad.yaml", "schema: bad\n")
            .unwrap_err()
            .contains("expected schema"));
        assert!(load_assay_manifest_text("bad.yaml", "schema: bad\n")
            .unwrap_err()
            .contains("expected schema"));
    }

    #[test]
    fn validate_manifest_path_collects_yaml_files_recursively_and_ignores_other_schemas() {
        let dir = temp_dir("collect");
        fs::write(dir.join("variant.yaml"), variant_yaml("rs1")).unwrap();
        fs::write(dir.join("panel.yml"), "schema: bioscript:panel:1.0\n").unwrap();
        fs::write(dir.join("notes.txt"), "ignored").unwrap();
        fs::create_dir_all(dir.join("nested")).unwrap();
        fs::write(dir.join("nested/missing-schema.yaml"), "name: missing\n").unwrap();

        let files = collect_yaml_files(&dir).unwrap();
        assert_eq!(files.len(), 3);

        let report = validate_variants_path(&dir).unwrap();
        assert_eq!(report.files_scanned, 3);
        assert_eq!(report.total_errors(), 1);
        assert!(report.render_text().contains("missing schema"));

        let panel_report = validate_panels_path(&dir).unwrap();
        assert_eq!(panel_report.files_scanned, 3);
        assert!(panel_report.total_errors() >= 1);

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn validate_file_helpers_skip_unrelated_schemas_and_handle_missing_schema() {
        let dir = temp_dir("files");
        let missing = dir.join("missing.yaml");
        fs::write(&missing, "name: missing\n").unwrap();
        assert_eq!(validate_variant_file(&missing).unwrap().issues[0].path, "schema");
        assert_eq!(validate_panel_file(&missing).unwrap().issues[0].path, "schema");
        assert_eq!(validate_assay_file(&missing).unwrap().issues[0].path, "schema");

        let panel = dir.join("panel.yaml");
        fs::write(&panel, "schema: bioscript:panel:1.0\n").unwrap();
        assert!(validate_variant_file(&panel).unwrap().issues.is_empty());

        let pgx = dir.join("pgx.yaml");
        fs::write(
            &pgx,
            r#"
schema: bioscript:pgx-findings:1.0
version: "1.0"
findings: []
"#,
        )
        .unwrap();
        assert!(!validate_variant_file(&pgx).unwrap().issues.is_empty());

        fs::remove_dir_all(dir).unwrap();
    }
}
