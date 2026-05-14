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
    if schema == "bioscript:variant-catalogue:1.0" {
        let mut issues = Vec::new();
        validate_variant_catalogue_root(&value, &mut issues);
        return Ok(FileReport {
            file: path.to_path_buf(),
            issues,
        });
    }
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
