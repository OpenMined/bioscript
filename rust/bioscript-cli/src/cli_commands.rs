fn run_prepare(args: Vec<String>) -> Result<(), String> {
    let mut root: Option<PathBuf> = None;
    let mut input_file: Option<String> = None;
    let mut reference_file: Option<String> = None;
    let mut input_format: Option<GenotypeSourceFormat> = None;
    let mut cache_dir: Option<PathBuf> = None;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--root" => {
                root = Some(PathBuf::from(
                    iter.next().ok_or("--root requires a directory")?,
                ));
            }
            "--input-file" => {
                input_file = Some(iter.next().ok_or("--input-file requires a path")?);
            }
            "--reference-file" => {
                reference_file = Some(iter.next().ok_or("--reference-file requires a path")?);
            }
            "--input-format" => {
                let value = iter.next().ok_or("--input-format requires a value")?;
                if !value.eq_ignore_ascii_case("auto") {
                    input_format = Some(
                        value
                            .parse::<GenotypeSourceFormat>()
                            .map_err(|err| format!("invalid --input-format: {err}"))?,
                    );
                }
            }
            "--cache-dir" => {
                cache_dir = Some(PathBuf::from(
                    iter.next().ok_or("--cache-dir requires a path")?,
                ));
            }
            other => {
                return Err(format!("unexpected argument: {other}"));
            }
        }
    }

    let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
    let effective_root = root.unwrap_or_else(|| cwd.clone());
    let effective_cache = cache_dir.unwrap_or_else(|| cwd.join(".bioscript-cache"));

    let request = PrepareRequest {
        root: effective_root,
        cwd,
        cache_dir: effective_cache,
        input_file,
        input_format,
        reference_file,
    };

    let prepared = prepare_indexes(&request)?;

    // print the flags that should be passed to a subsequent bioscript run
    let flags = shell_flags(&prepared);
    if flags.is_empty() {
        eprintln!("bioscript prepare: nothing to index");
    } else {
        println!("{flags}");
    }

    Ok(())
}

fn run_inspect(args: Vec<String>) -> Result<(), String> {
    let mut path: Option<PathBuf> = None;
    let mut options = InspectOptions::default();

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--input-index" => {
                options.input_index = Some(PathBuf::from(
                    iter.next().ok_or("--input-index requires a path")?,
                ));
            }
            "--reference-file" => {
                options.reference_file = Some(PathBuf::from(
                    iter.next().ok_or("--reference-file requires a path")?,
                ));
            }
            "--reference-index" => {
                options.reference_index = Some(PathBuf::from(
                    iter.next().ok_or("--reference-index requires a path")?,
                ));
            }
            "--detect-sex" => {
                options.detect_sex = true;
            }
            other if path.is_none() => {
                path = Some(PathBuf::from(other));
            }
            other => {
                return Err(format!("unexpected argument: {other}"));
            }
        }
    }

    let Some(path) = path else {
        return Err(
            "usage: bioscript inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>] [--detect-sex]"
                .to_owned(),
        );
    };

    let inspection = inspect_file(&path, &options).map_err(|err| err.to_string())?;
    println!("{}", inspection.render_text());
    Ok(())
}

fn run_validate_variants(args: Vec<String>) -> Result<(), String> {
    let mut path: Option<PathBuf> = None;
    let mut report_path: Option<PathBuf> = None;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        if arg == "--report" {
            let Some(value) = iter.next() else {
                return Err("--report requires a path".to_owned());
            };
            report_path = Some(PathBuf::from(value));
        } else if path.is_none() {
            path = Some(PathBuf::from(arg));
        } else {
            return Err(format!("unexpected argument: {arg}"));
        }
    }

    let Some(path) = path else {
        return Err("usage: bioscript validate-variants <path> [--report <file>]".to_owned());
    };

    let report = validate_variants_path(&path)?;
    let text = report.render_text();
    print!("{text}");

    if let Some(report_path) = report_path {
        if let Some(parent) = report_path.parent() {
            std::fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create report dir {}: {err}", parent.display())
            })?;
        }
        std::fs::write(&report_path, text)
            .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    }

    if report.has_errors() {
        return Err(format!(
            "validation found {} errors and {} warnings",
            report.total_errors(),
            report.total_warnings()
        ));
    }

    Ok(())
}

fn run_validate_panels(args: Vec<String>) -> Result<(), String> {
    let mut path: Option<PathBuf> = None;
    let mut report_path: Option<PathBuf> = None;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        if arg == "--report" {
            let Some(value) = iter.next() else {
                return Err("--report requires a path".to_owned());
            };
            report_path = Some(PathBuf::from(value));
        } else if path.is_none() {
            path = Some(PathBuf::from(arg));
        } else {
            return Err(format!("unexpected argument: {arg}"));
        }
    }

    let Some(path) = path else {
        return Err("usage: bioscript validate-panels <path> [--report <file>]".to_owned());
    };

    let report = validate_panels_path(&path)?;
    let text = report.render_text();
    print!("{text}");

    if let Some(report_path) = report_path {
        if let Some(parent) = report_path.parent() {
            std::fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create report dir {}: {err}", parent.display())
            })?;
        }
        std::fs::write(&report_path, text)
            .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    }

    if report.has_errors() {
        return Err(format!(
            "validation found {} errors and {} warnings",
            report.total_errors(),
            report.total_warnings()
        ));
    }

    Ok(())
}

fn is_yaml_manifest(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| matches!(ext, "yaml" | "yml"))
}

fn run_validate_assays(args: Vec<String>) -> Result<(), String> {
    let mut path: Option<PathBuf> = None;
    let mut report_path: Option<PathBuf> = None;

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        if arg == "--report" {
            let Some(value) = iter.next() else {
                return Err("--report requires a path".to_owned());
            };
            report_path = Some(PathBuf::from(value));
        } else if path.is_none() {
            path = Some(PathBuf::from(arg));
        } else {
            return Err(format!("unexpected argument: {arg}"));
        }
    }

    let Some(path) = path else {
        return Err("usage: bioscript validate-assays <path> [--report <file>]".to_owned());
    };

    let report = validate_assays_path(&path)?;
    let text = report.render_text();
    print!("{text}");

    if let Some(report_path) = report_path {
        if let Some(parent) = report_path.parent() {
            std::fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create report dir {}: {err}", parent.display())
            })?;
        }
        std::fs::write(&report_path, text)
            .map_err(|err| format!("failed to write {}: {err}", report_path.display()))?;
    }

    if report.has_errors() {
        return Err(format!(
            "validation found {} errors and {} warnings",
            report.total_errors(),
            report.total_warnings()
        ));
    }

    Ok(())
}

#[cfg(test)]
mod cli_command_tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-cli-command-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn valid_variant_yaml() -> &'static str {
        r#"
schema: bioscript:variant:1.0
version: "1.0"
name: Test variant
gene: ABC
identifiers:
  rsids: [rs1]
coordinates:
  grch38:
    chrom: "1"
    pos: 100
alleles:
  kind: snv
  ref: G
  alts: [A]
"#
    }

    #[test]
    fn prepare_and_inspect_commands_validate_arguments() {
        assert!(run_prepare(vec!["--root".to_owned()])
            .unwrap_err()
            .contains("--root requires"));
        assert!(run_prepare(vec!["--input-file".to_owned()])
            .unwrap_err()
            .contains("--input-file requires"));
        assert!(run_prepare(vec![
            "--input-format".to_owned(),
            "not-a-format".to_owned(),
        ])
        .unwrap_err()
        .contains("invalid --input-format"));
        assert!(run_prepare(vec!["--unexpected".to_owned()])
            .unwrap_err()
            .contains("unexpected argument"));

        assert!(run_inspect(Vec::new()).unwrap_err().contains("usage"));
        assert!(run_inspect(vec!["a.txt".to_owned(), "b.txt".to_owned()])
            .unwrap_err()
            .contains("unexpected argument"));
        assert!(run_inspect(vec!["sample.cram".to_owned(), "--input-index".to_owned()])
            .unwrap_err()
            .contains("--input-index requires"));
    }

    #[test]
    fn yaml_manifest_extension_matching_is_case_sensitive_by_contract() {
        assert!(is_yaml_manifest(Path::new("panel.yaml")));
        assert!(is_yaml_manifest(Path::new("panel.yml")));
        assert!(!is_yaml_manifest(Path::new("panel.YAML")));
        assert!(!is_yaml_manifest(Path::new("panel.json")));
    }

    #[test]
    fn validate_variants_writes_report_and_surfaces_errors() {
        let dir = temp_dir("variants");
        let valid = dir.join("variant.yaml");
        let report = dir.join("reports/variant.txt");
        fs::write(&valid, valid_variant_yaml()).unwrap();

        run_validate_variants(vec![
            valid.display().to_string(),
            "--report".to_owned(),
            report.display().to_string(),
        ])
        .unwrap();
        assert!(fs::read_to_string(&report).unwrap().contains("files_scanned"));

        let invalid = dir.join("invalid.yaml");
        fs::write(&invalid, "schema: bioscript:variant:1.0\n").unwrap();
        let err = run_validate_variants(vec![invalid.display().to_string()]).unwrap_err();
        assert!(err.contains("validation found"));

        assert!(run_validate_variants(Vec::new()).unwrap_err().contains("usage"));
        assert!(run_validate_variants(vec![valid.display().to_string(), "--report".to_owned()])
            .unwrap_err()
            .contains("--report requires"));
        assert!(run_validate_variants(vec![
            valid.display().to_string(),
            "extra".to_owned(),
        ])
        .unwrap_err()
        .contains("unexpected argument"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn validate_panels_and_assays_cover_report_and_error_paths() {
        let dir = temp_dir("panels-assays");
        let variant = dir.join("variant.yaml");
        fs::write(&variant, valid_variant_yaml()).unwrap();

        let panel = dir.join("panel.yaml");
        fs::write(
            &panel,
            r#"
schema: bioscript:panel:1.0
version: "1.0"
name: Test panel
members:
  - kind: variant
    path: variant.yaml
"#,
        )
        .unwrap();
        let panel_report = dir.join("reports/panel.txt");
        run_validate_panels(vec![
            panel.display().to_string(),
            "--report".to_owned(),
            panel_report.display().to_string(),
        ])
        .unwrap();
        assert!(panel_report.exists());

        let assay = dir.join("assay.yaml");
        fs::write(
            &assay,
            r#"
schema: bioscript:assay:1.0
version: "1.0"
name: Test assay
members:
  - kind: variant
    path: variant.yaml
"#,
        )
        .unwrap();
        let assay_report = dir.join("reports/assay.txt");
        run_validate_assays(vec![
            assay.display().to_string(),
            "--report".to_owned(),
            assay_report.display().to_string(),
        ])
        .unwrap();
        assert!(assay_report.exists());

        assert!(run_validate_panels(Vec::new()).unwrap_err().contains("usage"));
        assert!(run_validate_panels(vec![panel.display().to_string(), "--report".to_owned()])
            .unwrap_err()
            .contains("--report requires"));
        assert!(run_validate_panels(vec![panel.display().to_string(), "extra".to_owned()])
            .unwrap_err()
            .contains("unexpected argument"));
        assert!(run_validate_assays(Vec::new()).unwrap_err().contains("usage"));
        assert!(run_validate_assays(vec![assay.display().to_string(), "--report".to_owned()])
            .unwrap_err()
            .contains("--report requires"));
        assert!(run_validate_assays(vec![assay.display().to_string(), "extra".to_owned()])
            .unwrap_err()
            .contains("unexpected argument"));

        fs::remove_dir_all(dir).unwrap();
    }
}
