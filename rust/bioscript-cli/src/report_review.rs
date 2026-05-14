struct ReviewReportOptions {
    manifest_path: PathBuf,
    cases_path: PathBuf,
    output_dir: PathBuf,
    root: PathBuf,
    html: bool,
    filters: Vec<String>,
}

struct ReviewCase {
    id: String,
    label: String,
    genotypes: BTreeMap<String, Option<String>>,
}

fn run_review_report(args: Vec<String>) -> Result<(), String> {
    let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
    let mut manifest_path: Option<PathBuf> = None;
    let mut cases_path: Option<PathBuf> = None;
    let mut output_dir: Option<PathBuf> = None;
    let mut root: Option<PathBuf> = None;
    let mut html = false;
    let mut filters = Vec::new();

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--cases" => {
                cases_path = Some(PathBuf::from(iter.next().ok_or("--cases requires a path")?));
            }
            "--output-dir" => {
                output_dir = Some(PathBuf::from(
                    iter.next().ok_or("--output-dir requires a path")?,
                ));
            }
            "--root" => {
                root = Some(PathBuf::from(
                    iter.next().ok_or("--root requires a directory")?,
                ));
            }
            "--html" => html = true,
            "--filter" => filters.push(iter.next().ok_or("--filter requires key=value")?),
            value if value.starts_with('-') => return Err(format!("unexpected argument: {value}")),
            value => {
                if manifest_path.is_none() {
                    manifest_path = Some(PathBuf::from(value));
                } else {
                    return Err(format!("unexpected argument: {value}"));
                }
            }
        }
    }

    let Some(manifest_path) = manifest_path else {
        return Err("usage: bioscript review <manifest.yaml> --cases <cases.yaml> --output-dir <dir> [--html]".to_owned());
    };
    let cases_path = cases_path.ok_or("bioscript review requires --cases")?;
    let output_dir = output_dir.ok_or("bioscript review requires --output-dir")?;
    let root = root.unwrap_or(cwd);
    let manifest_path = if is_package_url(&manifest_path.to_string_lossy()) {
        prepare_package_entrypoint_from_arg(&root, &manifest_path)?
    } else {
        prepare_package_entrypoint_from_arg(&root, &absolutize(&root, &manifest_path))?
    };
    let options = ReviewReportOptions {
        manifest_path,
        cases_path: absolutize(&root, &cases_path),
        output_dir: absolutize(&root, &output_dir),
        root,
        html,
        filters,
    };
    generate_review_report(&options)
}

fn generate_review_report(options: &ReviewReportOptions) -> Result<(), String> {
    fs::create_dir_all(&options.output_dir).map_err(|err| {
        format!(
            "failed to create output dir {}: {err}",
            options.output_dir.display()
        )
    })?;

    let manifest_workspace = bioscript_reporting::FilesystemManifestWorkspace::new(&options.root);
    let manifest_path = options.manifest_path.display().to_string();
    let manifest_context =
        bioscript_reporting::load_report_manifest_context(&manifest_workspace, &manifest_path)?;
    let cases = load_review_cases(&options.cases_path)?;
    let mut observations = Vec::new();
    let mut analyses = Vec::new();
    let mut reports = Vec::new();

    for case in cases {
        let input_bytes = review_case_genotype_text(&case);
        let store = GenotypeStore::from_bytes(&format!("{}.txt", case.id), input_bytes.as_bytes())
            .map_err(|err| err.to_string())?;
        let input_observations = run_manifest_rows_with_store(
            &options.root,
            &options.manifest_path,
            &store,
            &case.id,
            &options.filters,
        )?
        .iter()
        .map(|row| {
            app_observation_from_manifest_row(
                &options.root,
                row,
                &manifest_context.assay_id,
                None,
                None,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
        observations.extend(input_observations.clone());

        let input_analyses = run_review_analyses(options, &case, &input_bytes)?;
        analyses.extend(input_analyses.clone());
        let synthetic_input = PathBuf::from(format!("review://{}", case.id));
        let synthetic_input_name = synthetic_input
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default();
        let synthetic_input_path = synthetic_input.display().to_string();
        let mut report = bioscript_reporting::app_input_report_json(
            bioscript_reporting::AppInputReportInput {
            assay_id: &manifest_context.assay_id,
            participant_id: &case.id,
            input_file_name: synthetic_input_name,
            input_file_path: &synthetic_input_path,
            observations: &input_observations,
            analyses: &input_analyses,
            findings: &manifest_context.findings,
            provenance: &manifest_context.provenance,
            input_inspection: None,
            manifest_metadata: &manifest_context.manifest_metadata,
            },
        );
        if let Some(object) = report.as_object_mut() {
            object.insert(
                "review_case".to_owned(),
                serde_json::json!({
                    "id": case.id,
                    "label": case.label,
                }),
            );
        }
        reports.push(report);
    }
    let review_temp_dir = options.output_dir.join(".review-temp");
    if review_temp_dir.exists() {
        fs::remove_dir_all(&review_temp_dir).map_err(|err| {
            format!(
                "failed to remove review temp dir {}: {err}",
                review_temp_dir.display()
            )
        })?;
    }

    write_app_observations(&options.output_dir, &observations, AppOutputFormat::Tsv)?;
    write_app_analyses(&options.output_dir, &analyses)?;
    write_app_reports(&options.output_dir, &reports, AppOutputFormat::Jsonl)?;
    if options.html {
        write_app_html(&options.output_dir, &observations, &reports)?;
    }
    println!(
        "review reports: {}",
        options.output_dir.join("reports.jsonl").display()
    );
    if options.html {
        println!("review html: {}", options.output_dir.join("index.html").display());
    }
    Ok(())
}

fn run_manifest_rows_with_store(
    runtime_root: &Path,
    manifest_path: &Path,
    store: &GenotypeStore,
    participant_id: &str,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    match manifest_schema(manifest_path)?.as_str() {
        "bioscript:variant:1.0" | "bioscript:variant" => {
            let manifest = load_variant_manifest(manifest_path)?;
            Ok(vec![run_variant_manifest_with_store(
                runtime_root,
                &manifest,
                store,
                Some(participant_id),
            )?])
        }
        "bioscript:panel:1.0" => {
            let manifest = load_panel_manifest(manifest_path)?;
            run_panel_manifest_with_store(runtime_root, &manifest, store, Some(participant_id), filters)
        }
        "bioscript:assay:1.0" => {
            let manifest = load_assay_manifest(manifest_path)?;
            run_assay_manifest_with_store(runtime_root, &manifest, store, Some(participant_id), filters)
        }
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

fn run_review_analyses(
    options: &ReviewReportOptions,
    case: &ReviewCase,
    input_bytes: &str,
) -> Result<Vec<serde_json::Value>, String> {
    let temp_dir = options.output_dir.join(".review-temp");
    fs::create_dir_all(&temp_dir).map_err(|err| {
        format!(
            "failed to create review temp dir {}: {err}",
            temp_dir.display()
        )
    })?;
    let temp_path = temp_dir.join(format!("{}.txt", case.id));
    fs::write(&temp_path, input_bytes)
        .map_err(|err| format!("failed to write review temp input {}: {err}", temp_path.display()))?;
    let loader = GenotypeLoadOptions {
        format: Some(GenotypeSourceFormat::Text),
        ..GenotypeLoadOptions::default()
    };
    let observation_rows = Vec::new();
    let analysis_options = ReportAnalysisOptions {
        runtime_root: &options.root,
        input_file: &temp_path,
        participant_id: &case.id,
        loader: &loader,
        output_dir: &options.output_dir,
        observation_rows: &observation_rows,
        filters: &options.filters,
        max_duration_ms: 1_000,
    };
    let result = run_manifest_analyses_for_report(&options.manifest_path, &analysis_options);
    let cleanup = fs::remove_file(&temp_path);
    if let Err(err) = cleanup {
        return Err(format!("failed to remove review temp input {}: {err}", temp_path.display()));
    }
    result
}

fn load_review_cases(path: &Path) -> Result<Vec<ReviewCase>, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read review cases {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse review cases {}: {err}", path.display()))?;
    let cases = value
        .get("cases")
        .and_then(serde_yaml::Value::as_sequence)
        .ok_or("review cases missing cases list")?;
    cases
        .iter()
        .map(review_case_from_yaml)
        .collect::<Result<Vec<_>, _>>()
}

fn review_case_from_yaml(value: &serde_yaml::Value) -> Result<ReviewCase, String> {
    let id = value
        .get("id")
        .and_then(serde_yaml::Value::as_str)
        .ok_or("review case missing id")?
        .to_owned();
    let label = value
        .get("label")
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or(&id)
        .to_owned();
    let genotypes_value = value
        .get("genotypes")
        .or_else(|| value.get("variants"))
        .and_then(serde_yaml::Value::as_mapping)
        .ok_or_else(|| format!("review case {id} missing genotypes"))?;
    let mut genotypes = BTreeMap::new();
    for (key, value) in genotypes_value {
        let Some(rsid) = key.as_str() else {
            return Err(format!("review case {id} has non-string genotype key"));
        };
        let genotype = if value.is_null() {
            None
        } else {
            Some(
                value
                    .as_str()
                    .ok_or_else(|| format!("review case {id} genotype {rsid} must be string or null"))?
                    .to_owned(),
            )
        };
        genotypes.insert(rsid.to_owned(), genotype);
    }
    Ok(ReviewCase {
        id,
        label,
        genotypes,
    })
}

fn review_case_genotype_text(case: &ReviewCase) -> String {
    let mut out = String::from("rsid\tgenotype\n");
    for (rsid, genotype) in &case.genotypes {
        let Some(genotype) = genotype else {
            continue;
        };
        let _ = writeln!(out, "{}\t{}", rsid.replace('\t', " "), genotype.replace('\t', " "));
    }
    out
}

#[cfg(test)]
mod review_report_tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-review-report-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn run_review_report_validates_required_arguments() {
        assert!(run_review_report(Vec::new()).unwrap_err().contains("usage"));
        assert!(run_review_report(vec!["manifest.yaml".to_owned()])
            .unwrap_err()
            .contains("--cases"));
        assert!(run_review_report(vec![
            "manifest.yaml".to_owned(),
            "--cases".to_owned(),
            "cases.yaml".to_owned(),
        ])
        .unwrap_err()
        .contains("--output-dir"));
        assert!(run_review_report(vec![
            "manifest.yaml".to_owned(),
            "--unknown".to_owned(),
        ])
        .unwrap_err()
        .contains("unexpected argument"));
        assert!(run_review_report(vec![
            "first.yaml".to_owned(),
            "second.yaml".to_owned(),
        ])
        .unwrap_err()
        .contains("unexpected argument"));
    }

    #[test]
    fn review_cases_load_labels_variants_and_null_genotypes() {
        let dir = temp_dir("cases");
        let cases_path = dir.join("cases.yaml");
        fs::write(
            &cases_path,
            r#"
cases:
  - id: c1
    label: First case
    genotypes:
      rs1: A/G
      rs2: null
  - id: c2
    variants:
      rs3: C/T
"#,
        )
        .unwrap();

        let cases = load_review_cases(&cases_path).unwrap();
        assert_eq!(cases.len(), 2);
        assert_eq!(cases[0].label, "First case");
        assert_eq!(cases[0].genotypes["rs1"], Some("A/G".to_owned()));
        assert_eq!(cases[0].genotypes["rs2"], None);
        assert_eq!(cases[1].label, "c2");
        assert_eq!(cases[1].genotypes["rs3"], Some("C/T".to_owned()));

        let text = review_case_genotype_text(&cases[0]);
        assert!(text.contains("rs1\tA/G"));
        assert!(!text.contains("rs2"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn review_case_parser_reports_shape_errors() {
        let missing_cases = serde_yaml::from_str::<serde_yaml::Value>("not_cases: []").unwrap();
        let path = temp_dir("errors").join("missing-cases.yaml");
        fs::write(&path, serde_yaml::to_string(&missing_cases).unwrap()).unwrap();
        assert!(review_cases_err(&path).contains("missing cases list"));
        let dir = path.parent().unwrap().to_path_buf();

        let missing_id = serde_yaml::from_str::<serde_yaml::Value>(
            r#"{label: no id, genotypes: {rs1: A/G}}"#,
        )
        .unwrap();
        assert!(review_case_err(&missing_id).contains("missing id"));

        let missing_genotypes =
            serde_yaml::from_str::<serde_yaml::Value>(r#"{id: c1}"#).unwrap();
        assert!(review_case_err(&missing_genotypes).contains("missing genotypes"));

        let bad_key =
            serde_yaml::from_str::<serde_yaml::Value>(r#"{id: c1, genotypes: {1: A/G}}"#)
                .unwrap();
        assert!(review_case_err(&bad_key).contains("non-string genotype key"));

        let bad_value =
            serde_yaml::from_str::<serde_yaml::Value>(r#"{id: c1, genotypes: {rs1: [A, G]}}"#)
                .unwrap();
        assert!(review_case_err(&bad_value).contains("must be string or null"));

        fs::remove_dir_all(dir).unwrap();
    }

    fn review_cases_err(path: &Path) -> String {
        match load_review_cases(path) {
            Ok(_) => panic!("expected review cases to fail"),
            Err(err) => err,
        }
    }

    fn review_case_err(value: &serde_yaml::Value) -> String {
        match review_case_from_yaml(value) {
            Ok(_) => panic!("expected review case to fail"),
            Err(err) => err,
        }
    }

    #[test]
    fn generate_review_report_writes_observations_reports_html_and_cleans_temp_input() {
        let dir = temp_dir("generate");
        let manifest = dir.join("variant.yaml");
        fs::write(
            &manifest,
            r#"
schema: bioscript:variant:1.0
version: "1.0"
name: rs1
gene: ABC
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
findings:
  - schema: bioscript:trait:1.0
    summary: Variant present
    binding:
      source: variant
      variant: variant.yaml
      key: outcome
      value: variant
provenance:
  sources:
    - kind: database
      label: Fixture
      url: https://example.test/rs1
"#,
        )
        .unwrap();
        let cases = dir.join("cases.yaml");
        fs::write(
            &cases,
            r#"
cases:
  - id: case1
    label: Case One
    genotypes:
      rs1: A/G
  - id: case2
    genotypes:
      rs1: null
"#,
        )
        .unwrap();
        let output = dir.join("out");
        let options = ReviewReportOptions {
            manifest_path: manifest,
            cases_path: cases,
            output_dir: output.clone(),
            root: dir.clone(),
            html: true,
            filters: Vec::new(),
        };

        generate_review_report(&options).unwrap();

        assert!(fs::read_to_string(output.join("observations.tsv"))
            .unwrap()
            .contains("case1"));
        let reports = fs::read_to_string(output.join("reports.jsonl")).unwrap();
        assert!(reports.contains("\"review_case\""));
        assert!(reports.contains("\"case1\""));
        assert!(fs::read_to_string(output.join("index.html"))
            .unwrap()
            .contains("<!doctype html>"));
        assert!(!output.join(".review-temp").exists());

        fs::remove_dir_all(dir).unwrap();
    }
}
