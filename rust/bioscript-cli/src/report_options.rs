#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AppOutputFormat {
    Tsv,
    Json,
    Jsonl,
    Both,
}

struct AppReportOptions {
    manifest_path: PathBuf,
    input_files: Vec<PathBuf>,
    output_dir: PathBuf,
    root: PathBuf,
    html: bool,
    observations_format: AppOutputFormat,
    reports_format: AppOutputFormat,
    loader: GenotypeLoadOptions,
    filters: Vec<String>,
}

fn run_app_report(args: Vec<String>) -> Result<(), String> {
    let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
    let mut manifest_path: Option<PathBuf> = None;
    let mut input_files: Vec<PathBuf> = Vec::new();
    let mut output_dir: Option<PathBuf> = None;
    let mut root: Option<PathBuf> = None;
    let mut html = false;
    let mut observations_format = AppOutputFormat::Tsv;
    let mut reports_format = AppOutputFormat::Jsonl;
    let mut filters = Vec::new();
    let mut loader = GenotypeLoadOptions::default();

    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--input-file" => input_files.push(PathBuf::from(
                iter.next().ok_or("--input-file requires a path")?,
            )),
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
            "--observations-format" => {
                observations_format = parse_app_output_format(
                    &iter
                        .next()
                        .ok_or("--observations-format requires a value")?,
                )?;
            }
            "--reports-format" => {
                reports_format = parse_app_output_format(
                    &iter.next().ok_or("--reports-format requires a value")?,
                )?;
            }
            "--input-format" => {
                let value = iter.next().ok_or("--input-format requires a value")?;
                if value.eq_ignore_ascii_case("auto") {
                    loader.format = None;
                } else {
                    loader.format =
                        Some(value.parse::<GenotypeSourceFormat>().map_err(|err| {
                            format!("invalid --input-format value {value}: {err}")
                        })?);
                }
            }
            "--input-index" => {
                loader.input_index = Some(PathBuf::from(
                    iter.next().ok_or("--input-index requires a path")?,
                ));
            }
            "--reference-file" => {
                loader.reference_file = Some(PathBuf::from(
                    iter.next().ok_or("--reference-file requires a path")?,
                ));
            }
            "--reference-index" => {
                loader.reference_index = Some(PathBuf::from(
                    iter.next().ok_or("--reference-index requires a path")?,
                ));
            }
            value if value.starts_with('-') => return Err(format!("unexpected argument: {value}")),
            value => {
                if manifest_path.is_none() {
                    manifest_path = Some(PathBuf::from(value));
                } else {
                    input_files.push(PathBuf::from(value));
                }
            }
        }
    }

    let Some(manifest_path) = manifest_path else {
        return Err("usage: bioscript report <manifest.yaml> --input-file <path> [--input-file <path>...] --output-dir <dir> [--html]".to_owned());
    };
    if input_files.is_empty() {
        return Err("bioscript report requires at least one --input-file".to_owned());
    }
    let output_dir = output_dir.ok_or("bioscript report requires --output-dir")?;
    let root = root.unwrap_or(cwd);
    normalize_loader_paths(&root, &mut loader);

    let options = AppReportOptions {
        manifest_path: absolutize(&root, &manifest_path),
        input_files: input_files
            .iter()
            .map(|path| absolutize(&root, path))
            .collect(),
        output_dir: absolutize(&root, &output_dir),
        root,
        html,
        observations_format,
        reports_format,
        loader,
        filters,
    };
    generate_app_report(&options)
}

fn parse_app_output_format(value: &str) -> Result<AppOutputFormat, String> {
    match value {
        "tsv" => Ok(AppOutputFormat::Tsv),
        "json" => Ok(AppOutputFormat::Json),
        "jsonl" => Ok(AppOutputFormat::Jsonl),
        "both" => Ok(AppOutputFormat::Both),
        other => Err(format!(
            "unsupported output format '{other}'; expected tsv, json, jsonl, or both"
        )),
    }
}

fn absolutize(root: &Path, path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        root.join(path)
    }
}

fn generate_app_report(options: &AppReportOptions) -> Result<(), String> {
    fs::create_dir_all(&options.output_dir).map_err(|err| {
        format!(
            "failed to create output dir {}: {err}",
            options.output_dir.display()
        )
    })?;

    let assay_id = app_assay_id(&options.manifest_path)?;
    let findings = load_manifest_findings(&options.root, &options.manifest_path)?;
    let provenance = load_manifest_provenance_links(&options.root, &options.manifest_path)?;
    let mut observations = Vec::new();
    let mut analyses = Vec::new();
    let mut reports = Vec::new();

    for input_file in &options.input_files {
        let participant_id = participant_id_from_path(input_file);
        let rows = run_manifest_rows_for_report(
            &options.root,
            &options.manifest_path,
            input_file,
            &participant_id,
            &options.loader,
            &options.filters,
        )?;
        let input_observations = rows
            .iter()
            .map(|row| app_observation_from_manifest_row(&options.root, row, &assay_id))
            .collect::<Result<Vec<_>, _>>()?;
        observations.extend(input_observations.clone());
        let input_analyses = run_manifest_analyses_for_report(
            &options.root,
            &options.manifest_path,
            input_file,
            &participant_id,
            &options.loader,
            &options.output_dir,
        )?;
        analyses.extend(input_analyses.clone());
        let matched_findings = match_app_findings(&findings, &input_observations, &input_analyses);
        reports.push(app_report_json(
            &assay_id,
            &participant_id,
            input_file,
            &input_observations,
            &input_analyses,
            &matched_findings,
            &provenance,
        ));
    }

    write_app_observations(
        &options.output_dir,
        &observations,
        options.observations_format,
    )?;
    write_app_analyses(&options.output_dir, &analyses)?;
    write_app_reports(&options.output_dir, &reports, options.reports_format)?;
    if options.html {
        write_app_html(&options.output_dir, &observations, &reports)?;
    }

    println!(
        "observations: {}",
        options.output_dir.join("observations.tsv").display()
    );
    println!(
        "analysis: {}",
        options.output_dir.join("analysis.jsonl").display()
    );
    println!(
        "reports: {}",
        options.output_dir.join("reports.jsonl").display()
    );
    if options.html {
        println!("html: {}", options.output_dir.join("index.html").display());
    }
    Ok(())
}
