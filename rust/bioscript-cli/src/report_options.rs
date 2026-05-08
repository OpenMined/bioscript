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
    open_report: bool,
    observations_format: AppOutputFormat,
    reports_format: AppOutputFormat,
    loader: GenotypeLoadOptions,
    filters: Vec<String>,
    analysis_max_duration_ms: u64,
    detect_sex: bool,
    sample_sex: Option<InferredSex>,
}

struct AppReportCliState {
    cwd: PathBuf,
    manifest_path: Option<PathBuf>,
    input_files: Vec<PathBuf>,
    output_dir: Option<PathBuf>,
    root: Option<PathBuf>,
    html: bool,
    open_report: bool,
    observations_format: AppOutputFormat,
    reports_format: AppOutputFormat,
    loader: GenotypeLoadOptions,
    filters: Vec<String>,
    analysis_max_duration_ms: u64,
    detect_sex: bool,
    sample_sex: Option<InferredSex>,
}

fn run_app_report(args: Vec<String>) -> Result<(), String> {
    let mut state = AppReportCliState::new()?;
    let mut iter = args.into_iter();
    while let Some(arg) = iter.next() {
        state.consume_arg(&arg, &mut iter)?;
    }
    generate_app_report(&state.finish()?)
}

impl AppReportCliState {
    fn new() -> Result<Self, String> {
        Ok(Self {
            cwd: env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?,
            manifest_path: None,
            input_files: Vec::new(),
            output_dir: None,
            root: None,
            html: false,
            open_report: false,
            observations_format: AppOutputFormat::Tsv,
            reports_format: AppOutputFormat::Jsonl,
            loader: GenotypeLoadOptions::default(),
            filters: Vec::new(),
            analysis_max_duration_ms: 1_000,
            detect_sex: false,
            sample_sex: None,
        })
    }

    fn consume_arg(
        &mut self,
        arg: &str,
        iter: &mut std::vec::IntoIter<String>,
    ) -> Result<(), String> {
        match arg {
            "--input-file" => self.input_files.push(PathBuf::from(next_arg(iter, "--input-file")?)),
            "--output-dir" => self.output_dir = Some(PathBuf::from(next_arg(iter, "--output-dir")?)),
            "--root" => self.root = Some(PathBuf::from(next_arg(iter, "--root")?)),
            "--html" => self.html = true,
            "--open" => {
                self.html = true;
                self.open_report = true;
            }
            "--filter" => self.filters.push(next_arg(iter, "--filter")?),
            "--detect-sex" => self.detect_sex = true,
            "--sample-sex" => {
                self.sample_sex = Some(parse_sample_sex(&next_arg(iter, "--sample-sex")?)?);
            }
            "--observations-format" => {
                self.observations_format =
                    parse_app_output_format(&next_arg(iter, "--observations-format")?)?;
            }
            "--reports-format" => {
                self.reports_format =
                    parse_app_output_format(&next_arg(iter, "--reports-format")?)?;
            }
            "--analysis-max-duration-ms" => self.parse_analysis_timeout(iter)?,
            "--input-format" => self.parse_input_format(iter)?,
            "--input-index" => {
                self.loader.input_index = Some(PathBuf::from(next_arg(iter, "--input-index")?));
            }
            "--reference-file" => {
                self.loader.reference_file =
                    Some(PathBuf::from(next_arg(iter, "--reference-file")?));
            }
            "--reference-index" => {
                self.loader.reference_index =
                    Some(PathBuf::from(next_arg(iter, "--reference-index")?));
            }
            "--allow-md5-mismatch" => self.loader.allow_reference_md5_mismatch = true,
            value if value.starts_with('-') => return Err(format!("unexpected argument: {value}")),
            value => self.consume_path(value),
        }
        Ok(())
    }

    fn parse_analysis_timeout(
        &mut self,
        iter: &mut std::vec::IntoIter<String>,
    ) -> Result<(), String> {
        let value = next_arg(iter, "--analysis-max-duration-ms")?;
        self.analysis_max_duration_ms = value
            .parse::<u64>()
            .map_err(|err| format!("invalid --analysis-max-duration-ms value {value}: {err}"))?;
        Ok(())
    }

    fn parse_input_format(&mut self, iter: &mut std::vec::IntoIter<String>) -> Result<(), String> {
        let value = next_arg(iter, "--input-format")?;
        self.loader.format = if value.eq_ignore_ascii_case("auto") {
            None
        } else {
            Some(value.parse::<GenotypeSourceFormat>().map_err(|err| {
                format!("invalid --input-format value {value}: {err}")
            })?)
        };
        Ok(())
    }

    fn consume_path(&mut self, value: &str) {
        if self.manifest_path.is_none() {
            self.manifest_path = Some(PathBuf::from(value));
        } else {
            self.input_files.push(PathBuf::from(value));
        }
    }

    fn finish(mut self) -> Result<AppReportOptions, String> {
        let Some(manifest_path) = self.manifest_path else {
            return Err("usage: bioscript report <manifest.yaml> --input-file <path> [--input-file <path>...] --output-dir <dir> [--html]".to_owned());
        };
        if self.input_files.is_empty() {
            return Err("bioscript report requires at least one --input-file".to_owned());
        }
        let output_dir = self.output_dir.ok_or("bioscript report requires --output-dir")?;
        let root = self.root.unwrap_or(self.cwd);
        normalize_loader_paths(&root, &mut self.loader);
        let manifest_path = if is_package_url(&manifest_path.to_string_lossy()) {
            prepare_package_entrypoint_from_arg(&root, &manifest_path)?
        } else {
            prepare_package_entrypoint_from_arg(&root, &absolutize(&root, &manifest_path))?
        };
        Ok(AppReportOptions {
            manifest_path,
            input_files: self
                .input_files
                .iter()
                .map(|path| absolutize(&root, path))
                .collect(),
            output_dir: absolutize(&root, &output_dir),
            root,
            html: self.html,
            open_report: self.open_report,
            observations_format: self.observations_format,
            reports_format: self.reports_format,
            loader: self.loader,
            filters: self.filters,
            analysis_max_duration_ms: self.analysis_max_duration_ms,
            detect_sex: self.detect_sex,
            sample_sex: self.sample_sex,
        })
    }
}

fn next_arg(iter: &mut std::vec::IntoIter<String>, flag: &str) -> Result<String, String> {
    iter.next().ok_or_else(|| format!("{flag} requires a value"))
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

fn parse_sample_sex(value: &str) -> Result<InferredSex, String> {
    match value.to_ascii_lowercase().as_str() {
        "male" | "m" => Ok(InferredSex::Male),
        "female" | "f" => Ok(InferredSex::Female),
        "unknown" | "u" => Ok(InferredSex::Unknown),
        other => Err(format!(
            "unsupported --sample-sex value '{other}'; expected male, female, or unknown"
        )),
    }
}

fn explicit_sample_sex_inference(sex: InferredSex) -> SexInference {
    SexInference {
        sex,
        confidence: SexDetectionConfidence::High,
        method: "explicit_sample_sex".to_owned(),
        evidence: vec!["source=sample_sex_cli".to_owned()],
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
    let manifest_metadata = report_manifest_metadata(&options.manifest_path)?;
    let findings = load_manifest_findings(&options.root, &options.manifest_path)?;
    let provenance = load_manifest_provenance_links(&options.root, &options.manifest_path)?;
    let mut observations = Vec::new();
    let mut analyses = Vec::new();
    let mut reports = Vec::new();

    for input_file in &options.input_files {
        let participant_id = participant_id_from_path(input_file);
        let inspect_options = InspectOptions {
            input_index: options.loader.input_index.clone(),
            reference_file: options.loader.reference_file.clone(),
            reference_index: options.loader.reference_index.clone(),
            detect_sex: options.detect_sex,
        };
        let mut input_inspection =
            inspect_file(input_file, &inspect_options).map_err(|err| err.to_string())?;
        if let Some(sample_sex) = options.sample_sex {
            input_inspection.inferred_sex = Some(explicit_sample_sex_inference(sample_sex));
        }
        let input_loader = loader_with_inspection(&options.loader, &input_inspection);
        let rows = run_manifest_rows_for_report(
            &options.root,
            &options.manifest_path,
            input_file,
            &participant_id,
            &input_loader,
            &options.filters,
        )?;
        let input_observations = rows
            .iter()
            .map(|row| {
                app_observation_from_manifest_row(
                    &options.root,
                    row,
                    &assay_id,
                    input_inspection.inferred_sex.as_ref(),
                    input_inspection.assembly,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        observations.extend(input_observations.clone());
        let analysis_options = ReportAnalysisOptions {
            runtime_root: &options.root,
            input_file,
            participant_id: &participant_id,
            loader: &input_loader,
            output_dir: &options.output_dir,
            filters: &options.filters,
            max_duration_ms: options.analysis_max_duration_ms,
        };
        let input_analyses =
            run_manifest_analyses_for_report(&options.manifest_path, &analysis_options)?;
        analyses.extend(input_analyses.clone());
        let matched_findings = match_app_findings(&findings, &input_observations, &input_analyses);
        reports.push(app_report_json(AppReportJsonInput {
            assay_id: &assay_id,
            participant_id: &participant_id,
            input_file,
            observations: &input_observations,
            analyses: &input_analyses,
            findings: &matched_findings,
            provenance: &provenance,
            input_inspection: Some(&input_inspection),
            manifest_metadata: &manifest_metadata,
        }));
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
    open_app_html_report_if_requested(options);
    print_app_report_paths(&options.output_dir, options.html);
    Ok(())
}

fn loader_with_inspection(
    base: &GenotypeLoadOptions,
    inspection: &bioscript_formats::FileInspection,
) -> GenotypeLoadOptions {
    let mut loader = base.clone();
    loader.assembly = inspection.assembly.or(loader.assembly);
    loader.inferred_sex = inspection
        .inferred_sex
        .as_ref()
        .map(|inference| inference.sex)
        .or(loader.inferred_sex);
    loader
}

fn open_app_html_report_if_requested(options: &AppReportOptions) {
    if options.open_report
        && let Err(err) = open_html_report(&options.output_dir.join("index.html"))
    {
        eprintln!("warning: {err}");
    }
}

fn print_app_report_paths(output_dir: &Path, include_html: bool) {
    println!("observations: {}", output_dir.join("observations.tsv").display());
    println!("analysis: {}", output_dir.join("analysis.jsonl").display());
    println!("reports: {}", output_dir.join("reports.jsonl").display());
    if include_html {
        println!("html: {}", output_dir.join("index.html").display());
    }
}

fn open_html_report(path: &Path) -> Result<(), String> {
    let opener = if cfg!(target_os = "macos") {
        "open"
    } else if cfg!(target_os = "windows") {
        "cmd"
    } else {
        "xdg-open"
    };
    let status = if cfg!(target_os = "windows") {
        let path_text = path.display().to_string();
        std::process::Command::new(opener)
            .args(["/C", "start", "", &path_text])
            .status()
    } else {
        std::process::Command::new(opener).arg(path).status()
    }
    .map_err(|err| format!("failed to open html report {}: {err}", path.display()))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!(
            "failed to open html report {}: opener exited with {status}",
            path.display()
        ))
    }
}
