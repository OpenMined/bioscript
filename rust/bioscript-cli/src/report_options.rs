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

    let manifest_workspace = bioscript_reporting::FilesystemManifestWorkspace::new(&options.root);
    let manifest_path = options.manifest_path.display().to_string();
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
        let store = match GenotypeStore::from_file_with_options(input_file, &input_loader) {
            Ok(store) => store,
            Err(err) => {
                // Aligned inputs (BAM, or a CRAM we can't pileup without a
                // reference) can't always be genotyped here: BAM variant
                // lookup is unimplemented and a reference-compressed CRAM
                // needs --reference-file. For an advanced assay the analysis
                // consumes the raw aligned reads directly (e.g. VNtyper runs
                // Kestrel over the MUC1 slice), so per-variant genotyping is
                // best-effort enrichment, not a hard gate. Degrade to an
                // empty store — variant members report as missing (a valid
                // partial result) — instead of aborting the whole report.
                let is_aligned = matches!(
                    input_inspection.detected_kind,
                    bioscript_formats::DetectedKind::AlignmentBam
                        | bioscript_formats::DetectedKind::AlignmentCram
                );
                if is_aligned {
                    eprintln!(
                        "bioscript: variant genotyping unavailable for aligned \
                         input {} ({err}); continuing with analyses only — \
                         variant members reported as missing",
                        input_file.display()
                    );
                    GenotypeStore::empty()
                } else {
                    return Err(err.to_string());
                }
            }
        };
        let analysis_runner = CliReportAnalysisRunner {
            runtime_root: &options.root,
            input_file,
            participant_id: &participant_id,
            loader: &input_loader,
            output_dir: &options.output_dir,
            max_duration_ms: options.analysis_max_duration_ms,
        };
        let input_file_name = input_file
            .file_name()
            .and_then(|value| value.to_str())
            .unwrap_or_default();
        let input_file_path = input_file.display().to_string();
        let run = bioscript_reporting::run_report(
            &manifest_workspace,
            &manifest_path,
            &store,
            &analysis_runner,
            bioscript_reporting::ReportInputContext {
                participant_id: &participant_id,
                input_file_name,
                input_file_path: &input_file_path,
                input_inspection: Some(&input_inspection),
            },
            bioscript_reporting::ReportRunOptions {
                filters: &options.filters,
            },
        )?;
        observations.extend(run.observations);
        analyses.extend(run.analyses);
        reports.push(run.report);
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

#[cfg(test)]
mod app_report_option_tests {
    use super::*;
    use bioscript_core::Assembly;
    use bioscript_formats::{
        DetectedKind, DetectionConfidence, FileContainer, FileInspection, InferredSex,
        SexDetectionConfidence, SexInference,
    };
    use std::time::{SystemTime, UNIX_EPOCH};

    fn args(items: &[&str]) -> std::vec::IntoIter<String> {
        items
            .iter()
            .map(|item| (*item).to_owned())
            .collect::<Vec<_>>()
            .into_iter()
    }

    fn inspection() -> FileInspection {
        FileInspection {
            path: PathBuf::from("sample.vcf"),
            container: FileContainer::Plain,
            detected_kind: DetectedKind::Vcf,
            confidence: DetectionConfidence::Authoritative,
            source: None,
            assembly: Some(Assembly::Grch37),
            phased: Some(false),
            selected_entry: None,
            has_index: Some(false),
            index_path: None,
            reference_matches: None,
            inferred_sex: Some(SexInference {
                sex: InferredSex::Male,
                confidence: SexDetectionConfidence::High,
                method: "fixture".to_owned(),
                evidence: vec!["test".to_owned()],
            }),
            evidence: Vec::new(),
            warnings: Vec::new(),
            duration_ms: 0,
        }
    }

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-app-report-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn output_format_and_sample_sex_parsers_accept_aliases_and_reject_unknowns() {
        assert_eq!(parse_app_output_format("tsv").unwrap(), AppOutputFormat::Tsv);
        assert_eq!(parse_app_output_format("json").unwrap(), AppOutputFormat::Json);
        assert_eq!(
            parse_app_output_format("jsonl").unwrap(),
            AppOutputFormat::Jsonl
        );
        assert_eq!(parse_app_output_format("both").unwrap(), AppOutputFormat::Both);
        assert!(parse_app_output_format("xml").unwrap_err().contains("unsupported"));

        assert_eq!(parse_sample_sex("male").unwrap(), InferredSex::Male);
        assert_eq!(parse_sample_sex("M").unwrap(), InferredSex::Male);
        assert_eq!(parse_sample_sex("female").unwrap(), InferredSex::Female);
        assert_eq!(parse_sample_sex("f").unwrap(), InferredSex::Female);
        assert_eq!(parse_sample_sex("unknown").unwrap(), InferredSex::Unknown);
        assert_eq!(parse_sample_sex("U").unwrap(), InferredSex::Unknown);
        assert!(parse_sample_sex("other").unwrap_err().contains("unsupported"));
    }

    #[test]
    fn cli_state_consumes_flags_and_finishes_with_normalized_paths() {
        let mut state = AppReportCliState::new().unwrap();
        let mut iter = args(&[
            "panel.yaml",
            "--input-file",
            "sample.vcf",
            "--output-dir",
            "out",
            "--root",
            ".",
            "--html",
            "--filter",
            "tag=pgx",
            "--detect-sex",
            "--sample-sex",
            "female",
            "--observations-format",
            "both",
            "--reports-format",
            "json",
            "--analysis-max-duration-ms",
            "2500",
            "--input-format",
            "vcf",
            "--input-index",
            "sample.vcf.tbi",
            "--reference-file",
            "ref.fa",
            "--reference-index",
            "ref.fa.fai",
            "--allow-md5-mismatch",
        ]);
        while let Some(arg) = iter.next() {
            state.consume_arg(&arg, &mut iter).unwrap();
        }

        let options = state.finish().unwrap();
        assert_eq!(options.manifest_path, PathBuf::from("./panel.yaml"));
        assert_eq!(options.input_files, vec![PathBuf::from("./sample.vcf")]);
        assert_eq!(options.output_dir, PathBuf::from("./out"));
        assert!(options.html);
        assert!(!options.open_report);
        assert_eq!(options.filters, vec!["tag=pgx"]);
        assert_eq!(options.observations_format, AppOutputFormat::Both);
        assert_eq!(options.reports_format, AppOutputFormat::Json);
        assert_eq!(options.analysis_max_duration_ms, 2500);
        assert!(options.detect_sex);
        assert_eq!(options.sample_sex, Some(InferredSex::Female));
        assert_eq!(options.loader.format, Some(GenotypeSourceFormat::Vcf));
        assert_eq!(
            options.loader.input_index,
            Some(PathBuf::from("./sample.vcf.tbi"))
        );
        assert!(options.loader.allow_reference_md5_mismatch);
    }

    #[test]
    fn cli_state_reports_required_argument_errors() {
        let missing_manifest = finish_err(AppReportCliState::new().unwrap());
        assert!(missing_manifest.contains("usage: bioscript report"));

        let mut state = AppReportCliState::new().unwrap();
        let mut iter = args(&["manifest.yaml"]);
        while let Some(arg) = iter.next() {
            state.consume_arg(&arg, &mut iter).unwrap();
        }
        assert!(finish_err(state).contains("at least one --input-file"));

        let mut state = AppReportCliState::new().unwrap();
        let mut iter = args(&["manifest.yaml", "--input-file", "sample.txt"]);
        while let Some(arg) = iter.next() {
            state.consume_arg(&arg, &mut iter).unwrap();
        }
        assert!(finish_err(state).contains("--output-dir"));

        let mut state = AppReportCliState::new().unwrap();
        let mut iter = args(&["--analysis-max-duration-ms", "not-a-number"]);
        let arg = iter.next().unwrap();
        assert!(state.consume_arg(&arg, &mut iter).unwrap_err().contains("invalid"));
    }

    #[test]
    fn loader_and_path_helpers_preserve_explicit_settings() {
        let base = GenotypeLoadOptions {
            assembly: Some(Assembly::Grch38),
            inferred_sex: None,
            ..GenotypeLoadOptions::default()
        };
        let loader = loader_with_inspection(&base, &inspection());
        assert_eq!(loader.assembly, Some(Assembly::Grch37));
        assert_eq!(loader.inferred_sex, Some(InferredSex::Male));

        let explicit = explicit_sample_sex_inference(InferredSex::Unknown);
        assert_eq!(explicit.sex, InferredSex::Unknown);
        assert_eq!(explicit.method, "explicit_sample_sex");

        let root = Path::new("/tmp/bioscript-root");
        assert_eq!(absolutize(root, Path::new("a/b")), root.join("a/b"));
        assert_eq!(
            absolutize(root, Path::new("/already/absolute")),
            PathBuf::from("/already/absolute")
        );
    }

    #[test]
    fn run_app_report_generates_outputs_for_text_input_and_explicit_sex() {
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
    chrom: "X"
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
"#,
        )
        .unwrap();
        let input = dir.join("sample.txt");
        fs::write(&input, "rsid\tgenotype\nrs1\tG\n").unwrap();
        let output = dir.join("out");

        run_app_report(vec![
            manifest.display().to_string(),
            "--input-file".to_owned(),
            input.display().to_string(),
            "--output-dir".to_owned(),
            output.display().to_string(),
            "--root".to_owned(),
            dir.display().to_string(),
            "--html".to_owned(),
            "--input-format".to_owned(),
            "text".to_owned(),
            "--sample-sex".to_owned(),
            "male".to_owned(),
            "--observations-format".to_owned(),
            "both".to_owned(),
            "--reports-format".to_owned(),
            "both".to_owned(),
        ])
        .unwrap();

        assert!(fs::read_to_string(output.join("observations.tsv"))
            .unwrap()
            .contains("sample"));
        assert!(fs::read_to_string(output.join("observations.jsonl"))
            .unwrap()
            .contains("detected_sex=male"));
        assert!(fs::read_to_string(output.join("reports.jsonl"))
            .unwrap()
            .contains("Variant present"));
        assert!(fs::read_to_string(output.join("reports.json"))
            .unwrap()
            .contains("report-set"));
        assert!(fs::read_to_string(output.join("index.html"))
            .unwrap()
            .contains("<!doctype html>"));

        fs::remove_dir_all(dir).unwrap();
    }

    fn finish_err(state: AppReportCliState) -> String {
        match state.finish() {
            Ok(_) => panic!("expected report CLI state to fail"),
            Err(err) => err,
        }
    }
}
