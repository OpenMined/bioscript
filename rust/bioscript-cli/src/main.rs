use std::{
    collections::BTreeMap,
    env,
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
    process::ExitCode,
    time::{Duration, Instant},
};

use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, InspectOptions, PrepareRequest,
    inspect_file, prepare_indexes, shell_flags,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig, StageTiming};
use bioscript_schema::{
    AssayManifest, PanelInterpretation, PanelManifest, VariantManifest, load_assay_manifest,
    load_panel_manifest, load_variant_manifest, validate_assays_path, validate_panels_path,
    validate_variants_path,
};
use monty::ResourceLimits;

fn main() -> ExitCode {
    match run_cli() {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("bioscript: {err}");
            ExitCode::FAILURE
        }
    }
}

#[allow(clippy::too_many_lines)]
fn run_cli() -> Result<(), String> {
    let mut args = env::args().skip(1);
    if let Some(first) = args.next() {
        if first == "report" {
            return run_app_report(args.collect());
        }
        if first == "validate-variants" {
            return run_validate_variants(args.collect());
        }
        if first == "validate-panels" {
            return run_validate_panels(args.collect());
        }
        if first == "validate-assays" {
            return run_validate_assays(args.collect());
        }
        if first == "prepare" {
            return run_prepare(args.collect());
        }
        if first == "inspect" {
            return run_inspect(args.collect());
        }
    }

    let mut args = env::args().skip(1);
    let mut script_path: Option<PathBuf> = None;
    let mut root: Option<PathBuf> = None;
    let mut input_file: Option<String> = None;
    let mut output_file: Option<String> = None;
    let mut participant_id: Option<String> = None;
    let mut trace_report: Option<PathBuf> = None;
    let mut timing_report: Option<PathBuf> = None;
    let mut filters: Vec<String> = Vec::new();
    let mut auto_index = false;
    let mut cache_dir: Option<PathBuf> = None;
    let mut loader = GenotypeLoadOptions::default();
    let mut limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(100))
        .max_memory(8 * 1024 * 1024)
        .max_allocations(200_000)
        .gc_interval(1000)
        .max_recursion_depth(Some(200));

    while let Some(arg) = args.next() {
        if arg == "--root" {
            let Some(value) = args.next() else {
                return Err("--root requires a directory".to_owned());
            };
            root = Some(PathBuf::from(value));
        } else if arg == "--input-file" {
            let Some(value) = args.next() else {
                return Err("--input-file requires a path".to_owned());
            };
            input_file = Some(value);
        } else if arg == "--output-file" {
            let Some(value) = args.next() else {
                return Err("--output-file requires a path".to_owned());
            };
            output_file = Some(value);
        } else if arg == "--participant-id" {
            let Some(value) = args.next() else {
                return Err("--participant-id requires a value".to_owned());
            };
            participant_id = Some(value);
        } else if arg == "--trace-report" {
            let Some(value) = args.next() else {
                return Err("--trace-report requires a path".to_owned());
            };
            trace_report = Some(PathBuf::from(value));
        } else if arg == "--timing-report" {
            let Some(value) = args.next() else {
                return Err("--timing-report requires a path".to_owned());
            };
            timing_report = Some(PathBuf::from(value));
        } else if arg == "--filter" {
            let Some(value) = args.next() else {
                return Err("--filter requires key=value".to_owned());
            };
            filters.push(value);
        } else if arg == "--input-format" {
            let Some(value) = args.next() else {
                return Err("--input-format requires a value".to_owned());
            };
            if value.eq_ignore_ascii_case("auto") {
                loader.format = None;
            } else {
                let parsed = value
                    .parse::<GenotypeSourceFormat>()
                    .map_err(|err| format!("invalid --input-format value {value}: {err}"))?;
                loader.format = Some(parsed);
            }
        } else if arg == "--input-index" {
            let Some(value) = args.next() else {
                return Err("--input-index requires a path".to_owned());
            };
            loader.input_index = Some(PathBuf::from(value));
        } else if arg == "--reference-file" {
            let Some(value) = args.next() else {
                return Err("--reference-file requires a path".to_owned());
            };
            loader.reference_file = Some(PathBuf::from(value));
        } else if arg == "--reference-index" {
            let Some(value) = args.next() else {
                return Err("--reference-index requires a path".to_owned());
            };
            loader.reference_index = Some(PathBuf::from(value));
        } else if arg == "--max-duration-ms" {
            let Some(value) = args.next() else {
                return Err("--max-duration-ms requires an integer".to_owned());
            };
            let parsed = value
                .parse::<u64>()
                .map_err(|err| format!("invalid --max-duration-ms value {value}: {err}"))?;
            limits = limits.max_duration(Duration::from_millis(parsed));
        } else if arg == "--max-memory-bytes" {
            let Some(value) = args.next() else {
                return Err("--max-memory-bytes requires an integer".to_owned());
            };
            let parsed = value
                .parse::<usize>()
                .map_err(|err| format!("invalid --max-memory-bytes value {value}: {err}"))?;
            limits = limits.max_memory(parsed);
        } else if arg == "--max-allocations" {
            let Some(value) = args.next() else {
                return Err("--max-allocations requires an integer".to_owned());
            };
            let parsed = value
                .parse::<usize>()
                .map_err(|err| format!("invalid --max-allocations value {value}: {err}"))?;
            limits = limits.max_allocations(parsed);
        } else if arg == "--auto-index" {
            auto_index = true;
        } else if arg == "--cache-dir" {
            let Some(value) = args.next() else {
                return Err("--cache-dir requires a path".to_owned());
            };
            cache_dir = Some(PathBuf::from(value));
        } else if arg == "--max-recursion-depth" {
            let Some(value) = args.next() else {
                return Err("--max-recursion-depth requires an integer".to_owned());
            };
            let parsed = value
                .parse::<usize>()
                .map_err(|err| format!("invalid --max-recursion-depth value {value}: {err}"))?;
            limits = limits.max_recursion_depth(Some(parsed));
        } else if script_path.is_none() {
            script_path = Some(PathBuf::from(arg));
        } else {
            return Err(format!("unexpected argument: {arg}"));
        }
    }

    let Some(script_path) = script_path else {
        return Err(
            "usage: bioscript <script.py|manifest.yaml> [--root <dir>] [--input-file <path>] [--output-file <path>] [--participant-id <id>] [--trace-report <path>] [--timing-report <path>] [--filter key=value] [--input-format auto|text|zip|vcf|cram] [--input-index <path>] [--reference-file <path>] [--reference-index <path>] [--auto-index] [--cache-dir <path>] [--max-duration-ms N] [--max-memory-bytes N] [--max-allocations N] [--max-recursion-depth N]\n       bioscript report <manifest.yaml> --input-file <path> [--input-file <path>...] --output-dir <dir> [--html] [--root <dir>] [--input-format auto|text|zip|vcf|cram]\n       bioscript validate-variants <path> [--report <file>]\n       bioscript validate-panels <path> [--report <file>]\n       bioscript validate-assays <path> [--report <file>]\n       bioscript prepare [--root <dir>] [--input-file <path>] [--reference-file <path>] [--input-format auto|text|zip|vcf|cram] [--cache-dir <path>]\n       bioscript inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>]"
                .to_owned(),
        );
    };

    let runtime_root = match root {
        Some(dir) => dir,
        None => {
            env::current_dir().map_err(|err| format!("failed to get current directory: {err}"))?
        }
    };
    normalize_loader_paths(&runtime_root, &mut loader);

    // auto-index: detect and build missing indexes for CRAM/BAM/FASTA
    let mut cli_timings: Vec<StageTiming> = Vec::new();
    if auto_index {
        let auto_index_started = Instant::now();
        let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
        let effective_cache = cache_dir
            .clone()
            .unwrap_or_else(|| cwd.join(".bioscript-cache"));
        let request = PrepareRequest {
            root: runtime_root.clone(),
            cwd: cwd.clone(),
            cache_dir: effective_cache,
            input_file: input_file.clone(),
            input_format: loader.format,
            reference_file: loader
                .reference_file
                .as_ref()
                .map(|p| p.to_string_lossy().to_string()),
        };
        let prepared = prepare_indexes(&request)?;
        if let Some(idx) = prepared.input_index
            && loader.input_index.is_none()
        {
            eprintln!("bioscript: auto-indexed input -> {}", idx.display());
            loader.input_index = Some(idx);
        }
        if let Some(ref_file) = prepared.reference_file {
            loader.reference_file = Some(ref_file);
        }
        if let Some(ref_idx) = prepared.reference_index
            && loader.reference_index.is_none()
        {
            eprintln!("bioscript: auto-indexed reference -> {}", ref_idx.display());
            loader.reference_index = Some(ref_idx);
        }
        cli_timings.push(StageTiming {
            stage: "auto_index".to_owned(),
            duration_ms: auto_index_started.elapsed().as_millis(),
            detail: "prepare_indexes".to_owned(),
        });
    }

    if is_yaml_manifest(&script_path) {
        let manifest_started = Instant::now();
        let manifest_options = ManifestRunOptions {
            input_file: input_file.as_deref(),
            output_file: output_file.as_deref(),
            participant_id: participant_id.as_deref(),
            trace_report: trace_report.as_deref(),
            loader: &loader,
            filters: &filters,
        };
        run_manifest(&runtime_root, &script_path, &manifest_options)?;
        cli_timings.push(StageTiming {
            stage: "manifest_run".to_owned(),
            duration_ms: manifest_started.elapsed().as_millis(),
            detail: script_path.display().to_string(),
        });
        if let Some(timing_path) = timing_report {
            write_timing_report(&timing_path, &cli_timings)?;
        }
        return Ok(());
    }

    let runtime = BioscriptRuntime::with_config(runtime_root, RuntimeConfig { limits, loader })
        .map_err(|err| err.to_string())?;
    let mut inputs = Vec::new();
    if let Some(input_file) = input_file {
        inputs.push(("input_file", monty::MontyObject::String(input_file)));
    }
    if let Some(output_file) = output_file {
        inputs.push(("output_file", monty::MontyObject::String(output_file)));
    }
    if let Some(participant_id) = participant_id {
        inputs.push(("participant_id", monty::MontyObject::String(participant_id)));
    }

    runtime
        .run_file(&script_path, trace_report.as_deref(), inputs)
        .map_err(|err| err.to_string())?;
    if let Some(timing_path) = timing_report {
        let mut all_timings = cli_timings;
        all_timings.extend(runtime.timing_snapshot());
        write_timing_report(&timing_path, &all_timings)?;
    }
    Ok(())
}

fn write_timing_report(path: &PathBuf, timings: &[StageTiming]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|err| {
            format!(
                "failed to create timing report dir {}: {err}",
                parent.display()
            )
        })?;
    }
    let mut output = String::from("stage\tduration_ms\tdetail\n");
    for timing in timings {
        let _ = writeln!(
            output,
            "{}\t{}\t{}",
            timing.stage,
            timing.duration_ms,
            timing.detail.replace('\t', " ")
        );
    }
    fs::write(path, output)
        .map_err(|err| format!("failed to write timing report {}: {err}", path.display()))
}

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
            "usage: bioscript inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>]"
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

#[allow(clippy::too_many_lines)]
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

fn run_manifest_rows_for_report(
    runtime_root: &Path,
    manifest_path: &Path,
    input_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_text = input_file.display().to_string();
    match manifest_schema(manifest_path)?.as_str() {
        "bioscript:variant:1.0" | "bioscript:variant" => {
            let manifest = load_variant_manifest(manifest_path)?;
            Ok(vec![run_variant_manifest(
                runtime_root,
                &manifest,
                Some(&input_text),
                Some(participant_id),
                loader,
            )?])
        }
        "bioscript:panel:1.0" => {
            let manifest = load_panel_manifest(manifest_path)?;
            run_panel_manifest(
                runtime_root,
                &manifest,
                Some(&input_text),
                Some(participant_id),
                loader,
                filters,
            )
        }
        "bioscript:assay:1.0" => {
            let manifest = load_assay_manifest(manifest_path)?;
            run_assay_manifest(
                runtime_root,
                &manifest,
                Some(&input_text),
                Some(participant_id),
                loader,
                filters,
            )
        }
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

fn run_manifest_analyses_for_report(
    runtime_root: &Path,
    manifest_path: &Path,
    input_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
    output_dir: &Path,
) -> Result<Vec<serde_json::Value>, String> {
    match manifest_schema(manifest_path)?.as_str() {
        "bioscript:panel:1.0" => {
            let manifest = load_panel_manifest(manifest_path)?;
            let mut analyses = Vec::new();
            analyses.extend(run_interpretations_for_report(
                runtime_root,
                &manifest.path,
                &manifest.name,
                &manifest.interpretations,
                input_file,
                participant_id,
                loader,
                output_dir,
            )?);
            for member in &manifest.members {
                if member.kind != "assay" {
                    continue;
                }
                let Some(path) = &member.path else {
                    continue;
                };
                let resolved = resolve_manifest_path(runtime_root, &manifest.path, path)?;
                analyses.extend(run_manifest_analyses_for_report(
                    runtime_root,
                    &resolved,
                    input_file,
                    participant_id,
                    loader,
                    output_dir,
                )?);
            }
            Ok(analyses)
        }
        "bioscript:assay:1.0" => {
            let manifest = load_assay_manifest(manifest_path)?;
            run_interpretations_for_report(
                runtime_root,
                &manifest.path,
                &manifest.name,
                &manifest.interpretations,
                input_file,
                participant_id,
                loader,
                output_dir,
            )
        }
        "bioscript:variant:1.0" | "bioscript:variant" => Ok(Vec::new()),
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

#[allow(clippy::too_many_arguments)]
fn run_interpretations_for_report(
    runtime_root: &Path,
    manifest_path: &Path,
    manifest_name: &str,
    interpretations: &[PanelInterpretation],
    input_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
    output_dir: &Path,
) -> Result<Vec<serde_json::Value>, String> {
    let mut outputs = Vec::new();
    for interpretation in interpretations {
        if interpretation.kind != "bioscript" {
            return Err(format!(
                "analysis '{}' uses unsupported kind '{}'",
                interpretation.id, interpretation.kind
            ));
        }
        let script_path = resolve_manifest_path(runtime_root, manifest_path, &interpretation.path)?;
        let format = interpretation
            .output_format
            .as_deref()
            .unwrap_or("json")
            .to_ascii_lowercase();
        let analysis_dir = output_dir.join("analysis").join(participant_id);
        fs::create_dir_all(&analysis_dir).map_err(|err| {
            format!(
                "failed to create analysis output dir {}: {err}",
                analysis_dir.display()
            )
        })?;
        let extension = match format.as_str() {
            "tsv" => "tsv",
            "json" => "json",
            "jsonl" => "jsonl",
            other => return Err(format!("unsupported analysis output_format '{other}'")),
        };
        let output_file = analysis_dir.join(format!("{}.{}", interpretation.id, extension));
        run_bioscript_analysis_script(
            runtime_root,
            &script_path,
            input_file,
            &output_file,
            participant_id,
            loader,
        )?;
        let rows = parse_analysis_output(&output_file, &format)?;
        outputs.push(serde_json::json!({
            "schema": "bioscript:analysis-output:1.0",
            "version": "1.0",
            "participant_id": participant_id,
            "assay_id": manifest_name,
            "analysis_id": interpretation.id,
            "kind": interpretation.kind,
            "output_format": format,
            "manifest_path": manifest_path.strip_prefix(runtime_root).unwrap_or(manifest_path).display().to_string(),
            "script_path": script_path.strip_prefix(runtime_root).unwrap_or(&script_path).display().to_string(),
            "output_file": output_file.strip_prefix(runtime_root).unwrap_or(&output_file).display().to_string(),
            "derived_from": interpretation.derived_from.clone(),
            "emits": interpretation.emits.iter().map(|emit| serde_json::json!({
                "key": emit.key.clone(),
                "label": emit.label.clone(),
                "value_type": emit.value_type.clone(),
                "format": emit.format.clone(),
            })).collect::<Vec<_>>(),
            "logic": interpretation.logic.as_ref().map(|logic| serde_json::json!({
                "description": logic.description.clone(),
                "source": logic.source.as_ref().map(|source| serde_json::json!({
                    "name": source.name.clone(),
                    "url": source.url.clone(),
                })),
            })),
            "rows": rows,
        }));
    }
    Ok(outputs)
}

fn run_bioscript_analysis_script(
    runtime_root: &Path,
    script_path: &Path,
    input_file: &Path,
    output_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
) -> Result<(), String> {
    let limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(1000))
        .max_memory(16 * 1024 * 1024)
        .max_allocations(400_000)
        .gc_interval(1000)
        .max_recursion_depth(Some(200));
    let runtime = BioscriptRuntime::with_config(
        runtime_root.to_path_buf(),
        RuntimeConfig {
            limits,
            loader: loader.clone(),
        },
    )
    .map_err(|err| err.to_string())?;
    runtime
        .run_file(
            script_path,
            None,
            vec![
                (
                    "input_file",
                    monty::MontyObject::String(runtime_path_string(runtime_root, input_file)),
                ),
                (
                    "output_file",
                    monty::MontyObject::String(runtime_path_string(runtime_root, output_file)),
                ),
                (
                    "participant_id",
                    monty::MontyObject::String(participant_id.to_owned()),
                ),
            ],
        )
        .map(|_| ())
        .map_err(|err| err.to_string())
}

fn runtime_path_string(runtime_root: &Path, path: &Path) -> String {
    path.strip_prefix(runtime_root)
        .unwrap_or(path)
        .display()
        .to_string()
}

fn parse_analysis_output(path: &Path, format: &str) -> Result<Vec<serde_json::Value>, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read analysis output {}: {err}", path.display()))?;
    match format {
        "tsv" => parse_analysis_tsv(&text),
        "json" => {
            let value: serde_json::Value = serde_json::from_str(&text).map_err(|err| {
                format!("failed to parse analysis JSON {}: {err}", path.display())
            })?;
            Ok(match value {
                serde_json::Value::Array(rows) => rows,
                serde_json::Value::Object(mut object) => object
                    .remove("rows")
                    .and_then(|rows| rows.as_array().cloned())
                    .unwrap_or_else(|| vec![serde_json::Value::Object(object)]),
                other => vec![other],
            })
        }
        "jsonl" => text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| serde_json::from_str(line).map_err(|err| err.to_string()))
            .collect(),
        other => Err(format!("unsupported analysis output_format '{other}'")),
    }
}

fn parse_analysis_tsv(text: &str) -> Result<Vec<serde_json::Value>, String> {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let Some(header_line) = lines.next() else {
        return Ok(Vec::new());
    };
    let headers: Vec<&str> = header_line.split('\t').collect();
    let mut rows = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split('\t').collect();
        let mut object = serde_json::Map::new();
        for (idx, header) in headers.iter().enumerate() {
            object.insert(
                (*header).to_owned(),
                serde_json::Value::String(values.get(idx).copied().unwrap_or_default().to_owned()),
            );
        }
        rows.push(serde_json::Value::Object(object));
    }
    Ok(rows)
}

fn app_assay_id(path: &Path) -> Result<String, String> {
    match manifest_schema(path)?.as_str() {
        "bioscript:panel:1.0" => Ok(load_panel_manifest(path)?.name),
        "bioscript:assay:1.0" => Ok(load_assay_manifest(path)?.name),
        "bioscript:variant:1.0" | "bioscript:variant" => Ok(load_variant_manifest(path)?.name),
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

fn participant_id_from_path(path: &Path) -> String {
    let file_name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("participant");
    file_name
        .trim_end_matches(".txt.zip")
        .trim_end_matches(".csv.zip")
        .trim_end_matches(".vcf.gz")
        .trim_end_matches(".cram")
        .trim_end_matches(".zip")
        .trim_end_matches(".txt")
        .trim_end_matches(".csv")
        .to_owned()
}

fn app_observation_from_manifest_row(
    runtime_root: &Path,
    row: &BTreeMap<String, String>,
    assay_id: &str,
) -> Result<serde_json::Value, String> {
    let row_path = row.get("path").cloned().unwrap_or_default();
    let manifest_path = if Path::new(&row_path).is_absolute() {
        PathBuf::from(&row_path)
    } else {
        runtime_root.join(&row_path)
    };
    let manifest = load_variant_manifest(&manifest_path)?;
    let ref_allele = manifest.spec.reference.clone().unwrap_or_default();
    let genotype_display = row.get("genotype").cloned().unwrap_or_default();
    let alt_alleles = variant_alt_alleles(&manifest_path)?;
    let alt_allele = observed_alt_allele(&genotype_display, &ref_allele, &alt_alleles)
        .or_else(|| manifest.spec.alternate.clone())
        .unwrap_or_default();
    let (genotype, zygosity) = normalize_app_genotype(&genotype_display, &ref_allele, &alt_allele);
    let depth = parse_optional_u32(row.get("depth"));
    let ref_count = parse_optional_u32(row.get("ref_count"));
    let alt_count = parse_optional_u32(row.get("alt_count"));
    let allele_balance = match (alt_count, depth) {
        (Some(alt_count), Some(depth)) if depth > 0 => {
            Some(f64::from(alt_count) / f64::from(depth))
        }
        _ => None,
    };
    let assembly = row.get("assembly").cloned().unwrap_or_default();
    let locus = if assembly.eq_ignore_ascii_case("grch37") {
        manifest.spec.grch37.as_ref()
    } else {
        manifest
            .spec
            .grch38
            .as_ref()
            .or(manifest.spec.grch37.as_ref())
    };
    let outcome = if genotype == "./." {
        "no_call"
    } else if zygosity == "hom_ref" {
        "reference"
    } else if zygosity == "het" || zygosity == "hom_alt" {
        "variant"
    } else {
        "unknown"
    };
    let evidence_raw = row.get("evidence").cloned().unwrap_or_default();
    Ok(serde_json::json!({
        "participant_id": row.get("participant_id").cloned().unwrap_or_default(),
        "assay_id": assay_id,
        "assay_version": "1.0",
        "variant_key": manifest.name,
        "variant_path": row_path,
        "rsid": row.get("matched_rsid").filter(|value| !value.is_empty()).cloned().or_else(|| manifest.spec.rsids.first().cloned()),
        "assembly": if assembly.is_empty() { serde_json::Value::Null } else { serde_json::Value::String(assembly.to_uppercase()) },
        "chrom": locus.map_or(String::new(), |locus| locus.chrom.clone()),
        "pos_start": locus.map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.start)),
        "pos_end": locus.map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.end)),
        "ref": ref_allele,
        "alt": alt_allele,
        "kind": manifest.spec.kind.map_or("unknown".to_owned(), |kind| format!("{kind:?}").to_lowercase()),
        "match_status": if row.get("matched_rsid").is_some_and(|value| !value.is_empty()) || !genotype_display.is_empty() { "found" } else { "not_found" },
        "coverage_status": depth.map_or("covered", |depth| if depth > 0 { "covered" } else { "not_covered" }),
        "call_status": if genotype == "./." { "no_call" } else { "called" },
        "genotype": genotype,
        "genotype_display": genotype_display,
        "zygosity": zygosity,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "depth": depth,
        "genotype_quality": serde_json::Value::Null,
        "allele_balance": allele_balance,
        "outcome": outcome,
        "evidence_type": if row.get("backend").is_some_and(|value| value == "cram") { "mpileup" } else { "genotype_file" },
        "evidence_raw": evidence_raw,
        "facets": serde_json::Value::Null,
    }))
}

fn variant_alt_alleles(path: &Path) -> Result<Vec<String>, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read variant YAML {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse variant YAML {}: {err}", path.display()))?;
    let Some(items) = value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alleles".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| {
            mapping
                .get(serde_yaml::Value::String("observed_alts".to_owned()))
                .or_else(|| mapping.get(serde_yaml::Value::String("alts".to_owned())))
        })
        .and_then(serde_yaml::Value::as_sequence)
    else {
        return Ok(Vec::new());
    };
    Ok(items
        .iter()
        .filter_map(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
        .collect())
}

fn observed_alt_allele(
    genotype_display: &str,
    ref_allele: &str,
    alts: &[String],
) -> Option<String> {
    if ref_allele.len() != 1 {
        return None;
    }
    let ref_ch = ref_allele.chars().next()?;
    genotype_display
        .chars()
        .filter(|ch| ch.is_ascii_alphabetic() && *ch != ref_ch)
        .find_map(|ch| {
            alts.iter()
                .find(|alt| alt.len() == 1 && alt.starts_with(ch))
                .cloned()
        })
}

fn normalize_app_genotype(display: &str, ref_allele: &str, alt_allele: &str) -> (String, String) {
    if display.is_empty() {
        return ("./.".to_owned(), "unknown".to_owned());
    }
    let alleles: Vec<char> = display
        .chars()
        .filter(|ch| ch.is_ascii_alphabetic())
        .collect();
    if alleles.len() != 2 || ref_allele.len() != 1 || alt_allele.len() != 1 {
        return (display.to_owned(), "unknown".to_owned());
    }
    let ref_ch = ref_allele.chars().next().unwrap_or_default();
    let alt_ch = alt_allele.chars().next().unwrap_or_default();
    let alt_count = alleles.iter().filter(|allele| **allele == alt_ch).count();
    let ref_count = alleles.iter().filter(|allele| **allele == ref_ch).count();
    match (ref_count, alt_count) {
        (2, 0) => ("0/0".to_owned(), "hom_ref".to_owned()),
        (1, 1) => ("0/1".to_owned(), "het".to_owned()),
        (0, 2) => ("1/1".to_owned(), "hom_alt".to_owned()),
        _ => (display.to_owned(), "unknown".to_owned()),
    }
}

fn parse_optional_u32(value: Option<&String>) -> Option<u32> {
    value.and_then(|value| value.parse::<u32>().ok())
}

fn load_manifest_findings(
    root: &Path,
    manifest_path: &Path,
) -> Result<Vec<serde_json::Value>, String> {
    let value = load_yaml_value(manifest_path)?;
    let schema = value
        .get("schema")
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or_default();
    let mut findings = Vec::new();

    if matches!(
        schema,
        "bioscript:variant:1.0"
            | "bioscript:variant"
            | "bioscript:assay:1.0"
            | "bioscript:panel:1.0"
            | "bioscript:pgx-findings:1.0"
    ) {
        if let Some(items) = value
            .get("findings")
            .and_then(serde_yaml::Value::as_sequence)
        {
            for item in items {
                let json_item = yaml_to_json(item.clone())?;
                let include = json_item
                    .get("include")
                    .and_then(serde_json::Value::as_str)
                    .map(str::to_owned);
                if let Some(include) = include {
                    let include_path = resolve_manifest_path(root, manifest_path, &include)?;
                    let mut included = load_manifest_findings(root, &include_path)?;
                    let inherited_binding = json_item.get("binding").cloned();
                    for included_item in &mut included {
                        if inherited_binding.is_some()
                            && included_item.get("binding").is_none()
                            && included_item.get("effects").is_none()
                        {
                            if let Some(object) = included_item.as_object_mut() {
                                object.insert(
                                    "binding".to_owned(),
                                    inherited_binding.clone().unwrap_or(serde_json::Value::Null),
                                );
                            }
                        }
                    }
                    findings.extend(included);
                    continue;
                }
                if json_item.get("include").is_none() {
                    findings.push(json_item);
                }
            }
        }
    }

    if matches!(schema, "bioscript:assay:1.0" | "bioscript:panel:1.0") {
        if let Some(items) = value
            .get("members")
            .and_then(serde_yaml::Value::as_sequence)
        {
            for member in items {
                let Some(kind) = member.get("kind").and_then(serde_yaml::Value::as_str) else {
                    continue;
                };
                if !matches!(kind, "variant" | "assay") {
                    continue;
                }
                let Some(path) = member.get("path").and_then(serde_yaml::Value::as_str) else {
                    continue;
                };
                let member_path = resolve_manifest_path(root, manifest_path, path)?;
                findings.extend(load_manifest_findings(root, &member_path)?);
            }
        }
    }

    Ok(findings)
}

fn load_yaml_value(path: &Path) -> Result<serde_yaml::Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read YAML {}: {err}", path.display()))?;
    serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))
}

fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, String> {
    serde_json::to_value(value).map_err(|err| format!("failed to convert YAML to JSON: {err}"))
}

fn load_manifest_provenance_links(
    root: &Path,
    manifest_path: &Path,
) -> Result<Vec<serde_json::Value>, String> {
    let value = load_yaml_value(manifest_path)?;
    let schema = value
        .get("schema")
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or_default();
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    collect_manifest_provenance_entries(&value, &mut links)?;

    if matches!(
        schema,
        "bioscript:variant:1.0"
            | "bioscript:variant"
            | "bioscript:assay:1.0"
            | "bioscript:panel:1.0"
            | "bioscript:pgx-findings:1.0"
    ) {
        if let Some(items) = value
            .get("findings")
            .and_then(serde_yaml::Value::as_sequence)
        {
            for item in items {
                let json_item = yaml_to_json(item.clone())?;
                let Some(include) = json_item.get("include").and_then(serde_json::Value::as_str)
                else {
                    continue;
                };
                let include_path = resolve_manifest_path(root, manifest_path, include)?;
                for item in load_manifest_provenance_links(root, &include_path)? {
                    if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                        links.entry(url.to_owned()).or_insert(item);
                    }
                }
            }
        }
    }

    if matches!(schema, "bioscript:assay:1.0" | "bioscript:panel:1.0") {
        if let Some(items) = value
            .get("members")
            .and_then(serde_yaml::Value::as_sequence)
        {
            for member in items {
                let Some(kind) = member.get("kind").and_then(serde_yaml::Value::as_str) else {
                    continue;
                };
                if !matches!(kind, "variant" | "assay") {
                    continue;
                }
                let Some(path) = member.get("path").and_then(serde_yaml::Value::as_str) else {
                    continue;
                };
                let member_path = resolve_manifest_path(root, manifest_path, path)?;
                for item in load_manifest_provenance_links(root, &member_path)? {
                    if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                        links.entry(url.to_owned()).or_insert(item);
                    }
                }
            }
        }
    }

    Ok(links.into_values().collect())
}

fn collect_manifest_provenance_entries(
    value: &serde_yaml::Value,
    links: &mut BTreeMap<String, serde_json::Value>,
) -> Result<(), String> {
    if let Some(sources) = value
        .get("provenance")
        .and_then(|provenance| provenance.get("sources"))
        .and_then(serde_yaml::Value::as_sequence)
    {
        for source in sources {
            let json = yaml_to_json(source.clone())?;
            if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
                links.entry(url.to_owned()).or_insert(json);
            }
        }
    }
    if let Some(source) = value.get("source") {
        let json = yaml_to_json(source.clone())?;
        if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
            links.entry(url.to_owned()).or_insert(json);
        }
    }
    Ok(())
}

fn match_app_findings(
    findings: &[serde_json::Value],
    observations: &[serde_json::Value],
    analyses: &[serde_json::Value],
) -> Vec<serde_json::Value> {
    let mut matched = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for finding in findings {
        if let Some(effects) = finding.get("effects").and_then(serde_json::Value::as_array) {
            for effect in effects {
                if let Some(observation) = app_finding_match_observation(effect, observations) {
                    let mut item = finding.clone();
                    if let Some(object) = item.as_object_mut() {
                        object.remove("effects");
                        object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                        object.insert("matched_effect".to_owned(), effect.clone());
                        object.insert(
                            "matched_observation".to_owned(),
                            app_finding_observation_context(observation),
                        );
                    }
                    let key = app_finding_dedupe_key(&item);
                    if seen.insert(key) {
                        matched.push(item);
                    }
                } else if let Some(analysis) = app_finding_match_analysis(effect, analyses) {
                    let mut item = finding.clone();
                    if let Some(object) = item.as_object_mut() {
                        object.remove("effects");
                        object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                        object.insert("matched_effect".to_owned(), effect.clone());
                        object.insert("matched_analysis".to_owned(), analysis);
                    }
                    let key = app_finding_dedupe_key(&item);
                    if seen.insert(key) {
                        matched.push(item);
                    }
                }
            }
        } else if let Some(observation) = app_finding_match_observation(finding, observations) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert(
                    "matched_observation".to_owned(),
                    app_finding_observation_context(observation),
                );
            }
            let key = app_finding_dedupe_key(&item);
            if seen.insert(key) {
                matched.push(item);
            }
        } else if let Some(analysis) = app_finding_match_analysis(finding, analyses) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert("matched_analysis".to_owned(), analysis);
            }
            let key = app_finding_dedupe_key(&item);
            if seen.insert(key) {
                matched.push(item);
            }
        }
    }
    matched
}

fn app_finding_match_observation<'a>(
    finding: &serde_json::Value,
    observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let Some(binding) = finding.get("binding") else {
        return None;
    };
    match binding.get("source").and_then(serde_json::Value::as_str) {
        Some("variant") => app_variant_binding_match_observation(binding, observations),
        _ => None,
    }
}

fn app_finding_match_analysis(
    finding: &serde_json::Value,
    analyses: &[serde_json::Value],
) -> Option<serde_json::Value> {
    let binding = finding.get("binding")?;
    if binding.get("source").and_then(serde_json::Value::as_str) != Some("analysis") {
        return None;
    }
    let analysis_id = binding
        .get("analysis_id")
        .or_else(|| binding.get("analysis"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let key = binding.get("key").and_then(serde_json::Value::as_str)?;
    for analysis in analyses {
        if !analysis_id.is_empty()
            && analysis
                .get("analysis_id")
                .and_then(serde_json::Value::as_str)
                != Some(analysis_id)
        {
            continue;
        }
        let Some(rows) = analysis.get("rows").and_then(serde_json::Value::as_array) else {
            continue;
        };
        for row in rows {
            if app_binding_matches_value(row.get(key), binding) {
                return Some(serde_json::json!({
                    "participant_id": analysis.get("participant_id").cloned().unwrap_or(serde_json::Value::Null),
                    "assay_id": analysis.get("assay_id").cloned().unwrap_or(serde_json::Value::Null),
                    "analysis_id": analysis.get("analysis_id").cloned().unwrap_or(serde_json::Value::Null),
                    "key": key,
                    "value": row.get(key).cloned().unwrap_or(serde_json::Value::Null),
                    "row": row,
                }));
            }
        }
    }
    None
}

fn app_variant_binding_match_observation<'a>(
    binding: &serde_json::Value,
    observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let operator = binding
        .get("operator")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("equals");
    if matches!(operator, "dosage_equals" | "dosage_in") {
        let allele = binding
            .get("allele")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        return observations
            .iter()
            .filter(|observation| !app_variant_ref_mismatch(binding, observation))
            .find(|observation| {
                let dosage = app_observation_allele_dosage(observation, allele);
                app_binding_matches_dosage(dosage, binding)
            });
    }

    let key = binding
        .get("key")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if key.is_empty() {
        return None;
    }
    observations
        .iter()
        .filter(|observation| !app_variant_ref_mismatch(binding, observation))
        .find(|observation| app_binding_matches_value(observation.get(key), binding))
}

fn app_finding_observation_context(observation: &serde_json::Value) -> serde_json::Value {
    serde_json::json!({
        "participant_id": observation.get("participant_id").cloned().unwrap_or(serde_json::Value::Null),
        "rsid": observation.get("rsid").cloned().unwrap_or(serde_json::Value::Null),
        "ref": observation.get("ref").cloned().unwrap_or(serde_json::Value::Null),
        "alt": observation.get("alt").cloned().unwrap_or(serde_json::Value::Null),
        "genotype_display": observation.get("genotype_display").cloned().unwrap_or(serde_json::Value::Null),
        "outcome": observation.get("outcome").cloned().unwrap_or(serde_json::Value::Null),
    })
}

fn app_variant_ref_mismatch(binding: &serde_json::Value, observation: &serde_json::Value) -> bool {
    let variant_ref = binding
        .get("variant")
        .or_else(|| binding.get("path"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if variant_ref.is_empty() {
        return false;
    }
    let basename = Path::new(variant_ref)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(variant_ref);
    let candidates = [
        observation
            .get("variant_key")
            .and_then(serde_json::Value::as_str),
        observation
            .get("variant_path")
            .and_then(serde_json::Value::as_str),
        observation.get("rsid").and_then(serde_json::Value::as_str),
    ];
    !candidates.into_iter().flatten().any(|candidate| {
        candidate == variant_ref
            || Path::new(candidate)
                .file_name()
                .and_then(|value| value.to_str())
                .is_some_and(|value| value == basename)
    })
}

fn app_observation_allele_dosage(observation: &serde_json::Value, allele: &str) -> Option<i64> {
    if allele.is_empty() {
        return None;
    }
    let ref_allele = observation
        .get("ref")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let alt_allele = observation
        .get("alt")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let zygosity = observation
        .get("zygosity")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if allele == ref_allele {
        return match zygosity {
            "hom_ref" => Some(2),
            "het" => Some(1),
            "hom_alt" => Some(0),
            _ => None,
        };
    }
    if allele == alt_allele {
        return match zygosity {
            "hom_ref" => Some(0),
            "het" => Some(1),
            "hom_alt" => Some(2),
            _ => None,
        };
    }
    let display = observation
        .get("genotype_display")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if allele.len() == 1 {
        let allele_ch = allele.chars().next()?.to_ascii_uppercase();
        return Some(
            display
                .chars()
                .filter(|ch| ch.to_ascii_uppercase() == allele_ch)
                .count()
                .try_into()
                .ok()?,
        );
    }
    None
}

fn app_binding_matches_value(
    actual: Option<&serde_json::Value>,
    binding: &serde_json::Value,
) -> bool {
    let actual = actual.and_then(value_as_string).unwrap_or_default();
    match binding
        .get("operator")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("equals")
    {
        "equals" => binding
            .get("value")
            .and_then(value_as_string)
            .is_some_and(|value| value == actual),
        "in" => binding
            .get("values")
            .and_then(serde_json::Value::as_array)
            .is_some_and(|values| {
                values
                    .iter()
                    .filter_map(value_as_string)
                    .any(|value| value == actual)
            }),
        _ => false,
    }
}

fn app_binding_matches_dosage(dosage: Option<i64>, binding: &serde_json::Value) -> bool {
    let Some(dosage) = dosage else {
        return false;
    };
    match binding
        .get("operator")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "dosage_equals" => binding
            .get("value")
            .and_then(serde_json::Value::as_i64)
            .is_some_and(|value| value == dosage),
        "dosage_in" => binding
            .get("values")
            .and_then(serde_json::Value::as_array)
            .is_some_and(|values| {
                values
                    .iter()
                    .filter_map(serde_json::Value::as_i64)
                    .any(|value| value == dosage)
            }),
        _ => false,
    }
}

fn value_as_string(value: &serde_json::Value) -> Option<String> {
    match value {
        serde_json::Value::String(value) => Some(value.clone()),
        serde_json::Value::Number(value) => Some(value.to_string()),
        serde_json::Value::Bool(value) => Some(value.to_string()),
        _ => None,
    }
}

fn app_finding_dedupe_key(finding: &serde_json::Value) -> String {
    let effect_key = finding
        .get("matched_effect")
        .and_then(|effect| {
            effect
                .get("id")
                .or_else(|| effect.get("label"))
                .or_else(|| effect.get("text"))
        })
        .and_then(value_as_string)
        .unwrap_or_default();
    if let Some(evidence) = finding.get("evidence") {
        let source = evidence
            .get("source")
            .and_then(value_as_string)
            .unwrap_or_default();
        let kind = evidence
            .get("kind")
            .and_then(value_as_string)
            .unwrap_or_default();
        let id = evidence
            .get("id")
            .and_then(value_as_string)
            .unwrap_or_default();
        if !source.is_empty() || !kind.is_empty() || !id.is_empty() {
            return format!("evidence|{source}|{kind}|{id}|{effect_key}");
        }
        if let Some(url) = evidence.get("url").and_then(value_as_string) {
            return format!("evidence_url|{url}|{effect_key}");
        }
    }
    if let Some(id) = finding.get("id").and_then(value_as_string) {
        return format!("id|{id}|{effect_key}");
    }
    format!(
        "content|{}|{}|{}|{}",
        finding
            .get("schema")
            .and_then(value_as_string)
            .unwrap_or_default(),
        finding
            .get("label")
            .and_then(value_as_string)
            .unwrap_or_default(),
        finding
            .get("notes")
            .and_then(value_as_string)
            .unwrap_or_default(),
        effect_key
    )
}

fn app_report_json(
    assay_id: &str,
    participant_id: &str,
    input_file: &Path,
    observations: &[serde_json::Value],
    analyses: &[serde_json::Value],
    findings: &[serde_json::Value],
    provenance: &[serde_json::Value],
) -> serde_json::Value {
    let called = observations
        .iter()
        .filter(|item| {
            item.get("call_status").and_then(serde_json::Value::as_str) == Some("called")
        })
        .count();
    serde_json::json!({
        "schema": "bioscript:report:1.0",
        "version": "1.0",
        "participant_id": participant_id,
        "assay_id": assay_id,
        "assay_version": "1.0",
        "input": {
            "file_name": input_file.file_name().and_then(|value| value.to_str()).unwrap_or_default(),
            "file_path": input_file.display().to_string(),
        },
        "report_status": if called == observations.len() { "complete" } else { "partial" },
        "derived_from": observations.iter().filter_map(|item| item.get("variant_key").cloned()).collect::<Vec<_>>(),
        "analyses": analyses,
        "findings": findings,
        "provenance": provenance,
        "metrics": {
            "n_sites_tested": observations.len(),
            "n_sites_called": called,
            "n_sites_missing": observations.len().saturating_sub(called),
            "n_analyses": analyses.len(),
            "n_findings_matched": findings.len(),
        }
    })
}

fn write_app_observations(
    output_dir: &Path,
    observations: &[serde_json::Value],
    format: AppOutputFormat,
) -> Result<(), String> {
    if matches!(format, AppOutputFormat::Tsv | AppOutputFormat::Both) {
        let mut out = bioscript_core::OBSERVATION_TSV_HEADERS.join("\t");
        out.push('\n');
        for observation in observations {
            let line = bioscript_core::OBSERVATION_TSV_HEADERS
                .iter()
                .map(|header| json_field_as_tsv(observation.get(*header)))
                .collect::<Vec<_>>()
                .join("\t");
            out.push_str(&line);
            out.push('\n');
        }
        fs::write(output_dir.join("observations.tsv"), out)
            .map_err(|err| format!("failed to write observations.tsv: {err}"))?;
    }
    if matches!(format, AppOutputFormat::Jsonl | AppOutputFormat::Both) {
        write_jsonl(&output_dir.join("observations.jsonl"), observations)?;
    }
    if matches!(format, AppOutputFormat::Json) {
        write_json_pretty(
            &output_dir.join("observations.json"),
            &serde_json::json!({"observations": observations}),
        )?;
    }
    Ok(())
}

fn write_app_analyses(output_dir: &Path, analyses: &[serde_json::Value]) -> Result<(), String> {
    write_jsonl(&output_dir.join("analysis.jsonl"), analyses)
}

fn write_app_reports(
    output_dir: &Path,
    reports: &[serde_json::Value],
    format: AppOutputFormat,
) -> Result<(), String> {
    if matches!(format, AppOutputFormat::Jsonl | AppOutputFormat::Both) {
        write_jsonl(&output_dir.join("reports.jsonl"), reports)?;
    }
    if matches!(format, AppOutputFormat::Json | AppOutputFormat::Both) {
        write_json_pretty(
            &output_dir.join("reports.json"),
            &serde_json::json!({
                "schema": "bioscript:report-set:1.0",
                "version": "1.0",
                "reports": reports,
            }),
        )?;
    }
    Ok(())
}

fn write_jsonl(path: &Path, rows: &[serde_json::Value]) -> Result<(), String> {
    let mut out = String::new();
    for row in rows {
        let line = serde_json::to_string(row).map_err(|err| err.to_string())?;
        out.push_str(&line);
        out.push('\n');
    }
    fs::write(path, out).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn write_json_pretty(path: &Path, value: &serde_json::Value) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|err| err.to_string())?;
    fs::write(path, text).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}

fn write_app_html(
    output_dir: &Path,
    observations: &[serde_json::Value],
    reports: &[serde_json::Value],
) -> Result<(), String> {
    let mut out = String::from(
        r##"<!doctype html><meta charset="utf-8"><title>BioScript report</title><style>body{font-family:system-ui,sans-serif;margin:0;background:#f7f8fa;color:#1f2933}.wrap{max-width:1440px;margin:0 auto;padding:24px}h1{margin:0 0 10px}h2{margin:32px 0 10px;scroll-margin-top:82px}.nav{position:sticky;top:0;z-index:20;display:flex;gap:8px;flex-wrap:wrap;align-items:center;margin:16px -24px 22px;padding:10px 24px;background:rgba(247,248,250,.96);border-block:1px solid #d8dee6;backdrop-filter:saturate(160%) blur(8px)}.nav a{border:1px solid #cbd5df;background:#fff;color:#1f2933;text-decoration:none;padding:7px 10px;border-radius:6px}.nav a:hover{background:#eef2f6}.table-tools,.level-filter{display:flex;justify-content:space-between;gap:12px;align-items:center;margin:6px 0}.table-tools input{width:min(420px,100%);border:1px solid #cbd5df;border-radius:6px;padding:7px 9px;font:inherit;background:#fff}.table-tools label,.level-filter label{display:flex;gap:5px;align-items:center}.level-filter{justify-content:flex-start;flex-wrap:wrap}.level-filter a{display:inline-grid;place-items:center;width:20px;height:20px;border:1px solid #cbd5df;border-radius:50%;text-decoration:none;color:#1f2933;background:#fff;font-weight:700}.filter-action{border:1px solid #cbd5df;background:#fff;color:#1f2933;border-radius:6px;padding:4px 7px;font:inherit;cursor:pointer}.filter-action:hover{background:#eef2f6}.table-wrap{overflow:auto;border:1px solid #d8dee6;background:white;border-radius:8px}table{border-collapse:collapse;width:100%;font-size:13px}td,th{border-bottom:1px solid #e5e9ef;padding:6px 8px;text-align:left;vertical-align:top}th{position:sticky;top:0;background:#eef2f6;z-index:1;white-space:nowrap;cursor:pointer;user-select:none}.sort-mark{font-size:10px;color:#667085;margin-left:3px}.row-variant td{background:#fff7cc}.row-reference td{background:#eaf7ee}.refs-hidden .row-reference{display:none}.genotype-hit{font-weight:700}.allele-hit{color:#075985;background:#dff4ff;border-radius:3px;padding:0 2px}.level-badge,.pgx-badge{display:inline-block;min-width:2.2em;text-align:center;border-radius:999px;padding:2px 8px;color:#fff;font-weight:700}.level-1,.pgx-informative{background:#0abc72}.level-2,.pgx-actionable{background:#2a74df}.level-3,.pgx-recommended{background:#ffc107;color:#1f2933}.level-4,.pgx-required{background:#c53b3b}.pgx-no-clinical,.pgx-criteria{background:#667085}.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}.muted{color:#667085}.effect{max-width:760px;min-width:360px}.logic-note{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 12px;margin:8px 0 10px}.logic-note p{margin:0 0 6px}.provenance-list{background:#fff;border:1px solid #d8dee6;border-radius:8px;padding:10px 18px}.provenance-list li{margin:8px 0}pre{white-space:pre-wrap;background:#fff;padding:12px;border:1px solid #d8dee6;border-radius:8px;max-height:520px;overflow:auto}</style><script>const sortState={};function cellText(row,i){return(row.cells[i]?.innerText||"").trim()}function cmp(a,b){const an=Number(a),bn=Number(b);if(a!==""&&b!==""&&!Number.isNaN(an)&&!Number.isNaN(bn))return an-bn;return a.localeCompare(b,undefined,{numeric:true,sensitivity:"base"})}function sortTable(id,col){const table=document.getElementById(id);const tbody=table.tBodies[0];const key=id+":"+col;const dir=sortState[key]==="asc"?"desc":"asc";sortState[key]=dir;table.querySelectorAll(".sort-mark").forEach(s=>s.textContent="");table.tHead.rows[0].cells[col].querySelector(".sort-mark").textContent=dir==="asc"?"^":"v";Array.from(tbody.rows).sort((a,b)=>{const v=cmp(cellText(a,col),cellText(b,col));return dir==="asc"?v:-v}).forEach(r=>tbody.appendChild(r))}function applyTableFilters(id){const q=(document.querySelector('[data-filter-for="'+id+'"]')?.value||"").toLowerCase();document.querySelectorAll("#"+id+" tbody tr").forEach(row=>{let ok=row.innerText.toLowerCase().includes(q);if(id==="summaries-table"&&row.dataset.level){ok=ok&&!!document.querySelector('[data-level-filter="'+row.dataset.level+'"]:checked')}if(id==="labels-table"&&row.dataset.pgxLevel){ok=ok&&!!document.querySelector('[data-pgx-level-filter="'+row.dataset.pgxLevel+'"]:checked')}row.style.display=ok?"":"none"})}function setFilterGroup(group,id,checked){document.querySelectorAll(group==="level"?"[data-level-filter]":"[data-pgx-level-filter]").forEach(input=>input.checked=checked);applyTableFilters(id)}function toggleRefs(show){document.getElementById("report-wrap").classList.toggle("refs-hidden",!show)}</script><div class="wrap refs-hidden" id="report-wrap"><h1>BioScript Report</h1>"##,
    );
    let label_findings = collect_report_findings(reports, "bioscript:pgx-label:1.0");
    let summary_findings = collect_report_findings(reports, "bioscript:pgx-summary:1.0");
    let analysis_outputs = collect_report_analyses(reports);
    let _ = write!(
        out,
        "<div class=\"muted\">{} observation(s), {} analysis output(s), {} PGx label finding(s), {} PGx summary finding(s)</div>",
        observations.len(),
        analysis_outputs.len(),
        label_findings.len(),
        summary_findings.len()
    );
    out.push_str("<nav class=\"nav\"><a href=\"#observations\">Observations</a><a href=\"#analysis\">Analysis</a><a href=\"#labels\">PGx Labels</a><a href=\"#summaries\">PGx Summaries</a><a href=\"#provenance\">Provenance</a><a href=\"#json\">Raw JSON</a></nav>");
    out.push_str("<section id=\"observations\"><h2>Observations</h2>");
    render_observation_table(&mut out, observations);
    out.push_str("</section>");
    out.push_str("<section id=\"analysis\"><h2>Analysis</h2>");
    render_analysis_tables(&mut out, &analysis_outputs);
    out.push_str("</section>");
    out.push_str("<section id=\"labels\"><h2>PGx Label Annotations</h2>");
    render_pgx_label_table(&mut out, &label_findings);
    out.push_str("</section>");
    out.push_str("<section id=\"summaries\"><h2>PGx Summary Annotations</h2>");
    render_pgx_summary_table(&mut out, &summary_findings);
    out.push_str("</section>");
    out.push_str("<section id=\"provenance\"><h2>Provenance</h2>");
    render_provenance_links(&mut out, reports);
    out.push_str("</section>");
    out.push_str("<section id=\"json\"><h2>Raw Reports JSON</h2>");
    for report in reports {
        let text = serde_json::to_string_pretty(report).map_err(|err| err.to_string())?;
        let _ = write!(out, "<pre>{}</pre>", html_escape(&text));
    }
    out.push_str("</section></div>");
    fs::write(output_dir.join("index.html"), out)
        .map_err(|err| format!("failed to write index.html: {err}"))
}

fn collect_report_analyses(reports: &[serde_json::Value]) -> Vec<serde_json::Value> {
    reports
        .iter()
        .filter_map(|report| report.get("analyses").and_then(serde_json::Value::as_array))
        .flat_map(|analyses| analyses.iter())
        .cloned()
        .collect()
}

fn collect_report_findings(reports: &[serde_json::Value], schema: &str) -> Vec<serde_json::Value> {
    reports
        .iter()
        .filter_map(|report| report.get("findings").and_then(serde_json::Value::as_array))
        .flat_map(|findings| findings.iter())
        .filter(|finding| finding.get("schema").and_then(serde_json::Value::as_str) == Some(schema))
        .cloned()
        .collect()
}

fn render_analysis_tables(out: &mut String, analyses: &[serde_json::Value]) {
    if analyses.is_empty() {
        out.push_str("<p class=\"muted\">No analysis outputs.</p>");
        return;
    }
    for (index, analysis) in analyses.iter().enumerate() {
        let table_id = format!("analysis-table-{index}");
        let title = format!(
            "{} / {}",
            value_str(analysis, "participant_id"),
            value_str(analysis, "analysis_id")
        );
        let _ = write!(out, "<h3>{}</h3>", html_escape(&title));
        render_analysis_logic(out, analysis);
        let rows = analysis
            .get("rows")
            .and_then(serde_json::Value::as_array)
            .cloned()
            .unwrap_or_default();
        if rows.is_empty() {
            out.push_str("<p class=\"muted\">No rows emitted.</p>");
            continue;
        }
        let headers = analysis_row_headers(&rows);
        let header_refs = headers.iter().map(String::as_str).collect::<Vec<_>>();
        render_table_start(out, &table_id, &header_refs);
        for row in rows {
            out.push_str("<tr>");
            for header in &headers {
                table_cell(out, &json_field_as_tsv(row.get(header)));
            }
            out.push_str("</tr>");
        }
        render_table_end(out);
    }
}

fn analysis_row_headers(rows: &[serde_json::Value]) -> Vec<String> {
    let mut headers = Vec::new();
    for row in rows {
        let Some(object) = row.as_object() else {
            continue;
        };
        for key in object.keys() {
            if !headers.contains(key) {
                headers.push(key.clone());
            }
        }
    }
    headers
}

fn render_analysis_logic(out: &mut String, analysis: &serde_json::Value) {
    let Some(logic) = analysis.get("logic") else {
        return;
    };
    if logic.is_null() {
        return;
    }
    let description = logic
        .get("description")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let source = logic.get("source").unwrap_or(&serde_json::Value::Null);
    let source_name = source
        .get("name")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("source");
    let source_url = source
        .get("url")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    out.push_str("<div class=\"logic-note\">");
    if !description.is_empty() {
        let _ = write!(out, "<p>{}</p>", html_escape(description));
    }
    if !source_url.is_empty() {
        let _ = write!(
            out,
            "<p class=\"muted\">Logic source: <a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">{}</a></p>",
            html_escape(source_url),
            html_escape(source_name)
        );
    }
    out.push_str("</div>");
}

fn render_provenance_links(out: &mut String, reports: &[serde_json::Value]) {
    let mut links = BTreeMap::<String, String>::new();
    for report in reports {
        collect_provenance_links_from_value(report, &mut links);
    }
    if links.is_empty() {
        out.push_str("<p class=\"muted\">No provenance links.</p>");
        return;
    }
    out.push_str("<ul class=\"provenance-list\">");
    for (url, label) in links {
        let display = if label.is_empty() { &url } else { &label };
        let _ = write!(
            out,
            "<li><a href=\"{}\" target=\"_blank\" rel=\"noopener noreferrer\">{}</a><div class=\"muted mono\">{}</div></li>",
            html_escape(&url),
            html_escape(display),
            html_escape(&url)
        );
    }
    out.push_str("</ul>");
}

fn collect_provenance_links_from_value(
    value: &serde_json::Value,
    links: &mut BTreeMap<String, String>,
) {
    match value {
        serde_json::Value::Object(object) => {
            if let Some(url) = object.get("url").and_then(serde_json::Value::as_str)
                && url.starts_with("http")
            {
                let label = object
                    .get("name")
                    .or_else(|| object.get("label"))
                    .or_else(|| object.get("source"))
                    .and_then(value_as_string)
                    .unwrap_or_default();
                links.entry(url.to_owned()).or_insert(label);
            }
            for child in object.values() {
                collect_provenance_links_from_value(child, links);
            }
        }
        serde_json::Value::Array(items) => {
            for item in items {
                collect_provenance_links_from_value(item, links);
            }
        }
        _ => {}
    }
}

fn render_observation_table(out: &mut String, observations: &[serde_json::Value]) {
    let headers = [
        "participant_id",
        "rsid",
        "ref",
        "alt",
        "genotype_display",
        "genotype",
        "zygosity",
        "outcome",
        "match_status",
        "coverage_status",
        "call_status",
        "assembly",
        "chrom",
        "pos_start",
        "pos_end",
        "kind",
        "ref_count",
        "alt_count",
        "depth",
        "genotype_quality",
        "allele_balance",
        "evidence_type",
        "evidence_raw",
        "facets",
        "assay_id",
        "assay_version",
        "variant_key",
    ];
    render_table_start(out, "observations-table", &headers);
    for observation in observations {
        let _ = write!(out, "<tr class=\"{}\">", observation_row_class(observation));
        for header in headers {
            render_observation_cell(out, observation, header);
        }
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></div>");
}

fn observation_row_class(observation: &serde_json::Value) -> &'static str {
    match observation
        .get("outcome")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "variant" => "row-variant",
        "reference" => "row-reference",
        _ => "",
    }
}

fn render_observation_cell(out: &mut String, observation: &serde_json::Value, header: &str) {
    if header == "genotype_display" {
        let outcome = observation
            .get("outcome")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let value = json_field_as_tsv(observation.get(header));
        if outcome == "variant" {
            let alt = observation
                .get("alt")
                .and_then(serde_json::Value::as_str)
                .unwrap_or_default();
            let _ = write!(
                out,
                "<td class=\"genotype-hit\">{}</td>",
                highlight_allele(&value, alt)
            );
            return;
        }
    }
    let _ = write!(
        out,
        "<td>{}</td>",
        html_escape(&json_field_as_tsv(observation.get(header)))
    );
}

fn highlight_allele(value: &str, allele: &str) -> String {
    if value.is_empty() || allele.is_empty() {
        return html_escape(value);
    }
    if allele.chars().count() == 1 {
        let target = allele
            .chars()
            .next()
            .unwrap_or_default()
            .to_ascii_uppercase();
        let mut out = String::new();
        for ch in value.chars() {
            let escaped = html_escape(&ch.to_string());
            if ch.to_ascii_uppercase() == target {
                let _ = write!(out, "<span class=\"allele-hit\">{escaped}</span>");
            } else {
                out.push_str(&escaped);
            }
        }
        return out;
    }
    let escaped_value = html_escape(value);
    let escaped_allele = html_escape(allele);
    escaped_value.replace(
        &escaped_allele,
        &format!("<span class=\"allele-hit\">{escaped_allele}</span>"),
    )
}

fn render_pgx_label_table(out: &mut String, findings: &[serde_json::Value]) {
    let headers = [
        "Variant",
        "Ref/Alt",
        "Genes",
        "Drug(s)",
        "Regulator",
        "Action",
        "Label",
        "Evidence",
    ];
    render_pgx_label_filters(out);
    render_table_start(out, "labels-table", &headers);
    for finding in findings {
        let evidence = finding.get("evidence");
        let url = evidence
            .and_then(|value| value.get("url"))
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let pgx_level = value_str(finding, "pgx_action_level");
        let _ = write!(
            out,
            "<tr data-pgx-level=\"{}\">",
            html_escape(&pgx_level_slug(pgx_level))
        );
        table_cell(out, value_str(finding, "variant"));
        class_cell(out, &matched_ref_alt(finding), "mono");
        table_cell(out, &join_string_array(finding.get("genes")));
        table_cell(out, &join_drugs(finding));
        table_cell(out, &join_string_array(finding.get("regulatory_sources")));
        pgx_level_cell(out, pgx_level);
        table_cell(out, value_str(finding, "label"));
        link_cell(out, url);
        out.push_str("</tr>");
    }
    render_table_end(out);
}

fn render_pgx_summary_table(out: &mut String, findings: &[serde_json::Value]) {
    let headers = [
        "Variant",
        "Ref/Alt",
        "Genotype",
        "Drug(s)",
        "Category",
        "Level",
        "Phenotype",
        "Effect",
        "Evidence",
    ];
    render_evidence_level_filters(out);
    render_table_start(out, "summaries-table", &headers);
    for finding in findings {
        let effect = finding
            .get("matched_effect")
            .unwrap_or(&serde_json::Value::Null);
        let evidence = finding.get("evidence");
        let url = evidence
            .and_then(|value| value.get("url"))
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        let evidence_level = value_str(finding, "evidence_level");
        let _ = write!(
            out,
            "<tr data-level=\"{}\">",
            html_escape(&evidence_level_group(evidence_level))
        );
        table_cell(out, value_str(finding, "variant"));
        class_cell(out, &matched_ref_alt(finding), "mono");
        table_cell(out, value_str(effect, "label"));
        table_cell(out, &join_drugs(finding));
        table_cell(out, &join_string_array(finding.get("phenotype_categories")));
        evidence_level_cell(out, evidence_level);
        table_cell(out, &join_string_array(finding.get("phenotypes")));
        class_cell(out, value_str(effect, "text"), "effect");
        link_cell(out, url);
        out.push_str("</tr>");
    }
    render_table_end(out);
}

fn render_evidence_level_filters(out: &mut String) {
    out.push_str("<div class=\"level-filter\"><span class=\"muted\">Evidence:</span>");
    for (level, label) in [
        ("1", "Level 1"),
        ("1a", "Level 1A"),
        ("1b", "Level 1B"),
        ("2", "Level 2"),
        ("2a", "Level 2A"),
        ("2b", "Level 2B"),
        ("3", "Level 3"),
        ("4", "Level 4"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-level-filter=\"{level}\" checked onchange=\"applyTableFilters('summaries-table')\"> {label}</label>"
        );
    }
    out.push_str("<button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('level','summaries-table',true)\">Show all</button><button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('level','summaries-table',false)\">Hide all</button>");
    out.push_str("<a href=\"https://www.clinpgx.org/page/clinAnnLevels\" title=\"ClinPGx levels of evidence\">i</a></div>");
}

fn render_pgx_label_filters(out: &mut String) {
    out.push_str("<div class=\"level-filter\"><span class=\"muted\">PGx level:</span>");
    for (level, label) in [
        ("required", "Testing Required"),
        ("recommended", "Testing Recommended"),
        ("actionable", "Actionable PGx"),
        ("informative", "Informative PGx"),
        ("no-clinical", "No Clinical PGx"),
        ("criteria", "Criteria Not Met"),
    ] {
        let _ = write!(
            out,
            "<label><input type=\"checkbox\" data-pgx-level-filter=\"{level}\" checked onchange=\"applyTableFilters('labels-table')\"> {label}</label>"
        );
    }
    out.push_str("<button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('pgx-level','labels-table',true)\">Show all</button><button class=\"filter-action\" type=\"button\" onclick=\"setFilterGroup('pgx-level','labels-table',false)\">Hide all</button>");
    out.push_str("<a href=\"https://www.clinpgx.org/page/drugLabelLegend#pgx-level\" title=\"ClinPGx drug label PGx levels\">i</a></div>");
}

fn matched_ref_alt(finding: &serde_json::Value) -> String {
    let Some(observation) = finding.get("matched_observation") else {
        return String::new();
    };
    let ref_allele = observation
        .get("ref")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let alt_allele = observation
        .get("alt")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if ref_allele.is_empty() && alt_allele.is_empty() {
        String::new()
    } else {
        let alt_display = alt_allele.replace(',', "/");
        format!("{ref_allele}->{alt_display}")
    }
}

fn evidence_level_group(level: &str) -> String {
    let normalized = level.trim().to_ascii_lowercase();
    if normalized.starts_with("1a") {
        "1a".to_owned()
    } else if normalized.starts_with("1b") {
        "1b".to_owned()
    } else if normalized.starts_with('1') {
        "1".to_owned()
    } else if normalized.starts_with("2a") {
        "2a".to_owned()
    } else if normalized.starts_with("2b") {
        "2b".to_owned()
    } else if normalized.starts_with('2') {
        "2".to_owned()
    } else if normalized.starts_with('3') {
        "3".to_owned()
    } else if normalized.starts_with('4') {
        "4".to_owned()
    } else {
        "unknown".to_owned()
    }
}

fn evidence_level_color_group(level: &str) -> String {
    level
        .chars()
        .find(|ch| ch.is_ascii_digit())
        .map(|ch| ch.to_string())
        .unwrap_or_else(|| "unknown".to_owned())
}

fn evidence_level_cell(out: &mut String, level: &str) {
    if level.is_empty() {
        out.push_str("<td></td>");
        return;
    }
    let group = evidence_level_color_group(level);
    let _ = write!(
        out,
        "<td><span class=\"level-badge level-{}\">{}</span></td>",
        html_escape(&group),
        html_escape(level)
    );
}

fn pgx_level_slug(level: &str) -> String {
    let normalized = level.to_ascii_lowercase();
    if normalized.contains("required") {
        "required".to_owned()
    } else if normalized.contains("recommended") {
        "recommended".to_owned()
    } else if normalized.contains("actionable") {
        "actionable".to_owned()
    } else if normalized.contains("informative") {
        "informative".to_owned()
    } else if normalized.contains("no clinical") {
        "no-clinical".to_owned()
    } else if normalized.contains("criteria") {
        "criteria".to_owned()
    } else {
        "unknown".to_owned()
    }
}

fn pgx_level_cell(out: &mut String, level: &str) {
    if level.is_empty() {
        out.push_str("<td></td>");
        return;
    }
    let slug = pgx_level_slug(level);
    let _ = write!(
        out,
        "<td><span class=\"pgx-badge pgx-{}\">{}</span></td>",
        html_escape(&slug),
        html_escape(level)
    );
}

fn render_table_start(out: &mut String, table_id: &str, headers: &[&str]) {
    let escaped_id = html_escape(table_id);
    let refs_control = if table_id == "observations-table" {
        "<label><input type=\"checkbox\" onchange=\"toggleRefs(this.checked)\"> Show refs</label>"
    } else {
        ""
    };
    let _ = write!(
        out,
        "<div class=\"table-tools\"><input type=\"search\" placeholder=\"Filter table\" data-filter-for=\"{escaped_id}\" oninput=\"applyTableFilters('{escaped_id}')\">{refs_control}</div><div class=\"table-wrap\"><table id=\"{escaped_id}\"><thead><tr>"
    );
    for (index, header) in headers.iter().enumerate() {
        let _ = write!(
            out,
            "<th onclick=\"sortTable('{}',{})\">{}<span class=\"sort-mark\"></span></th>",
            escaped_id,
            index,
            html_escape(header)
        );
    }
    out.push_str("</tr></thead><tbody>");
}

fn render_table_end(out: &mut String) {
    out.push_str("</tbody></table></div>");
}

fn table_cell(out: &mut String, value: &str) {
    class_cell(out, value, "");
}

fn class_cell(out: &mut String, value: &str, class_name: &str) {
    if class_name.is_empty() {
        let _ = write!(out, "<td>{}</td>", html_escape(value));
    } else {
        let _ = write!(
            out,
            "<td class=\"{}\">{}</td>",
            class_name,
            html_escape(value)
        );
    }
}

fn link_cell(out: &mut String, url: &str) {
    if url.is_empty() {
        out.push_str("<td></td>");
    } else {
        let escaped = html_escape(url);
        let _ = write!(
            out,
            "<td><a href=\"{escaped}\" target=\"_blank\" rel=\"noopener noreferrer\">source</a></td>"
        );
    }
}

fn value_str<'a>(value: &'a serde_json::Value, key: &str) -> &'a str {
    value
        .get(key)
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
}

fn join_string_array(value: Option<&serde_json::Value>) -> String {
    value
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_json::Value::as_str)
                .collect::<Vec<_>>()
                .join(", ")
        })
        .unwrap_or_default()
}

fn join_drugs(finding: &serde_json::Value) -> String {
    finding
        .get("drugs")
        .and_then(serde_json::Value::as_array)
        .map(|items| {
            items
                .iter()
                .filter_map(|drug| drug.get("name").and_then(serde_json::Value::as_str))
                .collect::<Vec<_>>()
                .join(", ")
        })
        .unwrap_or_default()
}

fn html_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
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

fn is_yaml_manifest(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| matches!(ext, "yaml" | "yml"))
}

struct ManifestRunOptions<'a> {
    input_file: Option<&'a str>,
    output_file: Option<&'a str>,
    participant_id: Option<&'a str>,
    trace_report: Option<&'a Path>,
    loader: &'a GenotypeLoadOptions,
    filters: &'a [String],
}

fn run_manifest(
    runtime_root: &Path,
    manifest_path: &Path,
    options: &ManifestRunOptions<'_>,
) -> Result<(), String> {
    let schema = manifest_schema(manifest_path)?;
    let resolved_input = options
        .input_file
        .map(|value| resolve_cli_path(runtime_root, value));
    let resolved_output = options
        .output_file
        .map(|value| resolve_cli_path_buf(runtime_root, Path::new(value)));
    let resolved_trace = options
        .trace_report
        .map(|value| resolve_cli_path_buf(runtime_root, value));
    match schema.as_str() {
        "bioscript:variant:1.0" | "bioscript:variant" => {
            let manifest = load_variant_manifest(manifest_path)?;
            let row = run_variant_manifest(
                runtime_root,
                &manifest,
                resolved_input.as_deref(),
                options.participant_id,
                options.loader,
            )?;
            write_manifest_outputs(
                std::slice::from_ref(&row),
                resolved_output.as_deref(),
                resolved_trace.as_deref(),
            )?;
            Ok(())
        }
        "bioscript:panel:1.0" => {
            let manifest = load_panel_manifest(manifest_path)?;
            let rows = run_panel_manifest(
                runtime_root,
                &manifest,
                resolved_input.as_deref(),
                options.participant_id,
                options.loader,
                options.filters,
            )?;
            write_manifest_outputs(&rows, resolved_output.as_deref(), resolved_trace.as_deref())?;
            Ok(())
        }
        "bioscript:assay:1.0" => {
            let manifest = load_assay_manifest(manifest_path)?;
            let rows = run_assay_manifest(
                runtime_root,
                &manifest,
                resolved_input.as_deref(),
                options.participant_id,
                options.loader,
                options.filters,
            )?;
            write_manifest_outputs(&rows, resolved_output.as_deref(), resolved_trace.as_deref())?;
            Ok(())
        }
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

fn run_variant_manifest(
    runtime_root: &Path,
    manifest: &VariantManifest,
    input_file: Option<&str>,
    participant_id: Option<&str>,
    loader: &GenotypeLoadOptions,
) -> Result<BTreeMap<String, String>, String> {
    let input_file = input_file.ok_or("manifest execution requires --input-file")?;
    let store = GenotypeStore::from_file_with_options(Path::new(input_file), loader)
        .map_err(|err| err.to_string())?;
    let observation = store
        .lookup_variant(&manifest.spec)
        .map_err(|err| err.to_string())?;
    Ok(variant_row(
        runtime_root,
        &manifest.path,
        &manifest.name,
        &manifest.tags,
        &observation,
        participant_id,
    ))
}

fn run_panel_manifest(
    runtime_root: &Path,
    panel: &PanelManifest,
    input_file: Option<&str>,
    participant_id: Option<&str>,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_file = input_file.ok_or("manifest execution requires --input-file")?;
    let store = GenotypeStore::from_file_with_options(Path::new(input_file), loader)
        .map_err(|err| err.to_string())?;
    let mut rows = Vec::new();

    for member in &panel.members {
        let Some(path) = &member.path else {
            return Err("remote panel members are not executable yet".to_owned());
        };
        let resolved = resolve_manifest_path(runtime_root, &panel.path, path)?;
        if member.kind == "variant" {
            let manifest = load_variant_manifest(&resolved)?;
            if !matches_filters(&manifest, &resolved, filters) {
                continue;
            }
            let observation = store
                .lookup_variant(&manifest.spec)
                .map_err(|err| err.to_string())?;
            rows.push(variant_row(
                runtime_root,
                &resolved,
                &manifest.name,
                &manifest.tags,
                &observation,
                participant_id,
            ));
        } else if member.kind == "assay" {
            let assay = load_assay_manifest(&resolved)?;
            rows.extend(run_assay_manifest_with_store(
                runtime_root,
                &assay,
                &store,
                participant_id,
                filters,
            )?);
        } else {
            return Err(format!(
                "panel member kind '{}' is not executable",
                member.kind
            ));
        }
    }

    Ok(rows)
}

fn run_assay_manifest(
    runtime_root: &Path,
    assay: &AssayManifest,
    input_file: Option<&str>,
    participant_id: Option<&str>,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_file = input_file.ok_or("manifest execution requires --input-file")?;
    let store = GenotypeStore::from_file_with_options(Path::new(input_file), loader)
        .map_err(|err| err.to_string())?;
    run_assay_manifest_with_store(runtime_root, assay, &store, participant_id, filters)
}

fn run_assay_manifest_with_store(
    runtime_root: &Path,
    assay: &AssayManifest,
    store: &GenotypeStore,
    participant_id: Option<&str>,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let mut rows = Vec::new();

    for member in &assay.members {
        if member.kind != "variant" {
            return Err(format!(
                "assay member kind '{}' is not executable",
                member.kind
            ));
        }
        let Some(path) = &member.path else {
            return Err("remote assay members are not executable yet".to_owned());
        };
        let resolved = resolve_manifest_path(runtime_root, &assay.path, path)?;
        let manifest = load_variant_manifest(&resolved)?;
        if !matches_filters(&manifest, &resolved, filters) {
            continue;
        }
        let observation = store
            .lookup_variant(&manifest.spec)
            .map_err(|err| err.to_string())?;
        rows.push(variant_row(
            runtime_root,
            &resolved,
            &manifest.name,
            &manifest.tags,
            &observation,
            participant_id,
        ));
    }

    Ok(rows)
}

fn variant_row(
    runtime_root: &Path,
    path: &Path,
    name: &str,
    tags: &[String],
    observation: &bioscript_core::VariantObservation,
    participant_id: Option<&str>,
) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("kind".to_owned(), "variant".to_owned());
    row.insert("name".to_owned(), name.to_owned());
    row.insert(
        "path".to_owned(),
        path.strip_prefix(runtime_root)
            .unwrap_or(path)
            .display()
            .to_string(),
    );
    row.insert("tags".to_owned(), tags.join(","));
    row.insert("backend".to_owned(), observation.backend.clone());
    row.insert(
        "participant_id".to_owned(),
        participant_id.unwrap_or_default().to_owned(),
    );
    row.insert(
        "matched_rsid".to_owned(),
        observation.matched_rsid.clone().unwrap_or_default(),
    );
    row.insert(
        "assembly".to_owned(),
        observation
            .assembly
            .map(|value| match value {
                bioscript_core::Assembly::Grch37 => "grch37".to_owned(),
                bioscript_core::Assembly::Grch38 => "grch38".to_owned(),
            })
            .unwrap_or_default(),
    );
    row.insert(
        "genotype".to_owned(),
        observation.genotype.clone().unwrap_or_default(),
    );
    row.insert(
        "ref_count".to_owned(),
        observation
            .ref_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "alt_count".to_owned(),
        observation
            .alt_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "depth".to_owned(),
        observation
            .depth
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert("evidence".to_owned(), observation.evidence.join(" | "));
    row
}

fn write_manifest_outputs(
    rows: &[BTreeMap<String, String>],
    output_file: Option<&Path>,
    trace_report: Option<&Path>,
) -> Result<(), String> {
    let text = render_rows_as_tsv(rows);
    if let Some(output_file) = output_file {
        if let Some(parent) = output_file.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create output dir {}: {err}", parent.display())
            })?;
        }
        fs::write(output_file, &text)
            .map_err(|err| format!("failed to write output {}: {err}", output_file.display()))?;
    } else {
        print!("{text}");
    }

    if let Some(trace_report) = trace_report {
        if let Some(parent) = trace_report.parent() {
            fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create trace dir {}: {err}", parent.display()))?;
        }
        let mut trace = String::from("step\tline\tcode\n");
        for (idx, row) in rows.iter().enumerate() {
            let _ = writeln!(
                trace,
                "{}\t{}\t{}",
                idx + 1,
                idx + 1,
                row.get("path").cloned().unwrap_or_default()
            );
        }
        fs::write(trace_report, trace)
            .map_err(|err| format!("failed to write trace {}: {err}", trace_report.display()))?;
    }

    Ok(())
}

fn resolve_cli_path(root: &Path, value: &str) -> String {
    resolve_cli_path_buf(root, Path::new(value))
        .display()
        .to_string()
}

fn resolve_cli_path_buf(root: &Path, value: &Path) -> PathBuf {
    if value.is_absolute() {
        value.to_path_buf()
    } else {
        root.join(value)
    }
}

fn render_rows_as_tsv(rows: &[BTreeMap<String, String>]) -> String {
    let headers = [
        "kind",
        "name",
        "path",
        "tags",
        "participant_id",
        "backend",
        "matched_rsid",
        "assembly",
        "genotype",
        "ref_count",
        "alt_count",
        "depth",
        "evidence",
    ];
    let mut out = headers.join("\t");
    out.push('\n');
    for row in rows {
        let line = headers
            .iter()
            .map(|header| {
                row.get(*header)
                    .cloned()
                    .unwrap_or_default()
                    .replace('\t', " ")
            })
            .collect::<Vec<_>>()
            .join("\t");
        out.push_str(&line);
        out.push('\n');
    }
    out
}

fn matches_filters(manifest: &VariantManifest, path: &Path, filters: &[String]) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("kind", value)) => value == "variant",
        Some(("name", value)) => manifest.name.contains(value),
        Some(("path", value)) => path.display().to_string().contains(value),
        Some(("tag", value)) => manifest.tags.iter().any(|tag| tag == value),
        Some(_) | None => false,
    })
}

fn resolve_manifest_path(
    runtime_root: &Path,
    manifest_path: &Path,
    relative: &str,
) -> Result<PathBuf, String> {
    let base_dir = manifest_path
        .parent()
        .ok_or_else(|| format!("manifest has no parent: {}", manifest_path.display()))?;
    let joined = base_dir.join(relative);
    let canonical_root = runtime_root
        .canonicalize()
        .map_err(|err| format!("failed to resolve root {}: {err}", runtime_root.display()))?;
    let canonical_base = base_dir.canonicalize().map_err(|err| {
        format!(
            "failed to resolve manifest dir {}: {err}",
            base_dir.display()
        )
    })?;
    let canonical_joined = joined
        .canonicalize()
        .map_err(|err| format!("failed to resolve {}: {err}", joined.display()))?;
    let boundary = if canonical_base.starts_with(&canonical_root) {
        &canonical_root
    } else {
        &canonical_base
    };
    if !canonical_joined.starts_with(boundary) {
        return Err(format!(
            "manifest member path escapes bioscript root: {}",
            canonical_joined.display()
        ));
    }
    Ok(canonical_joined)
}

fn manifest_schema(path: &Path) -> Result<String, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))?;
    value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("schema".to_owned())))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
        .ok_or_else(|| format!("{} is missing schema", path.display()))
}

fn normalize_loader_paths(root: &Path, loader: &mut GenotypeLoadOptions) {
    if let Some(path) = loader.input_index.take() {
        loader.input_index = Some(resolve_cli_path_buf(root, &path));
    }
    if let Some(path) = loader.reference_file.take() {
        loader.reference_file = Some(resolve_cli_path_buf(root, &path));
    }
    if let Some(path) = loader.reference_index.take() {
        loader.reference_index = Some(resolve_cli_path_buf(root, &path));
    }
}
