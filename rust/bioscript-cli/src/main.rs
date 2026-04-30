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
    PanelManifest, VariantManifest, load_panel_manifest, load_variant_manifest,
    validate_panels_path, validate_variants_path,
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
    run_cli_args(env::args().skip(1).collect())
}

#[allow(clippy::too_many_lines)]
fn run_cli_args(raw_args: Vec<String>) -> Result<(), String> {
    let mut args = raw_args.clone().into_iter();
    if let Some(first) = args.next() {
        if first == "validate-variants" {
            return run_validate_variants(args.collect());
        }
        if first == "validate-panels" {
            return run_validate_panels(args.collect());
        }
        if first == "prepare" {
            return run_prepare(args.collect());
        }
        if first == "inspect" {
            return run_inspect(args.collect());
        }
    }

    let mut args = raw_args.into_iter();
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
        } else if arg == "--allow-md5-mismatch" {
            loader.allow_reference_md5_mismatch = true;
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
            "usage: bioscript <script.py|manifest.yaml> [--root <dir>] [--input-file <path>] [--output-file <path>] [--participant-id <id>] [--trace-report <path>] [--timing-report <path>] [--filter key=value] [--input-format auto|text|zip|vcf|cram] [--input-index <path>] [--reference-file <path>] [--reference-index <path>] [--allow-md5-mismatch] [--auto-index] [--cache-dir <path>] [--max-duration-ms N] [--max-memory-bytes N] [--max-allocations N] [--max-recursion-depth N]\n       bioscript validate-variants <path> [--report <file>]\n       bioscript validate-panels <path> [--report <file>]\n       bioscript prepare [--root <dir>] [--input-file <path>] [--reference-file <path>] [--input-format auto|text|zip|vcf|cram] [--cache-dir <path>]\n       bioscript inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>]"
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
        if member.kind != "variant" {
            return Err(format!(
                "panel member kind '{}' is not executable yet; panel execution is currently variant-only",
                member.kind
            ));
        }
        let Some(path) = &member.path else {
            return Err("remote panel members are not executable yet".to_owned());
        };
        let resolved = resolve_manifest_path(runtime_root, &panel.path, path)?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use bioscript_core::{Assembly, VariantObservation};
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock drift")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "bioscript-cli-unit-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn cli_private_helpers_render_rows_filters_paths_and_loader_paths() {
        let root = temp_dir("helpers-root");
        let manifest_path = root.join("panels/panel.yaml");
        let member_dir = root.join("panels/members");
        fs::create_dir_all(&member_dir).unwrap();
        let variant_path = member_dir.join("apol1.yaml");
        fs::write(&manifest_path, "schema: bioscript:panel:1.0\n").unwrap();
        fs::write(&variant_path, "schema: bioscript:variant:1.0\n").unwrap();

        let manifest = VariantManifest {
            name: "APOL1 G1".to_owned(),
            path: variant_path.clone(),
            tags: vec!["kidney".to_owned(), "apol1".to_owned()],
            spec: bioscript_core::VariantSpec::default(),
        };

        assert!(matches_filters(
            &manifest,
            &variant_path,
            &[
                "kind=variant".to_owned(),
                "name=APOL1".to_owned(),
                "tag=kidney".to_owned(),
                "path=apol1".to_owned(),
            ],
        ));
        assert!(!matches_filters(
            &manifest,
            &variant_path,
            &["kind=panel".to_owned()]
        ));
        assert!(!matches_filters(
            &manifest,
            &variant_path,
            &["bad".to_owned()]
        ));

        assert_eq!(
            resolve_manifest_path(&root, &manifest_path, "members/apol1.yaml").unwrap(),
            variant_path.canonicalize().unwrap()
        );
        let outside = temp_dir("helpers-outside").join("escape.yaml");
        fs::write(&outside, "schema: bioscript:variant:1.0\n").unwrap();
        let err =
            resolve_manifest_path(&root, &manifest_path, &outside.to_string_lossy()).unwrap_err();
        assert!(err.contains("escapes bioscript root"), "{err}");

        let observation = VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid: Some("rs1".to_owned()),
            assembly: Some(Assembly::Grch38),
            genotype: Some("AG".to_owned()),
            ref_count: Some(7),
            alt_count: Some(3),
            depth: Some(10),
            evidence: vec!["one\twith tab".to_owned(), "two".to_owned()],
            ..VariantObservation::default()
        };
        let row = variant_row(
            &root,
            &variant_path,
            "APOL1 G1",
            &["kidney".to_owned()],
            &observation,
            Some("p1"),
        );
        let tsv = render_rows_as_tsv(&[row]);
        assert!(tsv.contains("participant_id\tbackend"), "{tsv}");
        assert!(tsv.contains("p1\tvcf\trs1\tgrch38\tAG\t7\t3\t10"), "{tsv}");
        assert!(tsv.contains("one with tab | two"), "{tsv}");

        assert_eq!(
            resolve_cli_path(&root, "sample.txt"),
            root.join("sample.txt").display().to_string()
        );
        assert_eq!(
            resolve_cli_path_buf(&root, Path::new("/tmp/abs")),
            PathBuf::from("/tmp/abs")
        );

        let mut loader = GenotypeLoadOptions {
            input_index: Some(PathBuf::from("input.crai")),
            reference_file: Some(PathBuf::from("ref.fa")),
            reference_index: Some(PathBuf::from("ref.fa.fai")),
            ..GenotypeLoadOptions::default()
        };
        normalize_loader_paths(&root, &mut loader);
        assert_eq!(
            loader.input_index.as_deref(),
            Some(root.join("input.crai").as_path())
        );
        assert_eq!(
            loader.reference_file.as_deref(),
            Some(root.join("ref.fa").as_path())
        );
        assert_eq!(
            loader.reference_index.as_deref(),
            Some(root.join("ref.fa.fai").as_path())
        );
    }

    #[test]
    fn cli_private_helpers_cover_manifest_schema_and_timing_errors() {
        let dir = temp_dir("schema-timing");
        let valid = dir.join("valid.yaml");
        let missing_schema = dir.join("missing.yaml");
        let invalid_yaml = dir.join("invalid.yaml");
        fs::write(&valid, "schema: bioscript:variant:1.0\n").unwrap();
        fs::write(&missing_schema, "name: no schema\n").unwrap();
        fs::write(&invalid_yaml, "schema: [").unwrap();

        assert_eq!(manifest_schema(&valid).unwrap(), "bioscript:variant:1.0");
        assert!(
            manifest_schema(&missing_schema)
                .unwrap_err()
                .contains("missing schema")
        );
        assert!(
            manifest_schema(&invalid_yaml)
                .unwrap_err()
                .contains("failed to parse YAML")
        );
        assert!(
            manifest_schema(&dir.join("absent.yaml"))
                .unwrap_err()
                .contains("failed to read")
        );

        let timing_path = dir.join("nested/timing.tsv");
        write_timing_report(
            &timing_path,
            &[
                StageTiming {
                    stage: "one".to_owned(),
                    duration_ms: 2,
                    detail: "contains\ttab".to_owned(),
                },
                StageTiming {
                    stage: "two".to_owned(),
                    duration_ms: 3,
                    detail: "plain".to_owned(),
                },
            ],
        )
        .unwrap();
        let report = fs::read_to_string(&timing_path).unwrap();
        assert!(report.contains("stage\tduration_ms\tdetail"));
        assert!(report.contains("one\t2\tcontains tab"));
    }

    #[test]
    fn cli_arg_parser_reports_missing_and_invalid_values_without_spawning() {
        for (flag, expected) in [
            ("--input-format", "--input-format requires a value"),
            ("--max-duration-ms", "--max-duration-ms requires an integer"),
            (
                "--max-memory-bytes",
                "--max-memory-bytes requires an integer",
            ),
            ("--max-allocations", "--max-allocations requires an integer"),
            (
                "--max-recursion-depth",
                "--max-recursion-depth requires an integer",
            ),
        ] {
            let err = run_cli_args(vec![flag.to_owned()]).unwrap_err();
            assert!(err.contains(expected), "{flag}: {err}");
        }

        for (flag, value, expected) in [
            (
                "--input-format",
                "unknown",
                "invalid --input-format value unknown",
            ),
            (
                "--max-duration-ms",
                "nan",
                "invalid --max-duration-ms value nan",
            ),
            (
                "--max-memory-bytes",
                "nan",
                "invalid --max-memory-bytes value nan",
            ),
            (
                "--max-allocations",
                "nan",
                "invalid --max-allocations value nan",
            ),
            (
                "--max-recursion-depth",
                "nan",
                "invalid --max-recursion-depth value nan",
            ),
        ] {
            let err = run_cli_args(vec![flag.to_owned(), value.to_owned()]).unwrap_err();
            assert!(err.contains(expected), "{flag}: {err}");
        }
    }
}
