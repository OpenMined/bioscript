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

fn run_cli() -> Result<(), String> {
    let args: Vec<String> = env::args().skip(1).collect();
    if dispatch_subcommand(&args)? {
        return Ok(());
    }

    let mut options = parse_cli_options(args)?;
    let script_path = options.script_path.clone().ok_or_else(|| USAGE.to_owned())?;
    let runtime_root = options
        .root
        .clone()
        .map_or_else(env::current_dir, Ok)
        .map_err(|err| format!("failed to get current directory: {err}"))?;
    normalize_loader_paths(&runtime_root, &mut options.loader);
    let mut cli_timings = prepare_cli_indexes(&runtime_root, &mut options)?;

    if is_yaml_manifest(&script_path) {
        run_cli_manifest(&runtime_root, &script_path, &options, &mut cli_timings)?;
    } else {
        run_cli_script(&script_path, options, cli_timings)?;
    }
    Ok(())
}

const USAGE: &str = "usage: bioscript <script.py|manifest.yaml> [--root <dir>] [--input-file <path>] [--output-file <path>] [--participant-id <id>] [--trace-report <path>] [--timing-report <path>] [--filter key=value] [--input-format auto|text|zip|vcf|cram] [--input-index <path>] [--reference-file <path>] [--reference-index <path>] [--auto-index] [--cache-dir <path>] [--max-duration-ms N] [--max-memory-bytes N] [--max-allocations N] [--max-recursion-depth N]\n       bioscript report <manifest.yaml> --input-file <path> [--input-file <path>...] --output-dir <dir> [--html] [--root <dir>] [--input-format auto|text|zip|vcf|cram]\n       bioscript validate-variants <path> [--report <file>]\n       bioscript validate-panels <path> [--report <file>]\n       bioscript validate-assays <path> [--report <file>]\n       bioscript prepare [--root <dir>] [--input-file <path>] [--reference-file <path>] [--input-format auto|text|zip|vcf|cram] [--cache-dir <path>]\n       bioscript inspect <path> [--input-index <path>] [--reference-file <path>] [--reference-index <path>]";

struct CliOptions {
    script_path: Option<PathBuf>,
    root: Option<PathBuf>,
    input_file: Option<String>,
    output_file: Option<String>,
    participant_id: Option<String>,
    trace_report: Option<PathBuf>,
    timing_report: Option<PathBuf>,
    filters: Vec<String>,
    auto_index: bool,
    cache_dir: Option<PathBuf>,
    loader: GenotypeLoadOptions,
    limits: ResourceLimits,
}

fn dispatch_subcommand(args: &[String]) -> Result<bool, String> {
    let Some((first, rest)) = args.split_first() else {
        return Ok(false);
    };
    let rest = rest.to_vec();
    match first.as_str() {
        "report" => run_app_report(rest).map(|()| true),
        "validate-variants" => run_validate_variants(rest).map(|()| true),
        "validate-panels" => run_validate_panels(rest).map(|()| true),
        "validate-assays" => run_validate_assays(rest).map(|()| true),
        "prepare" => run_prepare(rest).map(|()| true),
        "inspect" => run_inspect(rest).map(|()| true),
        _ => Ok(false),
    }
}

fn parse_cli_options(args: Vec<String>) -> Result<CliOptions, String> {
    let mut args = args.into_iter();
    let mut options = default_cli_options();
    while let Some(arg) = args.next() {
        parse_cli_arg(arg, &mut args, &mut options)?;
    }
    Ok(options)
}

fn default_cli_options() -> CliOptions {
    CliOptions {
        script_path: None,
        root: None,
        input_file: None,
        output_file: None,
        participant_id: None,
        trace_report: None,
        timing_report: None,
        filters: Vec::new(),
        auto_index: false,
        cache_dir: None,
        loader: GenotypeLoadOptions::default(),
        limits: ResourceLimits::new()
        .max_duration(Duration::from_millis(100))
        .max_memory(8 * 1024 * 1024)
        .max_allocations(200_000)
        .gc_interval(1000)
            .max_recursion_depth(Some(200)),
    }
}

fn parse_cli_arg(
    arg: String,
    args: &mut impl Iterator<Item = String>,
    options: &mut CliOptions,
) -> Result<(), String> {
    if parse_cli_path_arg(&arg, args, options)? {
        return Ok(());
    }
    if parse_cli_loader_arg(&arg, args, options)? {
        return Ok(());
    }
    if parse_cli_limit_arg(&arg, args, options)? {
        return Ok(());
    }
    if arg == "--auto-index" {
        options.auto_index = true;
    } else if options.script_path.is_none() {
        options.script_path = Some(PathBuf::from(arg));
    } else {
        return Err(format!("unexpected argument: {arg}"));
    }
    Ok(())
}

fn parse_cli_path_arg(
    arg: &str,
    args: &mut impl Iterator<Item = String>,
    options: &mut CliOptions,
) -> Result<bool, String> {
    if arg == "--root" {
            let Some(value) = args.next() else {
                return Err("--root requires a directory".to_owned());
            };
        options.root = Some(PathBuf::from(value));
    } else if arg == "--input-file" {
            let Some(value) = args.next() else {
                return Err("--input-file requires a path".to_owned());
            };
        options.input_file = Some(value);
    } else if arg == "--output-file" {
            let Some(value) = args.next() else {
                return Err("--output-file requires a path".to_owned());
            };
        options.output_file = Some(value);
    } else if arg == "--participant-id" {
            let Some(value) = args.next() else {
                return Err("--participant-id requires a value".to_owned());
            };
        options.participant_id = Some(value);
    } else if arg == "--trace-report" {
            let Some(value) = args.next() else {
                return Err("--trace-report requires a path".to_owned());
            };
        options.trace_report = Some(PathBuf::from(value));
    } else if arg == "--timing-report" {
            let Some(value) = args.next() else {
                return Err("--timing-report requires a path".to_owned());
            };
        options.timing_report = Some(PathBuf::from(value));
    } else if arg == "--filter" {
            let Some(value) = args.next() else {
                return Err("--filter requires key=value".to_owned());
            };
        options.filters.push(value);
    } else if arg == "--cache-dir" {
        let Some(value) = args.next() else {
            return Err("--cache-dir requires a path".to_owned());
        };
        options.cache_dir = Some(PathBuf::from(value));
    } else {
        return Ok(false);
    }
    Ok(true)
}

fn parse_cli_loader_arg(
    arg: &str,
    args: &mut impl Iterator<Item = String>,
    options: &mut CliOptions,
) -> Result<bool, String> {
    if arg == "--input-format" {
            let Some(value) = args.next() else {
                return Err("--input-format requires a value".to_owned());
            };
            if value.eq_ignore_ascii_case("auto") {
            options.loader.format = None;
            } else {
                let parsed = value
                    .parse::<GenotypeSourceFormat>()
                    .map_err(|err| format!("invalid --input-format value {value}: {err}"))?;
            options.loader.format = Some(parsed);
            }
    } else if arg == "--input-index" {
            let Some(value) = args.next() else {
                return Err("--input-index requires a path".to_owned());
            };
        options.loader.input_index = Some(PathBuf::from(value));
    } else if arg == "--reference-file" {
            let Some(value) = args.next() else {
                return Err("--reference-file requires a path".to_owned());
            };
        options.loader.reference_file = Some(PathBuf::from(value));
    } else if arg == "--reference-index" {
            let Some(value) = args.next() else {
                return Err("--reference-index requires a path".to_owned());
            };
        options.loader.reference_index = Some(PathBuf::from(value));
    } else {
        return Ok(false);
    }
    Ok(true)
}

fn parse_cli_limit_arg(
    arg: &str,
    args: &mut impl Iterator<Item = String>,
    options: &mut CliOptions,
) -> Result<bool, String> {
    if arg == "--max-duration-ms" {
            let Some(value) = args.next() else {
                return Err("--max-duration-ms requires an integer".to_owned());
            };
            let parsed = value
                .parse::<u64>()
                .map_err(|err| format!("invalid --max-duration-ms value {value}: {err}"))?;
        options.limits = options.limits.clone().max_duration(Duration::from_millis(parsed));
    } else if arg == "--max-memory-bytes" {
            let Some(value) = args.next() else {
                return Err("--max-memory-bytes requires an integer".to_owned());
            };
            let parsed = value
                .parse::<usize>()
                .map_err(|err| format!("invalid --max-memory-bytes value {value}: {err}"))?;
        options.limits = options.limits.clone().max_memory(parsed);
    } else if arg == "--max-allocations" {
            let Some(value) = args.next() else {
                return Err("--max-allocations requires an integer".to_owned());
            };
            let parsed = value
                .parse::<usize>()
                .map_err(|err| format!("invalid --max-allocations value {value}: {err}"))?;
        options.limits = options.limits.clone().max_allocations(parsed);
    } else if arg == "--max-recursion-depth" {
            let Some(value) = args.next() else {
                return Err("--max-recursion-depth requires an integer".to_owned());
            };
            let parsed = value
                .parse::<usize>()
                .map_err(|err| format!("invalid --max-recursion-depth value {value}: {err}"))?;
        options.limits = options.limits.clone().max_recursion_depth(Some(parsed));
    } else {
        return Ok(false);
    }
    Ok(true)
}

fn prepare_cli_indexes(
    runtime_root: &Path,
    options: &mut CliOptions,
) -> Result<Vec<StageTiming>, String> {
    let mut cli_timings: Vec<StageTiming> = Vec::new();
    if options.auto_index {
        let auto_index_started = Instant::now();
        let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
        let effective_cache = options
            .cache_dir
            .clone()
            .unwrap_or_else(|| cwd.join(".bioscript-cache"));
        let request = PrepareRequest {
            root: runtime_root.to_path_buf(),
            cwd,
            cache_dir: effective_cache,
            input_file: options.input_file.clone(),
            input_format: options.loader.format,
            reference_file: options
                .loader
                .reference_file
                .as_ref()
                .map(|p| p.to_string_lossy().to_string()),
        };
        let prepared = prepare_indexes(&request)?;
        if let Some(idx) = prepared.input_index
            && options.loader.input_index.is_none()
        {
            eprintln!("bioscript: auto-indexed input -> {}", idx.display());
            options.loader.input_index = Some(idx);
        }
        if let Some(ref_file) = prepared.reference_file {
            options.loader.reference_file = Some(ref_file);
        }
        if let Some(ref_idx) = prepared.reference_index
            && options.loader.reference_index.is_none()
        {
            eprintln!("bioscript: auto-indexed reference -> {}", ref_idx.display());
            options.loader.reference_index = Some(ref_idx);
        }
        cli_timings.push(StageTiming {
            stage: "auto_index".to_owned(),
            duration_ms: auto_index_started.elapsed().as_millis(),
            detail: "prepare_indexes".to_owned(),
        });
    }
    Ok(cli_timings)
}

fn run_cli_manifest(
    runtime_root: &Path,
    script_path: &Path,
    options: &CliOptions,
    cli_timings: &mut Vec<StageTiming>,
) -> Result<(), String> {
    let manifest_started = Instant::now();
    let manifest_options = ManifestRunOptions {
        input_file: options.input_file.as_deref(),
        output_file: options.output_file.as_deref(),
        participant_id: options.participant_id.as_deref(),
        trace_report: options.trace_report.as_deref(),
        loader: &options.loader,
        filters: &options.filters,
    };
    run_manifest(runtime_root, script_path, &manifest_options)?;
    cli_timings.push(StageTiming {
        stage: "manifest_run".to_owned(),
        duration_ms: manifest_started.elapsed().as_millis(),
        detail: script_path.display().to_string(),
    });
    if let Some(timing_path) = &options.timing_report {
        write_timing_report(timing_path, cli_timings)?;
    }
    Ok(())
}

fn run_cli_script(
    script_path: &Path,
    options: CliOptions,
    cli_timings: Vec<StageTiming>,
) -> Result<(), String> {
    let runtime_root = options
        .root
        .map_or_else(env::current_dir, Ok)
        .map_err(|err| format!("failed to get current directory: {err}"))?;
    let runtime = BioscriptRuntime::with_config(
        runtime_root,
        RuntimeConfig {
            limits: options.limits,
            loader: options.loader,
        },
    )
    .map_err(|err| err.to_string())?;
    let mut inputs = Vec::new();
    if let Some(input_file) = options.input_file {
        inputs.push(("input_file", monty::MontyObject::String(input_file)));
    }
    if let Some(output_file) = options.output_file {
        inputs.push(("output_file", monty::MontyObject::String(output_file)));
    }
    if let Some(participant_id) = options.participant_id {
        inputs.push(("participant_id", monty::MontyObject::String(participant_id)));
    }

    runtime
        .run_file(script_path, options.trace_report.as_deref(), inputs)
        .map_err(|err| err.to_string())?;
    if let Some(timing_path) = options.timing_report {
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

