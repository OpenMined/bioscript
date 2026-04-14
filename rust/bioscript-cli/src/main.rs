use std::{
    env,
    fmt::Write as _,
    fs,
    path::PathBuf,
    process::ExitCode,
    time::{Duration, Instant},
};

use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, PrepareRequest, prepare_indexes, shell_flags,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig, StageTiming};
use bioscript_schema::validate_variants_path;
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
        if first == "validate-variants" {
            return run_validate_variants(args.collect());
        }
        if first == "prepare" {
            return run_prepare(args.collect());
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
            "usage: bioscript <script.py> [--root <dir>] [--input-file <path>] [--output-file <path>] [--participant-id <id>] [--trace-report <path>] [--input-format auto|text|zip|vcf|cram] [--input-index <path>] [--reference-file <path>] [--reference-index <path>] [--auto-index] [--cache-dir <path>] [--max-duration-ms N] [--max-memory-bytes N] [--max-allocations N] [--max-recursion-depth N]\n       bioscript prepare [--root <dir>] [--input-file <path>] [--reference-file <path>] [--input-format auto|text|zip|vcf|cram] [--cache-dir <path>]"
                .to_owned(),
        );
    };

    let runtime_root = match root {
        Some(dir) => dir,
        None => {
            env::current_dir().map_err(|err| format!("failed to get current directory: {err}"))?
        }
    };

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
