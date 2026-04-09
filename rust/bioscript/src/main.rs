use std::{env, path::PathBuf, process::ExitCode, time::Duration};

use bioscript::{BioscriptRuntime, GenotypeLoadOptions, GenotypeSourceFormat, RuntimeConfig, validate_variants_path};
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
    let mut args = env::args().skip(1);
    if let Some(first) = args.next()
        && first == "validate-variants"
    {
        return run_validate_variants(args.collect());
    }

    let mut args = env::args().skip(1);
    let mut script_path: Option<PathBuf> = None;
    let mut root: Option<PathBuf> = None;
    let mut input_file: Option<String> = None;
    let mut output_file: Option<String> = None;
    let mut participant_id: Option<String> = None;
    let mut trace_report: Option<PathBuf> = None;
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
            "usage: bioscript <script.py> [--root <dir>] [--input-file <path>] [--output-file <path>] [--participant-id <id>] [--trace-report <path>] [--input-format auto|text|zip|vcf|cram] [--input-index <path>] [--reference-file <path>] [--reference-index <path>] [--max-duration-ms N] [--max-memory-bytes N] [--max-allocations N] [--max-recursion-depth N]"
                .to_owned(),
        );
    };

    let runtime_root = match root {
        Some(dir) => dir,
        None => env::current_dir().map_err(|err| format!("failed to get current directory: {err}"))?,
    };

    let runtime = BioscriptRuntime::with_config(
        runtime_root,
        RuntimeConfig {
            limits,
            loader,
        },
    )
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
            std::fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create report dir {}: {err}", parent.display()))?;
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
