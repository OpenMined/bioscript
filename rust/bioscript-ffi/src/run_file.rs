use std::{
    env,
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
    time::Instant,
};

use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, PrepareRequest, prepare_indexes,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig, StageTiming};
use monty::MontyObject;

use crate::{
    limits::build_limits,
    types::{RunFileRequest, RunFileResult},
};

/// Runs a bioscript file request described by a JSON-compatible Rust struct.
///
/// # Errors
///
/// Returns an error string when request parsing, optional index preparation,
/// runtime construction, script execution, or report writing fails.
pub fn run_file_request(request: RunFileRequest) -> Result<RunFileResult, String> {
    let script_path = PathBuf::from(&request.script_path);
    let runtime_root = runtime_root(&request)?;
    let mut loader = build_loader(&request)?;
    let limits = build_limits(&request);

    let mut ffi_timings: Vec<StageTiming> = Vec::new();
    if request.auto_index.unwrap_or(false) {
        run_auto_index(&request, &runtime_root, &mut loader, &mut ffi_timings)?;
    }

    let runtime = BioscriptRuntime::with_config(runtime_root, RuntimeConfig { limits, loader })
        .map_err(|err| err.to_string())?;
    let inputs = runtime_inputs(&request);

    runtime
        .run_file(
            &script_path,
            request
                .trace_report_path
                .as_deref()
                .map(std::path::Path::new),
            inputs,
        )
        .map_err(|err| err.to_string())?;

    if let Some(timing_path) = request.timing_report_path {
        let mut all_timings = ffi_timings;
        all_timings.extend(runtime.timing_snapshot());
        write_timing_report(&PathBuf::from(timing_path), &all_timings)?;
    }

    Ok(RunFileResult { ok: true })
}

fn runtime_root(request: &RunFileRequest) -> Result<PathBuf, String> {
    match request.root.as_deref() {
        Some(dir) => Ok(PathBuf::from(dir)),
        None => env::current_dir().map_err(|err| format!("failed to get current directory: {err}")),
    }
}

fn run_auto_index(
    request: &RunFileRequest,
    runtime_root: &Path,
    loader: &mut GenotypeLoadOptions,
    ffi_timings: &mut Vec<StageTiming>,
) -> Result<(), String> {
    let auto_index_started = Instant::now();
    let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
    let effective_cache = request.cache_dir.as_ref().map_or_else(
        || runtime_root.join(".bioscript-cache"),
        |dir| {
            let path = PathBuf::from(dir);
            if path.is_absolute() {
                path
            } else {
                runtime_root.join(path)
            }
        },
    );
    let prepare_request = PrepareRequest {
        root: runtime_root.to_path_buf(),
        cwd,
        cache_dir: effective_cache,
        input_file: request.input_file.clone(),
        input_format: loader.format,
        reference_file: loader
            .reference_file
            .as_ref()
            .map(|p| p.to_string_lossy().to_string()),
    };
    let prepared = prepare_indexes(&prepare_request)?;
    if let Some(idx) = prepared.input_index
        && loader.input_index.is_none()
    {
        loader.input_index = Some(idx);
    }
    if let Some(ref_file) = prepared.reference_file {
        loader.reference_file = Some(ref_file);
    }
    if let Some(ref_idx) = prepared.reference_index
        && loader.reference_index.is_none()
    {
        loader.reference_index = Some(ref_idx);
    }
    ffi_timings.push(StageTiming {
        stage: "auto_index".to_owned(),
        duration_ms: auto_index_started.elapsed().as_millis(),
        detail: "prepare_indexes".to_owned(),
    });
    Ok(())
}

fn runtime_inputs(request: &RunFileRequest) -> Vec<(&'static str, MontyObject)> {
    let mut inputs = Vec::new();
    if let Some(input_file) = request.input_file.clone() {
        inputs.push(("input_file", MontyObject::String(input_file)));
    }
    if let Some(output_file) = request.output_file.clone() {
        inputs.push(("output_file", MontyObject::String(output_file)));
    }
    if let Some(participant_id) = request.participant_id.clone() {
        inputs.push(("participant_id", MontyObject::String(participant_id)));
    }
    inputs
}

pub(crate) fn build_loader(request: &RunFileRequest) -> Result<GenotypeLoadOptions, String> {
    let mut loader = GenotypeLoadOptions::default();
    if let Some(value) = request.input_format.as_deref() {
        if value.eq_ignore_ascii_case("auto") {
            loader.format = None;
        } else {
            let parsed = value
                .parse::<GenotypeSourceFormat>()
                .map_err(|err| format!("invalid inputFormat value {value}: {err}"))?;
            loader.format = Some(parsed);
        }
    }
    loader.input_index = request.input_index.clone().map(PathBuf::from);
    loader.reference_file = request.reference_file.clone().map(PathBuf::from);
    loader.reference_index = request.reference_index.clone().map(PathBuf::from);
    loader.allow_reference_md5_mismatch = request.allow_md5_mismatch.unwrap_or(false);
    Ok(loader)
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
            "{}\t{}\t{}\n",
            timing.stage,
            timing.duration_ms,
            timing.detail.replace('\t', " ")
        );
    }
    fs::write(path, output)
        .map_err(|err| format!("failed to write timing report {}: {err}", path.display()))
}
