use std::{
    env,
    ffi::{CStr, CString},
    fmt::Write as _,
    fs,
    os::raw::c_char,
    path::PathBuf,
    time::{Duration, Instant},
};

use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, PrepareRequest, prepare_indexes,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig, StageTiming};
use monty::{MontyObject, ResourceLimits};
use serde::{Deserialize, Serialize};

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunFileRequest {
    pub script_path: String,
    pub root: Option<String>,
    pub input_file: Option<String>,
    pub output_file: Option<String>,
    pub participant_id: Option<String>,
    pub trace_report_path: Option<String>,
    pub timing_report_path: Option<String>,
    pub input_format: Option<String>,
    pub input_index: Option<String>,
    pub reference_file: Option<String>,
    pub reference_index: Option<String>,
    pub auto_index: Option<bool>,
    pub cache_dir: Option<String>,
    pub max_duration_ms: Option<u64>,
    pub max_memory_bytes: Option<usize>,
    pub max_allocations: Option<usize>,
    pub max_recursion_depth: Option<usize>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunFileResult {
    pub ok: bool,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct FfiResult<T> {
    ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    value: Option<T>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<String>,
}

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
            root: runtime_root.clone(),
            cwd: cwd.clone(),
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
    }

    let runtime = BioscriptRuntime::with_config(runtime_root, RuntimeConfig { limits, loader })
        .map_err(|err| err.to_string())?;

    let mut inputs = Vec::new();
    if let Some(input_file) = request.input_file {
        inputs.push(("input_file", MontyObject::String(input_file)));
    }
    if let Some(output_file) = request.output_file {
        inputs.push(("output_file", MontyObject::String(output_file)));
    }
    if let Some(participant_id) = request.participant_id {
        inputs.push(("participant_id", MontyObject::String(participant_id)));
    }

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

fn build_loader(request: &RunFileRequest) -> Result<GenotypeLoadOptions, String> {
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
    Ok(loader)
}

fn build_limits(request: &RunFileRequest) -> ResourceLimits {
    let mut limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(100))
        .max_memory(8 * 1024 * 1024)
        .max_allocations(200_000)
        .gc_interval(1000)
        .max_recursion_depth(Some(200));

    if let Some(value) = request.max_duration_ms {
        limits = limits.max_duration(Duration::from_millis(value));
    }
    if let Some(value) = request.max_memory_bytes {
        limits = limits.max_memory(value);
    }
    if let Some(value) = request.max_allocations {
        limits = limits.max_allocations(value);
    }
    if let Some(value) = request.max_recursion_depth {
        limits = limits.max_recursion_depth(Some(value));
    }

    limits
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

#[unsafe(no_mangle)]
/// Executes a bioscript run request encoded as a UTF-8 JSON C string.
///
/// # Safety
///
/// `request_json` must either be null or point to a valid, NUL-terminated C
/// string that remains alive for the duration of this call.
pub unsafe extern "C" fn bioscript_run_file_json(request_json: *const c_char) -> *mut c_char {
    let response = unsafe {
        if request_json.is_null() {
            FfiResult::<RunFileResult> {
                ok: false,
                value: None,
                error: Some("request_json was null".to_owned()),
            }
        } else {
            match CStr::from_ptr(request_json).to_str() {
                Ok(value) => match serde_json::from_str::<RunFileRequest>(value) {
                    Ok(request) => match run_file_request(request) {
                        Ok(result) => FfiResult {
                            ok: true,
                            value: Some(result),
                            error: None,
                        },
                        Err(error) => FfiResult::<RunFileResult> {
                            ok: false,
                            value: None,
                            error: Some(error),
                        },
                    },
                    Err(error) => FfiResult::<RunFileResult> {
                        ok: false,
                        value: None,
                        error: Some(format!("invalid request JSON: {error}")),
                    },
                },
                Err(error) => FfiResult::<RunFileResult> {
                    ok: false,
                    value: None,
                    error: Some(format!("request_json was not valid UTF-8: {error}")),
                },
            }
        }
    };

    match serde_json::to_string(&response) {
        Ok(json) => match CString::new(json) {
            Ok(value) => value.into_raw(),
            Err(_) => std::ptr::null_mut(),
        },
        Err(_) => std::ptr::null_mut(),
    }
}

#[unsafe(no_mangle)]
/// Frees a string previously returned by [`bioscript_run_file_json`].
///
/// # Safety
///
/// `ptr` must be null or a pointer returned by [`CString::into_raw`] from this
/// library, and it must not be freed more than once.
pub unsafe extern "C" fn bioscript_free_string(ptr: *mut c_char) {
    if !ptr.is_null() {
        unsafe {
            let _ = CString::from_raw(ptr);
        }
    }
}

#[cfg(target_os = "android")]
pub mod android {
    use crate::{RunFileRequest, run_file_request};
    use jni::JNIEnv;
    use jni::objects::{JClass, JString};

    #[unsafe(no_mangle)]
    pub unsafe extern "system" fn Java_expo_modules_bioscript_ExpoBioscriptNativeBridge_runFileNative<
        'local,
    >(
        mut env: JNIEnv<'local>,
        _class: JClass<'local>,
        request_json: JString<'local>,
    ) -> JString<'local> {
        let request_string: String = match env.get_string(&request_json) {
            Ok(value) => value.into(),
            Err(error) => {
                return env
                    .new_string(
                        serde_json::json!({
                            "ok": false,
                            "error": format!("failed to read request json from JVM: {error}")
                        })
                        .to_string(),
                    )
                    .expect("jni new_string should succeed");
            }
        };

        let response = match serde_json::from_str::<RunFileRequest>(&request_string) {
            Ok(request) => match run_file_request(request) {
                Ok(value) => serde_json::json!({ "ok": true, "value": value }).to_string(),
                Err(error) => serde_json::json!({ "ok": false, "error": error }).to_string(),
            },
            Err(error) => {
                serde_json::json!({ "ok": false, "error": format!("invalid request JSON: {error}") })
                    .to_string()
            }
        };

        env.new_string(response)
            .expect("jni new_string should succeed")
    }
}
