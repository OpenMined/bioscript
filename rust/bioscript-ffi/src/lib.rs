use std::{
    collections::HashMap,
    env,
    ffi::{CStr, CString},
    fs,
    os::raw::c_char,
    path::{Component, Path, PathBuf},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
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
    pub script_contents: Option<String>,
    pub root: Option<String>,
    pub input_file: Option<String>,
    pub input_contents: Option<String>,
    pub output_file: Option<String>,
    pub file_contents: Option<HashMap<String, String>>,
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

fn unique_temp_dir(label: &str) -> Result<PathBuf, String> {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(|err| format!("failed to get system time: {err}"))?
        .as_nanos();
    let dir = env::temp_dir().join(format!("bioscript-ffi-{label}-{}-{nanos}", std::process::id()));
    fs::create_dir_all(&dir)
        .map_err(|err| format!("failed to create temp dir {}: {err}", dir.display()))?;
    Ok(dir)
}

fn resolve_staged_path(root: &Path, relative: &str) -> Result<PathBuf, String> {
    let candidate = Path::new(relative);
    if candidate.is_absolute() {
        return Err(format!("staged file path must be relative: {relative}"));
    }

    let mut result = root.to_path_buf();
    for component in candidate.components() {
        match component {
            Component::Normal(segment) => result.push(segment),
            Component::CurDir => {}
            Component::ParentDir => {
                return Err(format!("staged file path must not traverse upward: {relative}"));
            }
            Component::Prefix(_) | Component::RootDir => {
                return Err(format!("staged file path must be relative: {relative}"));
            }
        }
    }

    Ok(result)
}

fn stage_runtime_files(root: &Path, request: &RunFileRequest) -> Result<(), String> {
    if let Some(files) = &request.file_contents {
        for (relative_path, contents) in files {
            let destination = resolve_staged_path(root, relative_path)?;
            if let Some(parent) = destination.parent() {
                fs::create_dir_all(parent).map_err(|err| {
                    format!("failed to create staged directory {}: {err}", parent.display())
                })?;
            }
            fs::write(&destination, contents)
                .map_err(|err| format!("failed to stage file {}: {err}", destination.display()))?;
        }
    }

    if let Some(input_contents) = &request.input_contents
        && let Some(input_file) = &request.input_file
    {
        let destination = resolve_staged_path(root, input_file)?;
        if let Some(parent) = destination.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create input directory {}: {err}", parent.display())
            })?;
        }
        fs::write(&destination, input_contents)
            .map_err(|err| format!("failed to stage input file {}: {err}", destination.display()))?;
    }

    if let Some(script_contents) = &request.script_contents {
        let destination = resolve_staged_path(root, &request.script_path)?;
        if let Some(parent) = destination.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create script directory {}: {err}", parent.display())
            })?;
        }
        fs::write(&destination, script_contents)
            .map_err(|err| format!("failed to stage script {}: {err}", destination.display()))?;
    }

    Ok(())
}

pub fn run_file_request(request: RunFileRequest) -> Result<RunFileResult, String> {
    let should_stage_runtime_files = request.script_contents.is_some()
        || request.input_contents.is_some()
        || request
            .file_contents
            .as_ref()
            .map(|files| !files.is_empty())
            .unwrap_or(false);

    let temp_runtime_root = if should_stage_runtime_files && request.root.is_none() {
        Some(unique_temp_dir("runtime")?)
    } else {
        None
    };

    let runtime_root = match request.root.as_deref() {
        Some(dir) => PathBuf::from(dir),
        None => match temp_runtime_root.as_ref() {
            Some(dir) => dir.clone(),
            None => env::current_dir().map_err(|err| format!("failed to get current directory: {err}"))?,
        },
    };
    fs::create_dir_all(&runtime_root)
        .map_err(|err| format!("failed to create runtime root {}: {err}", runtime_root.display()))?;

    if should_stage_runtime_files {
        stage_runtime_files(&runtime_root, &request)?;
    }

    let script_path = if should_stage_runtime_files {
        resolve_staged_path(&runtime_root, &request.script_path)?
    } else {
        PathBuf::from(&request.script_path)
    };

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
    loader.input_index = request.input_index.map(PathBuf::from);
    loader.reference_file = request.reference_file.map(PathBuf::from);
    loader.reference_index = request.reference_index.map(PathBuf::from);

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

    let mut ffi_timings: Vec<StageTiming> = Vec::new();
    if request.auto_index.unwrap_or(false) {
        let auto_index_started = Instant::now();
        let cwd = env::current_dir().map_err(|err| format!("failed to get cwd: {err}"))?;
        let effective_cache = request
            .cache_dir
            .as_ref()
            .map(PathBuf::from)
            .map(|path| {
                if path.is_absolute() {
                    path
                } else {
                    runtime_root.join(path)
                }
            })
            .unwrap_or_else(|| runtime_root.join(".bioscript-cache"));
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

    let runtime = BioscriptRuntime::with_config(
        runtime_root,
        RuntimeConfig { limits, loader },
    )
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
            request.trace_report_path.as_deref().map(std::path::Path::new),
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

fn write_timing_report(path: &PathBuf, timings: &[StageTiming]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .map_err(|err| format!("failed to create timing report dir {}: {err}", parent.display()))?;
    }
    let mut output = String::from("stage\tduration_ms\tdetail\n");
    for timing in timings {
        output.push_str(&format!(
            "{}\t{}\t{}\n",
            timing.stage,
            timing.duration_ms,
            timing.detail.replace('\t', " ")
        ));
    }
    fs::write(path, output).map_err(|err| format!("failed to write timing report {}: {err}", path.display()))
}

#[unsafe(no_mangle)]
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

        env.new_string(response).expect("jni new_string should succeed")
    }
}

#[cfg(test)]
mod tests {
    use super::{run_file_request, RunFileRequest};
    use std::{
        fs,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn unique_path(label: &str) -> String {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock drift")
            .as_nanos();
        std::env::temp_dir()
            .join(format!("bioscript-ffi-test-{label}-{}-{nanos}", std::process::id()))
            .to_string_lossy()
            .to_string()
    }

    #[test]
    fn run_file_request_stages_in_memory_assay_package() {
        let root = unique_path("staging");
        let output_rel = "outputs/result.tsv";

        let request = RunFileRequest {
            script_path: "assets/assays/herc2/herc2.py".to_owned(),
            script_contents: Some(
                r#"bioscript.write_tsv(output_file, [{"gene": "HERC2", "row_status": "matched"}])
"#
                .to_owned(),
            ),
            root: Some(root.clone()),
            input_file: None,
            input_contents: None,
            output_file: Some(output_rel.to_owned()),
            file_contents: Some(std::collections::HashMap::from([(
                "assets/assays/herc2/assay.yaml".to_owned(),
                "schema: bioscript:assay\n".to_owned(),
            )])),
            participant_id: None,
            trace_report_path: None,
            timing_report_path: None,
            input_format: None,
            input_index: None,
            reference_file: None,
            reference_index: None,
            auto_index: None,
            cache_dir: None,
            max_duration_ms: Some(5_000),
            max_memory_bytes: None,
            max_allocations: None,
            max_recursion_depth: None,
        };

        run_file_request(request).expect("run_file_request should succeed");

        let output = fs::read_to_string(format!("{root}/{output_rel}"))
            .expect("expected output file to be written");
        assert!(output.contains("HERC2"));
        assert!(output.contains("matched"));

        let staged_assay = fs::read_to_string(format!("{root}/assets/assays/herc2/assay.yaml"))
            .expect("expected assay file to be staged");
        assert!(staged_assay.contains("bioscript:assay"));
    }
}
