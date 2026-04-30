use std::{
    env,
    ffi::{CStr, CString},
    fmt::Write as _,
    fs,
    os::raw::c_char,
    path::PathBuf,
    time::{Duration, Instant},
};

use bioscript_core::{Assembly, VariantKind, VariantObservation, VariantSpec};
use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, PrepareRequest, prepare_indexes,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig, StageTiming};
use bioscript_schema::load_variant_manifest_text_for_lookup;
use monty::{MontyObject, ResourceLimits};
use serde::{Deserialize, Serialize};

const DEFAULT_MAX_DURATION_MS: u64 = 100;
const DEFAULT_MAX_MEMORY_BYTES: usize = 8 * 1024 * 1024;
const DEFAULT_MAX_ALLOCATIONS: usize = 200_000;
const DEFAULT_MAX_RECURSION_DEPTH: usize = 200;
const HARD_MAX_DURATION_MS: u64 = 60_000;
const HARD_MAX_MEMORY_BYTES: usize = 256 * 1024 * 1024;
const HARD_MAX_ALLOCATIONS: usize = 10_000_000;
const HARD_MAX_RECURSION_DEPTH: usize = 10_000;

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
    pub allow_md5_mismatch: Option<bool>,
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

#[derive(Debug, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunVariantYamlRequest {
    pub yaml_path: String,
    pub genome_path: String,
    pub input_format: Option<String>,
    pub input_index: Option<String>,
    pub reference_file: Option<String>,
    pub reference_index: Option<String>,
    pub allow_md5_mismatch: Option<bool>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct RunVariantYamlResult {
    pub observations: Vec<VariantObservationResult>,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct VariantObservationResult {
    pub name: String,
    pub backend: String,
    #[serde(rename = "ref", skip_serializing_if = "Option::is_none")]
    pub reference: Option<String>,
    #[serde(rename = "alt", skip_serializing_if = "Option::is_none")]
    pub alternate: Option<String>,
    #[serde(rename = "matchedRsid", skip_serializing_if = "Option::is_none")]
    pub matched_rsid: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub assembly: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genotype: Option<String>,
    #[serde(rename = "refCount", skip_serializing_if = "Option::is_none")]
    pub ref_count: Option<u32>,
    #[serde(rename = "altCount", skip_serializing_if = "Option::is_none")]
    pub alt_count: Option<u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub depth: Option<u32>,
    #[serde(rename = "rawCounts")]
    pub raw_counts: std::collections::BTreeMap<String, u32>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub decision: Option<String>,
    pub evidence: Vec<String>,
}

#[derive(Debug, Clone)]
struct NamedVariantSpec {
    name: String,
    spec: VariantSpec,
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

/// Runs a BioScript variant YAML assay against a supported genome file.
///
/// The native desktop/mobile path uses this instead of the web WASM exports.
/// It intentionally mirrors the web variant YAML flow: compile YAML through
/// `bioscript-schema`, choose the preferred assembly-specific variant, and
/// execute lookup through `bioscript-formats`.
pub fn run_variant_yaml_request(
    request: RunVariantYamlRequest,
) -> Result<RunVariantYamlResult, String> {
    let yaml_path = PathBuf::from(&request.yaml_path);
    let yaml_text = fs::read_to_string(&yaml_path)
        .map_err(|err| format!("failed to read YAML assay {}: {err}", yaml_path.display()))?;
    let variants = compile_variant_yaml_named(
        yaml_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("variant.yaml"),
        &yaml_text,
    )?;
    let selected = select_preferred_assembly_variants(&request.genome_path, variants);

    let mut loader = GenotypeLoadOptions::default();
    if let Some(value) = request.input_format.as_deref() {
        if value.eq_ignore_ascii_case("auto") {
            loader.format = None;
        } else {
            loader.format = Some(
                value
                    .parse::<GenotypeSourceFormat>()
                    .map_err(|err| format!("invalid inputFormat value {value}: {err}"))?,
            );
        }
    }
    loader.input_index = request.input_index.map(PathBuf::from);
    loader.reference_file = request.reference_file.map(PathBuf::from);
    loader.reference_index = request.reference_index.map(PathBuf::from);
    loader.allow_reference_md5_mismatch = request.allow_md5_mismatch.unwrap_or(false);

    let genome_path = PathBuf::from(&request.genome_path);
    let store = GenotypeStore::from_file_with_options(&genome_path, &loader)
        .map_err(|err| format!("failed to load genome {}: {err}", genome_path.display()))?;
    let specs = selected
        .iter()
        .map(|variant| variant.spec.clone())
        .collect::<Vec<_>>();
    let observations = store
        .lookup_variants(&specs)
        .map_err(|err| format!("variant lookup failed: {err}"))?;

    Ok(RunVariantYamlResult {
        observations: selected
            .into_iter()
            .zip(observations)
            .map(|(variant, observation)| observation_result(variant, observation))
            .collect(),
    })
}

fn runtime_root(request: &RunFileRequest) -> Result<PathBuf, String> {
    match request.root.as_deref() {
        Some(dir) => Ok(PathBuf::from(dir)),
        None => env::current_dir().map_err(|err| format!("failed to get current directory: {err}")),
    }
}

fn compile_variant_yaml_named(name: &str, text: &str) -> Result<Vec<NamedVariantSpec>, String> {
    let manifest = load_variant_manifest_text_for_lookup(name, text)
        .map_err(|err| format!("compile variant YAML failed: {err}"))?;
    let mut out = Vec::new();
    if let Some(locus) = manifest.spec.grch38.clone() {
        let mut spec = manifest.spec.clone();
        spec.grch37 = None;
        spec.grch38 = Some(locus);
        out.push(NamedVariantSpec {
            name: manifest.name.clone(),
            spec,
        });
    }
    if let Some(locus) = manifest.spec.grch37.clone() {
        let mut spec = manifest.spec;
        spec.grch37 = Some(locus);
        spec.grch38 = None;
        out.push(NamedVariantSpec {
            name: if out.is_empty() {
                manifest.name.clone()
            } else {
                format!("{}_grch37", manifest.name)
            },
            spec,
        });
    }
    if out.is_empty() {
        return Err(format!("variant {} has no coordinates", manifest.name));
    }
    Ok(out)
}

fn select_preferred_assembly_variants(
    genome_name: &str,
    variants: Vec<NamedVariantSpec>,
) -> Vec<NamedVariantSpec> {
    let target = infer_genome_assembly(genome_name).unwrap_or(Assembly::Grch38);
    let mut groups: Vec<(String, Vec<NamedVariantSpec>)> = Vec::new();
    for variant in variants {
        let key = variant_group_key(&variant);
        if let Some((_, values)) = groups.iter_mut().find(|(candidate, _)| candidate == &key) {
            values.push(variant);
        } else {
            groups.push((key, vec![variant]));
        }
    }

    groups
        .into_iter()
        .flat_map(|(_, group)| {
            let has_multiple_assemblies = group
                .iter()
                .filter_map(|variant| variant_assembly(&variant.spec))
                .fold(Vec::new(), |mut assemblies, assembly| {
                    if !assemblies.contains(&assembly) {
                        assemblies.push(assembly);
                    }
                    assemblies
                })
                .len()
                > 1;
            if group.len() <= 1 || !has_multiple_assemblies {
                return group;
            }
            let fallback = group[0].clone();
            vec![
                group
                    .iter()
                    .find(|variant| variant_assembly(&variant.spec) == Some(target))
                    .or_else(|| {
                        group.iter().find(|variant| {
                            variant_assembly(&variant.spec) == Some(Assembly::Grch38)
                        })
                    })
                    .cloned()
                    .unwrap_or(fallback),
            ]
        })
        .collect()
}

fn infer_genome_assembly(name: &str) -> Option<Assembly> {
    let lower = name.to_ascii_lowercase();
    if lower.contains("grch38") || lower.contains("hg38") {
        return Some(Assembly::Grch38);
    }
    if lower.contains("grch37") || lower.contains("hg19") {
        return Some(Assembly::Grch37);
    }
    None
}

fn variant_group_key(variant: &NamedVariantSpec) -> String {
    if let Some(rsid) = variant.spec.rsids.first() {
        return format!(
            "{}|{}|{}|{}",
            rsid.to_ascii_lowercase(),
            variant
                .spec
                .reference
                .as_deref()
                .unwrap_or("")
                .to_ascii_uppercase(),
            variant
                .spec
                .alternate
                .as_deref()
                .unwrap_or("")
                .to_ascii_uppercase(),
            variant_kind_label(variant.spec.kind)
        );
    }
    variant
        .name
        .trim_end_matches("_grch37")
        .trim_end_matches("_grch38")
        .to_ascii_lowercase()
}

fn variant_assembly(spec: &VariantSpec) -> Option<Assembly> {
    if spec.grch37.is_some() {
        return Some(Assembly::Grch37);
    }
    if spec.grch38.is_some() {
        return Some(Assembly::Grch38);
    }
    None
}

fn observation_result(
    variant: NamedVariantSpec,
    observation: VariantObservation,
) -> VariantObservationResult {
    VariantObservationResult {
        name: variant.name,
        backend: observation.backend,
        reference: variant.spec.reference,
        alternate: variant.spec.alternate,
        matched_rsid: observation.matched_rsid,
        assembly: observation.assembly.map(assembly_label),
        genotype: observation.genotype,
        ref_count: observation.ref_count,
        alt_count: observation.alt_count,
        depth: observation.depth,
        raw_counts: observation.raw_counts,
        decision: observation.decision,
        evidence: observation.evidence,
    }
}

fn assembly_label(assembly: Assembly) -> String {
    match assembly {
        Assembly::Grch37 => "grch37",
        Assembly::Grch38 => "grch38",
    }
    .to_owned()
}

fn variant_kind_label(kind: Option<VariantKind>) -> &'static str {
    match kind {
        Some(VariantKind::Snp) => "snv",
        Some(VariantKind::Insertion) => "insertion",
        Some(VariantKind::Deletion) => "deletion",
        Some(VariantKind::Indel) => "indel",
        Some(VariantKind::Other) => "other",
        None => "",
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
    loader.allow_reference_md5_mismatch = request.allow_md5_mismatch.unwrap_or(false);
    Ok(loader)
}

fn build_limits(request: &RunFileRequest) -> ResourceLimits {
    let mut limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(DEFAULT_MAX_DURATION_MS))
        .max_memory(DEFAULT_MAX_MEMORY_BYTES)
        .max_allocations(DEFAULT_MAX_ALLOCATIONS)
        .gc_interval(1000)
        .max_recursion_depth(Some(DEFAULT_MAX_RECURSION_DEPTH));

    if let Some(value) = request.max_duration_ms {
        limits = limits.max_duration(Duration::from_millis(value.min(HARD_MAX_DURATION_MS)));
    }
    if let Some(value) = request.max_memory_bytes {
        limits = limits.max_memory(value.min(HARD_MAX_MEMORY_BYTES));
    }
    if let Some(value) = request.max_allocations {
        limits = limits.max_allocations(value.min(HARD_MAX_ALLOCATIONS));
    }
    if let Some(value) = request.max_recursion_depth {
        limits = limits.max_recursion_depth(Some(value.min(HARD_MAX_RECURSION_DEPTH)));
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

#[cfg(test)]
mod tests {
    use super::*;

    fn request_with_limits() -> RunFileRequest {
        RunFileRequest {
            script_path: "script.py".to_owned(),
            root: None,
            input_file: None,
            output_file: None,
            participant_id: None,
            trace_report_path: None,
            timing_report_path: None,
            input_format: None,
            input_index: None,
            reference_file: None,
            reference_index: None,
            allow_md5_mismatch: None,
            auto_index: None,
            cache_dir: None,
            max_duration_ms: Some(u64::MAX),
            max_memory_bytes: Some(usize::MAX),
            max_allocations: Some(usize::MAX),
            max_recursion_depth: Some(usize::MAX),
        }
    }

    #[test]
    fn ffi_resource_limits_are_clamped_to_hard_ceilings() {
        let limits = build_limits(&request_with_limits());

        assert_eq!(
            limits.max_duration,
            Some(Duration::from_millis(HARD_MAX_DURATION_MS))
        );
        assert_eq!(limits.max_memory, Some(HARD_MAX_MEMORY_BYTES));
        assert_eq!(limits.max_allocations, Some(HARD_MAX_ALLOCATIONS));
        assert_eq!(limits.max_recursion_depth, Some(HARD_MAX_RECURSION_DEPTH));
    }
}
