use std::{
    env,
    path::PathBuf,
    process::ExitCode,
    time::{Duration, Instant},
};

use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, PrepareRequest, prepare_indexes,
};
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig, StageTiming};
use monty::ResourceLimits;

mod commands;
mod manifest;
mod paths;

use commands::{run_inspect, run_prepare, run_validate_panels, run_validate_variants};
use manifest::{ManifestRunOptions, is_yaml_manifest, run_manifest};
use paths::{normalize_loader_paths, write_timing_report};

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::manifest::{
        manifest_schema, matches_filters, render_rows_as_tsv, resolve_manifest_path, variant_row,
    };
    use crate::paths::{normalize_loader_paths, resolve_cli_path, resolve_cli_path_buf};
    use bioscript_core::{Assembly, VariantObservation};
    use bioscript_schema::VariantManifest;
    use std::fs;
    use std::path::Path;
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
