fn run_manifest_rows_for_report(
    runtime_root: &Path,
    manifest_path: &Path,
    input_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_text = input_file.display().to_string();
    let store = GenotypeStore::from_file_with_options(Path::new(&input_text), loader)
        .map_err(|err| err.to_string())?;
    let workspace = bioscript_reporting::FilesystemManifestWorkspace::new(runtime_root);
    let manifest_path_text = manifest_path.display().to_string();
    let tasks =
        bioscript_reporting::collect_variant_manifest_tasks(&workspace, &manifest_path_text, filters)?;
    let observations = store
        .lookup_variants(
            &tasks
                .iter()
                .map(|task| task.manifest.spec.clone())
                .collect::<Vec<_>>(),
        )
        .map_err(|err| err.to_string())?;
    Ok(tasks
        .into_iter()
        .zip(observations)
        .map(|(task, observation)| {
            let resolved = Path::new(&task.manifest_path);
            variant_row(
                runtime_root,
                resolved,
                &task.manifest.name,
                &task.manifest.tags,
                &observation,
                Some(participant_id),
            )
        })
        .collect())
}

struct ReportAnalysisOptions<'a> {
    runtime_root: &'a Path,
    input_file: &'a Path,
    participant_id: &'a str,
    loader: &'a GenotypeLoadOptions,
    output_dir: &'a Path,
    observation_rows: &'a [BTreeMap<String, String>],
    filters: &'a [String],
    max_duration_ms: u64,
}

fn run_manifest_analyses_for_report(
    manifest_path: &Path,
    options: &ReportAnalysisOptions<'_>,
) -> Result<Vec<serde_json::Value>, String> {
    let workspace = bioscript_reporting::FilesystemManifestWorkspace::new(options.runtime_root);
    let manifest_path_text = manifest_path.display().to_string();
    let mut analyses = Vec::new();
    for task in
        bioscript_reporting::collect_analysis_manifest_tasks(&workspace, &manifest_path_text, options.filters)?
    {
        analyses.extend(run_interpretations_for_report(
            Path::new(&task.manifest_path),
            &task.manifest_name,
            &task.interpretations,
            options,
        )?);
    }
    Ok(analyses)
}

fn run_interpretations_for_report(
    manifest_path: &Path,
    manifest_name: &str,
    interpretations: &[PanelInterpretation],
    options: &ReportAnalysisOptions<'_>,
) -> Result<Vec<serde_json::Value>, String> {
    let mut outputs = Vec::new();
    for interpretation in interpretations {
        bioscript_reporting::validate_bioscript_interpretation(interpretation)?;
        let script_path =
            resolve_manifest_path(options.runtime_root, manifest_path, &interpretation.path)?;
        let analysis_format =
            bioscript_reporting::analysis_output_format(interpretation.output_format.as_deref())?;
        let analysis_dir = options.output_dir.join("analysis").join(options.participant_id);
        fs::create_dir_all(&analysis_dir).map_err(|err| {
            format!(
                "failed to create analysis output dir {}: {err}",
                analysis_dir.display()
            )
        })?;
        let output_file = options.output_dir.join(
            bioscript_reporting::analysis_output_relative_file(
                options.participant_id,
                &interpretation.id,
                analysis_format.extension,
            ),
        );
        let observations_file = options.output_dir.join(
            bioscript_reporting::analysis_observations_relative_file(
                options.participant_id,
                &interpretation.id,
            ),
        );
        fs::write(
            &observations_file,
            bioscript_reporting::render_manifest_rows_tsv(options.observation_rows),
        )
        .map_err(|err| {
            format!(
                "failed to write analysis observations {}: {err}",
                observations_file.display()
            )
        })?;
        run_bioscript_analysis_script(&BioscriptAnalysisScriptInput {
            runtime_root: options.runtime_root,
            script_path: &script_path,
            input_file: options.input_file,
            output_file: &output_file,
            observations_file: &observations_file,
            participant_id: options.participant_id,
            loader: options.loader,
            analysis_max_duration_ms: options.max_duration_ms,
        })?;
        let (rows, row_headers) = parse_analysis_output(&output_file, analysis_format.format)?;
        let manifest_path_text = manifest_path
            .strip_prefix(options.runtime_root)
            .unwrap_or(manifest_path)
            .display()
            .to_string();
        let script_path_text = script_path
            .strip_prefix(options.runtime_root)
            .unwrap_or(&script_path)
            .display()
            .to_string();
        let output_file_text = output_file
            .strip_prefix(options.runtime_root)
            .unwrap_or(&output_file)
            .display()
            .to_string();
        let observations_file_text = observations_file
            .strip_prefix(options.runtime_root)
            .unwrap_or(&observations_file)
            .display()
            .to_string();
        outputs.push(bioscript_reporting::analysis_output_json(
            bioscript_reporting::AnalysisOutputJsonInput {
                participant_id: options.participant_id,
                assay_id: manifest_name,
                interpretation,
                output_format: analysis_format.format,
                manifest_path: &manifest_path_text,
                script_path: &script_path_text,
                output_file: &output_file_text,
                observations_file: Some(&observations_file_text),
                row_headers,
                rows,
            },
        ));
    }
    Ok(outputs)
}

struct BioscriptAnalysisScriptInput<'a> {
    runtime_root: &'a Path,
    script_path: &'a Path,
    input_file: &'a Path,
    output_file: &'a Path,
    observations_file: &'a Path,
    participant_id: &'a str,
    loader: &'a GenotypeLoadOptions,
    analysis_max_duration_ms: u64,
}

fn run_bioscript_analysis_script(input: &BioscriptAnalysisScriptInput<'_>) -> Result<(), String> {
    let limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(input.analysis_max_duration_ms))
        .max_memory(16 * 1024 * 1024)
        .max_allocations(400_000)
        .gc_interval(1000)
        .max_recursion_depth(Some(200));
    let runtime = BioscriptRuntime::with_config(
        input.runtime_root.to_path_buf(),
        RuntimeConfig {
            limits,
            loader: input.loader.clone(),
            ..RuntimeConfig::default()
        },
    )
    .map_err(|err| err.to_string())?;
    runtime
        .run_file(
            input.script_path,
            None,
            vec![
                (
                    "input_file",
                    monty::MontyObject::String(runtime_path_string(
                        input.runtime_root,
                        input.input_file,
                    )),
                ),
                (
                    "output_file",
                    monty::MontyObject::String(runtime_path_string(
                        input.runtime_root,
                        input.output_file,
                    )),
                ),
                (
                    "observations_file",
                    monty::MontyObject::String(runtime_path_string(
                        input.runtime_root,
                        input.observations_file,
                    )),
                ),
                (
                    "participant_id",
                    monty::MontyObject::String(input.participant_id.to_owned()),
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

fn parse_analysis_output(
    path: &Path,
    format: &str,
) -> Result<(Vec<serde_json::Value>, Vec<String>), String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read analysis output {}: {err}", path.display()))?;
    bioscript_reporting::parse_analysis_output_text(&text, format)
        .map_err(|err| format!("failed to parse analysis output {}: {err}", path.display()))
}

fn participant_id_from_path(path: &Path) -> String {
    bioscript_reporting::participant_id_from_path(path)
}

#[cfg(test)]
mod app_report_execution_tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-report-execution-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn path_and_analysis_output_helpers_normalize_values() {
        let root = Path::new("/tmp/runtime-root");
        let nested = root.join("analysis/p1/out.json");
        assert_eq!(runtime_path_string(root, &nested), "analysis/p1/out.json");
        assert_eq!(
            runtime_path_string(root, Path::new("/outside/file.txt")),
            "/outside/file.txt"
        );
        assert_eq!(
            participant_id_from_path(Path::new("/data/sample.vcf.gz")),
            "sample"
        );

        let dir = temp_dir("analysis-output");
        let json = dir.join("rows.json");
        fs::write(&json, r#"{"rows":[{"score":2,"label":"ok"}]}"#).unwrap();
        let (rows, headers) = parse_analysis_output(&json, "json").unwrap();
        assert_eq!(rows[0]["score"], 2);
        assert!(headers.contains(&"score".to_owned()));

        let bad = dir.join("bad.json");
        fs::write(&bad, "{").unwrap();
        assert!(parse_analysis_output(&bad, "json")
            .unwrap_err()
            .contains("failed to parse analysis output"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn run_interpretations_rejects_unsupported_kinds_before_runtime() {
        let dir = temp_dir("unsupported-kind");
        let manifest = dir.join("panel.yaml");
        let input = dir.join("input.txt");
        fs::write(&input, "rsid\tgenotype\nrs1\tA/G\n").unwrap();
        let interpretation = PanelInterpretation {
            id: "not-bioscript".to_owned(),
            label: Some("Not BioScript".to_owned()),
            kind: "python".to_owned(),
            path: "analysis.py".to_owned(),
            output_format: Some("json".to_owned()),
            derived_from: Vec::new(),
            emits: Vec::new(),
            logic: None,
        };
        let loader = GenotypeLoadOptions::default();
        let options = ReportAnalysisOptions {
            runtime_root: &dir,
            input_file: &input,
            participant_id: "p1",
            loader: &loader,
            output_dir: &dir,
            observation_rows: &[],
            filters: &[],
            max_duration_ms: 10,
        };

        let err =
            run_interpretations_for_report(&manifest, "assay", &[interpretation], &options)
                .unwrap_err();
        assert!(err.contains("unsupported kind"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn run_manifest_rows_for_report_reads_text_input_and_variant_manifest() {
        let dir = temp_dir("manifest-rows");
        let manifest = dir.join("variant.yaml");
        fs::write(
            &manifest,
            r#"
schema: bioscript:variant:1.0
version: "1.0"
name: rs1
gene: ABC
identifiers:
  rsids: [rs1]
coordinates:
  grch38:
    chrom: "1"
    pos: 100
alleles:
  kind: snv
  ref: A
  alts: [G]
"#,
        )
        .unwrap();
        let input = dir.join("sample.txt");
        fs::write(&input, "rsid\tgenotype\nrs1\tA/G\n").unwrap();
        let loader = GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Text),
            ..GenotypeLoadOptions::default()
        };

        let rows =
            run_manifest_rows_for_report(&dir, &manifest, &input, "p1", &loader, &[]).unwrap();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0]["participant_id"], "p1");
        assert_eq!(rows[0]["matched_rsid"], "rs1");
        assert_eq!(rows[0]["genotype"], "AG");

        let missing_input = dir.join("missing.txt");
        assert!(run_manifest_rows_for_report(
            &dir,
            &manifest,
            &missing_input,
            "p1",
            &loader,
            &[],
        )
        .unwrap_err()
        .contains("No such file"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn run_interpretations_executes_bioscript_analysis_and_builds_json_output() {
        let dir = temp_dir("analysis-success");
        let manifest = dir.join("assay.yaml");
        let script = dir.join("analysis.bs");
        let input = dir.join("sample.txt");
        let output = dir.join("out");
        fs::write(&input, "rsid\tgenotype\nrs1\tA/G\n").unwrap();
        fs::write(
            &script,
            r#"
def main():
    bioscript.write_tsv(output_file, [
        {"participant": participant_id, "score": 7, "source": input_file, "observations": observations_file}
    ])

if __name__ == "__main__":
    main()
"#,
        )
        .unwrap();
        let interpretation = PanelInterpretation {
            id: "score".to_owned(),
            label: Some("Score".to_owned()),
            kind: "bioscript".to_owned(),
            path: "analysis.bs".to_owned(),
            output_format: Some("tsv".to_owned()),
            derived_from: Vec::new(),
            emits: Vec::new(),
            logic: None,
        };
        let rows = [BTreeMap::from([
            ("participant_id".to_owned(), "sample".to_owned()),
            ("matched_rsid".to_owned(), "rs1".to_owned()),
            ("genotype".to_owned(), "AG".to_owned()),
        ])];
        let loader = GenotypeLoadOptions::default();
        let options = ReportAnalysisOptions {
            runtime_root: &dir,
            input_file: &input,
            participant_id: "sample",
            loader: &loader,
            output_dir: &output,
            observation_rows: &rows,
            filters: &[],
            max_duration_ms: 1000,
        };

        let outputs =
            run_interpretations_for_report(&manifest, "assay-one", &[interpretation], &options)
                .unwrap();
        assert_eq!(outputs.len(), 1);
        assert_eq!(outputs[0]["assay_id"], "assay-one");
        assert_eq!(outputs[0]["participant_id"], "sample");
        assert_eq!(outputs[0]["rows"][0]["score"], "7");
        let headers = outputs[0]["row_headers"].as_array().unwrap();
        assert!(headers.contains(&serde_json::Value::String("participant".to_owned())));
        assert!(output
            .join("analysis/sample/score.observations.tsv")
            .exists());

        fs::remove_dir_all(dir).unwrap();
    }
}
