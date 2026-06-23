struct ReportAnalysisOptions<'a> {
    runtime_root: &'a Path,
    input_file: &'a Path,
    participant_id: &'a str,
    loader: &'a GenotypeLoadOptions,
    output_dir: &'a Path,
    observation_rows: &'a [BTreeMap<String, String>],
    max_duration_ms: u64,
}

struct CliReportAnalysisRunner<'a> {
    runtime_root: &'a Path,
    input_file: &'a Path,
    participant_id: &'a str,
    loader: &'a GenotypeLoadOptions,
    output_dir: &'a Path,
    max_duration_ms: u64,
}

impl bioscript_reporting::ReportAnalysisRunner for CliReportAnalysisRunner<'_> {
    fn run_analysis_task(
        &self,
        task: &bioscript_reporting::AnalysisManifestTask,
        observation_rows: &[BTreeMap<String, String>],
        _variant_observations: &[bioscript_core::VariantObservation],
        _observations: &[serde_json::Value],
    ) -> Result<Vec<serde_json::Value>, String> {
        let options = ReportAnalysisOptions {
            runtime_root: self.runtime_root,
            input_file: self.input_file,
            participant_id: self.participant_id,
            loader: self.loader,
            output_dir: self.output_dir,
            observation_rows,
            max_duration_ms: self.max_duration_ms,
        };
        run_interpretations_for_report(
            Path::new(&task.manifest_path),
            &task.manifest_name,
            &task.interpretations,
            &options,
        )
    }
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
            manifest_path,
            interpretation,
            script_path: &script_path,
            input_file: options.input_file,
            output_file: &output_file,
            observations_file: &observations_file,
            participant_id: options.participant_id,
            loader: options.loader,
            analysis_max_duration_ms: options.max_duration_ms,
        })?;
        let (rows, row_headers) = parse_analysis_output(&output_file, analysis_format.format)?;
        let manifest_path_text = runtime_path_string(options.runtime_root, manifest_path);
        let script_path_text = runtime_path_string(options.runtime_root, &script_path);
        let output_file_text = runtime_path_string(options.runtime_root, &output_file);
        let observations_file_text =
            runtime_path_string(options.runtime_root, &observations_file);
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

fn participant_id_from_path(path: &Path) -> String {
    bioscript_reporting::participant_id_from_path(path)
}

fn runtime_path_string(runtime_root: &Path, path: &Path) -> String {
    path.strip_prefix(runtime_root)
        .unwrap_or(path)
        .display()
        .to_string()
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
            assets: Vec::new(),
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
            max_duration_ms: 10,
        };

        let err =
            run_interpretations_for_report(&manifest, "assay", &[interpretation], &options)
                .unwrap_err();
        assert!(err.contains("unsupported kind"));

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
        fs::write(&manifest, "schema: bioscript:assay:1.0\nname: assay-one\n").unwrap();
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
            assets: Vec::new(),
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
            max_duration_ms: 1000,
        };

        let outputs =
            run_interpretations_for_report(&manifest, "assay-one", &[interpretation], &options)
                .unwrap();
        assert_eq!(outputs.len(), 1);
        assert_eq!(outputs[0]["assay_id"], "assay-one");
        assert_eq!(outputs[0]["participant_id"], "sample");
        assert_eq!(outputs[0]["rows"][0]["score"], "7");
        assert_eq!(outputs[0]["rows"][0]["source"], "sample.txt");
        let headers = outputs[0]["row_headers"].as_array().unwrap();
        assert!(headers.contains(&serde_json::Value::String("participant".to_owned())));
        assert!(output
            .join("analysis/sample/score.observations.tsv")
            .exists());

        fs::remove_dir_all(dir).unwrap();
    }
}
