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
            output_dir: options.output_dir,
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

struct BioscriptAnalysisScriptInput<'a> {
    runtime_root: &'a Path,
    manifest_path: &'a Path,
    interpretation: &'a PanelInterpretation,
    script_path: &'a Path,
    input_file: &'a Path,
    output_file: &'a Path,
    output_dir: &'a Path,
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
    let input_bytes = fs::read(input.input_file)
        .map_err(|err| format!("failed to read analysis input {}: {err}", input.input_file.display()))?;
    let observations_text = fs::read_to_string(input.observations_file).map_err(|err| {
        format!(
            "failed to read analysis observations {}: {err}",
            input.observations_file.display()
        )
    })?;
    let script_text = fs::read_to_string(input.script_path)
        .map_err(|err| format!("failed to read analysis script {}: {err}", input.script_path.display()))?;
    let virtual_input_file = virtual_genotype_input_path(input.loader);
    let virtual_observations_file = "/work/observations.tsv".to_owned();
    let output_extension = input
        .output_file
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or("tsv");
    let virtual_output_file = format!("/output/results.{output_extension}");
    let script_virtual_path = virtual_pipeline_path(input.script_path, "analysis.py");
    let manifest_virtual_path = virtual_pipeline_path(input.manifest_path, "manifest.yaml");
    let virtual_input_index = input
        .loader
        .input_index
        .as_ref()
        .map(|_| virtual_genotype_index_path(input.loader));
    let virtual_reference_file = input
        .loader
        .reference_file
        .as_ref()
        .map(|_| "/input/reference.fa".to_owned());
    let virtual_reference_index = input
        .loader
        .reference_index
        .as_ref()
        .map(|_| "/input/reference.fa.fai".to_owned());
    let AnalysisVirtualTextFiles {
        text_files: virtual_text_files,
        asset_paths,
    } = collect_analysis_virtual_text_files(
        input,
        &script_virtual_path,
        script_text,
        &manifest_virtual_path,
        &virtual_observations_file,
        observations_text,
    )?;
    let mut virtual_binary_files = BTreeMap::new();
    virtual_binary_files.insert(virtual_input_file.clone(), input_bytes);
    if let (Some(real_path), Some(virtual_path)) =
        (&input.loader.input_index, &virtual_input_index)
    {
        virtual_binary_files.insert(
            virtual_path.clone(),
            fs::read(real_path).map_err(|err| {
                format!("failed to read analysis input index {}: {err}", real_path.display())
            })?,
        );
    }
    if let (Some(real_path), Some(virtual_path)) =
        (&input.loader.reference_file, &virtual_reference_file)
    {
        virtual_binary_files.insert(
            virtual_path.clone(),
            fs::read(real_path).map_err(|err| {
                format!("failed to read analysis reference {}: {err}", real_path.display())
            })?,
        );
    }
    if let (Some(real_path), Some(virtual_path)) =
        (&input.loader.reference_index, &virtual_reference_index)
    {
        virtual_binary_files.insert(
            virtual_path.clone(),
            fs::read(real_path).map_err(|err| {
                format!(
                    "failed to read analysis reference index {}: {err}",
                    real_path.display()
                )
            })?,
        );
    }
    let mut runtime_loader = input.loader.clone();
    runtime_loader.input_index = virtual_input_index.as_ref().map(PathBuf::from);
    runtime_loader.reference_file = virtual_reference_file.as_ref().map(PathBuf::from);
    runtime_loader.reference_index = virtual_reference_index.as_ref().map(PathBuf::from);
    let context = analysis_context(
        input.participant_id,
        &virtual_input_file,
        virtual_input_index.as_deref(),
        virtual_reference_file.as_deref(),
        virtual_reference_index.as_deref(),
        &script_virtual_path,
        &manifest_virtual_path,
        &asset_paths,
        &virtual_observations_file,
        &virtual_output_file,
    );
    let runtime = BioscriptRuntime::with_config(
        input.runtime_root.to_path_buf(),
        RuntimeConfig {
            limits,
            loader: runtime_loader,
            context,
            virtual_binary_files,
            virtual_text_files,
            ..RuntimeConfig::default()
        },
    )
    .map_err(|err| err.to_string())?;
    runtime
        .run_file(
            &script_virtual_path,
            None,
            vec![
                ("input_file", monty::MontyObject::String(virtual_input_file)),
                (
                    "input_index",
                    optional_string_object(virtual_input_index.as_deref()),
                ),
                (
                    "input_bai",
                    optional_string_object(virtual_input_index.as_deref()),
                ),
                (
                    "reference_fasta",
                    optional_string_object(virtual_reference_file.as_deref()),
                ),
                (
                    "alignment_reference_fasta",
                    optional_string_object(virtual_reference_file.as_deref()),
                ),
                (
                    "reference_index",
                    optional_string_object(virtual_reference_index.as_deref()),
                ),
                (
                    "alignment_reference_index",
                    optional_string_object(virtual_reference_index.as_deref()),
                ),
                ("output_file", monty::MontyObject::String(virtual_output_file.clone())),
                ("observations_file", monty::MontyObject::String(virtual_observations_file)),
                ("asset_paths", monty_string_dict(&asset_paths)),
                (
                    "participant_id",
                    monty::MontyObject::String(input.participant_id.to_owned()),
                ),
            ],
        )
        .map_err(|err| err.to_string())?;
    let written = runtime.virtual_written_text_files();
    let written_binary = runtime.virtual_written_binary_files();
    let output_text = written
        .get(&virtual_output_file)
        .ok_or_else(|| format!("analysis did not write {virtual_output_file}"))?;
    if let Some(parent) = input.output_file.parent() {
        fs::create_dir_all(parent)
            .map_err(|err| format!("failed to create analysis output dir {}: {err}", parent.display()))?;
    }
    fs::write(input.output_file, output_text).map_err(|err| {
        format!(
            "failed to write analysis output {}: {err}",
            input.output_file.display()
        )
    })?;
    persist_virtual_output_files(
        &written,
        &written_binary,
        &virtual_output_file,
        input.output_dir,
    )
}

struct AnalysisVirtualTextFiles {
    text_files: BTreeMap<String, String>,
    asset_paths: BTreeMap<String, String>,
}

fn collect_analysis_virtual_text_files(
    input: &BioscriptAnalysisScriptInput<'_>,
    script_virtual_path: &str,
    script_text: String,
    manifest_virtual_path: &str,
    virtual_observations_file: &str,
    observations_text: String,
) -> Result<AnalysisVirtualTextFiles, String> {
    let mut virtual_text_files = BTreeMap::new();
    virtual_text_files.insert(script_virtual_path.to_owned(), script_text);
    virtual_text_files.insert(
        manifest_virtual_path.to_owned(),
        fs::read_to_string(input.manifest_path).map_err(|err| {
            format!(
                "failed to read analysis manifest {}: {err}",
                input.manifest_path.display()
            )
        })?,
    );
    virtual_text_files.insert(virtual_observations_file.to_owned(), observations_text);
    let mut asset_paths = BTreeMap::new();
    for asset in &input.interpretation.assets {
        let asset_path = resolve_manifest_path(input.runtime_root, input.manifest_path, &asset.path)?;
        let virtual_asset_path = virtual_pipeline_path(&asset_path, &asset.path);
        let text = fs::read_to_string(&asset_path)
            .map_err(|err| format!("failed to read analysis asset {}: {err}", asset_path.display()))?;
        virtual_text_files.insert(virtual_asset_path.clone(), text);
        asset_paths.insert(asset.id.clone(), virtual_asset_path);
    }
    Ok(AnalysisVirtualTextFiles {
        text_files: virtual_text_files,
        asset_paths,
    })
}

fn persist_virtual_output_files(
    written: &BTreeMap<String, String>,
    written_binary: &BTreeMap<String, Vec<u8>>,
    primary_virtual_output_file: &str,
    output_dir: &Path,
) -> Result<(), String> {
    for (virtual_path, text) in written {
        if virtual_path == primary_virtual_output_file {
            continue;
        }
        let Some(output_path) = virtual_output_path(output_dir, virtual_path)? else {
            continue;
        };
        if let Some(parent) = output_path.parent() {
            fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create output dir {}: {err}", parent.display()))?;
        }
        fs::write(&output_path, text)
            .map_err(|err| format!("failed to write analysis output {}: {err}", output_path.display()))?;
    }
    for (virtual_path, bytes) in written_binary {
        let Some(output_path) = virtual_output_path(output_dir, virtual_path)? else {
            continue;
        };
        if let Some(parent) = output_path.parent() {
            fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create output dir {}: {err}", parent.display()))?;
        }
        fs::write(&output_path, bytes)
            .map_err(|err| format!("failed to write analysis output {}: {err}", output_path.display()))?;
    }
    Ok(())
}

fn virtual_output_path(output_dir: &Path, virtual_path: &str) -> Result<Option<PathBuf>, String> {
    let Some(relative) = virtual_path.strip_prefix("/output/") else {
        return Ok(None);
    };
    let relative_path = Path::new(relative);
    if relative_path
        .components()
        .any(|component| !matches!(component, std::path::Component::Normal(_)))
    {
        return Err(format!("analysis wrote invalid output path {virtual_path}"));
    }
    Ok(Some(output_dir.join(relative_path)))
}

fn virtual_pipeline_path(path: &Path, fallback: &str) -> String {
    let name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(fallback);
    format!("/input/pipeline/{name}")
}

fn virtual_genotype_input_path(loader: &GenotypeLoadOptions) -> String {
    match loader.format {
        Some(GenotypeSourceFormat::Bam) => "/input/genotypes.bam",
        Some(GenotypeSourceFormat::Cram) => "/input/genotypes.cram",
        _ => "/input/genotypes",
    }
    .to_owned()
}

fn virtual_genotype_index_path(loader: &GenotypeLoadOptions) -> String {
    match loader.format {
        Some(GenotypeSourceFormat::Bam) => "/input/genotypes.bam.bai",
        Some(GenotypeSourceFormat::Cram) => "/input/genotypes.cram.crai",
        _ => "/input/input.index",
    }
    .to_owned()
}

fn analysis_context(
    participant_id: &str,
    input_file: &str,
    input_index: Option<&str>,
    reference_file: Option<&str>,
    reference_index: Option<&str>,
    script_path: &str,
    manifest_path: &str,
    asset_paths: &BTreeMap<String, String>,
    observations_file: &str,
    output_file: &str,
) -> BTreeMap<String, monty::MontyObject> {
    BTreeMap::from([
        (
            "participant_id".to_owned(),
            monty::MontyObject::String(participant_id.to_owned()),
        ),
        (
            "input_files".to_owned(),
            monty_string_dict(&optional_input_files(
                input_file,
                input_index,
                reference_file,
                reference_index,
            )),
        ),
        (
            "pipeline_files".to_owned(),
            monty_string_dict(&BTreeMap::from([
                ("manifest".to_owned(), manifest_path.to_owned()),
                ("analysis".to_owned(), script_path.to_owned()),
            ])),
        ),
        ("assets".to_owned(), monty_string_dict(asset_paths)),
        (
            "observations_file".to_owned(),
            monty::MontyObject::String(observations_file.to_owned()),
        ),
        (
            "output_file".to_owned(),
            monty::MontyObject::String(output_file.to_owned()),
        ),
        (
            "input_index".to_owned(),
            optional_string_object(input_index),
        ),
        (
            "input_bai".to_owned(),
            optional_string_object(input_index),
        ),
        (
            "reference_fasta".to_owned(),
            optional_string_object(reference_file),
        ),
        (
            "alignment_reference_fasta".to_owned(),
            optional_string_object(reference_file),
        ),
        (
            "reference_index".to_owned(),
            optional_string_object(reference_index),
        ),
        (
            "alignment_reference_index".to_owned(),
            optional_string_object(reference_index),
        ),
    ])
}

fn optional_input_files(
    input_file: &str,
    input_index: Option<&str>,
    reference_file: Option<&str>,
    reference_index: Option<&str>,
) -> BTreeMap<String, String> {
    let mut files = BTreeMap::from([("genotypes".to_owned(), input_file.to_owned())]);
    if let Some(path) = input_index {
        files.insert("input_index".to_owned(), path.to_owned());
        files.insert("bai".to_owned(), path.to_owned());
    }
    if let Some(path) = reference_file {
        files.insert("reference_fasta".to_owned(), path.to_owned());
    }
    if let Some(path) = reference_index {
        files.insert("reference_index".to_owned(), path.to_owned());
    }
    files
}

fn optional_string_object(value: Option<&str>) -> monty::MontyObject {
    value.map_or(monty::MontyObject::None, |value| {
        monty::MontyObject::String(value.to_owned())
    })
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
        let headers = outputs[0]["row_headers"].as_array().unwrap();
        assert!(headers.contains(&serde_json::Value::String("participant".to_owned())));
        assert!(output
            .join("analysis/sample/score.observations.tsv")
            .exists());

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn run_interpretations_virtualizes_auxiliary_alignment_inputs() {
        let dir = temp_dir("analysis-auxiliary-inputs");
        let manifest = dir.join("assay.yaml");
        let script = dir.join("analysis.bs");
        let input = dir.join("sample.bam");
        let input_index = dir.join("sample.bam.bai");
        let reference = dir.join("reference.fa");
        let reference_index = dir.join("reference.fa.fai");
        let output = dir.join("out");
        fs::write(&input, b"alignment bytes").unwrap();
        fs::write(&input_index, b"index bytes").unwrap();
        fs::write(&reference, b">chr1\nACGT\n").unwrap();
        fs::write(&reference_index, b"chr1\t4\t6\t4\t5\n").unwrap();
        fs::write(&manifest, "schema: bioscript:assay:1.0\nname: assay-one\n").unwrap();
        fs::write(
            &script,
            r#"
def main():
    input_files = bioscript.context["input_files"]
    bioscript.write_tsv(output_file, [
        {
            "input_file": input_file,
            "input_index": input_index,
            "input_bai": input_bai,
            "context_input_index": bioscript.context["input_index"],
            "context_input_bai": bioscript.context["input_bai"],
            "input_files_index": input_files["input_index"],
            "reference_fasta": reference_fasta,
            "alignment_reference_fasta": alignment_reference_fasta,
            "reference_index": reference_index,
            "alignment_reference_index": alignment_reference_index,
            "input_files_reference": input_files["reference_fasta"],
            "input_files_reference_index": input_files["reference_index"],
        }
    ])

if __name__ == "__main__":
    main()
"#,
        )
        .unwrap();
        let interpretation = PanelInterpretation {
            id: "paths".to_owned(),
            label: Some("Paths".to_owned()),
            kind: "bioscript".to_owned(),
            path: "analysis.bs".to_owned(),
            output_format: Some("tsv".to_owned()),
            derived_from: Vec::new(),
            assets: Vec::new(),
            emits: Vec::new(),
            logic: None,
        };
        let loader = GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Bam),
            input_index: Some(input_index),
            reference_file: Some(reference),
            reference_index: Some(reference_index),
            ..GenotypeLoadOptions::default()
        };
        let options = ReportAnalysisOptions {
            runtime_root: &dir,
            input_file: &input,
            participant_id: "sample",
            loader: &loader,
            output_dir: &output,
            observation_rows: &[],
            max_duration_ms: 1000,
        };

        let outputs =
            run_interpretations_for_report(&manifest, "assay-one", &[interpretation], &options)
                .unwrap();
        let row = &outputs[0]["rows"][0];
        assert_eq!(row["input_file"], "/input/genotypes.bam");
        assert_eq!(row["input_index"], "/input/genotypes.bam.bai");
        assert_eq!(row["input_bai"], "/input/genotypes.bam.bai");
        assert_eq!(row["context_input_index"], "/input/genotypes.bam.bai");
        assert_eq!(row["context_input_bai"], "/input/genotypes.bam.bai");
        assert_eq!(row["input_files_index"], "/input/genotypes.bam.bai");
        assert_eq!(row["reference_fasta"], "/input/reference.fa");
        assert_eq!(row["alignment_reference_fasta"], "/input/reference.fa");
        assert_eq!(row["reference_index"], "/input/reference.fa.fai");
        assert_eq!(row["alignment_reference_index"], "/input/reference.fa.fai");
        assert_eq!(row["input_files_reference"], "/input/reference.fa");
        assert_eq!(row["input_files_reference_index"], "/input/reference.fa.fai");

        fs::remove_dir_all(dir).unwrap();
    }
}
