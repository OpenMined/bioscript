use super::*;

pub(crate) struct WasmReportAnalysisRunner<'a> {
    pub(crate) workspace: &'a PackageWorkspace,
    pub(crate) input_name: &'a str,
    pub(crate) input_bytes: &'a [u8],
    pub(crate) input_index_bytes: Option<&'a [u8]>,
    pub(crate) reference_bytes: Option<&'a [u8]>,
    pub(crate) reference_index_bytes: Option<&'a [u8]>,
    pub(crate) participant_id: &'a str,
    pub(crate) loader: &'a GenotypeLoadOptions,
    pub(crate) options: &'a ReportOptionsInput,
    pub(crate) extra_artifacts: std::cell::RefCell<Vec<ReportArtifactOutput>>,
}

impl bioscript_reporting::ReportAnalysisRunner for WasmReportAnalysisRunner<'_> {
    fn run_analysis_task(
        &self,
        task: &bioscript_reporting::AnalysisManifestTask,
        _observation_rows: &[BTreeMap<String, String>],
        variant_observations: &[VariantObservation],
        _observations: &[serde_json::Value],
    ) -> Result<Vec<serde_json::Value>, String> {
        self.workspace
            .run_interpretations(
                &task.manifest_path,
                &task.manifest_name,
                &task.interpretations,
                self.input_name,
                self.input_bytes,
                self.input_index_bytes,
                self.reference_bytes,
                self.reference_index_bytes,
                variant_observations,
                self.participant_id,
                self.loader,
                self.options,
                &self.extra_artifacts,
            )
            .map_err(|err| format!("{err:?}"))
    }
}

impl PackageWorkspace {
    #[allow(clippy::too_many_arguments)]
    fn run_interpretations(
        &self,
        manifest_path: &str,
        manifest_name: &str,
        interpretations: &[PanelInterpretation],
        _input_name: &str,
        input_bytes: &[u8],
        input_index_bytes: Option<&[u8]>,
        reference_bytes: Option<&[u8]>,
        reference_index_bytes: Option<&[u8]>,
        preloaded_observations: &[VariantObservation],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
        extra_artifacts: &std::cell::RefCell<Vec<ReportArtifactOutput>>,
    ) -> Result<Vec<serde_json::Value>, JsError> {
        let mut outputs = Vec::new();
        for interpretation in interpretations {
            bioscript_reporting::validate_bioscript_interpretation(interpretation)
                .map_err(|err| JsError::new(&err))?;
            let script_path = self.resolve(manifest_path, &interpretation.path)?;
            let analysis_format = bioscript_reporting::analysis_output_format(
                interpretation.output_format.as_deref(),
            )
            .map_err(|err| JsError::new(&err))?;
            let analysis_output_file = bioscript_reporting::analysis_output_relative_file(
                participant_id,
                &interpretation.id,
                analysis_format.extension,
            );
            let output_file = options
                .output_dir
                .as_deref()
                .filter(|dir| !dir.is_empty())
                .map(|dir| format!("{}/{}", dir.trim_end_matches('/'), analysis_output_file))
                .unwrap_or(analysis_output_file);
            let observations_output_file = bioscript_reporting::analysis_observations_relative_file(
                participant_id,
                &interpretation.id,
            );
            let _observations_file = options
                .output_dir
                .as_deref()
                .filter(|dir| !dir.is_empty())
                .map(|dir| format!("{}/{}", dir.trim_end_matches('/'), observations_output_file))
                .unwrap_or(observations_output_file);
            let output_extension = Path::new(&output_file)
                .extension()
                .and_then(|value| value.to_str())
                .unwrap_or(analysis_format.extension);
            let virtual_input_file = virtual_genotype_input_path(loader);
            let virtual_input_index =
                input_index_bytes.map(|_| virtual_genotype_index_path(loader));
            let virtual_reference_file = reference_bytes.map(|_| "/input/reference.fa".to_owned());
            let virtual_reference_index =
                reference_index_bytes.map(|_| "/input/reference.fa.fai".to_owned());
            let virtual_output_file = format!("/output/results.{output_extension}");
            let virtual_observations_file = "/work/observations.tsv".to_owned();
            let script_virtual_path = virtual_pipeline_path(&script_path, "analysis.py");
            let manifest_virtual_path = virtual_pipeline_path(manifest_path, "manifest.yaml");
            let mut virtual_text_files = BTreeMap::new();
            virtual_text_files.insert(
                script_virtual_path.clone(),
                self.text(&script_path)?.to_owned(),
            );
            virtual_text_files.insert(
                manifest_virtual_path.clone(),
                self.text(manifest_path)?.to_owned(),
            );
            virtual_text_files.insert(
                virtual_observations_file.clone(),
                bioscript_reporting::render_analysis_observations_tsv(preloaded_observations),
            );
            let mut asset_paths = BTreeMap::new();
            for asset in &interpretation.assets {
                let asset_path = self.resolve(manifest_path, &asset.path)?;
                let virtual_asset_path = virtual_pipeline_path(&asset_path, &asset.path);
                virtual_text_files.insert(
                    virtual_asset_path.clone(),
                    self.text(&asset_path)?.to_owned(),
                );
                asset_paths.insert(asset.id.clone(), virtual_asset_path);
            }
            let runtime_observations = if input_bytes.is_empty() {
                preloaded_observations.to_vec()
            } else {
                Vec::new()
            };
            let mut virtual_binary_files = BTreeMap::new();
            virtual_binary_files.insert(virtual_input_file.clone(), input_bytes.to_vec());
            if let (Some(path), Some(bytes)) = (&virtual_input_index, input_index_bytes) {
                virtual_binary_files.insert(path.clone(), bytes.to_vec());
            }
            if let (Some(path), Some(bytes)) = (&virtual_reference_file, reference_bytes) {
                virtual_binary_files.insert(path.clone(), bytes.to_vec());
            }
            if let (Some(path), Some(bytes)) = (&virtual_reference_index, reference_index_bytes) {
                virtual_binary_files.insert(path.clone(), bytes.to_vec());
            }
            let mut runtime_loader = loader.clone();
            runtime_loader.input_index = virtual_input_index.as_ref().map(PathBuf::from);
            runtime_loader.reference_file = virtual_reference_file.as_ref().map(PathBuf::from);
            runtime_loader.reference_index = virtual_reference_index.as_ref().map(PathBuf::from);
            let limits = ResourceLimits::new()
                .max_duration(Duration::from_millis(options.analysis_max_duration_ms))
                .max_memory(16 * 1024 * 1024)
                .max_allocations(400_000)
                .gc_interval(1000)
                .max_recursion_depth(Some(200));
            let runtime = BioscriptRuntime::with_config(
                PathBuf::new(),
                RuntimeConfig {
                    limits,
                    loader: runtime_loader,
                    context: analysis_context(
                        participant_id,
                        &virtual_input_file,
                        virtual_input_index.as_deref(),
                        virtual_reference_file.as_deref(),
                        virtual_reference_index.as_deref(),
                        &script_virtual_path,
                        &manifest_virtual_path,
                        &asset_paths,
                        &virtual_observations_file,
                        &virtual_output_file,
                    ),
                    virtual_binary_files,
                    virtual_text_files: std::mem::take(&mut virtual_text_files),
                    preloaded_observations: runtime_observations,
                },
            )
            .map_err(|err| JsError::new(&format!("create analysis runtime failed: {err:?}")))?;
            runtime
                .run_file(
                    &script_virtual_path,
                    None,
                    vec![
                        (
                            "input_file",
                            MontyObject::String(virtual_input_file.clone()),
                        ),
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
                        (
                            "output_file",
                            MontyObject::String(virtual_output_file.clone()),
                        ),
                        (
                            "observations_file",
                            MontyObject::String(virtual_observations_file.clone()),
                        ),
                        ("asset_paths", monty_string_dict(&asset_paths)),
                        (
                            "participant_id",
                            MontyObject::String(participant_id.to_owned()),
                        ),
                    ],
                )
                .map_err(|err| {
                    JsError::new(&format!("analysis {} failed: {err:?}", interpretation.id))
                })?;
            let written = runtime.virtual_written_text_files();
            let text = written.get(&virtual_output_file).ok_or_else(|| {
                JsError::new(&format!(
                    "analysis {} did not write {virtual_output_file}",
                    interpretation.id
                ))
            })?;
            extra_artifacts.borrow_mut().extend(output_artifacts(
                &written,
                &runtime.virtual_written_binary_files(),
                &virtual_output_file,
            )?);
            let (rows, row_headers) = parse_analysis_output_text(text, analysis_format.format)?;
            outputs.push(bioscript_reporting::analysis_output_json(
                bioscript_reporting::AnalysisOutputJsonInput {
                    participant_id,
                    assay_id: manifest_name,
                    interpretation,
                    output_format: analysis_format.format,
                    manifest_path,
                    script_path: &script_path,
                    output_file: &output_file,
                    observations_file: Some(&virtual_observations_file),
                    row_headers,
                    rows,
                },
            ));
        }
        Ok(outputs)
    }
}

fn output_artifacts(
    text_files: &BTreeMap<String, String>,
    binary_files: &BTreeMap<String, Vec<u8>>,
    primary_virtual_output_file: &str,
) -> Result<Vec<ReportArtifactOutput>, JsError> {
    let mut artifacts = Vec::new();
    for (path, text) in text_files {
        if path == primary_virtual_output_file {
            continue;
        }
        let Some(relative) = safe_output_relative_path(path)? else {
            continue;
        };
        artifacts.push(ReportArtifactOutput {
            name: artifact_name(&relative),
            path: relative.clone(),
            mime_type: mime_type_for_path(&relative).to_owned(),
            text: Some(text.clone()),
            bytes: None,
            primary: false,
        });
    }
    for (path, bytes) in binary_files {
        let Some(relative) = safe_output_relative_path(path)? else {
            continue;
        };
        artifacts.push(ReportArtifactOutput {
            name: artifact_name(&relative),
            path: relative.clone(),
            mime_type: mime_type_for_path(&relative).to_owned(),
            text: None,
            bytes: Some(bytes.clone()),
            primary: false,
        });
    }
    artifacts.sort_by(|a, b| a.path.cmp(&b.path));
    Ok(artifacts)
}

fn safe_output_relative_path(path: &str) -> Result<Option<String>, JsError> {
    let Some(relative) = path.strip_prefix("/output/") else {
        return Ok(None);
    };
    let relative_path = Path::new(relative);
    if relative_path
        .components()
        .any(|component| !matches!(component, std::path::Component::Normal(_)))
    {
        return Err(JsError::new(&format!(
            "analysis wrote invalid output path {path}"
        )));
    }
    Ok(Some(relative.replace('\\', "/")))
}

fn artifact_name(path: &str) -> String {
    Path::new(path)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(path)
        .to_owned()
}

fn mime_type_for_path(path: &str) -> &'static str {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".html") {
        "text/html"
    } else if lower.ends_with(".tsv") {
        "text/tab-separated-values"
    } else if lower.ends_with(".json") || lower.ends_with(".jsonl") {
        "application/json"
    } else if lower.ends_with(".vcf") || lower.ends_with(".bed") || lower.ends_with(".log") {
        "text/plain"
    } else if lower.ends_with(".gz") {
        "application/gzip"
    } else if lower.ends_with(".bam")
        || lower.ends_with(".bai")
        || lower.ends_with(".cram")
        || lower.ends_with(".crai")
        || lower.ends_with(".csi")
    {
        "application/octet-stream"
    } else {
        "text/plain"
    }
}

fn virtual_pipeline_path(path: &str, fallback: &str) -> String {
    let name = Path::new(path)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(fallback);
    format!("/input/pipeline/{name}")
}

fn virtual_genotype_input_path(loader: &GenotypeLoadOptions) -> String {
    match loader.format {
        Some(bioscript_formats::GenotypeSourceFormat::Bam) => "/input/genotypes.bam",
        Some(bioscript_formats::GenotypeSourceFormat::Cram) => "/input/genotypes.cram",
        _ => "/input/genotypes",
    }
    .to_owned()
}

fn virtual_genotype_index_path(loader: &GenotypeLoadOptions) -> String {
    match loader.format {
        Some(bioscript_formats::GenotypeSourceFormat::Bam) => "/input/genotypes.bam.bai",
        Some(bioscript_formats::GenotypeSourceFormat::Cram) => "/input/genotypes.cram.crai",
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
) -> BTreeMap<String, MontyObject> {
    BTreeMap::from([
        (
            "participant_id".to_owned(),
            MontyObject::String(participant_id.to_owned()),
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
            MontyObject::String(observations_file.to_owned()),
        ),
        (
            "output_file".to_owned(),
            MontyObject::String(output_file.to_owned()),
        ),
        (
            "input_index".to_owned(),
            optional_string_object(input_index),
        ),
        ("input_bai".to_owned(), optional_string_object(input_index)),
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

fn optional_string_object(value: Option<&str>) -> MontyObject {
    value.map_or(MontyObject::None, |value| {
        MontyObject::String(value.to_owned())
    })
}

fn monty_string_dict(values: &BTreeMap<String, String>) -> MontyObject {
    MontyObject::Dict(
        values
            .iter()
            .map(|(key, value)| {
                (
                    MontyObject::String(key.clone()),
                    MontyObject::String(value.clone()),
                )
            })
            .collect(),
    )
}
