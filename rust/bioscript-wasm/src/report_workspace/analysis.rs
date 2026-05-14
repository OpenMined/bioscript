use super::*;

impl PackageWorkspace {
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn run_manifest_analyses(
        &self,
        manifest_path: &str,
        input_name: &str,
        input_bytes: &[u8],
        preloaded_observations: &[VariantObservation],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
    ) -> Result<Vec<serde_json::Value>, JsError> {
        let tasks = bioscript_reporting::collect_analysis_manifest_tasks(
            self,
            manifest_path,
            &options.filters,
        )
        .map_err(|err| JsError::new(&err))?;
        let mut analyses = Vec::new();
        for task in tasks {
            analyses.extend(self.run_interpretations(
                &task.manifest_path,
                &task.manifest_name,
                &task.interpretations,
                input_name,
                input_bytes,
                preloaded_observations,
                participant_id,
                loader,
                options,
            )?);
        }
        Ok(analyses)
    }

    #[allow(clippy::too_many_arguments)]
    fn run_interpretations(
        &self,
        manifest_path: &str,
        manifest_name: &str,
        interpretations: &[PanelInterpretation],
        _input_name: &str,
        input_bytes: &[u8],
        preloaded_observations: &[VariantObservation],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
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
            let virtual_input_file = "/input/genotypes".to_owned();
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
            let mut virtual_binary_files = BTreeMap::new();
            virtual_binary_files.insert(virtual_input_file.clone(), input_bytes.to_vec());
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
                    loader: loader.clone(),
                    context: analysis_context(
                        participant_id,
                        &virtual_input_file,
                        &script_virtual_path,
                        &manifest_virtual_path,
                        &asset_paths,
                        &virtual_observations_file,
                        &virtual_output_file,
                    ),
                    virtual_binary_files,
                    virtual_text_files,
                    preloaded_observations: preloaded_observations.to_vec(),
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

fn virtual_pipeline_path(path: &str, fallback: &str) -> String {
    let name = Path::new(path)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(fallback);
    format!("/input/pipeline/{name}")
}

fn analysis_context(
    participant_id: &str,
    input_file: &str,
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
            monty_string_dict(&BTreeMap::from([(
                "genotypes".to_owned(),
                input_file.to_owned(),
            )])),
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
    ])
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
