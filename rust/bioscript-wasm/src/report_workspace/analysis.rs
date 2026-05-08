use super::*;

impl PackageWorkspace {
    pub(super) fn run_manifest_analyses(
        &self,
        manifest_path: &str,
        input_name: &str,
        input_bytes: &[u8],
        preloaded_observations: &[VariantObservation],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
    ) -> Result<Vec<serde_json::Value>, JsError> {
        match self.schema(manifest_path)?.as_str() {
            "bioscript:panel:1.0" => {
                let panel = self.load_panel(manifest_path)?;
                let mut analyses = self.run_interpretations(
                    manifest_path,
                    &panel.name,
                    &panel.interpretations,
                    input_name,
                    input_bytes,
                    preloaded_observations,
                    participant_id,
                    loader,
                    options,
                )?;
                for member in &panel.members {
                    if member.kind != "assay" {
                        continue;
                    }
                    let Some(path) = &member.path else {
                        continue;
                    };
                    let resolved = self.resolve(manifest_path, path)?;
                    analyses.extend(self.run_manifest_analyses(
                        &resolved,
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
            "bioscript:assay:1.0" => {
                let assay = self.load_assay(manifest_path)?;
                self.run_interpretations(
                    manifest_path,
                    &assay.name,
                    &assay.interpretations,
                    input_name,
                    input_bytes,
                    preloaded_observations,
                    participant_id,
                    loader,
                    options,
                )
            }
            _ => Ok(Vec::new()),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run_interpretations(
        &self,
        manifest_path: &str,
        manifest_name: &str,
        interpretations: &[PanelInterpretation],
        input_name: &str,
        input_bytes: &[u8],
        preloaded_observations: &[VariantObservation],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
    ) -> Result<Vec<serde_json::Value>, JsError> {
        let mut outputs = Vec::new();
        for interpretation in interpretations {
            if interpretation.kind != "bioscript" {
                return Err(JsError::new(&format!(
                    "analysis '{}' uses unsupported kind '{}'",
                    interpretation.id, interpretation.kind
                )));
            }
            let script_path = self.resolve(manifest_path, &interpretation.path)?;
            let output_file = format!(
                "analysis/{participant_id}/{}.{}",
                interpretation.id,
                interpretation.output_format.as_deref().unwrap_or("json")
            );
            let mut virtual_text_files = self.files.clone();
            let mut virtual_binary_files = BTreeMap::new();
            virtual_binary_files.insert(input_name.to_owned(), input_bytes.to_vec());
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
                    virtual_binary_files,
                    virtual_text_files: std::mem::take(&mut virtual_text_files),
                    preloaded_observations: preloaded_observations.to_vec(),
                },
            )
            .map_err(|err| JsError::new(&format!("create analysis runtime failed: {err:?}")))?;
            runtime
                .run_file(
                    &script_path,
                    None,
                    vec![
                        ("input_file", MontyObject::String(input_name.to_owned())),
                        ("output_file", MontyObject::String(output_file.clone())),
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
            let text = written.get(&output_file).ok_or_else(|| {
                JsError::new(&format!(
                    "analysis {} did not write {output_file}",
                    interpretation.id
                ))
            })?;
            let format = interpretation
                .output_format
                .as_deref()
                .unwrap_or("json")
                .to_ascii_lowercase();
            let (rows, row_headers) = parse_analysis_output_text(text, &format)?;
            outputs.push(serde_json::json!({
                "schema": "bioscript:analysis-output:1.0",
                "version": "1.0",
                "participant_id": participant_id,
                "assay_id": manifest_name,
                "analysis_id": interpretation.id,
                "analysis_label": interpretation.label,
                "kind": interpretation.kind,
                "output_format": format,
                "manifest_path": manifest_path,
                "script_path": script_path,
                "output_file": output_file,
                "derived_from": interpretation.derived_from,
                "emits": interpretation.emits.iter().map(|emit| serde_json::json!({
                    "key": emit.key,
                    "label": emit.label,
                    "value_type": emit.value_type,
                    "format": emit.format,
                })).collect::<Vec<_>>(),
                "logic": interpretation.logic.as_ref().map(|logic| serde_json::json!({
                    "description": logic.description,
                    "source": logic.source.as_ref().map(|source| serde_json::json!({
                        "name": source.name,
                        "url": source.url,
                    })),
                })),
                "row_headers": row_headers,
                "rows": rows,
            }));
        }
        Ok(outputs)
    }
}
