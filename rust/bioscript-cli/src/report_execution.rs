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
        run_bioscript_analysis_script(
            options.runtime_root,
            &script_path,
            options.input_file,
            &output_file,
            &observations_file,
            options.participant_id,
            options.loader,
            options.max_duration_ms,
        )?;
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

fn run_bioscript_analysis_script(
    runtime_root: &Path,
    script_path: &Path,
    input_file: &Path,
    output_file: &Path,
    observations_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
    analysis_max_duration_ms: u64,
) -> Result<(), String> {
    let limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(analysis_max_duration_ms))
        .max_memory(16 * 1024 * 1024)
        .max_allocations(400_000)
        .gc_interval(1000)
        .max_recursion_depth(Some(200));
    let runtime = BioscriptRuntime::with_config(
        runtime_root.to_path_buf(),
        RuntimeConfig {
            limits,
            loader: loader.clone(),
            ..RuntimeConfig::default()
        },
    )
    .map_err(|err| err.to_string())?;
    runtime
        .run_file(
            script_path,
            None,
            vec![
                (
                    "input_file",
                    monty::MontyObject::String(runtime_path_string(runtime_root, input_file)),
                ),
                (
                    "output_file",
                    monty::MontyObject::String(runtime_path_string(runtime_root, output_file)),
                ),
                (
                    "observations_file",
                    monty::MontyObject::String(runtime_path_string(
                        runtime_root,
                        observations_file,
                    )),
                ),
                (
                    "participant_id",
                    monty::MontyObject::String(participant_id.to_owned()),
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
