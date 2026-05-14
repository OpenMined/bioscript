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
    manifest_path: &'a Path,
    interpretation: &'a PanelInterpretation,
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
    let virtual_input_file = "/input/genotypes".to_owned();
    let virtual_observations_file = "/work/observations.tsv".to_owned();
    let output_extension = input
        .output_file
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or("tsv");
    let virtual_output_file = format!("/output/results.{output_extension}");
    let script_virtual_path = virtual_pipeline_path(input.script_path, "analysis.py");
    let manifest_virtual_path = virtual_pipeline_path(input.manifest_path, "manifest.yaml");
    let mut virtual_text_files = BTreeMap::new();
    virtual_text_files.insert(script_virtual_path.clone(), script_text);
    virtual_text_files.insert(
        manifest_virtual_path.clone(),
        fs::read_to_string(input.manifest_path).map_err(|err| {
            format!(
                "failed to read analysis manifest {}: {err}",
                input.manifest_path.display()
            )
        })?,
    );
    virtual_text_files.insert(virtual_observations_file.clone(), observations_text);
    let mut asset_paths = BTreeMap::new();
    for asset in &input.interpretation.assets {
        let asset_path = resolve_manifest_path(input.runtime_root, input.manifest_path, &asset.path)?;
        let virtual_asset_path = virtual_pipeline_path(&asset_path, &asset.path);
        let text = fs::read_to_string(&asset_path)
            .map_err(|err| format!("failed to read analysis asset {}: {err}", asset_path.display()))?;
        virtual_text_files.insert(virtual_asset_path.clone(), text);
        asset_paths.insert(asset.id.clone(), virtual_asset_path);
    }
    let mut virtual_binary_files = BTreeMap::new();
    virtual_binary_files.insert(virtual_input_file.clone(), input_bytes);
    let context = analysis_context(
        input.participant_id,
        &virtual_input_file,
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
            loader: input.loader.clone(),
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
    persist_virtual_output_files(&written, &virtual_output_file, input.output_file)
}

fn persist_virtual_output_files(
    written: &BTreeMap<String, String>,
    primary_virtual_output_file: &str,
    primary_output_file: &Path,
) -> Result<(), String> {
    let Some(output_dir) = primary_output_file.parent() else {
        return Ok(());
    };
    for (virtual_path, text) in written {
        if virtual_path == primary_virtual_output_file {
            continue;
        }
        let Some(relative) = virtual_path.strip_prefix("/output/") else {
            continue;
        };
        let relative_path = Path::new(relative);
        if relative_path
            .components()
            .any(|component| !matches!(component, std::path::Component::Normal(_)))
        {
            return Err(format!("analysis wrote invalid output path {virtual_path}"));
        }
        let output_path = output_dir.join(relative_path);
        if let Some(parent) = output_path.parent() {
            fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create output dir {}: {err}", parent.display()))?;
        }
        fs::write(&output_path, text)
            .map_err(|err| format!("failed to write analysis output {}: {err}", output_path.display()))?;
    }
    Ok(())
}

fn virtual_pipeline_path(path: &Path, fallback: &str) -> String {
    let name = path
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
) -> BTreeMap<String, monty::MontyObject> {
    BTreeMap::from([
        (
            "participant_id".to_owned(),
            monty::MontyObject::String(participant_id.to_owned()),
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
            monty::MontyObject::String(observations_file.to_owned()),
        ),
        (
            "output_file".to_owned(),
            monty::MontyObject::String(output_file.to_owned()),
        ),
    ])
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
