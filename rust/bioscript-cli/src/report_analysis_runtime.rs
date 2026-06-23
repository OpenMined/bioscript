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
    let observations_text = fs::read_to_string(input.observations_file).map_err(|err| {
        format!(
            "failed to read analysis observations {}: {err}",
            input.observations_file.display()
        )
    })?;
    let script_text = fs::read_to_string(input.script_path).map_err(|err| {
        format!(
            "failed to read analysis script {}: {err}",
            input.script_path.display()
        )
    })?;
    let virtual_input_file = "/input/genotypes".to_owned();
    let virtual_observations_file = "/work/observations.tsv".to_owned();
    let output_extension = input
        .output_file
        .extension()
        .and_then(|value| value.to_str())
        .unwrap_or("tsv");
    let virtual_output_file = format!("/output/results.{output_extension}");
    let (analysis_input_file, virtual_input_bytes) =
        analysis_input_file_arg(input.runtime_root, input.input_file)?;
    let script_virtual_path = virtual_pipeline_path(input.script_path, "analysis.py");
    let manifest_virtual_path = virtual_pipeline_path(input.manifest_path, "manifest.yaml");
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
    if let Some(input_bytes) = virtual_input_bytes {
        virtual_binary_files.insert(virtual_input_file.clone(), input_bytes);
    }
    let context = analysis_context(
        input.participant_id,
        &analysis_input_file,
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
                ("input_file", monty::MontyObject::String(analysis_input_file)),
                (
                    "output_file",
                    monty::MontyObject::String(virtual_output_file.clone()),
                ),
                (
                    "observations_file",
                    monty::MontyObject::String(virtual_observations_file),
                ),
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
        fs::create_dir_all(parent).map_err(|err| {
            format!(
                "failed to create analysis output dir {}: {err}",
                parent.display()
            )
        })?;
    }
    fs::write(input.output_file, output_text).map_err(|err| {
        format!(
            "failed to write analysis output {}: {err}",
            input.output_file.display()
        )
    })?;
    persist_virtual_output_files(&written, &virtual_output_file, input.output_file)
}

fn analysis_input_file_arg(
    runtime_root: &Path,
    input_file: &Path,
) -> Result<(String, Option<Vec<u8>>), String> {
    let canonical_root = runtime_root.canonicalize().map_err(|err| {
        format!(
            "failed to resolve runtime root {}: {err}",
            runtime_root.display()
        )
    })?;
    let canonical_input = input_file.canonicalize().map_err(|err| {
        format!(
            "failed to resolve analysis input {}: {err}",
            input_file.display()
        )
    })?;
    if let Ok(relative) = canonical_input.strip_prefix(&canonical_root) {
        let relative_text = relative.display().to_string();
        if !relative_text.is_empty() {
            return Ok((relative_text, None));
        }
    }

    let input_bytes = fs::read(input_file).map_err(|err| {
        format!(
            "failed to read analysis input {}: {err}",
            input_file.display()
        )
    })?;
    Ok(("/input/genotypes".to_owned(), Some(input_bytes)))
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
        let asset_path =
            resolve_manifest_path(input.runtime_root, input.manifest_path, &asset.path)?;
        let virtual_asset_path = virtual_pipeline_path(&asset_path, &asset.path);
        let text = fs::read_to_string(&asset_path).map_err(|err| {
            format!(
                "failed to read analysis asset {}: {err}",
                asset_path.display()
            )
        })?;
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
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create output dir {}: {err}", parent.display())
            })?;
        }
        fs::write(&output_path, text).map_err(|err| {
            format!(
                "failed to write analysis output {}: {err}",
                output_path.display()
            )
        })?;
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
