struct ManifestRunOptions<'a> {
    input_file: Option<&'a str>,
    output_file: Option<&'a str>,
    participant_id: Option<&'a str>,
    trace_report: Option<&'a Path>,
    loader: &'a GenotypeLoadOptions,
    filters: &'a [String],
}

fn run_manifest(
    runtime_root: &Path,
    manifest_path: &Path,
    options: &ManifestRunOptions<'_>,
) -> Result<(), String> {
    let workspace = bioscript_reporting::FilesystemManifestWorkspace::new(runtime_root);
    let schema = bioscript_reporting::report_manifest_schema(
        &workspace,
        &manifest_path.display().to_string(),
    )?;
    let resolved_input = options
        .input_file
        .map(|value| resolve_cli_path(runtime_root, value));
    let resolved_output = options
        .output_file
        .map(|value| resolve_cli_path_buf(runtime_root, Path::new(value)));
    let resolved_trace = options
        .trace_report
        .map(|value| resolve_cli_path_buf(runtime_root, value));
    match bioscript_reporting::report_manifest_kind(&schema)? {
        bioscript_reporting::ReportManifestKind::Variant => {
            let manifest = load_variant_manifest(manifest_path)?;
            let row = run_variant_manifest(
                runtime_root,
                &manifest,
                resolved_input.as_deref(),
                options.participant_id,
                options.loader,
            )?;
            write_manifest_outputs(
                std::slice::from_ref(&row),
                resolved_output.as_deref(),
                resolved_trace.as_deref(),
            )?;
            Ok(())
        }
        bioscript_reporting::ReportManifestKind::Panel => {
            let manifest = load_panel_manifest(manifest_path)?;
            let rows = run_panel_manifest(
                runtime_root,
                &manifest,
                resolved_input.as_deref(),
                options.participant_id,
                options.loader,
                options.filters,
            )?;
            write_manifest_outputs(&rows, resolved_output.as_deref(), resolved_trace.as_deref())?;
            Ok(())
        }
        bioscript_reporting::ReportManifestKind::Assay => {
            let manifest = load_assay_manifest(manifest_path)?;
            let rows = run_assay_manifest(
                runtime_root,
                &manifest,
                resolved_input.as_deref(),
                options.participant_id,
                options.loader,
                options.filters,
            )?;
            write_manifest_outputs(&rows, resolved_output.as_deref(), resolved_trace.as_deref())?;
            Ok(())
        }
    }
}

fn run_variant_manifest(
    runtime_root: &Path,
    manifest: &VariantManifest,
    input_file: Option<&str>,
    participant_id: Option<&str>,
    loader: &GenotypeLoadOptions,
) -> Result<BTreeMap<String, String>, String> {
    let input_file = input_file.ok_or("manifest execution requires --input-file")?;
    let store = GenotypeStore::from_file_with_options(Path::new(input_file), loader)
        .map_err(|err| err.to_string())?;
    run_variant_manifest_with_store(runtime_root, manifest, &store, participant_id)
}

fn run_variant_manifest_with_store(
    runtime_root: &Path,
    manifest: &VariantManifest,
    store: &GenotypeStore,
    participant_id: Option<&str>,
) -> Result<BTreeMap<String, String>, String> {
    let observation = store
        .lookup_variant(&manifest.spec)
        .map_err(|err| err.to_string())?;
    Ok(variant_row(
        runtime_root,
        &manifest.path,
        &manifest.name,
        &manifest.tags,
        &observation,
        participant_id,
    ))
}

fn run_panel_manifest(
    runtime_root: &Path,
    panel: &PanelManifest,
    input_file: Option<&str>,
    participant_id: Option<&str>,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_file = input_file.ok_or("manifest execution requires --input-file")?;
    let store = GenotypeStore::from_file_with_options(Path::new(input_file), loader)
        .map_err(|err| err.to_string())?;
    run_panel_manifest_with_store(runtime_root, panel, &store, participant_id, filters)
}

fn run_panel_manifest_with_store(
    runtime_root: &Path,
    panel: &PanelManifest,
    store: &GenotypeStore,
    participant_id: Option<&str>,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let mut rows_by_member: Vec<Vec<BTreeMap<String, String>>> = vec![Vec::new(); panel.members.len()];
    let mut variant_entries = Vec::new();

    for (member_index, member) in panel.members.iter().enumerate() {
        match bioscript_reporting::panel_executable_member(&member.kind, member.path.as_deref())? {
            bioscript_reporting::ExecutablePanelMember::Variant(path) => {
                let resolved = resolve_manifest_path(runtime_root, &panel.path, path)?;
                let manifest = load_variant_manifest(&resolved)?;
                if !matches_filters(&manifest, &resolved, filters) {
                    continue;
                }
                variant_entries.push((member_index, resolved, manifest));
            }
            bioscript_reporting::ExecutablePanelMember::Assay(path) => {
                let resolved = resolve_manifest_path(runtime_root, &panel.path, path)?;
                let assay = load_assay_manifest(&resolved)?;
                rows_by_member[member_index] = run_assay_manifest_with_store(
                    runtime_root,
                    &assay,
                    store,
                    participant_id,
                    filters,
                )?;
            }
        }
    }

    let observations = store
        .lookup_variants(
            &variant_entries
                .iter()
                .map(|(_, _, manifest)| manifest.spec.clone())
                .collect::<Vec<_>>(),
        )
        .map_err(|err| err.to_string())?;

    for ((member_index, resolved, manifest), observation) in
        variant_entries.into_iter().zip(observations)
    {
        rows_by_member[member_index].push(variant_row(
            runtime_root,
            &resolved,
            &manifest.name,
            &manifest.tags,
            &observation,
            participant_id,
        ));
    }

    let rows = rows_by_member.into_iter().flatten().collect();
    Ok(rows)
}

fn run_assay_manifest(
    runtime_root: &Path,
    assay: &AssayManifest,
    input_file: Option<&str>,
    participant_id: Option<&str>,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_file = input_file.ok_or("manifest execution requires --input-file")?;
    let store = GenotypeStore::from_file_with_options(Path::new(input_file), loader)
        .map_err(|err| err.to_string())?;
    run_assay_manifest_with_store(runtime_root, assay, &store, participant_id, filters)
}

fn run_assay_manifest_with_store(
    runtime_root: &Path,
    assay: &AssayManifest,
    store: &GenotypeStore,
    participant_id: Option<&str>,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let mut entries = Vec::new();

    for member in &assay.members {
        match bioscript_reporting::assay_executable_member(&member.kind, member.path.as_deref())? {
            bioscript_reporting::ExecutableAssayMember::Variant(path) => {
                let resolved = resolve_manifest_path(runtime_root, &assay.path, path)?;
                let manifest = load_variant_manifest(&resolved)?;
                if !matches_filters(&manifest, &resolved, filters) {
                    continue;
                }
                entries.push((resolved, manifest));
            }
        }
    }

    let observations = store
        .lookup_variants(
            &entries
                .iter()
                .map(|(_, manifest)| manifest.spec.clone())
                .collect::<Vec<_>>(),
        )
        .map_err(|err| err.to_string())?;

    let rows = entries
        .into_iter()
        .zip(observations)
        .map(|((resolved, manifest), observation)| {
            variant_row(
            runtime_root,
            &resolved,
            &manifest.name,
            &manifest.tags,
            &observation,
            participant_id,
            )
        })
        .collect();

    Ok(rows)
}

fn variant_row(
    runtime_root: &Path,
    path: &Path,
    name: &str,
    tags: &[String],
    observation: &bioscript_core::VariantObservation,
    participant_id: Option<&str>,
) -> BTreeMap<String, String> {
    let row_path = path
        .strip_prefix(runtime_root)
        .unwrap_or(path)
        .display()
        .to_string();
    bioscript_reporting::variant_row(
        &row_path,
        name,
        tags,
        observation,
        participant_id.unwrap_or_default(),
    )
}

fn write_manifest_outputs(
    rows: &[BTreeMap<String, String>],
    output_file: Option<&Path>,
    trace_report: Option<&Path>,
) -> Result<(), String> {
    let text = bioscript_reporting::render_manifest_rows_tsv(rows);
    if let Some(output_file) = output_file {
        if let Some(parent) = output_file.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                format!("failed to create output dir {}: {err}", parent.display())
            })?;
        }
        fs::write(output_file, &text)
            .map_err(|err| format!("failed to write output {}: {err}", output_file.display()))?;
    } else {
        print!("{text}");
    }

    if let Some(trace_report) = trace_report {
        if let Some(parent) = trace_report.parent() {
            fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create trace dir {}: {err}", parent.display()))?;
        }
        let trace = bioscript_reporting::render_manifest_trace_tsv(rows);
        fs::write(trace_report, trace)
            .map_err(|err| format!("failed to write trace {}: {err}", trace_report.display()))?;
    }

    Ok(())
}

fn resolve_cli_path(root: &Path, value: &str) -> String {
    resolve_cli_path_buf(root, Path::new(value))
        .display()
        .to_string()
}

fn resolve_cli_path_buf(root: &Path, value: &Path) -> PathBuf {
    if value.is_absolute() {
        value.to_path_buf()
    } else {
        root.join(value)
    }
}

fn matches_filters(manifest: &VariantManifest, path: &Path, filters: &[String]) -> bool {
    bioscript_reporting::matches_variant_manifest_filters(
        manifest,
        &path.display().to_string(),
        filters,
    )
}

fn resolve_manifest_path(
    runtime_root: &Path,
    manifest_path: &Path,
    relative: &str,
) -> Result<PathBuf, String> {
    bioscript_reporting::resolve_filesystem_manifest_path(runtime_root, manifest_path, relative)
}

fn manifest_schema(path: &Path) -> Result<String, String> {
    let root = path.parent().unwrap_or_else(|| Path::new("."));
    let workspace = bioscript_reporting::FilesystemManifestWorkspace::new(root);
    bioscript_reporting::report_manifest_schema(&workspace, &path.display().to_string())
}

fn normalize_loader_paths(root: &Path, loader: &mut GenotypeLoadOptions) {
    if let Some(path) = loader.input_index.take() {
        loader.input_index = Some(resolve_cli_path_buf(root, &path));
    }
    if let Some(path) = loader.reference_file.take() {
        loader.reference_file = Some(resolve_cli_path_buf(root, &path));
    }
    if let Some(path) = loader.reference_index.take() {
        loader.reference_index = Some(resolve_cli_path_buf(root, &path));
    }
}
