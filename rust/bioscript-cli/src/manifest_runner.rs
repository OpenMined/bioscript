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
    let schema = manifest_schema(manifest_path)?;
    let resolved_input = options
        .input_file
        .map(|value| resolve_cli_path(runtime_root, value));
    let resolved_output = options
        .output_file
        .map(|value| resolve_cli_path_buf(runtime_root, Path::new(value)));
    let resolved_trace = options
        .trace_report
        .map(|value| resolve_cli_path_buf(runtime_root, value));
    match schema.as_str() {
        "bioscript:variant:1.0" | "bioscript:variant" => {
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
        "bioscript:panel:1.0" => {
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
        "bioscript:assay:1.0" => {
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
        other => Err(format!("unsupported manifest schema '{other}'")),
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
        let Some(path) = &member.path else {
            return Err("remote panel members are not executable yet".to_owned());
        };
        let resolved = resolve_manifest_path(runtime_root, &panel.path, path)?;
        if member.kind == "variant" {
            let manifest = load_variant_manifest(&resolved)?;
            if !matches_filters(&manifest, &resolved, filters) {
                continue;
            }
            variant_entries.push((member_index, resolved, manifest));
        } else if member.kind == "assay" {
            let assay = load_assay_manifest(&resolved)?;
            rows_by_member[member_index] = run_assay_manifest_with_store(
                runtime_root,
                &assay,
                store,
                participant_id,
                filters,
            )?;
        } else {
            return Err(format!(
                "panel member kind '{}' is not executable",
                member.kind
            ));
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
        if member.kind != "variant" {
            return Err(format!(
                "assay member kind '{}' is not executable",
                member.kind
            ));
        }
        let Some(path) = &member.path else {
            return Err("remote assay members are not executable yet".to_owned());
        };
        let resolved = resolve_manifest_path(runtime_root, &assay.path, path)?;
        let manifest = load_variant_manifest(&resolved)?;
        if !matches_filters(&manifest, &resolved, filters) {
            continue;
        }
        entries.push((resolved, manifest));
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
    let mut row = BTreeMap::new();
    row.insert("kind".to_owned(), "variant".to_owned());
    row.insert("name".to_owned(), name.to_owned());
    row.insert(
        "path".to_owned(),
        path.strip_prefix(runtime_root)
            .unwrap_or(path)
            .display()
            .to_string(),
    );
    row.insert("tags".to_owned(), tags.join(","));
    row.insert("backend".to_owned(), observation.backend.clone());
    row.insert(
        "participant_id".to_owned(),
        participant_id.unwrap_or_default().to_owned(),
    );
    row.insert(
        "matched_rsid".to_owned(),
        observation.matched_rsid.clone().unwrap_or_default(),
    );
    row.insert(
        "assembly".to_owned(),
        observation
            .assembly
            .map(|value| match value {
                bioscript_core::Assembly::Grch37 => "grch37".to_owned(),
                bioscript_core::Assembly::Grch38 => "grch38".to_owned(),
            })
            .unwrap_or_default(),
    );
    row.insert(
        "genotype".to_owned(),
        observation.genotype.clone().unwrap_or_default(),
    );
    row.insert(
        "ref_count".to_owned(),
        observation
            .ref_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "alt_count".to_owned(),
        observation
            .alt_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "depth".to_owned(),
        observation
            .depth
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "raw_counts".to_owned(),
        serde_json::to_string(&observation.raw_counts).unwrap_or_default(),
    );
    row.insert("evidence".to_owned(), observation.evidence.join(" | "));
    row
}

fn write_manifest_outputs(
    rows: &[BTreeMap<String, String>],
    output_file: Option<&Path>,
    trace_report: Option<&Path>,
) -> Result<(), String> {
    let text = render_rows_as_tsv(rows);
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
        let mut trace = String::from("step\tline\tcode\n");
        for (idx, row) in rows.iter().enumerate() {
            let _ = writeln!(
                trace,
                "{}\t{}\t{}",
                idx + 1,
                idx + 1,
                row.get("path").cloned().unwrap_or_default()
            );
        }
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

fn render_rows_as_tsv(rows: &[BTreeMap<String, String>]) -> String {
    let headers = [
        "kind",
        "name",
        "path",
        "tags",
        "participant_id",
        "backend",
        "matched_rsid",
        "assembly",
        "genotype",
        "ref_count",
        "alt_count",
        "depth",
        "evidence",
    ];
    let mut out = headers.join("\t");
    out.push('\n');
    for row in rows {
        let line = headers
            .iter()
            .map(|header| {
                row.get(*header)
                    .cloned()
                    .unwrap_or_default()
                    .replace('\t', " ")
            })
            .collect::<Vec<_>>()
            .join("\t");
        out.push_str(&line);
        out.push('\n');
    }
    out
}

fn matches_filters(manifest: &VariantManifest, path: &Path, filters: &[String]) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("kind", value)) => value == "variant",
        Some(("name", value)) => manifest.name.contains(value),
        Some(("path", value)) => path.display().to_string().contains(value),
        Some(("tag", value)) => manifest.tags.iter().any(|tag| tag == value),
        Some(_) | None => false,
    })
}

fn resolve_manifest_path(
    runtime_root: &Path,
    manifest_path: &Path,
    relative: &str,
) -> Result<PathBuf, String> {
    let base_dir = manifest_path
        .parent()
        .ok_or_else(|| format!("manifest has no parent: {}", manifest_path.display()))?;
    let joined = base_dir.join(relative);
    let canonical_root = runtime_root
        .canonicalize()
        .map_err(|err| format!("failed to resolve root {}: {err}", runtime_root.display()))?;
    let canonical_base = base_dir.canonicalize().map_err(|err| {
        format!(
            "failed to resolve manifest dir {}: {err}",
            base_dir.display()
        )
    })?;
    let canonical_joined = joined
        .canonicalize()
        .map_err(|err| format!("failed to resolve {}: {err}", joined.display()))?;
    let boundary = if canonical_base.starts_with(&canonical_root) {
        &canonical_root
    } else {
        &canonical_base
    };
    if !canonical_joined.starts_with(boundary) {
        return Err(format!(
            "manifest member path escapes bioscript root: {}",
            canonical_joined.display()
        ));
    }
    Ok(canonical_joined)
}

fn manifest_schema(path: &Path) -> Result<String, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))?;
    value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("schema".to_owned())))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
        .ok_or_else(|| format!("{} is missing schema", path.display()))
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
