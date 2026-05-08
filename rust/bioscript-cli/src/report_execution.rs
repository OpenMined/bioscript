fn run_manifest_rows_for_report(
    runtime_root: &Path,
    manifest_path: &Path,
    input_file: &Path,
    participant_id: &str,
    loader: &GenotypeLoadOptions,
    filters: &[String],
) -> Result<Vec<BTreeMap<String, String>>, String> {
    let input_text = input_file.display().to_string();
    match manifest_schema(manifest_path)?.as_str() {
        "bioscript:variant:1.0" | "bioscript:variant" => {
            let manifest = load_variant_manifest(manifest_path)?;
            Ok(vec![run_variant_manifest(
                runtime_root,
                &manifest,
                Some(&input_text),
                Some(participant_id),
                loader,
            )?])
        }
        "bioscript:panel:1.0" => {
            let manifest = load_panel_manifest(manifest_path)?;
            run_panel_manifest(
                runtime_root,
                &manifest,
                Some(&input_text),
                Some(participant_id),
                loader,
                filters,
            )
        }
        "bioscript:assay:1.0" => {
            let manifest = load_assay_manifest(manifest_path)?;
            run_assay_manifest(
                runtime_root,
                &manifest,
                Some(&input_text),
                Some(participant_id),
                loader,
                filters,
            )
        }
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

struct ReportAnalysisOptions<'a> {
    runtime_root: &'a Path,
    input_file: &'a Path,
    participant_id: &'a str,
    loader: &'a GenotypeLoadOptions,
    output_dir: &'a Path,
    filters: &'a [String],
    max_duration_ms: u64,
}

fn run_manifest_analyses_for_report(
    manifest_path: &Path,
    options: &ReportAnalysisOptions<'_>,
) -> Result<Vec<serde_json::Value>, String> {
    match manifest_schema(manifest_path)?.as_str() {
        "bioscript:panel:1.0" => {
            let manifest = load_panel_manifest(manifest_path)?;
            let mut analyses = Vec::new();
            if options.filters.is_empty() {
                analyses.extend(run_interpretations_for_report(
                    &manifest.path,
                    &manifest.name,
                    &manifest.interpretations,
                    options,
                )?);
            }
            for member in &manifest.members {
                if member.kind != "assay" {
                    continue;
                }
                let Some(path) = &member.path else {
                    continue;
                };
                let resolved =
                    resolve_manifest_path(options.runtime_root, &manifest.path, path)?;
                if !analysis_path_matches_filters(&resolved, options.filters) {
                    continue;
                }
                analyses.extend(run_manifest_analyses_for_report(&resolved, options)?);
            }
            Ok(analyses)
        }
        "bioscript:assay:1.0" => {
            let manifest = load_assay_manifest(manifest_path)?;
            run_interpretations_for_report(
                &manifest.path,
                &manifest.name,
                &manifest.interpretations,
                options,
            )
        }
        "bioscript:variant:1.0" | "bioscript:variant" => Ok(Vec::new()),
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

fn analysis_path_matches_filters(path: &Path, filters: &[String]) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("path", value)) => path.display().to_string().contains(value),
        _ => false,
    })
}

fn run_interpretations_for_report(
    manifest_path: &Path,
    manifest_name: &str,
    interpretations: &[PanelInterpretation],
    options: &ReportAnalysisOptions<'_>,
) -> Result<Vec<serde_json::Value>, String> {
    let mut outputs = Vec::new();
    for interpretation in interpretations {
        if interpretation.kind != "bioscript" {
            return Err(format!(
                "analysis '{}' uses unsupported kind '{}'",
                interpretation.id, interpretation.kind
            ));
        }
        let script_path =
            resolve_manifest_path(options.runtime_root, manifest_path, &interpretation.path)?;
        let format = interpretation
            .output_format
            .as_deref()
            .unwrap_or("json")
            .to_ascii_lowercase();
        let analysis_dir = options.output_dir.join("analysis").join(options.participant_id);
        fs::create_dir_all(&analysis_dir).map_err(|err| {
            format!(
                "failed to create analysis output dir {}: {err}",
                analysis_dir.display()
            )
        })?;
        let extension = match format.as_str() {
            "tsv" => "tsv",
            "json" => "json",
            "jsonl" => "jsonl",
            other => return Err(format!("unsupported analysis output_format '{other}'")),
        };
        let output_file = analysis_dir.join(format!("{}.{}", interpretation.id, extension));
        run_bioscript_analysis_script(
            options.runtime_root,
            &script_path,
            options.input_file,
            &output_file,
            options.participant_id,
            options.loader,
            options.max_duration_ms,
        )?;
        let (rows, row_headers) = parse_analysis_output(&output_file, &format)?;
        outputs.push(serde_json::json!({
            "schema": "bioscript:analysis-output:1.0",
            "version": "1.0",
            "participant_id": options.participant_id,
            "assay_id": manifest_name,
            "analysis_id": interpretation.id,
            "analysis_label": interpretation.label.clone(),
            "kind": interpretation.kind,
            "output_format": format,
            "manifest_path": manifest_path.strip_prefix(options.runtime_root).unwrap_or(manifest_path).display().to_string(),
            "script_path": script_path.strip_prefix(options.runtime_root).unwrap_or(&script_path).display().to_string(),
            "output_file": output_file.strip_prefix(options.runtime_root).unwrap_or(&output_file).display().to_string(),
            "derived_from": interpretation.derived_from.clone(),
            "emits": interpretation.emits.iter().map(|emit| serde_json::json!({
                "key": emit.key.clone(),
                "label": emit.label.clone(),
                "value_type": emit.value_type.clone(),
                "format": emit.format.clone(),
            })).collect::<Vec<_>>(),
            "logic": interpretation.logic.as_ref().map(|logic| serde_json::json!({
                "description": logic.description.clone(),
                "source": logic.source.as_ref().map(|source| serde_json::json!({
                    "name": source.name.clone(),
                    "url": source.url.clone(),
                })),
            })),
            "row_headers": row_headers,
            "rows": rows,
        }));
    }
    Ok(outputs)
}

fn run_bioscript_analysis_script(
    runtime_root: &Path,
    script_path: &Path,
    input_file: &Path,
    output_file: &Path,
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
    match format {
        "tsv" => Ok(parse_analysis_tsv(&text)),
        "json" => {
            let value: serde_json::Value = serde_json::from_str(&text).map_err(|err| {
                format!("failed to parse analysis JSON {}: {err}", path.display())
            })?;
            let rows = match value {
                serde_json::Value::Array(rows) => rows,
                serde_json::Value::Object(mut object) => object
                    .remove("rows")
                    .and_then(|rows| rows.as_array().cloned())
                    .unwrap_or_else(|| vec![serde_json::Value::Object(object)]),
                other => vec![other],
            };
            let headers = analysis_headers_from_rows(&rows);
            Ok((rows, headers))
        }
        "jsonl" => {
            let rows = text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| serde_json::from_str(line).map_err(|err| err.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
            let headers = analysis_headers_from_rows(&rows);
            Ok((rows, headers))
        }
        other => Err(format!("unsupported analysis output_format '{other}'")),
    }
}

fn parse_analysis_tsv(text: &str) -> (Vec<serde_json::Value>, Vec<String>) {
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let Some(header_line) = lines.next() else {
        return (Vec::new(), Vec::new());
    };
    let headers: Vec<&str> = header_line.split('\t').collect();
    let mut rows = Vec::new();
    for line in lines {
        let values: Vec<&str> = line.split('\t').collect();
        let mut object = serde_json::Map::new();
        for (idx, header) in headers.iter().enumerate() {
            object.insert(
                (*header).to_owned(),
                serde_json::Value::String(values.get(idx).copied().unwrap_or_default().to_owned()),
            );
        }
        rows.push(serde_json::Value::Object(object));
    }
    (rows, headers.iter().map(|header| (*header).to_owned()).collect())
}

fn analysis_headers_from_rows(rows: &[serde_json::Value]) -> Vec<String> {
    let mut headers = Vec::new();
    for row in rows {
        let Some(object) = row.as_object() else {
            continue;
        };
        for key in object.keys() {
            if !headers.contains(key) {
                headers.push(key.clone());
            }
        }
    }
    headers
}

fn app_assay_id(path: &Path) -> Result<String, String> {
    match manifest_schema(path)?.as_str() {
        "bioscript:panel:1.0" => Ok(load_panel_manifest(path)?.name),
        "bioscript:assay:1.0" => Ok(load_assay_manifest(path)?.name),
        "bioscript:variant:1.0" | "bioscript:variant" => Ok(load_variant_manifest(path)?.name),
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

fn participant_id_from_path(path: &Path) -> String {
    let file_name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("participant");
    file_name
        .trim_end_matches(".txt.zip")
        .trim_end_matches(".csv.zip")
        .trim_end_matches(".vcf.gz")
        .trim_end_matches(".cram")
        .trim_end_matches(".zip")
        .trim_end_matches(".txt")
        .trim_end_matches(".csv")
        .to_owned()
}
