fn write_app_observations(
    output_dir: &Path,
    observations: &[serde_json::Value],
    format: AppOutputFormat,
) -> Result<(), String> {
    if matches!(format, AppOutputFormat::Tsv | AppOutputFormat::Both) {
        let out = bioscript_reporting::render_observations_tsv(observations);
        fs::write(output_dir.join("observations.tsv"), out)
            .map_err(|err| format!("failed to write observations.tsv: {err}"))?;
    }
    if matches!(format, AppOutputFormat::Jsonl | AppOutputFormat::Both) {
        write_jsonl(&output_dir.join("observations.jsonl"), observations)?;
    }
    if matches!(format, AppOutputFormat::Json) {
        write_json_pretty(
            &output_dir.join("observations.json"),
            &serde_json::json!({"observations": observations}),
        )?;
    }
    Ok(())
}

fn write_app_analyses(output_dir: &Path, analyses: &[serde_json::Value]) -> Result<(), String> {
    write_jsonl(&output_dir.join("analysis.jsonl"), analyses)
}

fn write_app_reports(
    output_dir: &Path,
    reports: &[serde_json::Value],
    format: AppOutputFormat,
) -> Result<(), String> {
    if matches!(format, AppOutputFormat::Jsonl | AppOutputFormat::Both) {
        write_jsonl(&output_dir.join("reports.jsonl"), reports)?;
    }
    if matches!(format, AppOutputFormat::Json | AppOutputFormat::Both) {
        write_json_pretty(
            &output_dir.join("reports.json"),
            &serde_json::json!({
                "schema": "bioscript:report-set:1.0",
                "version": "1.0",
                "reports": reports,
            }),
        )?;
    }
    Ok(())
}

fn write_jsonl(path: &Path, rows: &[serde_json::Value]) -> Result<(), String> {
    let out = bioscript_reporting::render_jsonl(rows)?;
    fs::write(path, out).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn write_json_pretty(path: &Path, value: &serde_json::Value) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value).map_err(|err| err.to_string())?;
    fs::write(path, text).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn write_app_html(
    output_dir: &Path,
    observations: &[serde_json::Value],
    reports: &[serde_json::Value],
) -> Result<(), String> {
    let out = bioscript_reporting::render_app_html_document(observations, reports)?;
    fs::write(output_dir.join("index.html"), out)
        .map_err(|err| format!("failed to write index.html: {err}"))
}
