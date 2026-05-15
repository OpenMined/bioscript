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

#[cfg(test)]
mod app_report_output_tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-report-output-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn writes_observation_outputs_for_each_format() {
        let rows = vec![serde_json::json!({
            "participant_id": "p1",
            "gene": "CYP2D6",
            "genotype": "A/G"
        })];

        let dir = temp_dir("observations");
        write_app_observations(&dir, &rows, AppOutputFormat::Both).unwrap();
        let tsv = fs::read_to_string(dir.join("observations.tsv")).unwrap();
        assert!(tsv.contains("participant_id"));
        assert!(tsv.contains("p1"));
        assert!(fs::read_to_string(dir.join("observations.jsonl"))
            .unwrap()
            .contains("\"participant_id\":\"p1\""));
        assert!(!dir.join("observations.json").exists());

        write_app_observations(&dir, &rows, AppOutputFormat::Json).unwrap();
        let json: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(dir.join("observations.json")).unwrap())
                .unwrap();
        assert_eq!(json["observations"][0]["genotype"], "A/G");

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn writes_analysis_reports_and_html_outputs() {
        let dir = temp_dir("reports");
        let analyses = vec![serde_json::json!({"id": "analysis-1"})];
        let reports = vec![serde_json::json!({
            "participant": {"id": "p1"},
            "observations": []
        })];
        let observations = vec![serde_json::json!({"participant_id": "p1"})];

        write_app_analyses(&dir, &analyses).unwrap();
        assert!(fs::read_to_string(dir.join("analysis.jsonl"))
            .unwrap()
            .contains("analysis-1"));

        write_app_reports(&dir, &reports, AppOutputFormat::Both).unwrap();
        assert!(fs::read_to_string(dir.join("reports.jsonl"))
            .unwrap()
            .contains("\"participant\""));
        let report_set: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(dir.join("reports.json")).unwrap()).unwrap();
        assert_eq!(report_set["schema"], "bioscript:report-set:1.0");

        write_app_html(&dir, &observations, &reports).unwrap();
        assert!(fs::read_to_string(dir.join("index.html"))
            .unwrap()
            .contains("<!doctype html>"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn write_helpers_report_filesystem_errors() {
        let missing_parent = env::temp_dir()
            .join(format!("bioscript-missing-parent-{}", std::process::id()))
            .join("out.jsonl");
        let err = write_jsonl(&missing_parent, &[]).unwrap_err();
        assert!(err.contains("failed to write"));
    }
}
