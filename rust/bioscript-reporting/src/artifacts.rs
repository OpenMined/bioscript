pub struct ReportArtifactTexts {
    pub observations_tsv: String,
    pub analysis_jsonl: String,
    pub reports_jsonl: String,
    pub html: String,
    pub text_output: String,
}

pub fn render_input_report_artifact_texts(
    input: crate::AppInputReportInput<'_>,
) -> Result<ReportArtifactTexts, String> {
    let report = crate::app_input_report_json(input);
    render_report_artifact_texts(input.observations, input.analyses, &[report])
}

pub fn render_report_artifact_texts(
    observations: &[serde_json::Value],
    analyses: &[serde_json::Value],
    reports: &[serde_json::Value],
) -> Result<ReportArtifactTexts, String> {
    Ok(ReportArtifactTexts {
        observations_tsv: render_observations_tsv(observations),
        analysis_jsonl: render_jsonl(analyses)?,
        reports_jsonl: render_jsonl(reports)?,
        html: crate::render_app_html_document(observations, reports)?,
        text_output: standard_text_output(),
    })
}

pub fn standard_text_output() -> String {
    "observations: observations.tsv\nanalysis: analysis.jsonl\nreports: reports.jsonl\nhtml: index.html\n"
        .to_owned()
}

pub fn render_observations_tsv(observations: &[serde_json::Value]) -> String {
    let mut out = bioscript_core::OBSERVATION_TSV_HEADERS.join("\t");
    out.push('\n');
    for observation in observations {
        let line = bioscript_core::OBSERVATION_TSV_HEADERS
            .iter()
            .map(|header| json_field_as_tsv(observation.get(*header)))
            .collect::<Vec<_>>()
            .join("\t");
        out.push_str(&line);
        out.push('\n');
    }
    out
}

pub fn render_jsonl(rows: &[serde_json::Value]) -> Result<String, String> {
    let mut out = String::new();
    for row in rows {
        out.push_str(&serde_json::to_string(row).map_err(|err| err.to_string())?);
        out.push('\n');
    }
    Ok(out)
}

pub fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}
