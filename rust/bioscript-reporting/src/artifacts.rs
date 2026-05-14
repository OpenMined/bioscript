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

#[cfg(test)]
mod tests {
    use serde_json::json;

    use crate::AppInputReportInput;

    use super::{
        json_field_as_tsv, render_input_report_artifact_texts, render_jsonl,
        render_observations_tsv, render_report_artifact_texts, standard_text_output,
    };

    fn observation() -> serde_json::Value {
        json!({
            "participant_id": "P001",
            "assay_id": "assay",
            "assay_version": "1.0",
            "variant_key": "rs123",
            "rsid": "rs123",
            "call_status": "called",
            "outcome": "variant"
        })
    }

    #[test]
    fn artifact_renderers_emit_standard_files() {
        let observations = vec![observation()];
        let analyses = vec![json!({"analysis_id": "a1", "rows": []})];
        let reports = vec![json!({
            "schema": "bioscript:report:1.0",
            "participant_id": "P001",
            "manifest": {"name": "panel"},
            "input": {"file_name": "sample.txt"},
            "analyses": analyses,
            "findings": [],
            "provenance": []
        })];

        let artifacts = render_report_artifact_texts(&observations, &analyses, &reports).unwrap();

        assert!(artifacts.observations_tsv.starts_with("participant_id\t"));
        assert!(artifacts.analysis_jsonl.contains("\"analysis_id\":\"a1\""));
        assert!(artifacts.reports_jsonl.contains("bioscript:report:1.0"));
        assert!(artifacts.html.contains("BioScript"));
        assert_eq!(artifacts.text_output, standard_text_output());
    }

    #[test]
    fn input_report_artifact_renderer_builds_report_json() {
        let observations = vec![observation()];
        let analyses = Vec::new();
        let findings = Vec::new();
        let provenance = Vec::new();
        let manifest = json!({"name": "panel"});

        let artifacts = render_input_report_artifact_texts(AppInputReportInput {
            assay_id: "assay",
            participant_id: "P001",
            input_file_name: "sample.txt",
            input_file_path: "/data/sample.txt",
            observations: &observations,
            analyses: &analyses,
            findings: &findings,
            provenance: &provenance,
            input_inspection: None,
            manifest_metadata: &manifest,
        })
        .unwrap();

        assert!(artifacts.reports_jsonl.contains("\"report_status\":\"complete\""));
        assert!(artifacts.html.contains("panel"));
    }

    #[test]
    fn low_level_serializers_escape_tsv_and_jsonl_rows() {
        let observations = vec![json!({
            "participant_id": "P\t001",
            "assay_id": "assay",
            "variant_key": "rs123",
            "call_status": "called",
            "facets": {"note": "a\nb"}
        })];

        let tsv = render_observations_tsv(&observations);
        assert!(tsv.contains("P 001"));
        assert!(!tsv.contains("P\t001"));
        assert_eq!(
            json_field_as_tsv(Some(&json!("line\tbreak\nvalue"))),
            "line break value"
        );
        assert_eq!(json_field_as_tsv(Some(&serde_json::Value::Null)), "");
        assert_eq!(json_field_as_tsv(None), "");

        let jsonl = render_jsonl(&[json!({"a": 1}), json!({"b": true})]).unwrap();
        assert_eq!(jsonl.lines().count(), 2);
    }
}
