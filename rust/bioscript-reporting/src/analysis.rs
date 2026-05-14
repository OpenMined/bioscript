use std::path::Path;

use bioscript_core::{Assembly, VariantObservation};
use bioscript_schema::PanelInterpretation;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnalysisOutputFormat {
    pub format: &'static str,
    pub extension: &'static str,
}

pub fn analysis_output_format(value: Option<&str>) -> Result<AnalysisOutputFormat, String> {
    match value.unwrap_or("json").to_ascii_lowercase().as_str() {
        "tsv" => Ok(AnalysisOutputFormat {
            format: "tsv",
            extension: "tsv",
        }),
        "json" => Ok(AnalysisOutputFormat {
            format: "json",
            extension: "json",
        }),
        "jsonl" => Ok(AnalysisOutputFormat {
            format: "jsonl",
            extension: "jsonl",
        }),
        other => Err(format!("unsupported analysis output_format '{other}'")),
    }
}

pub fn parse_analysis_output_text(
    text: &str,
    format: &str,
) -> Result<(Vec<serde_json::Value>, Vec<String>), String> {
    match format {
        "tsv" => Ok(parse_analysis_tsv(text)),
        "json" => {
            let value: serde_json::Value = serde_json::from_str(text)
                .map_err(|err| format!("failed to parse analysis JSON: {err}"))?;
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

pub fn validate_bioscript_interpretation(
    interpretation: &PanelInterpretation,
) -> Result<(), String> {
    if interpretation.kind == "bioscript" {
        return Ok(());
    }
    Err(format!(
        "analysis '{}' uses unsupported kind '{}'",
        interpretation.id, interpretation.kind
    ))
}

pub fn analysis_output_relative_file(
    participant_id: &str,
    interpretation_id: &str,
    extension: &str,
) -> String {
    format!("analysis/{participant_id}/{interpretation_id}.{extension}")
}

pub fn analysis_observations_relative_file(
    participant_id: &str,
    interpretation_id: &str,
) -> String {
    format!("analysis/{participant_id}/{interpretation_id}.observations.tsv")
}

pub fn render_analysis_observations_tsv(observations: &[VariantObservation]) -> String {
    let headers = [
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
    for observation in observations {
        let assembly = observation
            .assembly
            .map(|value| match value {
                Assembly::Grch37 => "grch37",
                Assembly::Grch38 => "grch38",
            })
            .unwrap_or_default();
        let values = [
            observation.matched_rsid.clone().unwrap_or_default(),
            assembly.to_owned(),
            observation.genotype.clone().unwrap_or_default(),
            observation
                .ref_count
                .map_or_else(String::new, |value| value.to_string()),
            observation
                .alt_count
                .map_or_else(String::new, |value| value.to_string()),
            observation
                .depth
                .map_or_else(String::new, |value| value.to_string()),
            observation.evidence.join(" | "),
        ];
        out.push_str(&values.join("\t"));
        out.push('\n');
    }
    out
}

pub struct AnalysisOutputJsonInput<'a> {
    pub participant_id: &'a str,
    pub assay_id: &'a str,
    pub interpretation: &'a PanelInterpretation,
    pub output_format: &'a str,
    pub manifest_path: &'a str,
    pub script_path: &'a str,
    pub output_file: &'a str,
    pub observations_file: Option<&'a str>,
    pub row_headers: Vec<String>,
    pub rows: Vec<serde_json::Value>,
}

pub fn analysis_output_json(input: AnalysisOutputJsonInput<'_>) -> serde_json::Value {
    let AnalysisOutputJsonInput {
        participant_id,
        assay_id,
        interpretation,
        output_format,
        manifest_path,
        script_path,
        output_file,
        observations_file,
        row_headers,
        rows,
    } = input;

    serde_json::json!({
        "schema": "bioscript:analysis-output:1.0",
        "version": "1.0",
        "participant_id": participant_id,
        "assay_id": assay_id,
        "analysis_id": interpretation.id,
        "analysis_label": interpretation.label,
        "kind": interpretation.kind,
        "output_format": output_format,
        "manifest_path": manifest_path,
        "script_path": script_path,
        "output_file": output_file,
        "observations_file": observations_file,
        "derived_from": interpretation.derived_from,
        "assets": [],
        "emits": interpretation.emits.iter().map(|emit| serde_json::json!({
            "key": emit.key,
            "label": emit.label,
            "value_type": emit.value_type,
            "format": emit.format,
        })).collect::<Vec<_>>(),
        "logic": interpretation.logic.as_ref().map(|logic| serde_json::json!({
            "description": logic.description,
            "source": logic.source.as_ref().map(|source| serde_json::json!({
                "name": source.name,
                "url": source.url,
            })),
        })),
        "row_headers": row_headers,
        "rows": rows,
    })
}

pub fn participant_id_from_path(path: &Path) -> String {
    let file_name = path
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("participant");
    participant_id_from_name(file_name)
}

pub fn participant_id_from_name(file_name: &str) -> String {
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
    (
        rows,
        headers.iter().map(|header| (*header).to_owned()).collect(),
    )
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

#[cfg(test)]
mod tests {
    use super::{
        AnalysisOutputFormat, AnalysisOutputJsonInput, analysis_observations_relative_file,
        analysis_output_format, analysis_output_json, analysis_output_relative_file,
        parse_analysis_output_text, participant_id_from_name, render_analysis_observations_tsv,
        validate_bioscript_interpretation,
    };
    use bioscript_core::{Assembly, VariantObservation};
    use bioscript_schema::{
        PanelInterpretation, PanelInterpretationLogic, PanelInterpretationLogicSource,
    };

    #[test]
    fn analysis_output_format_defaults_lowercases_and_rejects_unknown() {
        assert_eq!(
            analysis_output_format(None).unwrap(),
            AnalysisOutputFormat {
                format: "json",
                extension: "json"
            }
        );
        assert_eq!(
            analysis_output_format(Some("TSV")).unwrap(),
            AnalysisOutputFormat {
                format: "tsv",
                extension: "tsv"
            }
        );
        assert_eq!(
            analysis_output_format(Some("xml")).unwrap_err(),
            "unsupported analysis output_format 'xml'"
        );
    }

    #[test]
    fn participant_id_suffix_stripping_matches_cli_report_path() {
        assert_eq!(
            participant_id_from_name("NA06985.clean.vcf.gz"),
            "NA06985.clean"
        );
        assert_eq!(
            participant_id_from_name("genome_hu50B3F5_v5_Full.zip"),
            "genome_hu50B3F5_v5_Full"
        );
        assert_eq!(participant_id_from_name("sample.cram"), "sample");
        assert_eq!(participant_id_from_name("sample name.txt"), "sample name");
    }

    #[test]
    fn parses_tsv_analysis_output_with_missing_cells() {
        let (rows, headers) =
            parse_analysis_output_text("status\tnote\nnormal\nvariant\tflag\n", "tsv").unwrap();
        assert_eq!(headers, vec!["status", "note"]);
        assert_eq!(rows[0]["status"], "normal");
        assert_eq!(rows[0]["note"], "");
        assert_eq!(rows[1]["status"], "variant");
        assert_eq!(rows[1]["note"], "flag");
    }

    #[test]
    fn json_headers_include_keys_from_all_rows() {
        let (rows, headers) =
            parse_analysis_output_text(r#"[{"a":1},{"b":2,"a":3}]"#, "json").unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(headers, vec!["a", "b"]);
    }

    #[test]
    fn interpretation_validation_and_output_file_naming_are_shared() {
        let mut interpretation = PanelInterpretation {
            id: "apoe_epsilon".to_owned(),
            label: None,
            kind: "bioscript".to_owned(),
            path: "apoe.py".to_owned(),
            output_format: None,
            derived_from: Vec::new(),
            emits: Vec::new(),
            logic: None,
        };

        validate_bioscript_interpretation(&interpretation).unwrap();
        assert_eq!(
            analysis_output_relative_file("sample", &interpretation.id, "tsv"),
            "analysis/sample/apoe_epsilon.tsv"
        );
        assert_eq!(
            analysis_observations_relative_file("sample", &interpretation.id),
            "analysis/sample/apoe_epsilon.observations.tsv"
        );

        interpretation.kind = "shell".to_owned();
        assert_eq!(
            validate_bioscript_interpretation(&interpretation).unwrap_err(),
            "analysis 'apoe_epsilon' uses unsupported kind 'shell'"
        );
    }

    #[test]
    fn analysis_output_json_uses_shared_report_shape() {
        let interpretation = PanelInterpretation {
            id: "apoe_epsilon".to_owned(),
            label: Some("APOE epsilon".to_owned()),
            kind: "bioscript".to_owned(),
            path: "apoe.py".to_owned(),
            output_format: Some("tsv".to_owned()),
            derived_from: vec!["rs429358.yaml".to_owned()],
            emits: Vec::new(),
            logic: Some(PanelInterpretationLogic {
                description: Some("APOE logic".to_owned()),
                source: Some(PanelInterpretationLogicSource {
                    name: Some("ClinPGx".to_owned()),
                    url: Some("https://example.test".to_owned()),
                }),
            }),
        };
        let value = analysis_output_json(AnalysisOutputJsonInput {
            participant_id: "sample",
            assay_id: "pgx",
            interpretation: &interpretation,
            output_format: "tsv",
            manifest_path: "assets/APOE/assay.yaml",
            script_path: "assets/APOE/apoe.py",
            output_file: "analysis/sample/apoe_epsilon.tsv",
            observations_file: Some("analysis/sample/apoe_epsilon.observations.tsv"),
            row_headers: vec!["apoe_status".to_owned()],
            rows: vec![serde_json::json!({"apoe_status": "e3/e3"})],
        });

        assert_eq!(value["schema"], "bioscript:analysis-output:1.0");
        assert_eq!(value["analysis_id"], "apoe_epsilon");
        assert_eq!(value["output_format"], "tsv");
        assert_eq!(
            value["observations_file"],
            "analysis/sample/apoe_epsilon.observations.tsv"
        );
        assert_eq!(value["derived_from"][0], "rs429358.yaml");
        assert_eq!(value["assets"].as_array().unwrap().len(), 0);
        assert_eq!(value["emits"].as_array().unwrap().len(), 0);
        assert_eq!(value["logic"]["source"]["name"], "ClinPGx");
        assert_eq!(value["rows"][0]["apoe_status"], "e3/e3");
    }

    #[test]
    fn renders_analysis_observations_for_runtime_scripts() {
        let text = render_analysis_observations_tsv(&[VariantObservation {
            matched_rsid: Some("rs1".to_owned()),
            assembly: Some(Assembly::Grch38),
            genotype: Some("AG".to_owned()),
            ref_count: Some(8),
            alt_count: Some(4),
            depth: Some(12),
            evidence: vec!["source".to_owned()],
            ..VariantObservation::default()
        }]);

        assert!(text.starts_with("matched_rsid\tassembly\tgenotype\t"));
        assert!(text.contains("rs1\tgrch38\tAG\t8\t4\t12\tsource\n"));
    }
}
