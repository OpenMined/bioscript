use super::*;

pub(super) fn artifact(name: &str, mime_type: &str, text: String) -> ReportArtifactOutput {
    ReportArtifactOutput {
        name: name.to_owned(),
        path: name.to_owned(),
        mime_type: mime_type.to_owned(),
        text,
    }
}

pub(super) fn encode_report_run_output(
    started_ms: f64,
    artifacts: bioscript_reporting::ReportArtifactTexts,
) -> Result<String, JsError> {
    serde_json::to_string(&ReportRunOutput {
        artifacts: vec![
            artifact(
                "observations.tsv",
                "text/tab-separated-values",
                artifacts.observations_tsv,
            ),
            artifact(
                "analysis.jsonl",
                "application/jsonl",
                artifacts.analysis_jsonl,
            ),
            artifact(
                "reports.jsonl",
                "application/jsonl",
                artifacts.reports_jsonl,
            ),
            artifact("index.html", "text/html", artifacts.html),
        ],
        duration_ms: (js_sys::Date::now() - started_ms).max(0.0) as u128,
        text_output: artifacts.text_output,
    })
    .map_err(|err| JsError::new(&format!("failed to encode report output: {err}")))
}

pub(super) fn normalize_package_path(path: &str) -> Result<String, JsError> {
    let mut out = PathBuf::new();
    for component in Path::new(path).components() {
        match component {
            std::path::Component::Normal(value) => out.push(value),
            std::path::Component::CurDir => {}
            _ => return Err(JsError::new(&format!("unsafe package path: {path}"))),
        }
    }
    Ok(out.display().to_string().replace('\\', "/"))
}

pub(super) fn default_analysis_max_duration_ms() -> u64 {
    30_000
}

pub(super) fn participant_id_from_name(path: &str) -> String {
    bioscript_reporting::participant_id_from_path(Path::new(path))
}

pub(super) fn parse_analysis_output_text(
    text: &str,
    format: &str,
) -> Result<(Vec<serde_json::Value>, Vec<String>), JsError> {
    bioscript_reporting::parse_analysis_output_text(text, format).map_err(|err| JsError::new(&err))
}

pub(super) fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

pub(super) fn variant_primary_source_from_yaml(value: &serde_yaml::Value) -> serde_json::Value {
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    let _ = bioscript_reporting::collect_manifest_provenance_entries(value, &mut links);
    if let Some(source) = links
        .values()
        .find(|source| source_url_contains(source, "ncbi.nlm.nih.gov/snp/rs"))
    {
        return source.clone();
    }
    if let Some(rsid) = value
        .get("identifiers")
        .and_then(|identifiers| identifiers.get("rsids"))
        .and_then(serde_yaml::Value::as_sequence)
        .and_then(|items| items.iter().find_map(serde_yaml::Value::as_str))
    {
        return serde_json::json!({
            "kind": "database",
            "label": "dbSNP / NCBI SNP",
            "url": format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"),
            "fields": ["identifiers.rsids"],
        });
    }
    links
        .into_values()
        .next()
        .unwrap_or(serde_json::Value::Null)
}

fn source_url_contains(source: &serde_json::Value, needle: &str) -> bool {
    source
        .get("url")
        .and_then(serde_json::Value::as_str)
        .is_some_and(|url| url.contains(needle))
}

pub(super) fn variant_observed_alt_alleles_from_yaml(value: &serde_yaml::Value) -> Vec<String> {
    value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alleles".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("observed_alts".to_owned())))
        .and_then(serde_yaml::Value::as_sequence)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_yaml::Value::as_str)
                .map(ToOwned::to_owned)
                .collect()
        })
        .unwrap_or_default()
}

pub(super) fn variant_alt_alleles_from_yaml(value: &serde_yaml::Value) -> Vec<String> {
    value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alleles".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alts".to_owned())))
        .and_then(serde_yaml::Value::as_sequence)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_yaml::Value::as_str)
                .map(ToOwned::to_owned)
                .collect()
        })
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::participant_id_from_name;

    #[test]
    fn participant_id_suffix_stripping_matches_cli_report_path() {
        assert_eq!(
            participant_id_from_name("NA06985.clean.vcf.gz"),
            "NA06985.clean"
        );
        assert_eq!(
            participant_id_from_name("/tmp/genome_hu50B3F5_v5_Full.zip"),
            "genome_hu50B3F5_v5_Full"
        );
        assert_eq!(participant_id_from_name("sample.cram"), "sample");
        assert_eq!(participant_id_from_name("sample name.txt"), "sample name");
    }
}
