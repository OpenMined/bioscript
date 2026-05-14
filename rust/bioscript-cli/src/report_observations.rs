fn app_observation_from_manifest_row(
    runtime_root: &Path,
    row: &BTreeMap<String, String>,
    assay_id: &str,
    inferred_sex: Option<&SexInference>,
    fallback_assembly: Option<bioscript_core::Assembly>,
) -> Result<serde_json::Value, String> {
    let row_path = row.get("path").cloned().unwrap_or_default();
    let manifest_path = if Path::new(&row_path).is_absolute() {
        PathBuf::from(&row_path)
    } else {
        runtime_root.join(&row_path)
    };
    let manifest = load_variant_manifest(&manifest_path)?;
    let gene = variant_manifest_gene(&manifest_path)?;
    let observed_alt_alleles = variant_observed_alt_alleles(&manifest_path)?;
    let source = variant_primary_source(&manifest_path)?;
    Ok(bioscript_reporting::app_observation_from_manifest_row(
        bioscript_reporting::AppObservationInput {
            row,
            row_path: &row_path,
            assay_id,
            manifest,
            gene,
            source,
            observed_alt_alleles,
            inferred_sex,
            fallback_assembly,
        },
    ))
}

fn load_yaml_value(path: &Path) -> Result<serde_yaml::Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read YAML {}: {err}", path.display()))?;
    serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))
}

fn variant_primary_source(path: &Path) -> Result<serde_json::Value, String> {
    let value = load_yaml_value(path)?;
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    bioscript_reporting::collect_manifest_provenance_entries(&value, &mut links)?;
    if let Some(source) = links
        .values()
        .find(|source| source_url_contains(source, "ncbi.nlm.nih.gov/snp/rs"))
    {
        return Ok(source.clone());
    }
    if let Some(rsid) = value
        .get("identifiers")
        .and_then(|identifiers| identifiers.get("rsids"))
        .and_then(serde_yaml::Value::as_sequence)
        .and_then(|items| items.iter().find_map(serde_yaml::Value::as_str))
    {
        return Ok(serde_json::json!({
            "kind": "database",
            "label": "dbSNP / NCBI SNP",
            "url": format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"),
            "fields": ["identifiers.rsids"],
        }));
    }
    Ok(links.into_values().next().unwrap_or(serde_json::Value::Null))
}

fn source_url_contains(source: &serde_json::Value, needle: &str) -> bool {
    source
        .get("url")
        .and_then(serde_json::Value::as_str)
        .is_some_and(|url| url.contains(needle))
}

fn variant_manifest_gene(path: &Path) -> Result<String, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read variant YAML {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse variant YAML {}: {err}", path.display()))?;
    Ok(value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("gene".to_owned())))
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or_default()
        .to_owned())
}

fn variant_observed_alt_alleles(path: &Path) -> Result<Vec<String>, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read variant YAML {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse variant YAML {}: {err}", path.display()))?;
    let Some(items) = value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alleles".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| {
            mapping
                .get(serde_yaml::Value::String("observed_alts".to_owned()))
        })
        .and_then(serde_yaml::Value::as_sequence)
    else {
        return Ok(Vec::new());
    };
    Ok(items
        .iter()
        .filter_map(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
        .collect())
}

#[cfg(test)]
mod app_report_observation_tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = env::temp_dir().join(format!(
            "bioscript-report-observations-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_variant_yaml(dir: &Path, name: &str, body: &str) -> PathBuf {
        let path = dir.join(name);
        fs::write(&path, body).unwrap();
        path
    }

    #[test]
    fn manifest_helpers_extract_gene_alts_and_primary_source() {
        let dir = temp_dir("variant");
        let path = write_variant_yaml(
            &dir,
            "variant.yaml",
            r#"
schema: bioscript:variant:1.0
name: Test variant
gene: CYP2D6
identifiers:
  rsids: [rs123]
alleles:
  observed_alts: [A, T]
evidence:
  references:
    - label: Primary
      url: https://www.ncbi.nlm.nih.gov/snp/rs123
"#,
        );

        assert_eq!(variant_manifest_gene(&path).unwrap(), "CYP2D6");
        assert_eq!(
            variant_observed_alt_alleles(&path).unwrap(),
            vec!["A".to_owned(), "T".to_owned()]
        );
        let source = variant_primary_source(&path).unwrap();
        assert!(source["url"].as_str().unwrap().contains("rs123"));
        assert!(source_url_contains(&source, "ncbi.nlm.nih.gov/snp"));

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn primary_source_falls_back_to_rsid_when_no_provenance_link_exists() {
        let dir = temp_dir("rsid");
        let path = write_variant_yaml(
            &dir,
            "variant.yaml",
            r#"
schema: bioscript:variant:1.0
name: Test variant
identifiers:
  rsids:
    - rs4242
"#,
        );

        let source = variant_primary_source(&path).unwrap();
        assert_eq!(source["label"], "dbSNP / NCBI SNP");
        assert_eq!(source["url"], "https://www.ncbi.nlm.nih.gov/snp/rs4242");

        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn yaml_helpers_report_read_and_parse_errors() {
        let dir = temp_dir("errors");
        let missing = dir.join("missing.yaml");
        assert!(load_yaml_value(&missing).unwrap_err().contains("failed to read YAML"));

        let bad = write_variant_yaml(&dir, "bad.yaml", "name: [unterminated");
        assert!(load_yaml_value(&bad).unwrap_err().contains("failed to parse YAML"));
        assert!(variant_manifest_gene(&bad)
            .unwrap_err()
            .contains("failed to parse variant YAML"));

        let no_alts = write_variant_yaml(&dir, "no-alts.yaml", "gene: ABC\n");
        assert!(variant_observed_alt_alleles(&no_alts).unwrap().is_empty());

        fs::remove_dir_all(dir).unwrap();
    }
}
