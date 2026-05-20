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
    let workspace = bioscript_reporting::FilesystemManifestWorkspace::new(runtime_root);
    let task = bioscript_reporting::load_variant_manifest_task_by_path(
        &workspace,
        &manifest_path.display().to_string(),
    )?;
    let manifest = task.manifest;
    let (gene, alt_alleles, observed_alt_alleles, source) = if row_path.contains('#') {
        (
            String::new(),
            manifest
                .spec
                .alternate
                .clone()
                .into_iter()
                .collect::<Vec<_>>(),
            manifest.spec.observed_alternates.clone(),
            serde_json::Value::Null,
        )
    } else {
        (
            variant_manifest_gene(&manifest_path)?,
            variant_alt_alleles(&manifest_path)?,
            variant_observed_alt_alleles(&manifest_path)?,
            variant_primary_source(&manifest_path)?,
        )
    };
    Ok(bioscript_reporting::app_observation_from_manifest_row(
        bioscript_reporting::AppObservationInput {
            row,
            row_path: &row_path,
            assay_id,
            manifest,
            gene,
            source,
            alt_alleles,
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

fn variant_alt_alleles(path: &Path) -> Result<Vec<String>, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read variant YAML {}: {err}", path.display()))?;
    let value: serde_yaml::Value = serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse variant YAML {}: {err}", path.display()))?;
    let Some(items) = value
        .as_mapping()
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alleles".to_owned())))
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alts".to_owned())))
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
