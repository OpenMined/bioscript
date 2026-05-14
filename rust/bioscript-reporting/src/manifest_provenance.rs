use std::collections::BTreeMap;

use super::{ManifestWorkspace, traversable_manifest_member_paths, yaml_string, yaml_to_json};

pub fn load_manifest_provenance_links(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<Vec<serde_json::Value>, String> {
    let value = workspace.load_yaml(path)?;
    let schema = yaml_string(&value, "schema").unwrap_or_default();
    let mut links = BTreeMap::new();
    collect_manifest_provenance_entries(&value, &mut links)?;

    for member_path in traversable_manifest_member_paths(&schema, &value) {
        let resolved = workspace.resolve(path, member_path)?;
        for source in load_manifest_provenance_links(workspace, &resolved)? {
            let key = provenance_source_key(&source);
            links.entry(key).or_insert(source);
        }
    }

    Ok(links.into_values().collect())
}

pub fn collect_manifest_provenance_entries(
    value: &serde_yaml::Value,
    links: &mut BTreeMap<String, serde_json::Value>,
) -> Result<(), String> {
    let Some(sources) = value
        .get("provenance")
        .and_then(|provenance| provenance.get("sources"))
        .and_then(serde_yaml::Value::as_sequence)
    else {
        return Ok(());
    };

    for source in sources {
        let json_source = yaml_to_json(source.clone())?;
        let key = provenance_source_key(&json_source);
        links.entry(key).or_insert(json_source);
    }

    Ok(())
}

fn provenance_source_key(source: &serde_json::Value) -> String {
    source
        .get("url")
        .and_then(serde_json::Value::as_str)
        .or_else(|| source.get("label").and_then(serde_json::Value::as_str))
        .map(ToOwned::to_owned)
        .unwrap_or_else(|| source.to_string())
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::manifest::tests::MapWorkspace;

    #[test]
    fn collect_manifest_provenance_entries_deduplicates_sources() {
        let value: serde_yaml::Value = serde_yaml::from_str(
            r#"
schema: bioscript:variant:1.0
provenance:
  sources:
    - kind: database
      label: dbSNP
      url: https://www.ncbi.nlm.nih.gov/snp/rs1
    - kind: database
      label: Duplicate dbSNP
      url: https://www.ncbi.nlm.nih.gov/snp/rs1
"#,
        )
        .unwrap();

        let mut links = BTreeMap::new();
        collect_manifest_provenance_entries(&value, &mut links).unwrap();
        assert_eq!(links.len(), 1);
        assert_eq!(
            links.values().next().unwrap()["label"],
            serde_json::Value::String("dbSNP".to_owned())
        );
    }

    #[test]
    fn load_manifest_provenance_links_walks_members() {
        let workspace = MapWorkspace {
            files: BTreeMap::from([
                (
                    "panel.yaml".to_owned(),
                    r#"
schema: bioscript:panel:1.0
members:
  - kind: variant
    path: rs1.yaml
  - kind: assay
    path: assay/manifest.yaml
provenance:
  sources:
    - kind: database
      label: Panel
      url: https://example.test/panel
"#
                    .to_owned(),
                ),
                (
                    "rs1.yaml".to_owned(),
                    r#"
schema: bioscript:variant:1.0
provenance:
  sources:
    - kind: database
      label: Variant
      url: https://example.test/variant
"#
                    .to_owned(),
                ),
                (
                    "assay/manifest.yaml".to_owned(),
                    r#"
schema: bioscript:assay:1.0
members:
  - kind: variant
    path: rs2.yaml
provenance:
  sources:
    - kind: database
      label: Assay
      url: https://example.test/assay
"#
                    .to_owned(),
                ),
                (
                    "assay/rs2.yaml".to_owned(),
                    r#"
schema: bioscript:variant:1.0
provenance:
  sources:
    - kind: database
      label: Nested Variant
      url: https://example.test/nested
"#
                    .to_owned(),
                ),
            ]),
        };

        let links = load_manifest_provenance_links(&workspace, "panel.yaml").unwrap();
        let urls = links
            .iter()
            .filter_map(|link| link.get("url").and_then(serde_json::Value::as_str))
            .collect::<Vec<_>>();
        assert_eq!(
            urls,
            vec![
                "https://example.test/assay",
                "https://example.test/nested",
                "https://example.test/panel",
                "https://example.test/variant",
            ]
        );
    }
}
