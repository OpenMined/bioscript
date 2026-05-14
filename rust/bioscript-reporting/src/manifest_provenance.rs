use std::collections::BTreeMap;

use super::{
    ManifestWorkspace, manifest_supports_findings, traversable_manifest_member_paths, yaml_string,
    yaml_to_json,
};

pub fn load_manifest_provenance_links(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<Vec<serde_json::Value>, String> {
    let value = workspace.load_yaml(path)?;
    let schema = yaml_string(&value, "schema").unwrap_or_default();
    let mut links = BTreeMap::new();

    collect_manifest_provenance_entries(&value, &mut links)?;

    if manifest_supports_findings(&schema)
        && let Some(items) = value
            .get("findings")
            .and_then(serde_yaml::Value::as_sequence)
    {
        for item in items {
            if let Some(include) = item.get("include").and_then(serde_yaml::Value::as_str) {
                let include_path = workspace.resolve(path, include)?;
                let included = load_manifest_provenance_links(workspace, &include_path)?;
                for source in included {
                    insert_provenance_link(&mut links, source);
                }
            }
        }
    }

    for member_path in traversable_manifest_member_paths(&schema, &value) {
        let resolved = workspace.resolve(path, member_path)?;
        let member_links = load_manifest_provenance_links(workspace, &resolved)?;
        for source in member_links {
            insert_provenance_link(&mut links, source);
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
        let source = yaml_to_json(source.clone())?;
        insert_provenance_link(links, source);
    }
    Ok(())
}

fn insert_provenance_link(
    links: &mut BTreeMap<String, serde_json::Value>,
    source: serde_json::Value,
) {
    let key = source
        .get("url")
        .or_else(|| source.get("url_template"))
        .or_else(|| source.get("label"))
        .and_then(serde_json::Value::as_str)
        .map_or_else(|| source.to_string(), ToOwned::to_owned);
    links.entry(key).or_insert(source);
}
