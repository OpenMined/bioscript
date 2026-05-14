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
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    collect_manifest_provenance_entries(&value, &mut links)?;

    if manifest_supports_findings(&schema)
        && let Some(items) = value
            .get("findings")
            .and_then(serde_yaml::Value::as_sequence)
    {
        for item in items {
            let json_item = yaml_to_json(item.clone())?;
            let Some(include) = json_item.get("include").and_then(serde_json::Value::as_str) else {
                continue;
            };
            let include_path = workspace.resolve(path, include)?;
            for item in load_manifest_provenance_links(workspace, &include_path)? {
                if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                    links.entry(url.to_owned()).or_insert(item);
                }
            }
        }
    }

    for member_path in traversable_manifest_member_paths(&schema, &value) {
        let resolved = workspace.resolve(path, member_path)?;
        for item in load_manifest_provenance_links(workspace, &resolved)? {
            if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                links.entry(url.to_owned()).or_insert(item);
            }
        }
    }

    Ok(links.into_values().collect())
}

pub fn collect_manifest_provenance_entries(
    value: &serde_yaml::Value,
    links: &mut BTreeMap<String, serde_json::Value>,
) -> Result<(), String> {
    if let Some(sources) = value
        .get("provenance")
        .and_then(|provenance| provenance.get("sources"))
        .and_then(serde_yaml::Value::as_sequence)
    {
        for source in sources {
            let json = yaml_to_json(source.clone())?;
            if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
                links.entry(url.to_owned()).or_insert(json);
            }
        }
    }
    if let Some(source) = value.get("source") {
        let json = yaml_to_json(source.clone())?;
        if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
            links.entry(url.to_owned()).or_insert(json);
        }
    }
    Ok(())
}
