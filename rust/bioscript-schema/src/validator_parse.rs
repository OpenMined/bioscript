fn parse_downloads(root: &Value) -> Result<Vec<Download>, String> {
    let mut downloads = Vec::new();
    let Some(items) = value_at(root, &["downloads"]).and_then(Value::as_sequence) else {
        return Ok(downloads);
    };
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            return Err(format!("downloads[{idx}] must be a mapping"));
        };
        let id = mapping_required_string(mapping, "id", idx, "downloads")?;
        let url = mapping_required_string(mapping, "url", idx, "downloads")?;
        let sha256 = mapping_required_string(mapping, "sha256", idx, "downloads")?;
        let version = mapping_required_string(mapping, "version", idx, "downloads")?;
        let origin = normalize_download_url(&url)?;
        downloads.push(Download {
            id,
            url,
            origin,
            sha256,
            version,
        });
    }
    Ok(downloads)
}

fn parse_panel_members(root: &Value) -> Result<Vec<PanelMember>, String> {
    let mut members = Vec::new();
    let Some(items) = value_at(root, &["members"]).and_then(Value::as_sequence) else {
        return Ok(members);
    };
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            return Err(format!("members[{idx}] must be a mapping"));
        };
        members.push(PanelMember {
            kind: mapping_required_string(mapping, "kind", idx, "members")?,
            path: mapping
                .get(Value::String("path".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            download: mapping
                .get(Value::String("download".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            sha256: mapping
                .get(Value::String("sha256".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            version: mapping
                .get(Value::String("version".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
        });
    }
    Ok(members)
}

fn parse_panel_interpretations(root: &Value) -> Result<Vec<PanelInterpretation>, String> {
    let mut interpretations = Vec::new();
    let key = if value_at(root, &["analyses"]).is_some() {
        "analyses"
    } else {
        "interpretations"
    };
    let Some(items) = value_at(root, &[key]).and_then(Value::as_sequence) else {
        return Ok(interpretations);
    };
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            return Err(format!("{key}[{idx}] must be a mapping"));
        };
        interpretations.push(PanelInterpretation {
            id: mapping_required_string(mapping, "id", idx, key)?,
            label: mapping
                .get(Value::String("label".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            kind: mapping_required_string(mapping, "kind", idx, key)?,
            path: mapping_required_string(mapping, "path", idx, key)?,
            output_format: mapping
                .get(Value::String("output_format".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            derived_from: mapping_sequence_of_strings(mapping, "derived_from", idx, key)?,
            emits: parse_panel_interpretation_emits(mapping, idx)?,
            logic: parse_panel_interpretation_logic(mapping)?,
        });
    }
    Ok(interpretations)
}

fn parse_panel_interpretation_logic(
    mapping: &Mapping,
) -> Result<Option<PanelInterpretationLogic>, String> {
    let Some(logic) = mapping.get(Value::String("logic".to_owned())) else {
        return Ok(None);
    };
    let Some(logic_mapping) = logic.as_mapping() else {
        return Err("analysis logic must be a mapping".to_owned());
    };
    let source = match logic_mapping.get(Value::String("source".to_owned())) {
        Some(source) => {
            let Some(source_mapping) = source.as_mapping() else {
                return Err("analysis logic.source must be a mapping".to_owned());
            };
            Some(PanelInterpretationLogicSource {
                name: source_mapping
                    .get(Value::String("name".to_owned()))
                    .and_then(Value::as_str)
                    .map(ToOwned::to_owned),
                url: source_mapping
                    .get(Value::String("url".to_owned()))
                    .and_then(Value::as_str)
                    .map(ToOwned::to_owned),
            })
        }
        None => None,
    };
    Ok(Some(PanelInterpretationLogic {
        source,
        description: logic_mapping
            .get(Value::String("description".to_owned()))
            .and_then(Value::as_str)
            .map(ToOwned::to_owned),
    }))
}

fn parse_panel_interpretation_emits(
    mapping: &Mapping,
    interpretation_idx: usize,
) -> Result<Vec<PanelInterpretationEmit>, String> {
    let Some(items) = mapping
        .get(Value::String("emits".to_owned()))
        .and_then(Value::as_sequence)
    else {
        return Ok(Vec::new());
    };
    let mut emits = Vec::new();
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            return Err(format!(
                "interpretations[{interpretation_idx}].emits[{idx}] must be a mapping"
            ));
        };
        emits.push(PanelInterpretationEmit {
            key: mapping_required_string(
                mapping,
                "key",
                idx,
                &format!("interpretations[{interpretation_idx}].emits"),
            )?,
            label: mapping
                .get(Value::String("label".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            value_type: mapping
                .get(Value::String("value_type".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
            format: mapping
                .get(Value::String("format".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned),
        });
    }
    Ok(emits)
}

fn mapping_sequence_of_strings(
    mapping: &Mapping,
    field: &str,
    idx: usize,
    parent: &str,
) -> Result<Vec<String>, String> {
    let value = mapping
        .get(Value::String(field.to_owned()))
        .ok_or_else(|| format!("{parent}[{idx}].{field} is required"))?;
    let items = value
        .as_sequence()
        .ok_or_else(|| format!("{parent}[{idx}].{field} must be a sequence"))?;
    items
        .iter()
        .enumerate()
        .map(|(item_idx, item)| {
            item.as_str()
                .map(ToOwned::to_owned)
                .ok_or_else(|| format!("{parent}[{idx}].{field}[{item_idx}] must be a string"))
        })
        .collect()
}

fn mapping_required_string(
    mapping: &Mapping,
    field: &str,
    idx: usize,
    parent: &str,
) -> Result<String, String> {
    mapping
        .get(Value::String(field.to_owned()))
        .and_then(Value::as_str)
        .filter(|value| !value.trim().is_empty())
        .map(ToOwned::to_owned)
        .ok_or_else(|| format!("{parent}[{idx}].{field} missing or empty"))
}
