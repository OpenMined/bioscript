fn app_finding_match_observation<'a>(
    finding: &serde_json::Value,
    observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let binding = finding.get("binding")?;
    match binding.get("source").and_then(serde_json::Value::as_str) {
        Some("variant") => app_variant_binding_match_observation(binding, observations),
        _ => None,
    }
}

fn app_finding_match_analysis(
    finding: &serde_json::Value,
    analyses: &[serde_json::Value],
) -> Option<serde_json::Value> {
    let binding = finding.get("binding")?;
    if binding.get("source").and_then(serde_json::Value::as_str) != Some("analysis") {
        return None;
    }
    let analysis_id = binding
        .get("analysis_id")
        .or_else(|| binding.get("analysis"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let key = binding.get("key").and_then(serde_json::Value::as_str)?;
    for analysis in analyses {
        if !analysis_id.is_empty()
            && analysis
                .get("analysis_id")
                .and_then(serde_json::Value::as_str)
                != Some(analysis_id)
        {
            continue;
        }
        let Some(rows) = analysis.get("rows").and_then(serde_json::Value::as_array) else {
            continue;
        };
        for row in rows {
            if app_binding_matches_value(row.get(key), binding) {
                return Some(serde_json::json!({
                    "participant_id": analysis.get("participant_id").cloned().unwrap_or(serde_json::Value::Null),
                    "assay_id": analysis.get("assay_id").cloned().unwrap_or(serde_json::Value::Null),
                    "analysis_id": analysis.get("analysis_id").cloned().unwrap_or(serde_json::Value::Null),
                    "key": key,
                    "value": row.get(key).cloned().unwrap_or(serde_json::Value::Null),
                    "row": row,
                }));
            }
        }
    }
    None
}

fn app_variant_binding_match_observation<'a>(
    binding: &serde_json::Value,
    observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let operator = binding
        .get("operator")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("equals");
    if matches!(operator, "dosage_equals" | "dosage_in") {
        let allele = binding
            .get("allele")
            .and_then(serde_json::Value::as_str)
            .unwrap_or_default();
        return observations
            .iter()
            .filter(|observation| !app_variant_ref_mismatch(binding, observation))
            .find(|observation| {
                let dosage = app_observation_allele_dosage(observation, allele);
                app_binding_matches_dosage(dosage, binding)
            });
    }

    let key = binding
        .get("key")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if key.is_empty() {
        return None;
    }
    observations
        .iter()
        .filter(|observation| !app_variant_ref_mismatch(binding, observation))
        .find(|observation| app_binding_matches_value(observation.get(key), binding))
}

fn app_finding_observation_context(observation: &serde_json::Value) -> serde_json::Value {
    serde_json::json!({
        "participant_id": observation.get("participant_id").cloned().unwrap_or(serde_json::Value::Null),
        "rsid": observation.get("rsid").cloned().unwrap_or(serde_json::Value::Null),
        "ref": observation.get("ref").cloned().unwrap_or(serde_json::Value::Null),
        "alt": observation.get("alt").cloned().unwrap_or(serde_json::Value::Null),
        "genotype_display": observation.get("genotype_display").cloned().unwrap_or(serde_json::Value::Null),
        "outcome": observation.get("outcome").cloned().unwrap_or(serde_json::Value::Null),
    })
}

fn app_variant_ref_mismatch(binding: &serde_json::Value, observation: &serde_json::Value) -> bool {
    let variant_ref = binding
        .get("variant")
        .or_else(|| binding.get("path"))
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if variant_ref.is_empty() {
        return false;
    }
    let basename = Path::new(variant_ref)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(variant_ref);
    let candidates = [
        observation
            .get("variant_key")
            .and_then(serde_json::Value::as_str),
        observation
            .get("variant_path")
            .and_then(serde_json::Value::as_str),
        observation.get("rsid").and_then(serde_json::Value::as_str),
    ];
    !candidates.into_iter().flatten().any(|candidate| {
        candidate == variant_ref
            || Path::new(candidate)
                .file_name()
                .and_then(|value| value.to_str())
                .is_some_and(|value| value == basename)
    })
}

fn app_observation_allele_dosage(observation: &serde_json::Value, allele: &str) -> Option<i64> {
    if allele.is_empty() {
        return None;
    }
    let ref_allele = observation
        .get("ref")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let alt_allele = observation
        .get("alt")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    let zygosity = observation
        .get("zygosity")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if allele == ref_allele {
        return match zygosity {
            "hom_ref" => Some(2),
            "het" => Some(1),
            "hom_alt" => Some(0),
            _ => None,
        };
    }
    if allele == alt_allele {
        return match zygosity {
            "hom_ref" => Some(0),
            "het" => Some(1),
            "hom_alt" => Some(2),
            _ => None,
        };
    }
    let display = observation
        .get("genotype_display")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default();
    if allele.len() == 1 {
        let allele_ch = allele.chars().next()?.to_ascii_uppercase();
        return display
            .chars()
            .filter(|ch| ch.to_ascii_uppercase() == allele_ch)
            .count()
            .try_into()
            .ok();
    }
    None
}

fn app_binding_matches_value(
    actual: Option<&serde_json::Value>,
    binding: &serde_json::Value,
) -> bool {
    let actual = actual.and_then(value_as_string).unwrap_or_default();
    match binding
        .get("operator")
        .and_then(serde_json::Value::as_str)
        .unwrap_or("equals")
    {
        "equals" => binding
            .get("value")
            .and_then(value_as_string)
            .is_some_and(|value| value == actual),
        "in" => binding
            .get("values")
            .and_then(serde_json::Value::as_array)
            .is_some_and(|values| {
                values
                    .iter()
                    .filter_map(value_as_string)
                    .any(|value| value == actual)
            }),
        _ => false,
    }
}

fn app_binding_matches_dosage(dosage: Option<i64>, binding: &serde_json::Value) -> bool {
    let Some(dosage) = dosage else {
        return false;
    };
    match binding
        .get("operator")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "dosage_equals" => binding
            .get("value")
            .and_then(serde_json::Value::as_i64)
            .is_some_and(|value| value == dosage),
        "dosage_in" => binding
            .get("values")
            .and_then(serde_json::Value::as_array)
            .is_some_and(|values| {
                values
                    .iter()
                    .filter_map(serde_json::Value::as_i64)
                    .any(|value| value == dosage)
            }),
        _ => false,
    }
}

fn value_as_string(value: &serde_json::Value) -> Option<String> {
    match value {
        serde_json::Value::String(value) => Some(value.clone()),
        serde_json::Value::Number(value) => Some(value.to_string()),
        serde_json::Value::Bool(value) => Some(value.to_string()),
        _ => None,
    }
}

fn app_finding_dedupe_key(finding: &serde_json::Value) -> String {
    let effect_key = finding
        .get("matched_effect")
        .and_then(|effect| {
            effect
                .get("id")
                .or_else(|| effect.get("label"))
                .or_else(|| effect.get("text"))
        })
        .and_then(value_as_string)
        .unwrap_or_default();
    if let Some(evidence) = finding.get("evidence") {
        let source = evidence
            .get("source")
            .and_then(value_as_string)
            .unwrap_or_default();
        let kind = evidence
            .get("kind")
            .and_then(value_as_string)
            .unwrap_or_default();
        let id = evidence
            .get("id")
            .and_then(value_as_string)
            .unwrap_or_default();
        if !source.is_empty() || !kind.is_empty() || !id.is_empty() {
            return format!("evidence|{source}|{kind}|{id}|{effect_key}");
        }
        if let Some(url) = evidence.get("url").and_then(value_as_string) {
            return format!("evidence_url|{url}|{effect_key}");
        }
    }
    if let Some(id) = finding.get("id").and_then(value_as_string) {
        return format!("id|{id}|{effect_key}");
    }
    format!(
        "content|{}|{}|{}|{}",
        finding
            .get("schema")
            .and_then(value_as_string)
            .unwrap_or_default(),
        finding
            .get("label")
            .and_then(value_as_string)
            .unwrap_or_default(),
        finding
            .get("notes")
            .and_then(value_as_string)
            .unwrap_or_default(),
        effect_key
    )
}

