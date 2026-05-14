use std::path::Path;

pub fn match_app_findings(
    findings: &[serde_json::Value],
    observations: &[serde_json::Value],
    analyses: &[serde_json::Value],
) -> Vec<serde_json::Value> {
    let mut matched = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for finding in findings {
        if let Some(effects) = finding.get("effects").and_then(serde_json::Value::as_array) {
            for effect in effects {
                if let Some(observation) = app_finding_match_observation(effect, observations) {
                    let mut item = finding.clone();
                    if let Some(object) = item.as_object_mut() {
                        object.remove("effects");
                        object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                        object.insert("matched_effect".to_owned(), effect.clone());
                        object.insert(
                            "matched_observation".to_owned(),
                            app_finding_observation_context(observation),
                        );
                    }
                    if seen.insert(app_finding_dedupe_key(&item)) {
                        matched.push(item);
                    }
                } else if let Some(analysis) = app_finding_match_analysis(effect, analyses) {
                    let mut item = finding.clone();
                    if let Some(object) = item.as_object_mut() {
                        object.remove("effects");
                        object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                        object.insert("matched_effect".to_owned(), effect.clone());
                        object.insert("matched_analysis".to_owned(), analysis);
                    }
                    if seen.insert(app_finding_dedupe_key(&item)) {
                        matched.push(item);
                    }
                }
            }
        } else if let Some(observation) = app_finding_match_observation(finding, observations) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert(
                    "matched_observation".to_owned(),
                    app_finding_observation_context(observation),
                );
            }
            if seen.insert(app_finding_dedupe_key(&item)) {
                matched.push(item);
            }
        } else if let Some(analysis) = app_finding_match_analysis(finding, analyses) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert("matched_analysis".to_owned(), analysis);
            }
            if seen.insert(app_finding_dedupe_key(&item)) {
                matched.push(item);
            }
        }
    }
    matched
}

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
            .filter(|observation| app_binding_chromosome_count_matches(binding, observation))
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
    if key == "alt" {
        let expected_alleles = app_binding_expected_values(binding);
        return observations
            .iter()
            .filter(|observation| !app_variant_ref_mismatch(binding, observation))
            .filter(|observation| app_binding_chromosome_count_matches(binding, observation))
            .find(|observation| {
                app_binding_matches_value(observation.get(key), binding)
                    && expected_alleles.iter().any(|expected| {
                        app_observation_allele_dosage(observation, expected)
                            .is_some_and(|dosage| dosage > 0)
                    })
            });
    }
    observations
        .iter()
        .filter(|observation| !app_variant_ref_mismatch(binding, observation))
        .filter(|observation| app_binding_chromosome_count_matches(binding, observation))
        .find(|observation| app_binding_matches_value(observation.get(key), binding))
}

fn app_finding_observation_context(observation: &serde_json::Value) -> serde_json::Value {
    serde_json::json!({
        "participant_id": observation.get("participant_id").cloned().unwrap_or(serde_json::Value::Null),
        "rsid": observation.get("rsid").cloned().unwrap_or(serde_json::Value::Null),
        "gene": observation.get("gene").cloned().unwrap_or(serde_json::Value::Null),
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
            "hem_ref" | "het" => Some(1),
            "hom_alt" | "hem_alt" => Some(0),
            _ => None,
        };
    }
    if allele == alt_allele {
        return match zygosity {
            "hom_ref" | "hem_ref" => Some(0),
            "het" | "hem_alt" => Some(1),
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

fn app_binding_chromosome_count_matches(
    binding: &serde_json::Value,
    observation: &serde_json::Value,
) -> bool {
    let Some(expected) = binding
        .get("chromosome_count")
        .and_then(serde_json::Value::as_i64)
    else {
        return true;
    };
    app_observation_chromosome_count(observation).is_some_and(|actual| actual == expected)
}

fn app_observation_chromosome_count(observation: &serde_json::Value) -> Option<i64> {
    match observation
        .get("zygosity")
        .and_then(serde_json::Value::as_str)
        .unwrap_or_default()
    {
        "hem_ref" | "hem_alt" => Some(1),
        "hom_ref" | "het" | "hom_alt" => Some(2),
        _ => None,
    }
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

fn app_binding_expected_values(binding: &serde_json::Value) -> Vec<String> {
    let mut values = Vec::new();
    if let Some(value) = binding.get("value").and_then(value_as_string) {
        values.push(value.clone());
    }
    if let Some(array) = binding.get("values").and_then(serde_json::Value::as_array) {
        values.extend(array.iter().filter_map(value_as_string));
    }
    values
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

#[cfg(test)]
mod report_matching_tests {
    use super::*;

    #[test]
    fn alt_binding_requires_observed_allele_dosage() {
        let binding = serde_json::json!({
            "source": "variant",
            "key": "alt",
            "value": "G"
        });
        let observations = vec![
            serde_json::json!({
                "variant_path": "rs1.yaml",
                "ref": "A",
                "alt": "G",
                "genotype_display": "AA",
                "zygosity": "hom_ref"
            }),
            serde_json::json!({
                "variant_path": "rs2.yaml",
                "ref": "A",
                "alt": "G",
                "genotype_display": "AG",
                "zygosity": "het"
            }),
        ];

        let matched = app_variant_binding_match_observation(&binding, &observations)
            .expect("het alt observation should match");
        assert_eq!(
            matched
                .get("genotype_display")
                .and_then(serde_json::Value::as_str),
            Some("AG")
        );
    }

    #[test]
    fn alt_in_binding_requires_observed_allele_dosage() {
        let binding = serde_json::json!({
            "source": "variant",
            "key": "alt",
            "operator": "in",
            "values": ["G", "T"]
        });
        let observations = vec![serde_json::json!({
            "variant_path": "rs1.yaml",
            "ref": "A",
            "alt": "T",
            "genotype_display": "AT",
            "zygosity": "het"
        })];

        assert!(app_variant_binding_match_observation(&binding, &observations).is_some());
    }

    #[test]
    fn hemizygous_observations_count_as_single_allele_dosage() {
        let observations = vec![serde_json::json!({
            "variant_path": "rs3813929.yaml",
            "ref": "C",
            "alt": "T",
            "genotype": "1",
            "genotype_display": "T",
            "zygosity": "hem_alt"
        })];

        let include_binding = serde_json::json!({
            "source": "variant",
            "variant": "rs3813929.yaml",
            "key": "alt",
            "value": "T"
        });
        assert!(app_variant_binding_match_observation(&include_binding, &observations).is_some());

        let effect_binding = serde_json::json!({
            "source": "variant",
            "variant": "rs3813929.yaml",
            "allele": "T",
            "operator": "dosage_equals",
            "value": 1,
            "chromosome_count": 1
        });
        assert!(app_variant_binding_match_observation(&effect_binding, &observations).is_some());
    }

    #[test]
    fn chromosome_count_binding_separates_one_x_and_two_x_rows() {
        let observations = vec![serde_json::json!({
            "variant_path": "rs3813929.yaml",
            "ref": "C",
            "alt": "T",
            "genotype": "0/1",
            "genotype_display": "CT",
            "zygosity": "het"
        })];

        let one_x_binding = serde_json::json!({
            "source": "variant",
            "variant": "rs3813929.yaml",
            "allele": "T",
            "operator": "dosage_equals",
            "value": 1,
            "chromosome_count": 1
        });
        assert!(app_variant_binding_match_observation(&one_x_binding, &observations).is_none());

        let two_x_binding = serde_json::json!({
            "source": "variant",
            "variant": "rs3813929.yaml",
            "allele": "T",
            "operator": "dosage_equals",
            "value": 1,
            "chromosome_count": 2
        });
        assert!(app_variant_binding_match_observation(&two_x_binding, &observations).is_some());
    }

    #[test]
    fn match_app_findings_matches_variant_effects_and_deduplicates_evidence() {
        let findings = vec![
            serde_json::json!({
                "schema": "bioscript:finding:1.0",
                "label": "Repeated",
                "evidence": {"source": "db", "kind": "guideline", "id": "cpic-1"},
                "effects": [
                    {
                        "id": "effect-a",
                        "binding": {
                            "source": "variant",
                            "variant": "variants/rs1.yaml",
                            "key": "outcome",
                            "value": "variant"
                        }
                    },
                    {
                        "id": "effect-a",
                        "binding": {
                            "source": "variant",
                            "variant": "rs1.yaml",
                            "key": "outcome",
                            "value": "variant"
                        }
                    }
                ]
            })
        ];
        let observations = vec![serde_json::json!({
            "participant_id": "p1",
            "variant_key": "rs1",
            "variant_path": "variants/rs1.yaml",
            "rsid": "rs1",
            "gene": "ABC",
            "ref": "A",
            "alt": "G",
            "genotype_display": "AG",
            "zygosity": "het",
            "outcome": "variant"
        })];

        let matched = match_app_findings(&findings, &observations, &[]);
        assert_eq!(matched.len(), 1);
        assert_eq!(matched[0]["matched"], true);
        assert_eq!(matched[0]["matched_effect"]["id"], "effect-a");
        assert_eq!(matched[0]["matched_observation"]["participant_id"], "p1");
        assert!(matched[0].get("effects").is_none());
    }

    #[test]
    fn match_app_findings_matches_analysis_bindings_with_alias_and_value_types() {
        let findings = vec![
            serde_json::json!({
                "id": "analysis-finding",
                "binding": {
                    "source": "analysis",
                    "analysis": "star-allele",
                    "key": "score",
                    "operator": "in",
                    "values": [1, 2, 3]
                }
            }),
            serde_json::json!({
                "id": "non-match",
                "binding": {
                    "source": "analysis",
                    "analysis_id": "other",
                    "key": "score",
                    "value": 2
                }
            })
        ];
        let analyses = vec![
            serde_json::json!({
                "participant_id": "p1",
                "assay_id": "assay",
                "analysis_id": "star-allele",
                "rows": [
                    {"score": 2, "label": "ok"},
                    {"score": 5, "label": "skip"}
                ]
            }),
            serde_json::json!({
                "participant_id": "p1",
                "assay_id": "assay",
                "analysis_id": "other",
                "rows": "not rows"
            })
        ];

        let matched = match_app_findings(&findings, &[], &analyses);
        assert_eq!(matched.len(), 1);
        assert_eq!(matched[0]["matched_analysis"]["analysis_id"], "star-allele");
        assert_eq!(matched[0]["matched_analysis"]["key"], "score");
        assert_eq!(matched[0]["matched_analysis"]["value"], 2);
    }

    #[test]
    fn binding_helpers_cover_value_dosage_reference_and_dedupe_edges() {
        let observation = serde_json::json!({
            "variant_key": "rs1",
            "variant_path": "nested/rs1.yaml",
            "rsid": "rs1",
            "ref": "A",
            "alt": "G",
            "genotype_display": "AG",
            "zygosity": "het"
        });
        assert!(!app_variant_ref_mismatch(
            &serde_json::json!({"variant": "rs1.yaml"}),
            &observation
        ));
        assert!(app_variant_ref_mismatch(
            &serde_json::json!({"path": "rs2.yaml"}),
            &observation
        ));

        assert_eq!(app_observation_allele_dosage(&observation, "A"), Some(1));
        assert_eq!(app_observation_allele_dosage(&observation, "G"), Some(1));
        assert_eq!(app_observation_allele_dosage(&observation, "T"), Some(0));
        assert_eq!(app_observation_allele_dosage(&observation, "DEL"), None);
        assert_eq!(app_observation_chromosome_count(&observation), Some(2));

        assert!(app_binding_matches_value(
            Some(&serde_json::json!(true)),
            &serde_json::json!({"value": true})
        ));
        assert!(app_binding_matches_value(
            Some(&serde_json::json!(42)),
            &serde_json::json!({"operator": "in", "values": ["41", 42]})
        ));
        assert!(!app_binding_matches_value(
            Some(&serde_json::json!({"object": true})),
            &serde_json::json!({"value": "true"})
        ));
        assert_eq!(
            app_binding_expected_values(&serde_json::json!({
                "value": "A",
                "values": ["B", 3, false, {"ignored": true}]
            })),
            vec!["A", "B", "3", "false"]
        );

        assert!(app_binding_matches_dosage(
            Some(2),
            &serde_json::json!({"operator": "dosage_in", "values": [1, 2]})
        ));
        assert!(!app_binding_matches_dosage(
            None,
            &serde_json::json!({"operator": "dosage_equals", "value": 0})
        ));

        assert_eq!(
            app_finding_dedupe_key(&serde_json::json!({
                "evidence": {"url": "https://example.test/evidence"},
                "matched_effect": {"label": "effect"}
            })),
            "evidence_url|https://example.test/evidence|effect"
        );
        assert_eq!(
            app_finding_dedupe_key(&serde_json::json!({
                "schema": "s",
                "label": "l",
                "notes": "n",
                "matched_effect": {"text": "t"}
            })),
            "content|s|l|n|t"
        );
    }

    #[test]
    fn variant_binding_rejects_missing_keys_and_unsupported_operators() {
        let observations = vec![serde_json::json!({
            "variant_path": "rs1.yaml",
            "ref": "A",
            "alt": "G",
            "genotype_display": "AG",
            "zygosity": "het"
        })];
        assert!(app_variant_binding_match_observation(
            &serde_json::json!({"source": "variant"}),
            &observations
        )
        .is_none());
        assert!(app_variant_binding_match_observation(
            &serde_json::json!({
                "source": "variant",
                "key": "alt",
                "operator": "unknown",
                "value": "G"
            }),
            &observations
        )
        .is_none());
        assert!(app_variant_binding_match_observation(
            &serde_json::json!({
                "source": "variant",
                "allele": "",
                "operator": "dosage_equals",
                "value": 1
            }),
            &observations
        )
        .is_none());
    }
}
