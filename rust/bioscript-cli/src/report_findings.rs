fn load_manifest_findings(
    root: &Path,
    manifest_path: &Path,
) -> Result<Vec<serde_json::Value>, String> {
    let value = load_yaml_value(manifest_path)?;
    let schema = value
        .get("schema")
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or_default();
    let mut findings = Vec::new();

    if matches!(
        schema,
        "bioscript:variant:1.0"
            | "bioscript:variant"
            | "bioscript:assay:1.0"
            | "bioscript:panel:1.0"
            | "bioscript:pgx-findings:1.0"
    ) && let Some(items) = value
        .get("findings")
        .and_then(serde_yaml::Value::as_sequence)
    {
        for item in items {
            let json_item = yaml_to_json(item.clone())?;
            let include = json_item
                .get("include")
                .and_then(serde_json::Value::as_str)
                .map(str::to_owned);
            if let Some(include) = include {
                let include_path = resolve_manifest_path(root, manifest_path, &include)?;
                let mut included = load_manifest_findings(root, &include_path)?;
                let inherited_binding = json_item.get("binding").cloned();
                for included_item in &mut included {
                    if inherited_binding.is_some()
                        && included_item.get("binding").is_none()
                        && included_item.get("effects").is_none()
                        && let Some(object) = included_item.as_object_mut()
                    {
                        object.insert(
                            "binding".to_owned(),
                            inherited_binding.clone().unwrap_or(serde_json::Value::Null),
                        );
                    }
                }
                findings.extend(included);
                continue;
            }
            if json_item.get("include").is_none() {
                findings.push(json_item);
            }
        }
    }

    if matches!(schema, "bioscript:assay:1.0" | "bioscript:panel:1.0")
        && let Some(items) = value
            .get("members")
            .and_then(serde_yaml::Value::as_sequence)
    {
        for member in items {
            let Some(kind) = member.get("kind").and_then(serde_yaml::Value::as_str) else {
                continue;
            };
            if !matches!(kind, "variant" | "assay") {
                continue;
            }
            let Some(path) = member.get("path").and_then(serde_yaml::Value::as_str) else {
                continue;
            };
            let member_path = resolve_manifest_path(root, manifest_path, path)?;
            findings.extend(load_manifest_findings(root, &member_path)?);
        }
    }

    Ok(findings)
}

fn load_yaml_value(path: &Path) -> Result<serde_yaml::Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|err| format!("failed to read YAML {}: {err}", path.display()))?;
    serde_yaml::from_str(&text)
        .map_err(|err| format!("failed to parse YAML {}: {err}", path.display()))
}

fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, String> {
    serde_json::to_value(value).map_err(|err| format!("failed to convert YAML to JSON: {err}"))
}

fn load_manifest_provenance_links(
    root: &Path,
    manifest_path: &Path,
) -> Result<Vec<serde_json::Value>, String> {
    let value = load_yaml_value(manifest_path)?;
    let schema = value
        .get("schema")
        .and_then(serde_yaml::Value::as_str)
        .unwrap_or_default();
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    collect_manifest_provenance_entries(&value, &mut links)?;

    if matches!(
        schema,
        "bioscript:variant:1.0"
            | "bioscript:variant"
            | "bioscript:assay:1.0"
            | "bioscript:panel:1.0"
            | "bioscript:pgx-findings:1.0"
    ) && let Some(items) = value
        .get("findings")
        .and_then(serde_yaml::Value::as_sequence)
    {
        for item in items {
            let json_item = yaml_to_json(item.clone())?;
            let Some(include) = json_item.get("include").and_then(serde_json::Value::as_str) else {
                continue;
            };
            let include_path = resolve_manifest_path(root, manifest_path, include)?;
            for item in load_manifest_provenance_links(root, &include_path)? {
                if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                    links.entry(url.to_owned()).or_insert(item);
                }
            }
        }
    }

    if matches!(schema, "bioscript:assay:1.0" | "bioscript:panel:1.0")
        && let Some(items) = value
            .get("members")
            .and_then(serde_yaml::Value::as_sequence)
    {
        for member in items {
            let Some(kind) = member.get("kind").and_then(serde_yaml::Value::as_str) else {
                continue;
            };
            if !matches!(kind, "variant" | "assay") {
                continue;
            }
            let Some(path) = member.get("path").and_then(serde_yaml::Value::as_str) else {
                continue;
            };
            let member_path = resolve_manifest_path(root, manifest_path, path)?;
            for item in load_manifest_provenance_links(root, &member_path)? {
                if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                    links.entry(url.to_owned()).or_insert(item);
                }
            }
        }
    }

    Ok(links.into_values().collect())
}

fn collect_manifest_provenance_entries(
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

fn match_app_findings(
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
                    let key = app_finding_dedupe_key(&item);
                    if seen.insert(key) {
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
                    let key = app_finding_dedupe_key(&item);
                    if seen.insert(key) {
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
            let key = app_finding_dedupe_key(&item);
            if seen.insert(key) {
                matched.push(item);
            }
        } else if let Some(analysis) = app_finding_match_analysis(finding, analyses) {
            let mut item = finding.clone();
            if let Some(object) = item.as_object_mut() {
                object.insert("matched".to_owned(), serde_json::Value::Bool(true));
                object.insert("matched_analysis".to_owned(), analysis);
            }
            let key = app_finding_dedupe_key(&item);
            if seen.insert(key) {
                matched.push(item);
            }
        }
    }
    matched
}

#[cfg(test)]
mod report_observations_tests {
    use super::*;

    #[test]
    fn raw_counts_can_fill_display_for_homozygous_and_heterozygous_observations() {
        assert_eq!(
            genotype_display_from_raw_counts(r#"{"T": 24}"#).as_deref(),
            Some("TT")
        );
        assert_eq!(
            genotype_display_from_raw_counts(r#"{"C": 12, "T": 10}"#).as_deref(),
            Some("CT")
        );
    }

    #[test]
    fn non_reportable_alleles_are_classified_as_observed_or_unknown() {
        assert_eq!(
            classify_non_reportable_alleles("TT", "C", "G", &["T".to_owned()]),
            Some("observed_alt")
        );
        assert_eq!(
            classify_non_reportable_alleles("AT", "C", "G", &["T".to_owned()]),
            Some("unknown_alt")
        );
        assert_eq!(
            classify_non_reportable_alleles("CG", "C", "G", &["T".to_owned()]),
            None
        );
    }

    #[test]
    fn single_allele_sex_chromosome_calls_are_treated_as_hemizygous() {
        assert_eq!(
            normalize_app_genotype("G", "C", "G", "X", None),
            ("1".to_owned(), "hem_alt".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("C", "C", "G", "chrX", None),
            ("0".to_owned(), "hem_ref".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("G", "C", "G", "1", None),
            ("G".to_owned(), "unknown".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("GG", "C", "G", "X", None),
            ("1/1".to_owned(), "hom_alt".to_owned())
        );
    }

    #[test]
    fn confident_male_sex_chromosome_duplicate_calls_are_hemizygous() {
        let inferred_sex = SexInference {
            sex: InferredSex::Male,
            confidence: SexDetectionConfidence::High,
            method: "vcf_non_par_x_gt".to_owned(),
            evidence: vec!["called_y_snps=1200".to_owned()],
        };
        assert_eq!(
            normalize_app_genotype("GG", "C", "G", "X", Some(&inferred_sex)),
            ("1".to_owned(), "hem_alt".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("CC", "C", "G", "chrX", Some(&inferred_sex)),
            ("0".to_owned(), "hem_ref".to_owned())
        );
    }
}
