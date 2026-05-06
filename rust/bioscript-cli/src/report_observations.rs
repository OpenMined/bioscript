fn app_observation_from_manifest_row(
    runtime_root: &Path,
    row: &BTreeMap<String, String>,
    assay_id: &str,
) -> Result<serde_json::Value, String> {
    let row_path = row.get("path").cloned().unwrap_or_default();
    let manifest_path = if Path::new(&row_path).is_absolute() {
        PathBuf::from(&row_path)
    } else {
        runtime_root.join(&row_path)
    };
    let manifest = load_variant_manifest(&manifest_path)?;
    let ref_allele = manifest.spec.reference.clone().unwrap_or_default();
    let genotype_display = row.get("genotype").cloned().unwrap_or_default();
    let alt_alleles = variant_alt_alleles(&manifest_path)?;
    let alt_allele = observed_alt_allele(&genotype_display, &ref_allele, &alt_alleles)
        .or_else(|| manifest.spec.alternate.clone())
        .unwrap_or_default();
    let (genotype, zygosity) = normalize_app_genotype(&genotype_display, &ref_allele, &alt_allele);
    let depth = parse_optional_u32(row.get("depth"));
    let ref_count = parse_optional_u32(row.get("ref_count"));
    let alt_count = parse_optional_u32(row.get("alt_count"));
    let allele_balance = match (alt_count, depth) {
        (Some(alt_count), Some(depth)) if depth > 0 => {
            Some(f64::from(alt_count) / f64::from(depth))
        }
        _ => None,
    };
    let assembly = row.get("assembly").cloned().unwrap_or_default();
    let locus = if assembly.eq_ignore_ascii_case("grch37") {
        manifest.spec.grch37.as_ref()
    } else {
        manifest
            .spec
            .grch38
            .as_ref()
            .or(manifest.spec.grch37.as_ref())
    };
    let outcome = if genotype == "./." {
        "no_call"
    } else if zygosity == "hom_ref" {
        "reference"
    } else if zygosity == "het" || zygosity == "hom_alt" {
        "variant"
    } else {
        "unknown"
    };
    let evidence_raw = row.get("evidence").cloned().unwrap_or_default();
    Ok(serde_json::json!({
        "participant_id": row.get("participant_id").cloned().unwrap_or_default(),
        "assay_id": assay_id,
        "assay_version": "1.0",
        "variant_key": manifest.name,
        "variant_path": row_path,
        "rsid": row.get("matched_rsid").filter(|value| !value.is_empty()).cloned().or_else(|| manifest.spec.rsids.first().cloned()),
        "assembly": if assembly.is_empty() { serde_json::Value::Null } else { serde_json::Value::String(assembly.to_uppercase()) },
        "chrom": locus.map_or(String::new(), |locus| locus.chrom.clone()),
        "pos_start": locus.map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.start)),
        "pos_end": locus.map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.end)),
        "ref": ref_allele,
        "alt": alt_allele,
        "kind": manifest.spec.kind.map_or("unknown".to_owned(), |kind| format!("{kind:?}").to_lowercase()),
        "match_status": if row.get("matched_rsid").is_some_and(|value| !value.is_empty()) || !genotype_display.is_empty() { "found" } else { "not_found" },
        "coverage_status": depth.map_or("covered", |depth| if depth > 0 { "covered" } else { "not_covered" }),
        "call_status": if genotype == "./." { "no_call" } else { "called" },
        "genotype": genotype,
        "genotype_display": genotype_display,
        "zygosity": zygosity,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "depth": depth,
        "genotype_quality": serde_json::Value::Null,
        "allele_balance": allele_balance,
        "outcome": outcome,
        "evidence_type": if row.get("backend").is_some_and(|value| value == "cram") { "mpileup" } else { "genotype_file" },
        "evidence_raw": evidence_raw,
        "facets": serde_json::Value::Null,
    }))
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
        .and_then(|mapping| {
            mapping
                .get(serde_yaml::Value::String("observed_alts".to_owned()))
                .or_else(|| mapping.get(serde_yaml::Value::String("alts".to_owned())))
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

fn observed_alt_allele(
    genotype_display: &str,
    ref_allele: &str,
    alts: &[String],
) -> Option<String> {
    if ref_allele.len() != 1 {
        return None;
    }
    let ref_ch = ref_allele.chars().next()?;
    genotype_display
        .chars()
        .filter(|ch| ch.is_ascii_alphabetic() && *ch != ref_ch)
        .find_map(|ch| {
            alts.iter()
                .find(|alt| alt.len() == 1 && alt.starts_with(ch))
                .cloned()
        })
}

fn normalize_app_genotype(display: &str, ref_allele: &str, alt_allele: &str) -> (String, String) {
    if display.is_empty() {
        return ("./.".to_owned(), "unknown".to_owned());
    }
    let alleles: Vec<char> = display.chars().filter(char::is_ascii_alphabetic).collect();
    if alleles.len() != 2 || ref_allele.len() != 1 || alt_allele.len() != 1 {
        return (display.to_owned(), "unknown".to_owned());
    }
    let ref_ch = ref_allele.chars().next().unwrap_or_default();
    let alt_ch = alt_allele.chars().next().unwrap_or_default();
    let alt_count = alleles.iter().filter(|allele| **allele == alt_ch).count();
    let ref_count = alleles.iter().filter(|allele| **allele == ref_ch).count();
    match (ref_count, alt_count) {
        (2, 0) => ("0/0".to_owned(), "hom_ref".to_owned()),
        (1, 1) => ("0/1".to_owned(), "het".to_owned()),
        (0, 2) => ("1/1".to_owned(), "hom_alt".to_owned()),
        _ => (display.to_owned(), "unknown".to_owned()),
    }
}

fn parse_optional_u32(value: Option<&String>) -> Option<u32> {
    value.and_then(|value| value.parse::<u32>().ok())
}

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

