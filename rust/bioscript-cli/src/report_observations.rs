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
    let ref_allele = manifest.spec.reference.clone().unwrap_or_default();
    let reportable_alt = manifest.spec.alternate.clone().unwrap_or_default();
    let observed_alt_alleles = variant_observed_alt_alleles(&manifest_path)?;
    let mut genotype_display = row
        .get("genotype")
        .filter(|value| !value.is_empty())
        .cloned()
        .or_else(|| genotype_display_from_raw_counts(row.get("raw_counts")?))
        .unwrap_or_default();
    let depth = parse_optional_u32(row.get("depth"));
    let ref_count = parse_optional_u32(row.get("ref_count"));
    let alt_count = parse_optional_u32(row.get("alt_count"));
    if let Some(normalized_display) =
        deletion_copy_number_display(row, &manifest, depth, alt_count)
    {
        genotype_display = normalized_display;
    }
    let weak_indel_match = is_weak_delimited_indel_match(row, &manifest, &genotype_display);
    let allele_balance = match (alt_count, depth) {
        (Some(alt_count), Some(depth)) if depth > 0 => {
            Some(f64::from(alt_count) / f64::from(depth))
        }
        _ => None,
    };
    let assembly = row
        .get("assembly")
        .filter(|value| !value.is_empty())
        .cloned()
        .or_else(|| fallback_assembly.map(assembly_row_value))
        .unwrap_or_default();
    let locus = if assembly.eq_ignore_ascii_case("grch37") {
        manifest.spec.grch37.as_ref()
    } else {
        manifest
            .spec
            .grch38
            .as_ref()
            .or(manifest.spec.grch37.as_ref())
    };
    let chrom = locus.map_or(String::new(), |locus| locus.chrom.clone());
    let (genotype, zygosity) = normalize_app_genotype(
        &genotype_display,
        &ref_allele,
        &reportable_alt,
        manifest.spec.kind,
        &chrom,
        inferred_sex,
    );
    let non_reportable_status =
        classify_non_reportable_alleles(&genotype_display, &ref_allele, &reportable_alt, &observed_alt_alleles);
    let call = observation_call_values(
        depth,
        non_reportable_status,
        &genotype,
        &zygosity,
        &genotype_display,
    );
    let evidence_raw = observation_evidence_raw(row, &chrom, inferred_sex);
    let source = variant_primary_source(&manifest_path)?;
    let kind = manifest
        .spec
        .kind
        .map_or("unknown".to_owned(), |kind| format!("{kind:?}").to_lowercase());
    Ok(render_app_observation_json(AppObservationJson {
        allele_balance,
        alt_count,
        assay_id: assay_id.to_owned(),
        assembly,
        call,
        chrom,
        depth,
        evidence_raw,
        gene,
        genotype,
        genotype_display,
        kind,
        locus: locus.cloned(),
        manifest,
        non_reportable_status,
        observed_alt_alleles,
        ref_allele,
        ref_count,
        reportable_alt,
        row: row.clone(),
        row_path,
        source,
        weak_indel_match,
        zygosity,
    }))
}

fn assembly_row_value(assembly: bioscript_core::Assembly) -> String {
    match assembly {
        bioscript_core::Assembly::Grch37 => "grch37".to_owned(),
        bioscript_core::Assembly::Grch38 => "grch38".to_owned(),
    }
}

fn hemizygous_display_genotype(display: &str) -> String {
    display
        .chars()
        .find(char::is_ascii_alphabetic)
        .map_or_else(|| display.to_owned(), |allele| allele.to_string())
}

fn deletion_copy_number_display(
    row: &BTreeMap<String, String>,
    manifest: &bioscript_schema::VariantManifest,
    depth: Option<u32>,
    alt_count: Option<u32>,
) -> Option<String> {
    if !matches!(manifest.spec.kind, Some(bioscript_core::VariantKind::Deletion)) {
        return None;
    }
    if row.get("backend").map(String::as_str) != Some("cram") {
        return None;
    }
    if manifest.spec.reference.as_deref().unwrap_or_default().len() <= 1 {
        return None;
    }
    let depth = depth?;
    if depth == 0 {
        return None;
    }
    let alt_fraction = f64::from(alt_count.unwrap_or(0)) / f64::from(depth);
    if alt_fraction >= 0.8 {
        Some("DD".to_owned())
    } else if alt_fraction <= 0.2 {
        Some("II".to_owned())
    } else {
        Some("DI".to_owned())
    }
}

fn variant_primary_source(path: &Path) -> Result<serde_json::Value, String> {
    let value = load_yaml_value(path)?;
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    collect_manifest_provenance_entries(&value, &mut links)?;
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

fn normalize_app_genotype(
    display: &str,
    ref_allele: &str,
    alt_allele: &str,
    kind: Option<bioscript_core::VariantKind>,
    chrom: &str,
    inferred_sex: Option<&SexInference>,
) -> (String, String) {
    if display.is_empty() {
        return ("./.".to_owned(), "unknown".to_owned());
    }
    if matches!(kind, Some(bioscript_core::VariantKind::Deletion))
        && ref_allele.len() != 1
        && display
            .chars()
            .filter(char::is_ascii_alphabetic)
            .all(|allele| matches!(allele.to_ascii_uppercase(), 'I' | 'D'))
    {
        return normalize_app_genotype(display, "I", "D", None, chrom, inferred_sex);
    }
    let alleles: Vec<char> = display.chars().filter(char::is_ascii_alphabetic).collect();
    if ref_allele.len() != 1 || alt_allele.len() != 1 {
        return (display.to_owned(), "unknown".to_owned());
    }
    let ref_ch = ref_allele.chars().next().unwrap_or_default();
    let alt_ch = alt_allele.chars().next().unwrap_or_default();
    if alleles.len() == 1 && is_haploid_sex_chromosome(chrom) {
        let allele = alleles[0];
        if allele == ref_ch {
            return ("0".to_owned(), "hem_ref".to_owned());
        }
        if allele == alt_ch {
            return ("1".to_owned(), "hem_alt".to_owned());
        }
        return (display.to_owned(), "unknown".to_owned());
    }
    if alleles.len() != 2 {
        return (display.to_owned(), "unknown".to_owned());
    }
    if is_confident_male_sex_chromosome(chrom, inferred_sex) && alleles[0] == alleles[1] {
        let allele = alleles[0];
        if allele == ref_ch {
            return ("0".to_owned(), "hem_ref".to_owned());
        }
        if allele == alt_ch {
            return ("1".to_owned(), "hem_alt".to_owned());
        }
        return (display.to_owned(), "unknown".to_owned());
    }
    let alt_count = alleles.iter().filter(|allele| **allele == alt_ch).count();
    let ref_count = alleles.iter().filter(|allele| **allele == ref_ch).count();
    match (ref_count, alt_count) {
        (2, 0) => ("0/0".to_owned(), "hom_ref".to_owned()),
        (1, 1) => ("0/1".to_owned(), "het".to_owned()),
        (0, 2) => ("1/1".to_owned(), "hom_alt".to_owned()),
        _ => (display.to_owned(), "unknown".to_owned()),
    }
}

fn is_confident_male_sex_chromosome(chrom: &str, inferred_sex: Option<&SexInference>) -> bool {
    is_haploid_sex_chromosome(chrom)
        && inferred_sex.is_some_and(|sex| {
            sex.sex == InferredSex::Male
                && matches!(
                    sex.confidence,
                    SexDetectionConfidence::High | SexDetectionConfidence::Medium
                )
        })
}

fn is_haploid_sex_chromosome(chrom: &str) -> bool {
    matches!(
        chrom
            .trim()
            .trim_start_matches("chr")
            .trim_start_matches("CHR")
            .to_ascii_uppercase()
            .as_str(),
        "X" | "Y" | "23" | "24"
    )
}

fn observation_evidence_raw(
    row: &BTreeMap<String, String>,
    chrom: &str,
    inferred_sex: Option<&SexInference>,
) -> String {
    let mut evidence_raw = row.get("evidence").cloned().unwrap_or_default();
    if !is_haploid_sex_chromosome(chrom) {
        return evidence_raw;
    }
    let Some(inferred_sex) = inferred_sex else {
        return evidence_raw;
    };
    let sex_evidence = sex_inference_evidence_raw(inferred_sex);
    if sex_evidence.is_empty() {
        return evidence_raw;
    }
    if evidence_raw.is_empty() {
        evidence_raw = sex_evidence;
    } else {
        evidence_raw.push_str(" | ");
        evidence_raw.push_str(&sex_evidence);
    }
    evidence_raw
}

fn sex_inference_evidence_raw(inferred_sex: &SexInference) -> String {
    let sex = match inferred_sex.sex {
        InferredSex::Male => "male",
        InferredSex::Female => "female",
        InferredSex::Unknown => "unknown",
    };
    let confidence = match inferred_sex.confidence {
        SexDetectionConfidence::High => "high",
        SexDetectionConfidence::Medium => "medium",
        SexDetectionConfidence::Low => "low",
    };
    let mut fields = vec![
        format!("detected_sex={sex}"),
        format!("sex_confidence={confidence}"),
        format!("sex_method={}", inferred_sex.method),
    ];
    fields.extend(
        inferred_sex
            .evidence
            .iter()
            .map(|item| format!("sex_{item}")),
    );
    fields.join(" ")
}

fn genotype_display_from_raw_counts(raw_counts: &str) -> Option<String> {
    let counts: serde_json::Map<String, serde_json::Value> = serde_json::from_str(raw_counts).ok()?;
    let mut items = counts
        .into_iter()
        .filter_map(|(base, count)| {
            let base = base.chars().next()?.to_ascii_uppercase();
            let count = count.as_u64()?;
            if matches!(base, 'A' | 'C' | 'G' | 'T') && count > 0 {
                Some((base, count))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    if items.is_empty() {
        return None;
    }
    items.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    let total = items.iter().map(|(_, count)| *count).sum::<u64>();
    let (top_base, top_count) = items[0];
    if total == 0 || items.len() == 1 || top_count.saturating_mul(10) >= total.saturating_mul(8) {
        return Some(format!("{top_base}{top_base}"));
    }
    Some(format!("{}{}", top_base, items[1].0))
}

fn classify_non_reportable_alleles(
    display: &str,
    ref_allele: &str,
    reportable_alt: &str,
    observed_alts: &[String],
) -> Option<&'static str> {
    if display.is_empty() || ref_allele.len() != 1 || reportable_alt.len() != 1 {
        return None;
    }
    let ref_ch = ref_allele.chars().next()?.to_ascii_uppercase();
    let alt_ch = reportable_alt.chars().next()?.to_ascii_uppercase();
    let non_reportable = display
        .chars()
        .filter(char::is_ascii_alphabetic)
        .map(|ch| ch.to_ascii_uppercase())
        .filter(|ch| *ch != ref_ch && *ch != alt_ch)
        .collect::<Vec<_>>();
    if non_reportable.is_empty() {
        return None;
    }
    if non_reportable.iter().all(|ch| {
        observed_alts.iter().any(|alt| {
            alt.len() == 1
                && alt
                    .chars()
                    .next()
                    .is_some_and(|alt_ch| alt_ch.to_ascii_uppercase() == *ch)
        })
    }) {
        Some("observed_alt")
    } else {
        Some("unknown_alt")
    }
}

fn is_weak_delimited_indel_match(
    row: &BTreeMap<String, String>,
    manifest: &bioscript_schema::VariantManifest,
    genotype_display: &str,
) -> bool {
    if !matches!(manifest.spec.kind, Some(bioscript_core::VariantKind::Deletion)) {
        return false;
    }
    if !matches!(row.get("backend").map(String::as_str), Some("text" | "zip")) {
        return false;
    }
    if manifest.spec.reference.as_deref().unwrap_or_default().len() <= 1 {
        return false;
    }
    genotype_display
        .chars()
        .filter(char::is_ascii_alphabetic)
        .all(|allele| matches!(allele.to_ascii_uppercase(), 'I' | 'D'))
}

fn observation_facets(
    non_reportable_status: Option<&str>,
    observed_alts: &[String],
) -> serde_json::Value {
    let mut facets = Vec::new();
    if let Some(status) = non_reportable_status {
        facets.push(status.to_owned());
        if status == "observed_alt" && !observed_alts.is_empty() {
            facets.push(format!("known_observed_alts={}", observed_alts.join(",")));
        }
    }
    if facets.is_empty() {
        serde_json::Value::Null
    } else {
        serde_json::Value::String(facets.join(";"))
    }
}

fn parse_optional_u32(value: Option<&String>) -> Option<u32> {
    value.and_then(|value| value.parse::<u32>().ok())
}
