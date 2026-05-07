use super::*;

pub(super) fn artifact(name: &str, mime_type: &str, text: String) -> ReportArtifactOutput {
    ReportArtifactOutput {
        name: name.to_owned(),
        path: name.to_owned(),
        mime_type: mime_type.to_owned(),
        text,
    }
}

pub(super) fn variant_row(
    path: &str,
    name: &str,
    tags: &[String],
    observation: &VariantObservation,
    participant_id: &str,
) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("kind".to_owned(), "variant".to_owned());
    row.insert("name".to_owned(), name.to_owned());
    row.insert("path".to_owned(), path.to_owned());
    row.insert("tags".to_owned(), tags.join(","));
    row.insert("backend".to_owned(), observation.backend.clone());
    row.insert("participant_id".to_owned(), participant_id.to_owned());
    row.insert("matched_rsid".to_owned(), observation.matched_rsid.clone().unwrap_or_default());
    row.insert("assembly".to_owned(), observation.assembly.map(assembly_row_value).unwrap_or_default());
    row.insert("genotype".to_owned(), observation.genotype.clone().unwrap_or_default());
    row.insert("ref_count".to_owned(), observation.ref_count.map_or_else(String::new, |value| value.to_string()));
    row.insert("alt_count".to_owned(), observation.alt_count.map_or_else(String::new, |value| value.to_string()));
    row.insert("depth".to_owned(), observation.depth.map_or_else(String::new, |value| value.to_string()));
    row.insert("raw_counts".to_owned(), serde_json::to_string(&observation.raw_counts).unwrap_or_default());
    row.insert("evidence".to_owned(), observation.evidence.join(" | "));
    row
}

pub(super) fn render_app_observations_tsv(observations: &[serde_json::Value]) -> Result<String, JsError> {
    let mut out = OBSERVATION_TSV_HEADERS.join("\t");
    out.push('\n');
    for observation in observations {
        let line = OBSERVATION_TSV_HEADERS
            .iter()
            .map(|header| json_field_as_tsv(observation.get(*header)))
            .collect::<Vec<_>>()
            .join("\t");
        out.push_str(&line);
        out.push('\n');
    }
    Ok(out)
}

pub(super) fn render_jsonl(rows: &[serde_json::Value]) -> Result<String, JsError> {
    let mut out = String::new();
    for row in rows {
        out.push_str(&serde_json::to_string(row).map_err(|err| JsError::new(&err.to_string()))?);
        out.push('\n');
    }
    Ok(out)
}

pub(super) fn json_field_as_tsv(value: Option<&serde_json::Value>) -> String {
    match value {
        Some(serde_json::Value::Null) | None => String::new(),
        Some(serde_json::Value::String(value)) => value.replace(['\t', '\n'], " "),
        Some(value) => value.to_string().replace(['\t', '\n'], " "),
    }
}

pub(super) fn normalize_package_path(path: &str) -> Result<String, JsError> {
    let mut out = PathBuf::new();
    for component in Path::new(path).components() {
        match component {
            std::path::Component::Normal(value) => out.push(value),
            std::path::Component::CurDir => {}
            _ => return Err(JsError::new(&format!("unsafe package path: {path}"))),
        }
    }
    Ok(out.display().to_string().replace('\\', "/"))
}

pub(super) fn default_analysis_max_duration_ms() -> u64 {
    30_000
}

pub(super) fn participant_id_from_name(path: &str) -> String {
    Path::new(path)
        .file_stem()
        .and_then(|value| value.to_str())
        .unwrap_or(path)
        .replace([' ', '\t', '\n'], "_")
}

pub(super) fn app_assay_id(path: &Path) -> Result<String, JsError> {
    path.file_stem()
        .and_then(|value| value.to_str())
        .map(ToOwned::to_owned)
        .ok_or_else(|| JsError::new(&format!("failed to derive assay id from {}", path.display())))
}

pub(super) fn matches_filters(manifest: &VariantManifest, path: &str, filters: &[String]) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("kind", value)) => value == "variant",
        Some(("name", value)) => manifest.name.contains(value),
        Some(("path", value)) => path.contains(value),
        Some(("tag", value)) => manifest.tags.iter().any(|tag| tag == value),
        Some(_) | None => false,
    })
}

pub(super) fn parse_analysis_output_text(
    text: &str,
    format: &str,
) -> Result<(Vec<serde_json::Value>, Vec<String>), JsError> {
    match format {
        "tsv" => Ok(parse_analysis_tsv(text)),
        "json" => {
            let value: serde_json::Value = serde_json::from_str(text)
                .map_err(|err| JsError::new(&format!("failed to parse analysis JSON: {err}")))?;
            let rows = match value {
                serde_json::Value::Array(rows) => rows,
                serde_json::Value::Object(mut object) => object
                    .remove("rows")
                    .and_then(|rows| rows.as_array().cloned())
                    .unwrap_or_else(|| vec![serde_json::Value::Object(object)]),
                other => vec![other],
            };
            let row_headers = rows
                .iter()
                .find_map(|row| row.as_object())
                .map(|object| object.keys().cloned().collect())
                .unwrap_or_default();
            Ok((rows, row_headers))
        }
        "jsonl" => {
            let mut rows: Vec<serde_json::Value> = Vec::new();
            for line in text.lines().filter(|line| !line.trim().is_empty()) {
                rows.push(serde_json::from_str(line)
                    .map_err(|err| JsError::new(&format!("failed to parse analysis JSONL: {err}")))?);
            }
            let row_headers = rows
                .iter()
                .find_map(|row| row.as_object())
                .map(|object| object.keys().cloned().collect())
                .unwrap_or_default();
            Ok((rows, row_headers))
        }
        other => Err(JsError::new(&format!("unsupported analysis output_format '{other}'"))),
    }
}

fn parse_analysis_tsv(text: &str) -> (Vec<serde_json::Value>, Vec<String>) {
    let mut lines = text.lines();
    let headers = lines
        .next()
        .map(|line| line.split('\t').map(ToOwned::to_owned).collect::<Vec<_>>())
        .unwrap_or_default();
    let rows = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let object = headers
                .iter()
                .enumerate()
                .map(|(index, header)| {
                    (
                        header.clone(),
                        serde_json::Value::String(fields.get(index).copied().unwrap_or_default().to_owned()),
                    )
                })
                .collect();
            serde_json::Value::Object(object)
        })
        .collect();
    (rows, headers)
}

pub(super) fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, JsError> {
    serde_json::to_value(value).map_err(|err| JsError::new(&format!("failed to convert YAML to JSON: {err}")))
}

pub(super) fn collect_manifest_provenance_entries(
    value: &serde_yaml::Value,
    links: &mut BTreeMap<String, serde_json::Value>,
) -> Result<(), JsError> {
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

pub(super) fn input_inspection_json(inspection: &bioscript_formats::FileInspection) -> serde_json::Value {
    serde_json::json!({
        "container": match inspection.container {
            bioscript_formats::FileContainer::Plain => "plain",
            bioscript_formats::FileContainer::Zip => "zip",
        },
        "format": match inspection.detected_kind {
            bioscript_formats::DetectedKind::GenotypeText => "genotype_text",
            bioscript_formats::DetectedKind::Vcf => "vcf",
            bioscript_formats::DetectedKind::AlignmentCram => "alignment_cram",
            bioscript_formats::DetectedKind::AlignmentBam => "alignment_bam",
            bioscript_formats::DetectedKind::ReferenceFasta => "reference_fasta",
            bioscript_formats::DetectedKind::Unknown => "unknown",
        },
        "format_confidence": match inspection.confidence {
            bioscript_formats::DetectionConfidence::Authoritative => "authoritative",
            bioscript_formats::DetectionConfidence::StrongHeuristic => "strong_heuristic",
            bioscript_formats::DetectionConfidence::WeakHeuristic => "weak_heuristic",
            bioscript_formats::DetectionConfidence::Unknown => "unknown",
        },
        "assembly": inspection.assembly.map(|assembly| match assembly {
            Assembly::Grch37 => "grch37",
            Assembly::Grch38 => "grch38",
        }),
        "selected_entry": inspection.selected_entry,
        "source": inspection.source.as_ref().map(|source| serde_json::json!({
            "vendor": source.vendor,
            "platform_version": source.platform_version,
            "evidence": source.evidence,
        })),
        "inferred_sex": inspection.inferred_sex.as_ref().map(|sex| serde_json::json!({
            "sex": inferred_sex_name(sex.sex),
            "confidence": sex_detection_confidence_name(sex.confidence),
            "method": sex.method,
            "evidence": sex.evidence,
        })),
        "evidence": inspection.evidence,
        "warnings": inspection.warnings,
        "duration_ms": inspection.duration_ms,
    })
}

pub(super) fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value.get(key).and_then(serde_yaml::Value::as_str).map(ToOwned::to_owned)
}

pub(super) fn yaml_string_sequence(value: &serde_yaml::Value, key: &str) -> Vec<serde_json::Value> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_sequence)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_yaml::Value::as_str)
                .map(serde_json::Value::from)
                .collect()
        })
        .unwrap_or_default()
}

pub(super) fn yaml_mapping_string(mapping: &serde_yaml::Mapping, key: &str) -> Option<String> {
    mapping
        .get(serde_yaml::Value::String(key.to_owned()))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

pub(super) fn variant_primary_source_from_yaml(value: &serde_yaml::Value) -> serde_json::Value {
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    let _ = collect_manifest_provenance_entries(value, &mut links);
    if let Some(rsid) = value
        .get("identifiers")
        .and_then(|identifiers| identifiers.get("rsids"))
        .and_then(serde_yaml::Value::as_sequence)
        .and_then(|items| items.iter().find_map(serde_yaml::Value::as_str))
    {
        return serde_json::json!({
            "kind": "database",
            "label": "dbSNP / NCBI SNP",
            "url": format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"),
            "fields": ["identifiers.rsids"],
        });
    }
    links.into_values().next().unwrap_or(serde_json::Value::Null)
}

pub(super) fn normalize_app_genotype(
    display: &str,
    ref_allele: &str,
    alt_allele: &str,
    chrom: &str,
    inferred_sex: Option<&SexInference>,
) -> (String, String) {
    if display.is_empty() {
        return ("./.".to_owned(), "unknown".to_owned());
    }
    let alleles: Vec<char> = display.chars().filter(char::is_ascii_alphabetic).collect();
    if ref_allele.len() != 1 || alt_allele.len() != 1 {
        return (display.to_owned(), "unknown".to_owned());
    }
    let ref_ch = ref_allele.chars().next().unwrap_or_default();
    let alt_ch = alt_allele.chars().next().unwrap_or_default();
    if is_confident_male_sex_chromosome(chrom, inferred_sex) && alleles.len() == 2 && alleles[0] == alleles[1] {
        let allele = alleles[0];
        if allele == ref_ch {
            return ("0".to_owned(), "hem_ref".to_owned());
        }
        if allele == alt_ch {
            return ("1".to_owned(), "hem_alt".to_owned());
        }
    }
    if alleles.len() != 2 {
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
    matches!(
        chrom.trim().trim_start_matches("chr").to_ascii_uppercase().as_str(),
        "X" | "Y" | "23" | "24"
    ) && inferred_sex.is_some_and(|sex| {
        sex.sex == InferredSex::Male
            && matches!(sex.confidence, SexDetectionConfidence::High | SexDetectionConfidence::Medium)
    })
}

pub(super) fn assembly_row_value(assembly: Assembly) -> String {
    match assembly {
        Assembly::Grch37 => "grch37".to_owned(),
        Assembly::Grch38 => "grch38".to_owned(),
    }
}

fn inferred_sex_name(value: InferredSex) -> &'static str {
    match value {
        InferredSex::Male => "male",
        InferredSex::Female => "female",
        InferredSex::Unknown => "unknown",
    }
}

fn sex_detection_confidence_name(value: SexDetectionConfidence) -> &'static str {
    match value {
        SexDetectionConfidence::High => "high",
        SexDetectionConfidence::Medium => "medium",
        SexDetectionConfidence::Low => "low",
    }
}

pub(super) fn parse_optional_u32(value: Option<&String>) -> Option<u32> {
    value.and_then(|value| value.parse::<u32>().ok())
}
