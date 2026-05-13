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
    row.insert(
        "matched_rsid".to_owned(),
        observation.matched_rsid.clone().unwrap_or_default(),
    );
    row.insert(
        "assembly".to_owned(),
        observation
            .assembly
            .map(assembly_row_value)
            .unwrap_or_default(),
    );
    row.insert(
        "genotype".to_owned(),
        observation.genotype.clone().unwrap_or_default(),
    );
    row.insert(
        "ref_count".to_owned(),
        observation
            .ref_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "alt_count".to_owned(),
        observation
            .alt_count
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "depth".to_owned(),
        observation
            .depth
            .map_or_else(String::new, |value| value.to_string()),
    );
    row.insert(
        "raw_counts".to_owned(),
        serde_json::to_string(&observation.raw_counts).unwrap_or_default(),
    );
    row.insert("evidence".to_owned(), observation.evidence.join(" | "));
    row
}

pub(super) fn render_app_observations_tsv(
    observations: &[serde_json::Value],
) -> Result<String, JsError> {
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

/// Derive the assay id from a manifest path — matches the CLI's
/// `bioscript-cli::report_execution::app_assay_id`, which loads the manifest
/// and returns its `name:` field (panels / assays / variants all carry one).
/// This function operates on a `PackageWorkspace` so it can find files in the
/// in-memory map without touching disk.
///
/// Previously the wasm derived the id from the manifest filename stem (e.g.
/// `manifest.yaml` -> `manifest`), which diverged from the CLI's `pgx-1`
/// (panel `name:` field) and cascaded into the HTML report's
/// `participant_id × assay_id` keys.
pub(super) fn app_assay_id_from_workspace(
    workspace: &PackageWorkspace,
    manifest_path: &str,
) -> Result<String, JsError> {
    match workspace.schema(manifest_path)?.as_str() {
        "bioscript:panel:1.0" => Ok(workspace.load_panel(manifest_path)?.name),
        "bioscript:assay:1.0" => Ok(workspace.load_assay(manifest_path)?.name),
        "bioscript:variant:1.0" | "bioscript:variant" => {
            Ok(workspace.load_variant(manifest_path)?.name)
        }
        other => Err(JsError::new(&format!(
            "unsupported manifest schema '{other}'"
        ))),
    }
}

pub(super) fn app_assay_id(path: &Path) -> Result<String, JsError> {
    path.file_stem()
        .and_then(|value| value.to_str())
        .map(ToOwned::to_owned)
        .ok_or_else(|| {
            JsError::new(&format!(
                "failed to derive assay id from {}",
                path.display()
            ))
        })
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
                rows.push(serde_json::from_str(line).map_err(|err| {
                    JsError::new(&format!("failed to parse analysis JSONL: {err}"))
                })?);
            }
            let row_headers = rows
                .iter()
                .find_map(|row| row.as_object())
                .map(|object| object.keys().cloned().collect())
                .unwrap_or_default();
            Ok((rows, row_headers))
        }
        other => Err(JsError::new(&format!(
            "unsupported analysis output_format '{other}'"
        ))),
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
                        serde_json::Value::String(
                            fields.get(index).copied().unwrap_or_default().to_owned(),
                        ),
                    )
                })
                .collect();
            serde_json::Value::Object(object)
        })
        .collect();
    (rows, headers)
}

pub(super) fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, JsError> {
    serde_json::to_value(value)
        .map_err(|err| JsError::new(&format!("failed to convert YAML to JSON: {err}")))
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

pub(super) fn input_inspection_json(
    inspection: &bioscript_formats::FileInspection,
) -> serde_json::Value {
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
        "format_confidence": detection_confidence_name(inspection.confidence),
        "assembly": inspection.assembly.map(|assembly| match assembly {
            Assembly::Grch37 => "grch37",
            Assembly::Grch38 => "grch38",
        }),
        "phased": inspection.phased,
        "selected_entry": inspection.selected_entry,
        "has_index": inspection.has_index,
        "index_path": inspection.index_path.as_ref().map(|path| path.display().to_string()),
        "reference_matches": inspection.reference_matches,
        "source": inspection.source.as_ref().map(|source| serde_json::json!({
            "vendor": source.vendor,
            "platform_version": source.platform_version,
            "confidence": detection_confidence_name(source.confidence),
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

fn detection_confidence_name(value: bioscript_formats::DetectionConfidence) -> &'static str {
    match value {
        bioscript_formats::DetectionConfidence::Authoritative => "authoritative",
        bioscript_formats::DetectionConfidence::StrongHeuristic => "strong_heuristic",
        bioscript_formats::DetectionConfidence::WeakHeuristic => "weak_heuristic",
        bioscript_formats::DetectionConfidence::Unknown => "unknown",
    }
}

pub(super) fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
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
    links
        .into_values()
        .next()
        .unwrap_or(serde_json::Value::Null)
}

pub(super) fn normalize_app_genotype(
    display: &str,
    ref_allele: &str,
    alt_allele: &str,
    kind: Option<VariantKind>,
    chrom: &str,
    inferred_sex: Option<&SexInference>,
) -> (String, String) {
    if display.is_empty() {
        return ("./.".to_owned(), "unknown".to_owned());
    }
    if matches!(kind, Some(VariantKind::Deletion))
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
    if is_confident_male_sex_chromosome(chrom, inferred_sex)
        && alleles.len() == 2
        && alleles[0] == alleles[1]
    {
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

pub(super) fn deletion_copy_number_display(
    row: &BTreeMap<String, String>,
    manifest: &VariantManifest,
    depth: Option<u32>,
    alt_count: Option<u32>,
) -> Option<String> {
    if !matches!(manifest.spec.kind, Some(VariantKind::Deletion)) {
        return None;
    }
    if !matches!(
        row.get("backend").map(String::as_str),
        Some("cram" | "bam")
    ) {
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

fn is_confident_male_sex_chromosome(chrom: &str, inferred_sex: Option<&SexInference>) -> bool {
    matches!(
        chrom
            .trim()
            .trim_start_matches("chr")
            .to_ascii_uppercase()
            .as_str(),
        "X" | "Y" | "23" | "24"
    ) && inferred_sex.is_some_and(|sex| {
        sex.sex == InferredSex::Male
            && matches!(
                sex.confidence,
                SexDetectionConfidence::High | SexDetectionConfidence::Medium
            )
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes_long_deletion_reference_tokens_as_insertion_deletion_copy_number() {
        assert_eq!(
            normalize_app_genotype(
                "II",
                "TTATAA",
                "<DEL:6>",
                Some(VariantKind::Deletion),
                "22",
                None,
            ),
            ("0/0".to_owned(), "hom_ref".to_owned())
        );
        assert_eq!(
            normalize_app_genotype(
                "ID",
                "TTATAA",
                "<DEL:6>",
                Some(VariantKind::Deletion),
                "22",
                None,
            ),
            ("0/1".to_owned(), "het".to_owned())
        );
    }

    #[test]
    fn displays_cram_long_deletion_copy_number_as_insertion_deletion_tokens() {
        let manifest = VariantManifest {
            path: PathBuf::from("rs71785313.yaml"),
            name: "APOL1_G2".to_owned(),
            tags: Vec::new(),
            spec: VariantSpec {
                reference: Some("TTATAA".to_owned()),
                alternate: Some("<DEL:6>".to_owned()),
                kind: Some(VariantKind::Deletion),
                ..VariantSpec::default()
            },
        };
        let mut row = BTreeMap::new();
        row.insert("backend".to_owned(), "cram".to_owned());

        assert_eq!(
            deletion_copy_number_display(&row, &manifest, Some(39), Some(0)).as_deref(),
            Some("II")
        );
        assert_eq!(
            deletion_copy_number_display(&row, &manifest, Some(39), Some(39)).as_deref(),
            Some("DD")
        );
        assert_eq!(
            deletion_copy_number_display(&row, &manifest, Some(40), Some(20)).as_deref(),
            Some("DI")
        );
    }
}
