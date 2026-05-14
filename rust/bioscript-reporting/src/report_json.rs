#[derive(Clone, Copy)]
pub struct AppReportJsonInput<'a> {
    pub assay_id: &'a str,
    pub participant_id: &'a str,
    pub input_file_name: &'a str,
    pub input_file_path: &'a str,
    pub observations: &'a [serde_json::Value],
    pub analyses: &'a [serde_json::Value],
    pub findings: &'a [serde_json::Value],
    pub provenance: &'a [serde_json::Value],
    pub input_inspection: Option<&'a bioscript_formats::FileInspection>,
    pub manifest_metadata: &'a serde_json::Value,
}

#[derive(Clone, Copy)]
pub struct AppInputReportInput<'a> {
    pub assay_id: &'a str,
    pub participant_id: &'a str,
    pub input_file_name: &'a str,
    pub input_file_path: &'a str,
    pub observations: &'a [serde_json::Value],
    pub analyses: &'a [serde_json::Value],
    pub findings: &'a [serde_json::Value],
    pub provenance: &'a [serde_json::Value],
    pub input_inspection: Option<&'a bioscript_formats::FileInspection>,
    pub manifest_metadata: &'a serde_json::Value,
}

pub fn app_input_report_json(input: AppInputReportInput<'_>) -> serde_json::Value {
    let matched_findings =
        crate::match_app_findings(input.findings, input.observations, input.analyses);
    app_report_json(AppReportJsonInput {
        assay_id: input.assay_id,
        participant_id: input.participant_id,
        input_file_name: input.input_file_name,
        input_file_path: input.input_file_path,
        observations: input.observations,
        analyses: input.analyses,
        findings: &matched_findings,
        provenance: input.provenance,
        input_inspection: input.input_inspection,
        manifest_metadata: input.manifest_metadata,
    })
}

pub fn app_report_json(input: AppReportJsonInput<'_>) -> serde_json::Value {
    let called = input
        .observations
        .iter()
        .filter(|item| {
            item.get("call_status").and_then(serde_json::Value::as_str) == Some("called")
        })
        .count();
    let input_debug = input.input_inspection.map(|inspection| {
        let mut value = input_inspection_json(inspection);
        if observations_have_imputed_vcf_references(input.observations)
            && let Some(object) = value.as_object_mut()
        {
            object.insert(
                "vcf_missing_reference_imputation".to_owned(),
                serde_json::Value::Bool(true),
            );
        }
        value
    });
    serde_json::json!({
        "schema": "bioscript:report:1.0",
        "version": "1.0",
        "participant_id": input.participant_id,
        "assay_id": input.assay_id,
        "assay_version": "1.0",
        "manifest": input.manifest_metadata,
        "input": {
            "file_name": input.input_file_name,
            "file_path": input.input_file_path,
            "debug": input_debug,
        },
        "report_status": if called == input.observations.len() { "complete" } else { "partial" },
        "derived_from": input.observations.iter().filter_map(|item| item.get("variant_key").cloned()).collect::<Vec<_>>(),
        "analyses": input.analyses,
        "findings": input.findings,
        "provenance": input.provenance,
        "metrics": {
            "n_sites_tested": input.observations.len(),
            "n_sites_called": called,
            "n_sites_missing": input.observations.len().saturating_sub(called),
            "n_analyses": input.analyses.len(),
            "n_findings_matched": input.findings.len(),
        }
    })
}

fn observations_have_imputed_vcf_references(observations: &[serde_json::Value]) -> bool {
    observations.iter().any(|observation| {
        observation
            .get("evidence_raw")
            .and_then(serde_json::Value::as_str)
            .is_some_and(|evidence| {
                evidence.contains("imputed reference genotype from absent variant-only VCF record")
            })
    })
}

fn input_inspection_json(inspection: &bioscript_formats::FileInspection) -> serde_json::Value {
    serde_json::json!({
        "container": file_container_name(inspection.container),
        "format": detected_kind_name(inspection.detected_kind),
        "format_confidence": detection_confidence_name(inspection.confidence),
        "assembly": inspection.assembly.map(assembly_name),
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

fn file_container_name(value: bioscript_formats::FileContainer) -> &'static str {
    match value {
        bioscript_formats::FileContainer::Plain => "plain",
        bioscript_formats::FileContainer::Zip => "zip",
    }
}

fn detected_kind_name(value: bioscript_formats::DetectedKind) -> &'static str {
    match value {
        bioscript_formats::DetectedKind::GenotypeText => "genotype_text",
        bioscript_formats::DetectedKind::Vcf => "vcf",
        bioscript_formats::DetectedKind::AlignmentCram => "alignment_cram",
        bioscript_formats::DetectedKind::AlignmentBam => "alignment_bam",
        bioscript_formats::DetectedKind::ReferenceFasta => "reference_fasta",
        bioscript_formats::DetectedKind::Unknown => "unknown",
    }
}

fn detection_confidence_name(value: bioscript_formats::DetectionConfidence) -> &'static str {
    match value {
        bioscript_formats::DetectionConfidence::Authoritative => "authoritative",
        bioscript_formats::DetectionConfidence::StrongHeuristic => "strong_heuristic",
        bioscript_formats::DetectionConfidence::WeakHeuristic => "weak_heuristic",
        bioscript_formats::DetectionConfidence::Unknown => "unknown",
    }
}

fn assembly_name(value: bioscript_core::Assembly) -> &'static str {
    match value {
        bioscript_core::Assembly::Grch37 => "grch37",
        bioscript_core::Assembly::Grch38 => "grch38",
    }
}

fn inferred_sex_name(value: bioscript_formats::InferredSex) -> &'static str {
    match value {
        bioscript_formats::InferredSex::Male => "male",
        bioscript_formats::InferredSex::Female => "female",
        bioscript_formats::InferredSex::Unknown => "unknown",
    }
}

fn sex_detection_confidence_name(value: bioscript_formats::SexDetectionConfidence) -> &'static str {
    match value {
        bioscript_formats::SexDetectionConfidence::High => "high",
        bioscript_formats::SexDetectionConfidence::Medium => "medium",
        bioscript_formats::SexDetectionConfidence::Low => "low",
    }
}
