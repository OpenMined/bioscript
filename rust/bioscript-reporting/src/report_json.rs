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

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use bioscript_core::Assembly;
    use bioscript_formats::{
        DetectedKind, DetectionConfidence, FileContainer, FileInspection, InferredSex,
        SexDetectionConfidence, SexInference, SourceMetadata,
    };
    use serde_json::json;

    use super::{AppInputReportInput, AppReportJsonInput, app_input_report_json, app_report_json};

    fn inspection() -> FileInspection {
        FileInspection {
            path: PathBuf::from("sample.vcf"),
            container: FileContainer::Zip,
            detected_kind: DetectedKind::Vcf,
            confidence: DetectionConfidence::Authoritative,
            source: Some(SourceMetadata {
                vendor: Some("ExampleVendor".to_owned()),
                platform_version: Some("v1".to_owned()),
                confidence: DetectionConfidence::StrongHeuristic,
                evidence: vec!["filename".to_owned()],
            }),
            assembly: Some(Assembly::Grch38),
            phased: Some(true),
            selected_entry: Some("sample.vcf".to_owned()),
            has_index: Some(true),
            index_path: Some(PathBuf::from("sample.vcf.tbi")),
            reference_matches: Some(true),
            inferred_sex: Some(SexInference {
                sex: InferredSex::Female,
                confidence: SexDetectionConfidence::High,
                method: "x_het".to_owned(),
                evidence: vec!["x heterozygosity".to_owned()],
            }),
            evidence: vec!["vcf header".to_owned()],
            warnings: vec!["demo warning".to_owned()],
            duration_ms: 42,
        }
    }

    fn observation(call_status: &str, evidence_raw: &str) -> serde_json::Value {
        json!({
            "variant_key": "rs123",
            "call_status": call_status,
            "outcome": if call_status == "called" { "variant" } else { "unknown" },
            "evidence_raw": evidence_raw,
            "rsid": "rs123",
            "gene": "CYP2C19",
            "ref": "G",
            "alt": "A",
            "genotype_display": "G/A",
            "zygosity": "het"
        })
    }

    #[test]
    fn app_report_json_counts_sites_and_serializes_input_debug() {
        let observations = vec![
            observation("called", "observed genotype"),
            observation(
                "missing",
                "imputed reference genotype from absent variant-only VCF record",
            ),
        ];
        let analyses = vec![json!({"analysis_id": "a1"})];
        let findings = vec![json!({"schema": "bioscript:pgx-summary:1.0"})];
        let provenance = vec![json!({"label": "dbSNP", "url": "https://example.test"})];
        let manifest = json!({"name": "panel", "label": "Panel"});

        let report = app_report_json(AppReportJsonInput {
            assay_id: "assay",
            participant_id: "P001",
            input_file_name: "sample.vcf",
            input_file_path: "/data/sample.vcf",
            observations: &observations,
            analyses: &analyses,
            findings: &findings,
            provenance: &provenance,
            input_inspection: Some(&inspection()),
            manifest_metadata: &manifest,
        });

        assert_eq!(report["report_status"], "partial");
        assert_eq!(report["metrics"]["n_sites_tested"], 2);
        assert_eq!(report["metrics"]["n_sites_called"], 1);
        assert_eq!(report["input"]["debug"]["container"], "zip");
        assert_eq!(report["input"]["debug"]["format"], "vcf");
        assert_eq!(
            report["input"]["debug"]["vcf_missing_reference_imputation"],
            true
        );
        assert_eq!(report["input"]["debug"]["inferred_sex"]["sex"], "female");
    }

    #[test]
    fn app_input_report_json_matches_findings_before_building_report() {
        let observations = vec![observation("called", "observed genotype")];
        let analyses = vec![json!({"analysis_id": "a1", "rows": []})];
        let findings = vec![json!({
            "schema": "bioscript:pgx-summary:1.0",
            "effects": [{
                "binding": {
                    "source": "variant",
                    "variant": "rs123",
                    "key": "alt",
                    "value": "A"
                },
                "text": "A allele observed"
            }]
        })];
        let provenance = Vec::new();
        let manifest = json!({"name": "panel"});

        let report = app_input_report_json(AppInputReportInput {
            assay_id: "assay",
            participant_id: "P001",
            input_file_name: "sample.txt",
            input_file_path: "/data/sample.txt",
            observations: &observations,
            analyses: &analyses,
            findings: &findings,
            provenance: &provenance,
            input_inspection: None,
            manifest_metadata: &manifest,
        });

        assert_eq!(report["report_status"], "complete");
        assert_eq!(report["metrics"]["n_findings_matched"], 1);
        assert_eq!(
            report["findings"][0]["matched_effect"]["text"],
            "A allele observed"
        );
    }

    #[test]
    fn enum_name_helpers_cover_all_known_values() {
        assert_eq!(super::file_container_name(FileContainer::Plain), "plain");
        assert_eq!(super::file_container_name(FileContainer::Zip), "zip");
        assert_eq!(super::detected_kind_name(DetectedKind::GenotypeText), "genotype_text");
        assert_eq!(super::detected_kind_name(DetectedKind::Vcf), "vcf");
        assert_eq!(
            super::detected_kind_name(DetectedKind::AlignmentCram),
            "alignment_cram"
        );
        assert_eq!(
            super::detected_kind_name(DetectedKind::AlignmentBam),
            "alignment_bam"
        );
        assert_eq!(
            super::detected_kind_name(DetectedKind::ReferenceFasta),
            "reference_fasta"
        );
        assert_eq!(super::detected_kind_name(DetectedKind::Unknown), "unknown");
        assert_eq!(
            super::detection_confidence_name(DetectionConfidence::WeakHeuristic),
            "weak_heuristic"
        );
        assert_eq!(
            super::detection_confidence_name(DetectionConfidence::Unknown),
            "unknown"
        );
        assert_eq!(super::assembly_name(Assembly::Grch37), "grch37");
        assert_eq!(super::assembly_name(Assembly::Grch38), "grch38");
        assert_eq!(super::inferred_sex_name(InferredSex::Male), "male");
        assert_eq!(super::inferred_sex_name(InferredSex::Female), "female");
        assert_eq!(super::inferred_sex_name(InferredSex::Unknown), "unknown");
        assert_eq!(
            super::sex_detection_confidence_name(SexDetectionConfidence::Medium),
            "medium"
        );
        assert_eq!(
            super::sex_detection_confidence_name(SexDetectionConfidence::Low),
            "low"
        );
    }
}
