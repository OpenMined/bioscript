use std::collections::BTreeMap;

use bioscript_core::{Assembly, GenomicLocus};
use bioscript_formats::SexInference;
use bioscript_schema::VariantManifest;

mod facets;
mod genotype_display;

use facets::{
    classify_non_reportable_alleles, is_weak_delimited_indel_match, observation_facets,
    parse_optional_u32,
};
use genotype_display::{
    assembly_row_value, deletion_copy_number_display, genotype_display_from_raw_counts,
    hemizygous_display_genotype, normalize_app_genotype, observation_evidence_raw,
};

pub struct AppObservationInput<'a> {
    pub row: &'a BTreeMap<String, String>,
    pub row_path: &'a str,
    pub assay_id: &'a str,
    pub manifest: VariantManifest,
    pub gene: String,
    pub source: serde_json::Value,
    pub observed_alt_alleles: Vec<String>,
    pub inferred_sex: Option<&'a SexInference>,
    pub fallback_assembly: Option<Assembly>,
}

struct AppObservationJson {
    allele_balance: Option<f64>,
    alt_count: Option<u32>,
    assay_id: String,
    assembly: String,
    call: ObservationCallValues,
    chrom: String,
    depth: Option<u32>,
    evidence_raw: String,
    gene: String,
    genotype: String,
    genotype_display: String,
    kind: String,
    locus: Option<GenomicLocus>,
    manifest: VariantManifest,
    non_reportable_status: Option<&'static str>,
    observed_alt_alleles: Vec<String>,
    ref_allele: String,
    ref_count: Option<u32>,
    reportable_alt: String,
    row: BTreeMap<String, String>,
    row_path: String,
    source: serde_json::Value,
    weak_indel_match: bool,
    zygosity: String,
}

struct ObservationCallValues {
    outcome: &'static str,
    status: &'static str,
    reported_genotype_display: String,
}

pub fn app_observation_from_manifest_row(input: AppObservationInput<'_>) -> serde_json::Value {
    let AppObservationInput {
        row,
        row_path,
        assay_id,
        manifest,
        gene,
        source,
        observed_alt_alleles,
        inferred_sex,
        fallback_assembly,
    } = input;
    let ref_allele = manifest.spec.reference.clone().unwrap_or_default();
    let reportable_alt = manifest.spec.alternate.clone().unwrap_or_default();
    let mut genotype_display = row
        .get("genotype")
        .filter(|value| !value.is_empty())
        .cloned()
        .or_else(|| genotype_display_from_raw_counts(row.get("raw_counts")?))
        .unwrap_or_default();
    let depth = parse_optional_u32(row.get("depth"));
    let ref_count = parse_optional_u32(row.get("ref_count"));
    let alt_count = parse_optional_u32(row.get("alt_count"));
    if let Some(normalized_display) = deletion_copy_number_display(row, &manifest, depth, alt_count)
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
    let non_reportable_status = classify_non_reportable_alleles(
        &genotype_display,
        &ref_allele,
        &reportable_alt,
        &observed_alt_alleles,
    );
    let call = observation_call_values(
        depth,
        non_reportable_status,
        &genotype,
        &zygosity,
        &genotype_display,
    );
    let evidence_raw = observation_evidence_raw(row, &chrom, inferred_sex);
    let kind = manifest.spec.kind.map_or("unknown".to_owned(), |kind| {
        format!("{kind:?}").to_lowercase()
    });
    render_app_observation_json(AppObservationJson {
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
        row_path: row_path.to_owned(),
        source,
        weak_indel_match,
        zygosity,
    })
}

fn observation_call_values(
    depth: Option<u32>,
    non_reportable_status: Option<&'static str>,
    genotype: &str,
    zygosity: &str,
    genotype_display: &str,
) -> ObservationCallValues {
    let outcome = if depth == Some(0) {
        "not_covered"
    } else if non_reportable_status == Some("observed_alt") {
        "observed_alt"
    } else if non_reportable_status == Some("unknown_alt") {
        "unknown_alt"
    } else if genotype == "./." {
        "no_call"
    } else if zygosity == "hom_ref" || zygosity == "hem_ref" {
        "reference"
    } else if zygosity == "het" || zygosity == "hom_alt" || zygosity == "hem_alt" {
        "variant"
    } else {
        "unknown"
    };
    let status = if matches!(outcome, "observed_alt" | "unknown_alt") {
        outcome
    } else if genotype == "./." {
        "no_call"
    } else {
        "called"
    };
    let reported_genotype_display = if matches!(zygosity, "hem_ref" | "hem_alt") {
        hemizygous_display_genotype(genotype_display)
    } else if genotype_display.is_empty() && matches!(outcome, "no_call" | "not_covered") {
        "??".to_owned()
    } else {
        genotype_display.to_owned()
    };
    ObservationCallValues {
        outcome,
        status,
        reported_genotype_display,
    }
}

fn render_app_observation_json(input: AppObservationJson) -> serde_json::Value {
    let AppObservationJson {
        allele_balance,
        alt_count,
        assay_id,
        assembly,
        call,
        chrom,
        depth,
        evidence_raw,
        gene,
        genotype,
        genotype_display,
        kind,
        locus,
        manifest,
        non_reportable_status,
        observed_alt_alleles,
        ref_allele,
        ref_count,
        reportable_alt,
        row,
        row_path,
        source,
        weak_indel_match,
        zygosity,
    } = input;
    let gene = if gene.is_empty() {
        manifest_gene_from_tags(&manifest).unwrap_or_default()
    } else {
        gene
    };
    let source = if source.is_null() {
        manifest_default_source(&row, &manifest)
    } else {
        source
    };
    serde_json::json!({
        "participant_id": row.get("participant_id").cloned().unwrap_or_default(),
        "assay_id": assay_id,
        "assay_version": "1.0",
        "variant_key": manifest.name,
        "variant_path": row_path,
        "rsid": row.get("matched_rsid").filter(|value| !value.is_empty()).cloned().or_else(|| manifest.spec.rsids.first().cloned()),
        "gene": gene,
        "assembly": if assembly.is_empty() { serde_json::Value::Null } else { serde_json::Value::String(assembly.to_uppercase()) },
        "chrom": chrom,
        "pos_start": locus.as_ref().map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.start)),
        "pos_end": locus.as_ref().map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.end)),
        "ref": ref_allele,
        "alt": reportable_alt,
        "kind": kind,
        "match_status": if row.get("matched_rsid").is_some_and(|value| !value.is_empty()) || !genotype_display.is_empty() { "found" } else { "not_found" },
        "coverage_status": depth.map_or("covered", |depth| if depth > 0 { "covered" } else { "not_covered" }),
        "call_status": call.status,
        "genotype": genotype,
        "genotype_display": call.reported_genotype_display,
        "zygosity": zygosity,
        "ref_count": ref_count,
        "alt_count": alt_count,
        "depth": depth,
        "genotype_quality": serde_json::Value::Null,
        "allele_balance": allele_balance,
        "outcome": call.outcome,
        "evidence_type": if row.get("backend").is_some_and(|value| value == "cram") { "mpileup" } else { "genotype_file" },
        "evidence_raw": evidence_raw,
        "source": source,
        "match_quality": if weak_indel_match { serde_json::Value::String("weak".to_owned()) } else { serde_json::Value::Null },
        "match_notes": if weak_indel_match {
            serde_json::Value::String("consumer genotype file reported an insertion/deletion token at the marker, not sequence-resolved evidence for the exact deletion allele".to_owned())
        } else {
            serde_json::Value::Null
        },
        "facets": observation_facets(non_reportable_status, &observed_alt_alleles),
    })
}

fn manifest_gene_from_tags(manifest: &VariantManifest) -> Option<String> {
    manifest.tags.iter().find_map(|tag| {
        tag.strip_prefix("gene:")
            .filter(|gene| !gene.is_empty())
            .map(ToOwned::to_owned)
    })
}

fn manifest_default_source(
    row: &BTreeMap<String, String>,
    manifest: &VariantManifest,
) -> serde_json::Value {
    let rsid = row
        .get("matched_rsid")
        .filter(|rsid| !rsid.is_empty())
        .or_else(|| row.get("rsid").filter(|rsid| !rsid.is_empty()))
        .or_else(|| manifest.spec.rsids.first().filter(|rsid| !rsid.is_empty()));
    let Some(rsid) = rsid else {
        return serde_json::Value::Null;
    };
    serde_json::json!({
        "kind": "database",
        "label": "dbSNP / NCBI SNP",
        "url": format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"),
        "fields": ["identifiers.rsids"],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioscript_core::{VariantKind, VariantSpec};
    use bioscript_formats::{InferredSex, SexDetectionConfidence};
    use std::path::PathBuf;

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
            normalize_app_genotype("G", "C", "G", None, "X", None),
            ("1".to_owned(), "hem_alt".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("C", "C", "G", None, "chrX", None),
            ("0".to_owned(), "hem_ref".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("G", "C", "G", None, "1", None),
            ("G".to_owned(), "unknown".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("GG", "C", "G", None, "X", None),
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
            normalize_app_genotype("GG", "C", "G", None, "X", Some(&inferred_sex)),
            ("1".to_owned(), "hem_alt".to_owned())
        );
        assert_eq!(
            normalize_app_genotype("CC", "C", "G", None, "chrX", Some(&inferred_sex)),
            ("0".to_owned(), "hem_ref".to_owned())
        );
    }

    fn manifest(
        kind: VariantKind,
        chrom: &str,
        reference: &str,
        alternate: &str,
    ) -> VariantManifest {
        VariantManifest {
            path: PathBuf::from("variants/rs1.yaml"),
            name: "rs1".to_owned(),
            tags: vec!["tag:test".to_owned()],
            spec: VariantSpec {
                rsids: vec!["rs1".to_owned()],
                grch37: Some(bioscript_core::GenomicLocus {
                    chrom: chrom.to_owned(),
                    start: 10,
                    end: 10,
                }),
                grch38: Some(bioscript_core::GenomicLocus {
                    chrom: chrom.to_owned(),
                    start: 20,
                    end: 20,
                }),
                reference: Some(reference.to_owned()),
                alternate: Some(alternate.to_owned()),
                kind: Some(kind),
                ..VariantSpec::default()
            },
        }
    }

    fn base_row() -> BTreeMap<String, String> {
        BTreeMap::from([
            ("participant_id".to_owned(), "p1".to_owned()),
            ("matched_rsid".to_owned(), "rs1".to_owned()),
            ("backend".to_owned(), "text".to_owned()),
            ("assembly".to_owned(), "grch38".to_owned()),
            ("genotype".to_owned(), "AG".to_owned()),
            ("ref_count".to_owned(), "8".to_owned()),
            ("alt_count".to_owned(), "7".to_owned()),
            ("depth".to_owned(), "15".to_owned()),
            ("evidence".to_owned(), "fixture".to_owned()),
        ])
    }

    #[test]
    fn app_observation_json_covers_called_variant_fields() {
        let row = base_row();
        let observation = app_observation_from_manifest_row(AppObservationInput {
            row: &row,
            row_path: "variants/rs1.yaml",
            assay_id: "assay",
            manifest: manifest(VariantKind::Snp, "1", "A", "G"),
            gene: "ABC".to_owned(),
            source: serde_json::json!({"kind": "database"}),
            observed_alt_alleles: Vec::new(),
            inferred_sex: None,
            fallback_assembly: None,
        });

        assert_eq!(observation["participant_id"], "p1");
        assert_eq!(observation["assay_id"], "assay");
        assert_eq!(observation["variant_key"], "rs1");
        assert_eq!(observation["assembly"], "GRCH38");
        assert_eq!(observation["pos_start"], 20);
        assert_eq!(observation["genotype"], "0/1");
        assert_eq!(observation["zygosity"], "het");
        assert_eq!(observation["outcome"], "variant");
        assert_eq!(observation["allele_balance"], serde_json::json!(7.0 / 15.0));
        assert_eq!(observation["facets"], serde_json::Value::Null);
    }

    #[test]
    fn app_observation_json_covers_no_call_not_covered_and_fallback_assembly() {
        let mut no_call = BTreeMap::from([
            ("participant_id".to_owned(), "p2".to_owned()),
            ("backend".to_owned(), "text".to_owned()),
            ("depth".to_owned(), "0".to_owned()),
        ]);
        let observation = app_observation_from_manifest_row(AppObservationInput {
            row: &no_call,
            row_path: "variants/rs1.yaml",
            assay_id: "assay",
            manifest: manifest(VariantKind::Snp, "1", "A", "G"),
            gene: "ABC".to_owned(),
            source: serde_json::Value::Null,
            observed_alt_alleles: Vec::new(),
            inferred_sex: None,
            fallback_assembly: Some(Assembly::Grch37),
        });
        assert_eq!(observation["assembly"], "GRCH37");
        assert_eq!(observation["pos_start"], 10);
        assert_eq!(observation["match_status"], "not_found");
        assert_eq!(observation["coverage_status"], "not_covered");
        assert_eq!(observation["call_status"], "no_call");
        assert_eq!(observation["genotype_display"], "??");
        assert_eq!(observation["outcome"], "not_covered");

        no_call.insert("depth".to_owned(), "12".to_owned());
        let observation = app_observation_from_manifest_row(AppObservationInput {
            row: &no_call,
            row_path: "variants/rs1.yaml",
            assay_id: "assay",
            manifest: manifest(VariantKind::Snp, "1", "A", "G"),
            gene: "ABC".to_owned(),
            source: serde_json::Value::Null,
            observed_alt_alleles: Vec::new(),
            inferred_sex: None,
            fallback_assembly: Some(Assembly::Grch38),
        });
        assert_eq!(observation["outcome"], "no_call");
    }

    #[test]
    fn app_observation_json_covers_non_reportable_and_sex_evidence_paths() {
        let inferred_sex = SexInference {
            sex: InferredSex::Male,
            confidence: SexDetectionConfidence::Medium,
            method: "fixture".to_owned(),
            evidence: vec!["signal=present".to_owned()],
        };
        let mut row = base_row();
        row.insert("genotype".to_owned(), "TT".to_owned());
        let observation = app_observation_from_manifest_row(AppObservationInput {
            row: &row,
            row_path: "variants/rsx.yaml",
            assay_id: "assay",
            manifest: manifest(VariantKind::Snp, "X", "A", "G"),
            gene: "ABC".to_owned(),
            source: serde_json::Value::Null,
            observed_alt_alleles: vec!["T".to_owned()],
            inferred_sex: Some(&inferred_sex),
            fallback_assembly: None,
        });

        assert_eq!(observation["outcome"], "observed_alt");
        assert_eq!(observation["call_status"], "observed_alt");
        assert_eq!(observation["facets"], "observed_alt;known_observed_alts=T");
        assert!(
            observation["evidence_raw"]
                .as_str()
                .unwrap()
                .contains("detected_sex=male")
        );
    }

    #[test]
    fn app_observation_json_covers_raw_counts_and_weak_indel_match() {
        let mut row = BTreeMap::from([
            ("participant_id".to_owned(), "p3".to_owned()),
            ("matched_rsid".to_owned(), "rs1".to_owned()),
            ("backend".to_owned(), "zip".to_owned()),
            ("raw_counts".to_owned(), r#"{"D": 8, "I": 6}"#.to_owned()),
            ("depth".to_owned(), "14".to_owned()),
        ]);
        row.insert("genotype".to_owned(), "ID".to_owned());
        let observation = app_observation_from_manifest_row(AppObservationInput {
            row: &row,
            row_path: "variants/rs1.yaml",
            assay_id: "assay",
            manifest: manifest(VariantKind::Deletion, "22", "TTATAA", "<DEL:6>"),
            gene: "APOL1".to_owned(),
            source: serde_json::Value::Null,
            observed_alt_alleles: Vec::new(),
            inferred_sex: None,
            fallback_assembly: Some(Assembly::Grch38),
        });
        assert_eq!(observation["kind"], "deletion");
        assert_eq!(observation["genotype"], "0/1");
        assert_eq!(observation["match_quality"], "weak");
        assert!(
            observation["match_notes"]
                .as_str()
                .unwrap()
                .contains("insertion/deletion token")
        );
    }
}
