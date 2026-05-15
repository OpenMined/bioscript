use std::collections::BTreeMap;

use bioscript_core::{Assembly, VariantKind};
use bioscript_formats::{InferredSex, SexDetectionConfidence, SexInference};
use bioscript_schema::VariantManifest;

pub(super) fn assembly_row_value(assembly: Assembly) -> String {
    match assembly {
        Assembly::Grch37 => "grch37".to_owned(),
        Assembly::Grch38 => "grch38".to_owned(),
    }
}

pub(super) fn hemizygous_display_genotype(display: &str) -> String {
    display
        .chars()
        .find(char::is_ascii_alphabetic)
        .map_or_else(|| display.to_owned(), |allele| allele.to_string())
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
    if !matches!(row.get("backend").map(String::as_str), Some("cram" | "bam")) {
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

pub(super) fn observation_evidence_raw(
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

pub(super) fn genotype_display_from_raw_counts(raw_counts: &str) -> Option<String> {
    let counts: serde_json::Map<String, serde_json::Value> =
        serde_json::from_str(raw_counts).ok()?;
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
