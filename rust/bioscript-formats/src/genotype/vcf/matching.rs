use bioscript_core::{Assembly, GenomicLocus, VariantKind, VariantObservation, VariantSpec};

use crate::inspect::InferredSex;

use super::ParsedVcfRow;

pub fn choose_variant_locus_for_assembly(
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> Option<GenomicLocus> {
    match assembly {
        Some(Assembly::Grch37) => variant.grch37.clone().or_else(|| variant.grch38.clone()),
        Some(Assembly::Grch38) => variant.grch38.clone().or_else(|| variant.grch37.clone()),
        None => variant.grch37.clone().or_else(|| variant.grch38.clone()),
    }
}

pub(super) fn first_single_base_allele(value: Option<&str>) -> Option<char> {
    let value = value?;
    let mut chars = value.chars();
    let base = chars.next()?;
    chars.next().is_none().then_some(base)
}

pub fn imputed_reference_observation(
    backend_name: &str,
    label: &str,
    variant: &VariantSpec,
    locus: &GenomicLocus,
    assembly: Option<Assembly>,
    inferred_sex: Option<InferredSex>,
    missing_evidence: &str,
) -> Option<VariantObservation> {
    let genotype = match variant.kind {
        None | Some(VariantKind::Snp) => {
            let reference = reference_for_assembly(variant, assembly)?;
            first_single_base_allele(variant.alternate.as_deref())?;
            reference_genotype_for_locus(reference, locus, inferred_sex)
        }
        Some(VariantKind::Deletion) => "II".to_owned(),
        Some(VariantKind::Insertion | VariantKind::Indel) => {
            let reference = variant.reference.as_deref()?;
            let alternate = variant.alternate.as_deref()?;
            let alternates = [alternate];
            let token = super::super::vcf_tokens::vcf_reference_token(reference, &alternates);
            format!("{token}{token}")
        }
        _ => return None,
    };
    let evidence_prefix = if missing_evidence.contains(label) {
        missing_evidence.to_owned()
    } else {
        format!("{label}: {missing_evidence}")
    };
    Some(VariantObservation {
        backend: backend_name.to_owned(),
        matched_rsid: variant.rsids.first().cloned(),
        assembly,
        genotype: Some(genotype),
        evidence: vec![format!(
            "{evidence_prefix} | imputed reference genotype from absent variant-only VCF record"
        )],
        ..VariantObservation::default()
    })
}

fn reference_genotype_for_locus(
    reference: char,
    locus: &GenomicLocus,
    inferred_sex: Option<InferredSex>,
) -> String {
    let reference = reference.to_ascii_uppercase();
    if inferred_sex == Some(InferredSex::Male) && is_sex_chromosome(&locus.chrom) {
        reference.to_string()
    } else {
        format!("{reference}{reference}")
    }
}

pub(crate) fn reference_for_assembly(
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> Option<char> {
    let value = variant
        .assembly_reference(assembly)
        .or(variant.reference.as_deref())?;
    first_single_base_allele(Some(value))
}

fn is_sex_chromosome(chrom: &str) -> bool {
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

pub(crate) fn normalize_chromosome_name(value: &str) -> String {
    value.trim().trim_start_matches("chr").to_ascii_lowercase()
}

pub(crate) fn vcf_row_matches_variant(
    row: &ParsedVcfRow,
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> bool {
    let Some(locus) = choose_variant_locus_for_assembly(variant, assembly) else {
        return false;
    };

    if normalize_chromosome_name(&row.chrom) != normalize_chromosome_name(&locus.chrom) {
        return false;
    }

    match variant.kind.unwrap_or(VariantKind::Other) {
        VariantKind::Snp => {
            let expected_reference = variant
                .assembly_reference(assembly)
                .or(variant.reference.as_deref());
            row.position == locus.start
                && expected_reference
                    .is_none_or(|reference| reference.eq_ignore_ascii_case(&row.reference))
                && snp_row_has_catalog_allele(row, variant)
        }
        VariantKind::Deletion => {
            let expected_len = variant.deletion_length.unwrap_or(0);
            row.position == locus.start.saturating_sub(1)
                && row.alternates.iter().any(|alternate| {
                    let actual_len = row.reference.len().saturating_sub(alternate.len());
                    (expected_len == 0 || actual_len == expected_len)
                        && alternate.len() < row.reference.len()
                })
        }
        VariantKind::Insertion => insertion_row_matches_variant(row, variant, &locus),
        VariantKind::Indel => indel_row_matches_variant(row, variant, &locus),
        VariantKind::Other => row.position == locus.start,
    }
}

pub(crate) fn vcf_row_genotype_for_variant(row: &ParsedVcfRow, variant: &VariantSpec) -> String {
    if !matches!(
        variant.kind,
        Some(VariantKind::Deletion | VariantKind::Insertion | VariantKind::Indel)
    ) {
        return row.genotype.clone();
    }

    let parts: Vec<&str> = row.genotype.split('/').collect();
    if parts.len() <= 1 {
        return row.genotype.clone();
    }

    let alternates: Vec<&str> = row.alternates.iter().map(String::as_str).collect();
    let mut tokens = Vec::with_capacity(parts.len());
    for part in parts {
        if part.eq_ignore_ascii_case(&row.reference) {
            tokens.push(super::super::vcf_tokens::vcf_reference_token(
                &row.reference,
                &alternates,
            ));
        } else if let Some(alternate) = row
            .alternates
            .iter()
            .find(|alternate| alternate.eq_ignore_ascii_case(part))
        {
            tokens.push(super::super::vcf_tokens::vcf_alt_token(
                &row.reference,
                alternate,
            ));
        } else {
            return row.genotype.clone();
        }
    }

    if tokens
        .iter()
        .all(|token| token.chars().count() == 1 && token != "--")
    {
        return super::super::normalize_genotype(&tokens.join(""));
    }

    row.genotype.clone()
}

fn snp_row_has_catalog_allele(row: &ParsedVcfRow, variant: &VariantSpec) -> bool {
    let Some(alternate) = variant.alternate.as_ref() else {
        return true;
    };
    if row
        .alternates
        .iter()
        .any(|candidate| candidate.eq_ignore_ascii_case(alternate))
    {
        return true;
    }
    variant.reference.as_ref().is_some_and(|reference| {
        row.alternates
            .iter()
            .any(|candidate| candidate.eq_ignore_ascii_case(reference))
    })
}

fn insertion_row_matches_variant(
    row: &ParsedVcfRow,
    variant: &VariantSpec,
    locus: &GenomicLocus,
) -> bool {
    indel_position_matches(row.position, locus)
        && row.alternates.iter().any(|alternate| {
            alternate.len() > row.reference.len()
                && row_matches_catalog_alleles(row, alternate, variant)
        })
}

fn indel_row_matches_variant(
    row: &ParsedVcfRow,
    variant: &VariantSpec,
    locus: &GenomicLocus,
) -> bool {
    indel_position_matches(row.position, locus)
        && row.alternates.iter().any(|alternate| {
            alternate.len() != row.reference.len()
                && row_matches_catalog_alleles(row, alternate, variant)
        })
}

fn indel_position_matches(row_position: i64, locus: &GenomicLocus) -> bool {
    row_position == locus.start || row_position == locus.start.saturating_sub(1)
}

fn row_matches_catalog_alleles(row: &ParsedVcfRow, alternate: &str, variant: &VariantSpec) -> bool {
    if variant
        .reference
        .as_ref()
        .is_some_and(|reference| !reference.eq_ignore_ascii_case(&row.reference))
    {
        return false;
    }
    variant
        .alternate
        .as_ref()
        .is_none_or(|expected| expected.eq_ignore_ascii_case(alternate))
}
