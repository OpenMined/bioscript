pub(super) fn normalize_chrom(value: &str) -> String {
    let normalized = value
        .trim()
        .trim_start_matches("chr")
        .trim_start_matches("CHR")
        .to_ascii_uppercase();
    match normalized.as_str() {
        "23" => "X".to_owned(),
        "24" => "Y".to_owned(),
        "25" => "XY".to_owned(),
        "26" | "M" => "MT".to_owned(),
        _ => normalized,
    }
}

pub(super) fn is_called_genotype_text(value: &str) -> bool {
    let value = value.trim();
    if value.is_empty() || matches!(value, "--" | "00" | "." | "./." | ".|.") {
        return false;
    }
    value
        .chars()
        .all(|ch| matches!(ch.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'))
}

pub(super) fn genotype_allele_count(value: &str) -> usize {
    value
        .chars()
        .filter(|ch| matches!(ch.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'))
        .count()
}

pub(super) fn is_genotype_text_het(value: &str) -> bool {
    let alleles: Vec<char> = value
        .chars()
        .filter(|ch| matches!(ch.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'))
        .map(|ch| ch.to_ascii_uppercase())
        .collect();
    alleles.len() == 2 && alleles[0] != alleles[1]
}

pub(super) fn is_called_vcf_gt(value: &str) -> bool {
    let value = value.trim();
    !value.is_empty()
        && !value.contains('.')
        && (value != "0" || matches!(value, "0" | "1" | "2" | "3"))
}

pub(super) fn vcf_gt_allele_count(gt: &str) -> usize {
    gt.split(['/', '|'])
        .filter(|part| !part.is_empty() && *part != ".")
        .count()
}

pub(super) fn is_vcf_gt_het(gt: &str) -> bool {
    let alleles: Vec<&str> = gt
        .split(['/', '|'])
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    alleles.len() == 2 && alleles[0] != alleles[1]
}

pub(super) fn is_non_par_x(pos: u32) -> bool {
    // Human GRCh38 non-PAR X used by bcftools +guess-ploidy.
    // GRCh37 differs slightly, but these bounds cover the common non-PAR body
    // and avoid both pseudoautosomal ends for this QC heuristic.
    (2_781_480..=154_931_043).contains(&pos)
}
