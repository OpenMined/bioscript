use bioscript_core::VariantSpec;

pub(crate) fn describe_query(variant: &VariantSpec) -> &'static str {
    if variant.has_coordinates() {
        "variant_by_locus"
    } else {
        "variant_by_rsid"
    }
}

pub(crate) fn variant_sort_key(variant: &VariantSpec) -> (u8, String, i64, i64, String) {
    if let Some(locus) = &variant.grch38 {
        return (
            0,
            chrom_sort_key(&locus.chrom),
            locus.start,
            locus.end,
            variant.rsids.first().cloned().unwrap_or_default(),
        );
    }
    if let Some(locus) = &variant.grch37 {
        return (
            1,
            chrom_sort_key(&locus.chrom),
            locus.start,
            locus.end,
            variant.rsids.first().cloned().unwrap_or_default(),
        );
    }
    (
        2,
        "~".to_owned(),
        i64::MAX,
        i64::MAX,
        variant.rsids.first().cloned().unwrap_or_default(),
    )
}

pub(crate) fn chrom_sort_key(raw: &str) -> String {
    let chrom = raw.trim().strip_prefix("chr").unwrap_or(raw.trim());
    if let Ok(value) = chrom.parse::<u32>() {
        return format!("{value:03}");
    }
    match chrom.to_ascii_uppercase().as_str() {
        "X" => "023".to_owned(),
        "Y" => "024".to_owned(),
        "M" | "MT" => "025".to_owned(),
        other => format!("999-{other}"),
    }
}

pub(crate) fn normalize_genotype(value: &str) -> String {
    let cleaned = value.trim().replace(' ', "").to_ascii_uppercase();
    if cleaned.is_empty() || matches!(cleaned.as_str(), "NA" | "N/A" | "#N/A" | "NONE") {
        return "--".to_owned();
    }
    if cleaned.contains('/') {
        let parts: Vec<&str> = cleaned.split('/').collect();
        if parts.iter().any(|part| part.is_empty() || *part == "-") {
            return "ID".to_owned();
        }
        return parts.concat();
    }
    cleaned
}
