use std::collections::BTreeMap;

use bioscript_core::VariantKind;
use bioscript_schema::VariantManifest;

pub(super) fn classify_non_reportable_alleles(
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

pub(super) fn is_weak_delimited_indel_match(
    row: &BTreeMap<String, String>,
    manifest: &VariantManifest,
    genotype_display: &str,
) -> bool {
    if !matches!(manifest.spec.kind, Some(VariantKind::Deletion)) {
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

pub(super) fn observation_facets(
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

pub(super) fn parse_optional_u32(value: Option<&String>) -> Option<u32> {
    value.and_then(|value| value.parse::<u32>().ok())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bioscript_core::VariantSpec;
    use std::path::PathBuf;

    fn deletion_manifest(reference: &str) -> VariantManifest {
        VariantManifest {
            path: PathBuf::from("variant.yaml"),
            name: "Deletion".to_owned(),
            tags: Vec::new(),
            spec: VariantSpec {
                kind: Some(VariantKind::Deletion),
                reference: Some(reference.to_owned()),
                alternate: Some("<DEL>".to_owned()),
                ..VariantSpec::default()
            },
        }
    }

    #[test]
    fn non_reportable_allele_classifier_distinguishes_known_and_unknown_alts() {
        assert_eq!(
            classify_non_reportable_alleles("A/T", "A", "G", &["T".to_owned()]),
            Some("observed_alt")
        );
        assert_eq!(
            classify_non_reportable_alleles("A/C", "A", "G", &["T".to_owned()]),
            Some("unknown_alt")
        );
        assert_eq!(
            classify_non_reportable_alleles("", "A", "G", &["T".to_owned()]),
            None
        );
        assert_eq!(
            classify_non_reportable_alleles("A/T", "AT", "G", &["T".to_owned()]),
            None
        );
        assert_eq!(
            classify_non_reportable_alleles("A/T", "A", "GT", &["T".to_owned()]),
            None
        );
    }

    #[test]
    fn weak_delimited_indel_match_requires_text_deletion_shape() {
        let manifest = deletion_manifest("TTATAA");
        let mut row = BTreeMap::new();
        row.insert("backend".to_owned(), "text".to_owned());
        assert!(is_weak_delimited_indel_match(&row, &manifest, "ID"));
        assert!(is_weak_delimited_indel_match(&row, &manifest, "D/D"));
        assert!(!is_weak_delimited_indel_match(&row, &manifest, "AG"));

        row.insert("backend".to_owned(), "cram".to_owned());
        assert!(!is_weak_delimited_indel_match(&row, &manifest, "ID"));

        row.insert("backend".to_owned(), "zip".to_owned());
        assert!(!is_weak_delimited_indel_match(
            &row,
            &deletion_manifest("A"),
            "ID"
        ));
        let mut snv = deletion_manifest("TTATAA");
        snv.spec.kind = Some(VariantKind::Snp);
        assert!(!is_weak_delimited_indel_match(&row, &snv, "ID"));
    }

    #[test]
    fn observation_facets_and_optional_integer_parsing_cover_edges() {
        assert_eq!(observation_facets(None, &[]), serde_json::Value::Null);
        assert_eq!(
            observation_facets(Some("unknown_alt"), &[]),
            serde_json::Value::String("unknown_alt".to_owned())
        );
        assert_eq!(
            observation_facets(Some("observed_alt"), &["T".to_owned(), "C".to_owned()]),
            serde_json::Value::String("observed_alt;known_observed_alts=T,C".to_owned())
        );

        let good = "42".to_owned();
        let bad = "not-a-number".to_owned();
        assert_eq!(parse_optional_u32(Some(&good)), Some(42));
        assert_eq!(parse_optional_u32(Some(&bad)), None);
        assert_eq!(parse_optional_u32(None), None);
    }
}
