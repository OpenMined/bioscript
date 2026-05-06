use bioscript_core::{GenomicLocus, VariantKind, VariantSpec};
use serde_yaml::Value;

use super::common::{i64_at_mapping, mapping_at, scalar_at, seq_of_strings, value_at};

pub(crate) fn variant_spec_from_root(root: &Value) -> Result<VariantSpec, String> {
    let rsids = seq_of_strings(root, &["identifiers", "rsids"]).unwrap_or_default();
    let grch37 = locus_from_root(root, "grch37")?;
    let grch38 = locus_from_root(root, "grch38")?;
    let reference = scalar_at(root, &["alleles", "ref"]);
    let alternate = preferred_alternate_from_root(root);
    let deletion_length = value_at(root, &["alleles", "deletion_length"])
        .and_then(Value::as_u64)
        .and_then(|value| usize::try_from(value).ok());
    let motifs = seq_of_strings(root, &["alleles", "motifs"]).unwrap_or_default();
    let kind = scalar_at(root, &["alleles", "kind"]).map(|kind| match kind.as_str() {
        "snv" => VariantKind::Snp,
        "deletion" => VariantKind::Deletion,
        "insertion" => VariantKind::Insertion,
        "indel" => VariantKind::Indel,
        _ => VariantKind::Other,
    });

    Ok(VariantSpec {
        rsids,
        grch37,
        grch38,
        reference,
        alternate,
        kind,
        deletion_length,
        motifs,
    })
}

fn preferred_alternate_from_root(root: &Value) -> Option<String> {
    let alts = seq_of_strings(root, &["alleles", "alts"])?;
    if let Some(finding_alt) = first_specific_finding_alt(root)
        && alts.iter().any(|alt| alt == &finding_alt)
    {
        return Some(finding_alt);
    }
    alts.first().cloned()
}

fn first_specific_finding_alt(root: &Value) -> Option<String> {
    let findings = value_at(root, &["findings"])?.as_sequence()?;
    for finding in findings {
        let Some(alt) = finding
            .as_mapping()
            .and_then(|mapping| mapping.get(Value::String("alt".to_owned())))
            .and_then(Value::as_str)
            .map(str::trim)
        else {
            continue;
        };
        if !alt.is_empty() && alt != "*" {
            return Some(alt.to_owned());
        }
    }
    None
}

fn locus_from_root(root: &Value, assembly: &str) -> Result<Option<GenomicLocus>, String> {
    let Some(mapping) = mapping_at(root, &["coordinates", assembly]) else {
        return Ok(None);
    };
    let chrom = mapping
        .get(Value::String("chrom".to_owned()))
        .and_then(Value::as_str)
        .ok_or_else(|| format!("coordinates.{assembly}.chrom missing"))?;
    let (start, end) = if let Some(pos) = i64_at_mapping(mapping, "pos") {
        (pos, pos)
    } else {
        let start = i64_at_mapping(mapping, "start")
            .ok_or_else(|| format!("coordinates.{assembly}.start missing"))?;
        let end = i64_at_mapping(mapping, "end")
            .ok_or_else(|| format!("coordinates.{assembly}.end missing"))?;
        (start, end)
    };
    Ok(Some(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    }))
}
