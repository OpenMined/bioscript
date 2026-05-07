use std::fs::File;

use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::csi::{self, BinningIndex};
use noodles::tabix;

use bioscript_core::{Assembly, RuntimeError, VariantKind, VariantObservation, VariantSpec};

use super::{describe_query, parse_vcf_record, vcf_row_matches_variant};

pub(super) fn observe_indexed_vcf_variant(
    indexed: &mut csi::io::IndexedReader<bgzf::io::Reader<File>, tabix::Index>,
    label: &str,
    variant: &VariantSpec,
    locus: &bioscript_core::GenomicLocus,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    let query_pos = if matches!(
        variant.kind,
        Some(VariantKind::Deletion | VariantKind::Insertion | VariantKind::Indel)
    ) {
        locus.start.saturating_sub(1)
    } else {
        locus.start
    };
    let locus_label = format!("{}:{query_pos}", locus.chrom);
    let Some(seq_name) = resolve_vcf_chrom_name(indexed.index(), &locus.chrom) else {
        return Ok(VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly,
            evidence: vec![format!(
                "{label}: tabix index has no contig matching {} (tried chr-prefixed and bare forms)",
                locus.chrom
            )],
            ..VariantObservation::default()
        });
    };
    let pos_usize = usize::try_from(query_pos).map_err(|err| {
        RuntimeError::Io(format!(
            "{label}: invalid VCF position {query_pos} for {locus_label}: {err}"
        ))
    })?;
    let position = Position::try_from(pos_usize).map_err(|err| {
        RuntimeError::Io(format!(
            "{label}: invalid VCF position {query_pos} for {locus_label}: {err}"
        ))
    })?;
    let region = Region::new(seq_name.as_str(), position..=position);
    let query = indexed.query(&region).map_err(|err| {
        RuntimeError::Io(format!("{label}: tabix query for {locus_label}: {err}"))
    })?;

    let mut saw_any = false;
    for record_result in query {
        let record = record_result
            .map_err(|err| RuntimeError::Io(format!("{label}: tabix record iter: {err}")))?;
        let line: &str = record.as_ref();
        let Some(row) = parse_vcf_record(line)? else {
            continue;
        };
        saw_any = true;
        if vcf_row_matches_variant(&row, variant, assembly) {
            return Ok(VariantObservation {
                backend: "vcf".to_owned(),
                matched_rsid: variant.rsids.first().cloned().or_else(|| row.rsid.clone()),
                assembly,
                genotype: Some(row.genotype.clone()),
                evidence: vec![format!("{label}: resolved by indexed locus {locus_label}")],
                ..VariantObservation::default()
            });
        }
    }

    let evidence = if saw_any {
        vec![format!(
            "{label}: {locus_label} present but no record matched {}",
            describe_query(variant)
        )]
    } else {
        vec![format!("{label}: no VCF record at {locus_label}")]
    };
    Ok(VariantObservation {
        backend: "vcf".to_owned(),
        matched_rsid: variant.rsids.first().cloned(),
        assembly,
        evidence,
        ..VariantObservation::default()
    })
}

fn resolve_vcf_chrom_name(index: &tabix::Index, user_chrom: &str) -> Option<String> {
    let header = index.header()?;
    let names = header.reference_sequence_names();

    let trimmed = user_chrom.trim();
    let stripped = trimmed.strip_prefix("chr").unwrap_or(trimmed);

    let candidates = [
        trimmed.to_owned(),
        stripped.to_owned(),
        format!("chr{stripped}"),
    ];
    for cand in &candidates {
        if names.contains(cand.as_bytes()) {
            return Some(cand.clone());
        }
    }

    let target = stripped.to_ascii_lowercase();
    for name in names {
        let as_str = std::str::from_utf8(name.as_ref()).ok()?;
        let as_stripped = as_str.strip_prefix("chr").unwrap_or(as_str);
        if as_stripped.eq_ignore_ascii_case(&target) {
            return Some(as_str.to_owned());
        }
    }
    None
}
