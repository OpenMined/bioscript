use std::io::{Read, Seek};

use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::csi::{self, BinningIndex};
use noodles::tabix;

use bioscript_core::{Assembly, GenomicLocus, RuntimeError, VariantObservation};

use super::parse_vcf_record;

/// Observe a SNP at `locus` over an already-built tabix-indexed bgzipped VCF
/// reader. Caller builds `csi::io::IndexedReader::new(reader, tabix_index)`
/// once and calls this per variant.
pub fn observe_vcf_snp_with_reader<R>(
    indexed: &mut csi::io::IndexedReader<bgzf::io::Reader<R>, tabix::Index>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError>
where
    R: Read + Seek,
{
    let locus_label = format!("{}:{}", locus.chrom, locus.start);

    let Some(seq_name) = resolve_vcf_chrom_name(indexed.index(), &locus.chrom) else {
        return Ok(VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid,
            assembly,
            evidence: vec![format!(
                "{label}: tabix index has no contig matching {} (tried chr-prefixed and bare forms)",
                locus.chrom
            )],
            ..VariantObservation::default()
        });
    };

    let pos_usize = usize::try_from(locus.start).map_err(|err| {
        RuntimeError::Io(format!(
            "{label}: invalid VCF position {} for {locus_label}: {err}",
            locus.start
        ))
    })?;
    let position = Position::try_from(pos_usize).map_err(|err| {
        RuntimeError::Io(format!(
            "{label}: invalid VCF position {} for {locus_label}: {err}",
            locus.start
        ))
    })?;
    let region = Region::new(seq_name.as_str(), position..=position);

    let query = indexed.query(&region).map_err(|err| {
        RuntimeError::Io(format!("{label}: tabix query for {locus_label}: {err}"))
    })?;

    let reference_str = reference.to_ascii_uppercase().to_string();
    let alternate_str = alternate.to_ascii_uppercase().to_string();

    let mut saw_any = false;
    for record_result in query {
        let record = record_result
            .map_err(|err| RuntimeError::Io(format!("{label}: tabix record iter: {err}")))?;
        let line: &str = record.as_ref();
        let Some(row) = parse_vcf_record(line)? else {
            continue;
        };
        if row.position != locus.start {
            continue;
        }
        saw_any = true;
        if !row.reference.eq_ignore_ascii_case(&reference_str) {
            continue;
        }
        if !row
            .alternates
            .iter()
            .any(|alt| alt.eq_ignore_ascii_case(&alternate_str))
        {
            continue;
        }

        return Ok(VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid: matched_rsid.or_else(|| row.rsid.clone()),
            assembly,
            genotype: Some(row.genotype.clone()),
            evidence: vec![format!("{label}: resolved by locus {locus_label}")],
            ..VariantObservation::default()
        });
    }

    let evidence = if saw_any {
        vec![format!(
            "{label}: {locus_label} present but ref={reference}/alt={alternate} did not match any record"
        )]
    } else {
        vec![format!("{label}: no VCF record at {locus_label}")]
    };
    Ok(VariantObservation {
        backend: "vcf".to_owned(),
        matched_rsid,
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
