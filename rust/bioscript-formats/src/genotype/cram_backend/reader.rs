use std::{
    collections::{BTreeMap, BTreeSet},
    io::{Read, Seek},
};

use noodles::cram;

use bioscript_core::{Assembly, GenomicLocus, RuntimeError, VariantObservation};

use crate::alignment;

use super::{
    classify_expected_indel, describe_copy_number_decision_rule, describe_locus,
    describe_snp_decision_rule, infer_copy_number_genotype, infer_snp_genotype,
    record_overlaps_locus, snp_pileup_with_reader,
};

/// Observe a SNP at `locus` over an already-built CRAM `IndexedReader` and
/// reference repository (held by the reader). Mirrors the internal
/// `CramBackend::observe_snp` but reader-based, so non-filesystem callers
/// (e.g. wasm with a JS-backed reader) don't need a `GenotypeStore` or paths.
///
/// `matched_rsid` and `assembly` are passed through to the returned
/// observation unchanged — callers that already know them (e.g. from
/// compiling a YAML variant) should supply them; otherwise `None`.
pub fn observe_cram_snp_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    let pileup = snp_pileup_with_reader(reader, label, locus, reference, alternate, false)?;
    let ref_count = pileup.filtered_ref_count;
    let alt_count = pileup.filtered_alt_count;
    let depth = pileup.filtered_depth;
    let evidence = pileup.evidence_lines(&describe_locus(locus), locus.start);

    Ok(VariantObservation {
        backend: "cram".to_owned(),
        matched_rsid,
        assembly,
        genotype: infer_snp_genotype(reference, alternate, ref_count, alt_count, depth),
        ref_count: Some(ref_count),
        alt_count: Some(alt_count),
        depth: Some(depth),
        raw_counts: pileup.raw_base_counts,
        decision: Some(describe_snp_decision_rule(
            reference, alternate, ref_count, alt_count, depth,
        )),
        evidence,
    })
}

/// Observe an insertion/indel-like variant at `locus` over an already-built
/// CRAM `IndexedReader`.
pub fn observe_cram_indel_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    reference: &str,
    alternate: &str,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    let mut alt_count = 0u32;
    let mut ref_count = 0u32;
    let mut depth = 0u32;
    let mut matching_alt_lengths = BTreeSet::new();

    alignment::for_each_cram_record_with_reader(reader, label, locus, |record| {
        if record.is_unmapped || !record_overlaps_locus(&record, locus) {
            return Ok(true);
        }
        let classification = classify_expected_indel(&record, locus, reference.len(), alternate)?;
        if !classification.covering {
            return Ok(true);
        }
        depth += 1;
        if classification.matches_alt {
            alt_count += 1;
            matching_alt_lengths.insert(classification.observed_len);
        } else if classification.reference_like {
            ref_count += 1;
        }
        Ok(true)
    })?;

    let evidence_label = if matching_alt_lengths.is_empty() {
        "none".to_owned()
    } else {
        matching_alt_lengths
            .into_iter()
            .map(|len| len.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };

    Ok(VariantObservation {
        backend: "cram".to_owned(),
        matched_rsid,
        assembly,
        genotype: infer_copy_number_genotype(reference, alternate, ref_count, alt_count, depth),
        ref_count: Some(ref_count),
        alt_count: Some(alt_count),
        depth: Some(depth),
        raw_counts: BTreeMap::new(),
        decision: Some(describe_copy_number_decision_rule(
            reference, alternate, ref_count, alt_count, depth,
        )),
        evidence: vec![format!(
            "observed indel at {} depth={} ref_count={} alt_count={} matching_alt_lengths={}",
            describe_locus(locus),
            depth,
            ref_count,
            alt_count,
            evidence_label
        )],
    })
}
