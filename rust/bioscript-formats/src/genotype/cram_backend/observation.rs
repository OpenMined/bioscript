use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};
use noodles::cram;

use crate::alignment::{self, AlignmentOpKind};
use crate::genotype::types::CramBackend;

use super::{
    anchor_window, classify_expected_indel, describe_copy_number_decision_rule, describe_locus,
    describe_snp_decision_rule, first_base, indel_at_anchor, infer_copy_number_genotype,
    infer_snp_genotype, observe_snp_pileup, record_overlaps_locus, snp_pileup_with_reader,
    spans_position,
};

impl CramBackend {
    pub(super) fn observe_with_reader(
        &self,
        reader: &mut cram::io::indexed_reader::IndexedReader<std::fs::File>,
        label: &str,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
    ) -> Result<VariantObservation, RuntimeError> {
        match variant.kind.unwrap_or(VariantKind::Other) {
            VariantKind::Snp => {
                self.observe_snp_with_reader(reader, label, variant, assembly, locus)
            }
            VariantKind::Deletion => {
                self.observe_deletion_with_reader(reader, label, variant, assembly, locus)
            }
            VariantKind::Insertion | VariantKind::Indel => {
                self.observe_indel_with_reader(reader, label, variant, assembly, locus)
            }
            VariantKind::Other => Err(RuntimeError::Unsupported(format!(
                "backend '{}' does not yet support {:?} observation for {}",
                self.backend_name(),
                variant.kind.unwrap_or(VariantKind::Other),
                self.path.display()
            ))),
        }
    }

    pub(super) fn observe_snp(
        &self,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
        reference_file: &Path,
    ) -> Result<VariantObservation, RuntimeError> {
        let reference = variant
            .reference
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| {
                RuntimeError::InvalidArguments("SNP variant requires ref/reference".to_owned())
            })?;
        let alternate = variant
            .alternate
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| {
                RuntimeError::InvalidArguments("SNP variant requires alt/alternate".to_owned())
            })?;

        let target_pos = locus.start;
        let pileup = observe_snp_pileup(
            &self.path,
            &self.options,
            reference_file,
            locus,
            reference,
            alternate,
        )?;
        let ref_count = pileup.filtered_ref_count;
        let alt_count = pileup.filtered_alt_count;
        let depth = pileup.filtered_depth;

        let evidence = pileup.evidence_lines(&describe_locus(locus), target_pos);

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
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

    fn observe_snp_with_reader(
        &self,
        reader: &mut cram::io::indexed_reader::IndexedReader<std::fs::File>,
        label: &str,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
    ) -> Result<VariantObservation, RuntimeError> {
        let reference = variant
            .reference
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| {
                RuntimeError::InvalidArguments("SNP variant requires ref/reference".to_owned())
            })?;
        let alternate = variant
            .alternate
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| {
                RuntimeError::InvalidArguments("SNP variant requires alt/alternate".to_owned())
            })?;

        let pileup = snp_pileup_with_reader(
            reader,
            label,
            locus,
            reference,
            alternate,
            self.options.allow_reference_md5_mismatch,
        )?;
        let ref_count = pileup.filtered_ref_count;
        let alt_count = pileup.filtered_alt_count;
        let depth = pileup.filtered_depth;
        let evidence = pileup.evidence_lines(&describe_locus(locus), locus.start);

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
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

    pub(super) fn observe_deletion(
        &self,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
        reference_file: &Path,
    ) -> Result<VariantObservation, RuntimeError> {
        let deletion_length = variant.deletion_length.ok_or_else(|| {
            RuntimeError::InvalidArguments("deletion variant requires deletion_length".to_owned())
        })?;
        let reference = variant.reference.clone().unwrap_or_else(|| "I".to_owned());
        let alternate = variant.alternate.clone().unwrap_or_else(|| "D".to_owned());
        let anchor_pos = locus.start.saturating_sub(1);

        let mut alt_count = 0u32;
        let mut ref_count = 0u32;
        let mut depth = 0u32;

        alignment::for_each_cram_record(
            &self.path,
            &self.options,
            reference_file,
            &anchor_window(locus),
            |record| {
                if record.is_unmapped || !spans_position(&record, anchor_pos) {
                    return Ok(true);
                }
                depth += 1;
                match indel_at_anchor(&record, anchor_pos) {
                    Some((AlignmentOpKind::Deletion, len)) if len == deletion_length => {
                        alt_count += 1;
                    }
                    _ => ref_count += 1,
                }
                Ok(true)
            },
        )?;

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_copy_number_genotype(
                &reference, &alternate, ref_count, alt_count, depth,
            ),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            raw_counts: BTreeMap::new(),
            decision: Some(describe_copy_number_decision_rule(
                &reference, &alternate, ref_count, alt_count, depth,
            )),
            evidence: vec![format!(
                "observed deletion anchor {}:{} len={} depth={} ref_count={} alt_count={}",
                locus.chrom, anchor_pos, deletion_length, depth, ref_count, alt_count
            )],
        })
    }

    fn observe_deletion_with_reader(
        &self,
        reader: &mut cram::io::indexed_reader::IndexedReader<std::fs::File>,
        label: &str,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
    ) -> Result<VariantObservation, RuntimeError> {
        let deletion_length = variant.deletion_length.ok_or_else(|| {
            RuntimeError::InvalidArguments("deletion variant requires deletion_length".to_owned())
        })?;
        let reference = variant.reference.clone().unwrap_or_else(|| "I".to_owned());
        let alternate = variant.alternate.clone().unwrap_or_else(|| "D".to_owned());
        let anchor_pos = locus.start.saturating_sub(1);

        let mut alt_count = 0u32;
        let mut ref_count = 0u32;
        let mut depth = 0u32;

        alignment::for_each_cram_record_with_reader(
            reader,
            label,
            &anchor_window(locus),
            |record| {
                if record.is_unmapped || !spans_position(&record, anchor_pos) {
                    return Ok(true);
                }
                depth += 1;
                match indel_at_anchor(&record, anchor_pos) {
                    Some((AlignmentOpKind::Deletion, len)) if len == deletion_length => {
                        alt_count += 1;
                    }
                    _ => ref_count += 1,
                }
                Ok(true)
            },
        )?;

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_copy_number_genotype(
                &reference, &alternate, ref_count, alt_count, depth,
            ),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            raw_counts: BTreeMap::new(),
            decision: Some(describe_copy_number_decision_rule(
                &reference, &alternate, ref_count, alt_count, depth,
            )),
            evidence: vec![format!(
                "observed deletion anchor {}:{} len={} depth={} ref_count={} alt_count={}",
                locus.chrom, anchor_pos, deletion_length, depth, ref_count, alt_count
            )],
        })
    }

    pub(super) fn observe_indel(
        &self,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
        reference_file: &Path,
    ) -> Result<VariantObservation, RuntimeError> {
        let reference = variant.reference.clone().ok_or_else(|| {
            RuntimeError::InvalidArguments("indel variant requires ref/reference".to_owned())
        })?;
        let alternate = variant.alternate.clone().ok_or_else(|| {
            RuntimeError::InvalidArguments("indel variant requires alt/alternate".to_owned())
        })?;
        let records =
            alignment::query_cram_records(&self.path, &self.options, reference_file, locus)?;

        let mut alt_count = 0u32;
        let mut ref_count = 0u32;
        let mut depth = 0u32;
        let mut matching_alt_lengths = BTreeSet::new();

        for record in records {
            if record.is_unmapped {
                continue;
            }
            if !record_overlaps_locus(&record, locus) {
                continue;
            }
            let classification =
                classify_expected_indel(&record, locus, reference.len(), &alternate)?;
            if !classification.covering {
                continue;
            }
            depth += 1;
            if classification.matches_alt {
                alt_count += 1;
                matching_alt_lengths.insert(classification.observed_len);
            } else if classification.reference_like {
                ref_count += 1;
            }
        }

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
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_copy_number_genotype(
                &reference, &alternate, ref_count, alt_count, depth,
            ),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            raw_counts: BTreeMap::new(),
            decision: Some(describe_copy_number_decision_rule(
                &reference, &alternate, ref_count, alt_count, depth,
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

    fn observe_indel_with_reader(
        &self,
        reader: &mut cram::io::indexed_reader::IndexedReader<std::fs::File>,
        label: &str,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
    ) -> Result<VariantObservation, RuntimeError> {
        let reference = variant.reference.clone().ok_or_else(|| {
            RuntimeError::InvalidArguments("indel variant requires ref/reference".to_owned())
        })?;
        let alternate = variant.alternate.clone().ok_or_else(|| {
            RuntimeError::InvalidArguments("indel variant requires alt/alternate".to_owned())
        })?;

        super::observe_cram_indel_with_reader(
            reader,
            label,
            locus,
            &reference,
            &alternate,
            variant.rsids.first().cloned(),
            Some(assembly),
        )
    }
}
