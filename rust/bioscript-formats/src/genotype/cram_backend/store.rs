use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Write as _,
    path::Path,
};

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};

use crate::alignment::{self, AlignmentOpKind};

use super::{
    anchor_window, choose_variant_locus, classify_expected_indel,
    describe_copy_number_decision_rule, describe_locus, describe_snp_decision_rule, first_base,
    indel_at_anchor, infer_copy_number_genotype, infer_snp_genotype, observe_snp_pileup,
    record_overlaps_locus, spans_position,
};
use crate::genotype::{describe_query, types::CramBackend};

impl CramBackend {
    pub(crate) fn backend_name(&self) -> &'static str {
        "cram"
    }

    pub(crate) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let Some(reference_file) = self.options.reference_file.as_ref() else {
            return Err(RuntimeError::Unsupported(format!(
                "backend '{}' cannot satisfy query '{}' for {} without --reference-file",
                self.backend_name(),
                describe_query(variant),
                self.path.display()
            )));
        };

        let Some((assembly, locus)) = choose_variant_locus(variant, reference_file) else {
            let mut detail = format!(
                "backend '{}' cannot satisfy query '{}' for {} using reference {}",
                self.backend_name(),
                describe_query(variant),
                self.path.display(),
                reference_file.display()
            );
            detail.push_str(". This backend needs GRCh37/GRCh38 coordinates, not only rsIDs");
            if let Some(reference_index) = self.options.reference_index.as_ref() {
                let _ = write!(detail, " (reference index {})", reference_index.display());
            }
            if let Some(input_index) = self.options.input_index.as_ref() {
                let _ = write!(detail, " (input index {})", input_index.display());
            }
            return Err(RuntimeError::Unsupported(detail));
        };

        let observation = match variant.kind.unwrap_or(VariantKind::Other) {
            VariantKind::Snp => self.observe_snp(variant, assembly, &locus, reference_file)?,
            VariantKind::Deletion => {
                self.observe_deletion(variant, assembly, &locus, reference_file)?
            }
            VariantKind::Insertion | VariantKind::Indel => {
                self.observe_indel(variant, assembly, &locus, reference_file)?
            }
            VariantKind::Other => {
                return Err(RuntimeError::Unsupported(format!(
                    "backend '{}' does not yet support {:?} observation for {}",
                    self.backend_name(),
                    variant.kind.unwrap_or(VariantKind::Other),
                    self.path.display()
                )));
            }
        };

        Ok(observation)
    }

    fn observe_snp(
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

    fn observe_deletion(
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

    fn observe_indel(
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
}
