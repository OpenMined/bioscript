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
    anchor_window, classify_expected_indel_lengths, describe_copy_number_decision_rule,
    describe_locus, describe_snp_decision_rule, first_base, indel_at_anchor,
    infer_copy_number_genotype, infer_snp_genotype, observe_snp_pileup, record_overlaps_locus,
    recount_snp_pileup_counts, select_observed_snp_alternate, snp_pileup_with_reader,
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
        let mut pileup = observe_snp_pileup(
            &self.path,
            &self.options,
            reference_file,
            locus,
            reference,
            alternate,
        )?;
        let alternate = select_observed_snp_alternate(
            reference,
            alternate,
            &variant.observed_alternates,
            &pileup.filtered_base_counts,
            &pileup.raw_base_counts,
        );
        recount_snp_pileup_counts(&mut pileup, reference, alternate);
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

        let mut pileup = snp_pileup_with_reader(
            reader,
            label,
            locus,
            reference,
            alternate,
            self.options.allow_reference_md5_mismatch,
        )?;
        let alternate = select_observed_snp_alternate(
            reference,
            alternate,
            &variant.observed_alternates,
            &pileup.filtered_base_counts,
            &pileup.raw_base_counts,
        );
        recount_snp_pileup_counts(&mut pileup, reference, alternate);
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
        let alternate_lengths = indel_alternate_lengths(variant, &alternate);
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
            let classification = classify_expected_indel_lengths(
                &record,
                locus,
                reference.len(),
                &alternate_lengths,
            )?;
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
        let alternate_lengths = indel_alternate_lengths(variant, &alternate);

        super::observe_cram_indel_with_reader(
            reader,
            label,
            locus,
            &reference,
            &alternate,
            &alternate_lengths,
            variant.rsids.first().cloned(),
            Some(assembly),
        )
    }
}

fn indel_alternate_lengths(variant: &VariantSpec, fallback_alternate: &str) -> Vec<usize> {
    let mut lengths = variant
        .observed_alternates
        .iter()
        .map(String::len)
        .filter(|len| *len > 0)
        .collect::<Vec<_>>();
    if lengths.is_empty() {
        lengths.push(fallback_alternate.len());
    }
    lengths.sort_unstable();
    lengths.dedup();
    lengths
}

#[cfg(test)]
mod tests {
    use std::{fs, path::PathBuf};

    use super::*;
    use crate::genotype::GenotypeLoadOptions;

    fn fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
    }

    fn locus() -> GenomicLocus {
        GenomicLocus {
            chrom: "chr_test".to_owned(),
            start: 1000,
            end: 1000,
        }
    }

    fn backend() -> CramBackend {
        let dir = fixtures_dir();
        CramBackend {
            path: dir.join("mini.cram"),
            options: GenotypeLoadOptions {
                input_index: Some(dir.join("mini.cram.crai")),
                reference_file: Some(dir.join("mini.fa")),
                ..GenotypeLoadOptions::default()
            },
        }
    }

    fn open_reader() -> cram::io::indexed_reader::IndexedReader<std::fs::File> {
        let dir = fixtures_dir();
        let reference = dir.join("mini.fa");
        let repository = alignment::build_reference_repository(&reference).unwrap();
        let index =
            alignment::parse_crai_bytes(&fs::read(dir.join("mini.cram.crai")).unwrap()).unwrap();
        alignment::build_cram_indexed_reader_from_reader(
            fs::File::open(dir.join("mini.cram")).unwrap(),
            index,
            repository,
        )
        .unwrap()
    }

    #[test]
    fn observe_with_reader_dispatches_snp_deletion_indel_and_unsupported_kind() {
        let backend = backend();
        let locus = locus();

        let mut reader = open_reader();
        let snp = VariantSpec {
            rsids: vec!["mini_snp".to_owned()],
            reference: Some("A".to_owned()),
            alternate: Some("C".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };
        let observation = backend
            .observe_with_reader(&mut reader, "mini.cram", &snp, Assembly::Grch38, &locus)
            .unwrap();
        assert_eq!(observation.backend, "cram");
        assert_eq!(observation.matched_rsid.as_deref(), Some("mini_snp"));
        assert_eq!(observation.genotype.as_deref(), Some("AA"));
        assert_eq!(observation.ref_count, Some(50));
        assert_eq!(observation.alt_count, Some(0));
        assert_eq!(observation.depth, Some(50));

        let mut reader = open_reader();
        let deletion = VariantSpec {
            rsids: vec!["mini_del".to_owned()],
            reference: Some("I".to_owned()),
            alternate: Some("D".to_owned()),
            kind: Some(VariantKind::Deletion),
            deletion_length: Some(1),
            ..VariantSpec::default()
        };
        let observation = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &deletion,
                Assembly::Grch38,
                &locus,
            )
            .unwrap();
        assert_eq!(observation.genotype.as_deref(), Some("II"));
        assert_eq!(observation.ref_count, Some(50));
        assert_eq!(observation.alt_count, Some(0));
        assert!(observation.evidence[0].contains("observed deletion anchor"));

        let mut reader = open_reader();
        let indel = VariantSpec {
            rsids: vec!["mini_indel".to_owned()],
            reference: Some("A".to_owned()),
            alternate: Some("AT".to_owned()),
            kind: Some(VariantKind::Insertion),
            ..VariantSpec::default()
        };
        let observation = backend
            .observe_with_reader(&mut reader, "mini.cram", &indel, Assembly::Grch38, &locus)
            .unwrap();
        assert_eq!(observation.matched_rsid.as_deref(), Some("mini_indel"));
        assert_eq!(observation.genotype.as_deref(), Some("AA"));

        let mut reader = open_reader();
        let err = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &VariantSpec {
                    kind: Some(VariantKind::Other),
                    ..VariantSpec::default()
                },
                Assembly::Grch38,
                &locus,
            )
            .unwrap_err();
        assert!(err.to_string().contains("does not yet support"));
    }

    #[test]
    fn observe_with_reader_reports_required_variant_fields() {
        let backend = backend();
        let locus = locus();

        let mut reader = open_reader();
        let err = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &VariantSpec {
                    alternate: Some("C".to_owned()),
                    kind: Some(VariantKind::Snp),
                    ..VariantSpec::default()
                },
                Assembly::Grch38,
                &locus,
            )
            .unwrap_err();
        assert!(err.to_string().contains("SNP variant requires ref"));

        let mut reader = open_reader();
        let err = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &VariantSpec {
                    reference: Some("A".to_owned()),
                    kind: Some(VariantKind::Snp),
                    ..VariantSpec::default()
                },
                Assembly::Grch38,
                &locus,
            )
            .unwrap_err();
        assert!(err.to_string().contains("SNP variant requires alt"));

        let mut reader = open_reader();
        let err = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &VariantSpec {
                    kind: Some(VariantKind::Deletion),
                    ..VariantSpec::default()
                },
                Assembly::Grch38,
                &locus,
            )
            .unwrap_err();
        assert!(err.to_string().contains("deletion_length"));

        let mut reader = open_reader();
        let err = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &VariantSpec {
                    alternate: Some("AT".to_owned()),
                    kind: Some(VariantKind::Indel),
                    ..VariantSpec::default()
                },
                Assembly::Grch38,
                &locus,
            )
            .unwrap_err();
        assert!(err.to_string().contains("indel variant requires ref"));

        let mut reader = open_reader();
        let err = backend
            .observe_with_reader(
                &mut reader,
                "mini.cram",
                &VariantSpec {
                    reference: Some("A".to_owned()),
                    kind: Some(VariantKind::Indel),
                    ..VariantSpec::default()
                },
                Assembly::Grch38,
                &locus,
            )
            .unwrap_err();
        assert!(err.to_string().contains("indel variant requires alt"));
    }
}
