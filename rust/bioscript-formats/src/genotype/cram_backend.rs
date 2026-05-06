use std::{
    collections::BTreeMap,
    io::{Read, Seek},
    path::Path,
};

use noodles::core::Position;
use noodles::cram;
use noodles::sam::alignment::{
    Record as _,
    record::{Cigar as _, QualityScores as _, Sequence as _, cigar::op::Kind as CigarOpKind},
};

use bioscript_core::{Assembly, GenomicLocus, RuntimeError, VariantSpec};

use crate::alignment;

use super::GenotypeLoadOptions;

mod indel;
mod reader;
mod store;

#[cfg(test)]
pub(crate) use indel::len_as_i64;
pub(crate) use indel::{
    classify_expected_indel, indel_at_anchor, record_overlaps_locus, spans_position,
};
pub use reader::{observe_cram_indel_with_reader, observe_cram_snp_with_reader};

const DEFAULT_MPILEUP_MIN_BASE_QUALITY: u8 = 13;
const DEFAULT_MPILEUP_MIN_MAPPING_QUALITY: u8 = 0;

pub(crate) fn choose_variant_locus(
    variant: &VariantSpec,
    reference_file: &Path,
) -> Option<(Assembly, GenomicLocus)> {
    match detect_reference_assembly(reference_file) {
        Some(Assembly::Grch38) => variant
            .grch38
            .clone()
            .map(|locus| (Assembly::Grch38, locus))
            .or_else(|| {
                variant
                    .grch37
                    .clone()
                    .map(|locus| (Assembly::Grch37, locus))
            }),
        Some(Assembly::Grch37) => variant
            .grch37
            .clone()
            .map(|locus| (Assembly::Grch37, locus))
            .or_else(|| {
                variant
                    .grch38
                    .clone()
                    .map(|locus| (Assembly::Grch38, locus))
            }),
        None => variant
            .grch38
            .clone()
            .map(|locus| (Assembly::Grch38, locus))
            .or_else(|| {
                variant
                    .grch37
                    .clone()
                    .map(|locus| (Assembly::Grch37, locus))
            }),
    }
}

pub(crate) fn detect_reference_assembly(reference_file: &Path) -> Option<Assembly> {
    let lower = reference_file.to_string_lossy().to_ascii_lowercase();
    if lower.contains("grch38") || lower.contains("hg38") || lower.contains("assembly38") {
        Some(Assembly::Grch38)
    } else if lower.contains("grch37") || lower.contains("hg19") || lower.contains("assembly37") {
        Some(Assembly::Grch37)
    } else {
        None
    }
}

pub(crate) fn describe_locus(locus: &GenomicLocus) -> String {
    format!("{}:{}-{}", locus.chrom, locus.start, locus.end)
}

pub(crate) fn anchor_window(locus: &GenomicLocus) -> GenomicLocus {
    let anchor = locus.start.saturating_sub(1);
    GenomicLocus {
        chrom: locus.chrom.clone(),
        start: anchor,
        end: anchor,
    }
}

pub(crate) fn first_base(value: &str) -> Option<char> {
    value
        .trim()
        .chars()
        .next()
        .map(|ch| ch.to_ascii_uppercase())
}

pub(crate) fn infer_snp_genotype(
    reference: char,
    alternate: char,
    ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> Option<String> {
    if depth == 0 || ref_count + alt_count == 0 {
        return None;
    }
    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    if alt_fraction >= 0.8 {
        Some(format!("{alternate}{alternate}"))
    } else if alt_fraction <= 0.2 {
        Some(format!("{reference}{reference}"))
    } else {
        Some(format!("{reference}{alternate}"))
    }
}

pub(crate) fn describe_snp_decision_rule(
    reference: char,
    alternate: char,
    ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> String {
    if depth == 0 {
        return format!(
            "no covering reads for SNP; genotype unresolved (ref={reference}, alt={alternate})"
        );
    }
    if ref_count + alt_count == 0 {
        return format!(
            "no reads matched the declared SNP alleles; genotype unresolved; counts ref={ref_count} alt={alt_count} depth={depth} for {reference}>{alternate}"
        );
    }

    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    format!(
        "SNP genotype rule: alt_fraction={alt_fraction:.3} with thresholds ref<=0.200, het=(0.200,0.800), alt>=0.800; counts ref={ref_count} alt={alt_count} depth={depth} for {reference}>{alternate}"
    )
}

pub(crate) fn infer_copy_number_genotype(
    reference: &str,
    alternate: &str,
    _ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> Option<String> {
    if depth == 0 {
        return None;
    }
    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    if alt_fraction >= 0.8 {
        Some(format!("{alternate}{alternate}"))
    } else if alt_fraction <= 0.2 {
        Some(format!("{reference}{reference}"))
    } else {
        Some(format!("{reference}{alternate}"))
    }
}

pub(crate) fn describe_copy_number_decision_rule(
    reference: &str,
    alternate: &str,
    ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> String {
    if depth == 0 {
        return format!(
            "no covering reads for copy-number style variant; genotype unresolved (ref={reference}, alt={alternate})"
        );
    }

    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    format!(
        "copy-number genotype rule: alt_fraction={alt_fraction:.3} with thresholds ref<=0.200, het=(0.200,0.800), alt>=0.800; counts ref={ref_count} alt={alt_count} depth={depth} for {reference}->{alternate}"
    )
}

#[derive(Debug, Clone, Default)]
pub(crate) struct SnpPileupCounts {
    pub(crate) filtered_depth: u32,
    pub(crate) filtered_ref_count: u32,
    pub(crate) filtered_alt_count: u32,
    pub(crate) filtered_base_counts: BTreeMap<String, u32>,
    pub(crate) raw_depth: u32,
    pub(crate) raw_ref_count: u32,
    pub(crate) raw_alt_count: u32,
    pub(crate) raw_base_counts: BTreeMap<String, u32>,
    pub(crate) filtered_low_base_quality: u32,
    pub(crate) filtered_low_mapping_quality: u32,
    pub(crate) filtered_non_acgt: u32,
    pub(crate) filtered_unmapped: u32,
    pub(crate) filtered_secondary: u32,
    pub(crate) filtered_qc_fail: u32,
    pub(crate) filtered_duplicate: u32,
    pub(crate) filtered_improper_pair: u32,
    pub(crate) raw_forward_counts: BTreeMap<String, u32>,
    pub(crate) raw_reverse_counts: BTreeMap<String, u32>,
}

impl SnpPileupCounts {
    pub(crate) fn evidence_lines(&self, locus: &str, target_pos: i64) -> Vec<String> {
        vec![
            format!(
                "observed SNP pileup at {locus} target_pos={target_pos} filtered_depth={} ref_count={} alt_count={}",
                self.filtered_depth, self.filtered_ref_count, self.filtered_alt_count
            ),
            format!(
                "raw pileup depth={} ref_count={} alt_count={} raw_counts={:?}",
                self.raw_depth, self.raw_ref_count, self.raw_alt_count, self.raw_base_counts
            ),
            format!(
                "raw strand counts: forward={:?} reverse={:?}",
                self.raw_forward_counts, self.raw_reverse_counts
            ),
            format!(
                "filters applied: min_base_quality={} min_mapping_quality={} filtered_low_base_quality={} filtered_low_mapping_quality={} filtered_non_acgt={} filtered_unmapped={} filtered_secondary={} filtered_qc_fail={} filtered_duplicate={} filtered_improper_pair={}",
                DEFAULT_MPILEUP_MIN_BASE_QUALITY,
                DEFAULT_MPILEUP_MIN_MAPPING_QUALITY,
                self.filtered_low_base_quality,
                self.filtered_low_mapping_quality,
                self.filtered_non_acgt,
                self.filtered_unmapped,
                self.filtered_secondary,
                self.filtered_qc_fail,
                self.filtered_duplicate,
                self.filtered_improper_pair
            ),
        ]
    }
}

fn observe_snp_pileup(
    cram_path: &Path,
    options: &GenotypeLoadOptions,
    reference_file: &Path,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
) -> Result<SnpPileupCounts, RuntimeError> {
    let repository = alignment::build_reference_repository(reference_file)?;
    let mut reader =
        alignment::build_cram_indexed_reader_from_path(cram_path, options, repository)?;
    let label = cram_path.display().to_string();
    snp_pileup_with_reader(
        &mut reader,
        &label,
        locus,
        reference,
        alternate,
        options.allow_reference_md5_mismatch,
    )
}

fn snp_pileup_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    allow_reference_md5_mismatch: bool,
) -> Result<SnpPileupCounts, RuntimeError> {
    let mut counts = SnpPileupCounts::default();
    let target_position = Position::try_from(usize::try_from(locus.start).map_err(|_| {
        RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned())
    })?)
    .map_err(|_| RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned()))?;
    let reference_base = reference as u8;

    alignment::for_each_raw_cram_record_with_reader_inner(
        reader,
        label,
        locus,
        allow_reference_md5_mismatch,
        |record| {
            let flags = record
                .flags()
                .map_err(|err| RuntimeError::Io(format!("failed to read CRAM flags: {err}")))?;
            if flags.is_unmapped() {
                counts.filtered_unmapped += 1;
                return Ok(true);
            }
            if flags.is_secondary() {
                counts.filtered_secondary += 1;
                return Ok(true);
            }
            if flags.is_qc_fail() {
                counts.filtered_qc_fail += 1;
                return Ok(true);
            }
            if flags.is_duplicate() {
                counts.filtered_duplicate += 1;
                return Ok(true);
            }
            if flags.is_segmented() && !flags.is_properly_segmented() {
                counts.filtered_improper_pair += 1;
                return Ok(true);
            }

            let Some((base, base_quality)) =
                cram_base_quality_at_reference_position(&record, target_position, reference_base)?
            else {
                return Ok(true);
            };

            let normalized_base = normalize_pileup_base(base);
            record.mapping_quality().transpose().map_err(|err| {
                RuntimeError::Io(format!("failed to read CRAM mapping quality: {err}"))
            })?;
            let is_reverse = flags.is_reverse_complemented();
            if let Some(base) = normalized_base {
                counts.raw_depth += 1;
                *counts.raw_base_counts.entry(base.to_string()).or_insert(0) += 1;
                let strand_counts = if is_reverse {
                    &mut counts.raw_reverse_counts
                } else {
                    &mut counts.raw_forward_counts
                };
                *strand_counts.entry(base.to_string()).or_insert(0) += 1;
                if base == reference {
                    counts.raw_ref_count += 1;
                } else if base == alternate {
                    counts.raw_alt_count += 1;
                }
            }

            if base_quality < DEFAULT_MPILEUP_MIN_BASE_QUALITY {
                counts.filtered_low_base_quality += 1;
                return Ok(true);
            }

            let Some(base) = normalized_base else {
                counts.filtered_non_acgt += 1;
                return Ok(true);
            };

            counts.filtered_depth += 1;
            *counts
                .filtered_base_counts
                .entry(base.to_string())
                .or_insert(0) += 1;
            if base == reference {
                counts.filtered_ref_count += 1;
            } else if base == alternate {
                counts.filtered_alt_count += 1;
            }
            Ok(true)
        },
    )?;

    Ok(counts)
}

fn cram_base_quality_at_reference_position(
    record: &cram::Record<'_>,
    target_position: Position,
    reference_base: u8,
) -> Result<Option<(u8, u8)>, RuntimeError> {
    let Some(alignment_start) = record.alignment_start() else {
        return Ok(None);
    };
    let alignment_start = alignment_start
        .map_err(|err| RuntimeError::Io(format!("failed to read CRAM alignment start: {err}")))?;
    let mut reference_position = usize::from(alignment_start);
    let target = usize::from(target_position);
    let mut read_position = 0usize;
    let sequence = record.sequence();
    let qualities = record.quality_scores();

    for op in record.cigar().iter() {
        let op = op.map_err(|err| RuntimeError::Io(format!("failed to read CRAM CIGAR: {err}")))?;
        match op.kind() {
            CigarOpKind::Match | CigarOpKind::SequenceMatch | CigarOpKind::SequenceMismatch => {
                for offset in 0..op.len() {
                    if reference_position + offset == target {
                        let base = sequence
                            .get(read_position + offset)
                            .unwrap_or(reference_base);
                        let quality = qualities
                            .iter()
                            .nth(read_position + offset)
                            .transpose()
                            .map_err(|err| {
                                RuntimeError::Io(format!("failed to read CRAM base quality: {err}"))
                            })?
                            .unwrap_or(0);
                        return Ok(Some((base, quality)));
                    }
                }
                reference_position += op.len();
                read_position += op.len();
            }
            CigarOpKind::Insertion | CigarOpKind::SoftClip => {
                read_position += op.len();
            }
            CigarOpKind::Deletion | CigarOpKind::Skip => {
                if target >= reference_position && target < reference_position + op.len() {
                    return Ok(None);
                }
                reference_position += op.len();
            }
            CigarOpKind::HardClip | CigarOpKind::Pad => {}
        }
    }

    Ok(None)
}

pub(crate) fn normalize_pileup_base(base: u8) -> Option<char> {
    match (base as char).to_ascii_uppercase() {
        'A' | 'C' | 'G' | 'T' => Some((base as char).to_ascii_uppercase()),
        _ => None,
    }
}
