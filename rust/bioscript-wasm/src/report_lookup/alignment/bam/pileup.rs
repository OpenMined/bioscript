use super::*;

pub(super) fn bam_base_quality_at_reference_position(
    record: &noodles::bam::Record,
    target_position: noodles::core::Position,
) -> Result<Option<(u8, u8)>, RuntimeError> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let Some(alignment_start) = record
        .alignment_start()
        .transpose()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM alignment start: {err}")))?
    else {
        return Ok(None);
    };
    let mut reference_position = usize::from(alignment_start);
    let target = usize::from(target_position);
    let mut read_position = 0usize;
    let sequence = record.sequence();
    let qualities = record.quality_scores();

    for result in record.cigar().iter() {
        let op =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM CIGAR: {err}")))?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if target >= reference_position && target < reference_position + len {
                    let offset = target - reference_position;
                    let read_index = read_position + offset;
                    let Some(base) = sequence.get(read_index) else {
                        return Ok(None);
                    };
                    let quality = qualities
                        .as_ref()
                        .get(read_index)
                        .copied()
                        .unwrap_or(u8::MAX);
                    return Ok(Some((base, quality)));
                }
                reference_position += len;
                read_position += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                read_position += len;
            }
            Kind::Deletion | Kind::Skip => {
                if target >= reference_position && target < reference_position + len {
                    return Ok(None);
                }
                reference_position += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(None)
}

pub(super) fn bam_alignment_record(
    label: &str,
    record: &noodles::bam::Record,
) -> Result<bioscript_formats::alignment::AlignmentRecord, RuntimeError> {
    use noodles::sam::alignment::Record as _;

    let flags = record.flags();
    let is_unmapped = flags.is_unmapped();
    let start = record
        .alignment_start()
        .transpose()
        .map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read BAM alignment start from {label}: {err}"
            ))
        })?
        .map(|pos| i64::try_from(usize::from(pos)))
        .transpose()
        .map_err(|_| {
            RuntimeError::Unsupported(format!(
                "record alignment start exceeds i64 range in {label}"
            ))
        })?
        .unwrap_or(0);
    let end = record
        .alignment_end()
        .transpose()
        .map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read BAM alignment end from {label}: {err}"
            ))
        })?
        .map(|pos| i64::try_from(usize::from(pos)))
        .transpose()
        .map_err(|_| {
            RuntimeError::Unsupported(format!("record alignment end exceeds i64 range in {label}"))
        })?
        .unwrap_or(start);
    let cigar = record
        .cigar()
        .iter()
        .map(|result| {
            result.map(map_bam_op).map_err(|err| {
                RuntimeError::Io(format!("failed to read BAM CIGAR from {label}: {err}"))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(bioscript_formats::alignment::AlignmentRecord {
        start,
        end,
        is_unmapped,
        cigar,
    })
}

fn map_bam_op(
    op: noodles::sam::alignment::record::cigar::Op,
) -> bioscript_formats::alignment::AlignmentOp {
    use bioscript_formats::alignment::{AlignmentOp, AlignmentOpKind};
    use noodles::sam::alignment::record::cigar::op::Kind;

    let kind = match op.kind() {
        Kind::Match => AlignmentOpKind::Match,
        Kind::Insertion => AlignmentOpKind::Insertion,
        Kind::Deletion => AlignmentOpKind::Deletion,
        Kind::Skip => AlignmentOpKind::Skip,
        Kind::SoftClip => AlignmentOpKind::SoftClip,
        Kind::HardClip => AlignmentOpKind::HardClip,
        Kind::Pad => AlignmentOpKind::Pad,
        Kind::SequenceMatch => AlignmentOpKind::SequenceMatch,
        Kind::SequenceMismatch => AlignmentOpKind::SequenceMismatch,
    };

    AlignmentOp {
        kind,
        len: op.len(),
    }
}

pub(super) fn normalize_pileup_base(base: u8) -> Option<char> {
    match (base as char).to_ascii_uppercase() {
        'A' | 'C' | 'G' | 'T' => Some((base as char).to_ascii_uppercase()),
        _ => None,
    }
}

pub(super) struct IndelClassification {
    pub(super) covering: bool,
    pub(super) reference_like: bool,
    pub(super) matches_alt: bool,
    pub(super) observed_len: usize,
}

pub(super) fn record_overlaps_locus(
    record: &bioscript_formats::alignment::AlignmentRecord,
    locus: &GenomicLocus,
) -> bool {
    record.end >= locus.start && record.start <= locus.end
}

pub(super) fn spans_position(
    record: &bioscript_formats::alignment::AlignmentRecord,
    pos: i64,
) -> bool {
    pos >= record.start.saturating_sub(1) && pos <= record.end
}

pub(super) fn indel_at_anchor(
    record: &bioscript_formats::alignment::AlignmentRecord,
    anchor_pos: i64,
) -> Option<(bioscript_formats::alignment::AlignmentOpKind, usize)> {
    let mut ref_pos = record.start;

    for op in &record.cigar {
        match op.kind {
            bioscript_formats::alignment::AlignmentOpKind::Match
            | bioscript_formats::alignment::AlignmentOpKind::SequenceMatch
            | bioscript_formats::alignment::AlignmentOpKind::SequenceMismatch
            | bioscript_formats::alignment::AlignmentOpKind::Skip => {
                ref_pos += i64::try_from(op.len).ok()?;
            }
            bioscript_formats::alignment::AlignmentOpKind::Insertion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((
                        bioscript_formats::alignment::AlignmentOpKind::Insertion,
                        op.len,
                    ));
                }
            }
            bioscript_formats::alignment::AlignmentOpKind::Deletion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((
                        bioscript_formats::alignment::AlignmentOpKind::Deletion,
                        op.len,
                    ));
                }
                ref_pos += i64::try_from(op.len).ok()?;
            }
            bioscript_formats::alignment::AlignmentOpKind::SoftClip
            | bioscript_formats::alignment::AlignmentOpKind::HardClip
            | bioscript_formats::alignment::AlignmentOpKind::Pad => {}
        }
    }

    None
}

#[allow(dead_code)]
pub(super) fn classify_expected_indel(
    record: &bioscript_formats::alignment::AlignmentRecord,
    locus: &GenomicLocus,
    reference_len: usize,
    alternate: &str,
) -> Result<IndelClassification, RuntimeError> {
    classify_expected_indel_lengths(record, locus, reference_len, &[alternate.len()])
}

pub(super) fn classify_expected_indel_lengths(
    record: &bioscript_formats::alignment::AlignmentRecord,
    locus: &GenomicLocus,
    reference_len: usize,
    alternate_lengths: &[usize],
) -> Result<IndelClassification, RuntimeError> {
    let anchor_start = locus.start.saturating_sub(1);
    let anchor_end = locus.end;

    let covering = record.start <= locus.start && record.end >= locus.end;
    if !covering {
        return Ok(IndelClassification {
            covering: false,
            reference_like: false,
            matches_alt: false,
            observed_len: reference_len,
        });
    }

    let mut observed_len = reference_len;

    for anchor in anchor_start..=anchor_end {
        if let Some((kind, len)) = indel_at_anchor(record, anchor) {
            observed_len = match kind {
                bioscript_formats::alignment::AlignmentOpKind::Insertion => reference_len + len,
                bioscript_formats::alignment::AlignmentOpKind::Deletion => {
                    reference_len.saturating_sub(len)
                }
                _ => reference_len,
            };

            return Ok(IndelClassification {
                covering: true,
                reference_like: false,
                matches_alt: alternate_lengths.contains(&observed_len),
                observed_len,
            });
        }
    }

    Ok(IndelClassification {
        covering: true,
        reference_like: true,
        matches_alt: false,
        observed_len,
    })
}

#[derive(Default)]
pub(super) struct BamSnpPileupCounts {
    pub(super) filtered_depth: u32,
    pub(super) filtered_ref_count: u32,
    pub(super) filtered_alt_count: u32,
    pub(super) filtered_base_counts: BTreeMap<String, u32>,
    pub(super) raw_depth: u32,
    pub(super) raw_ref_count: u32,
    pub(super) raw_alt_count: u32,
    pub(super) raw_base_counts: BTreeMap<String, u32>,
    pub(super) filtered_low_base_quality: u32,
    pub(super) filtered_low_mapping_quality: u32,
    pub(super) filtered_non_acgt: u32,
    pub(super) filtered_unmapped: u32,
    pub(super) filtered_secondary: u32,
    pub(super) filtered_qc_fail: u32,
    pub(super) filtered_duplicate: u32,
    pub(super) filtered_improper_pair: u32,
    pub(super) raw_forward_counts: BTreeMap<String, u32>,
    pub(super) raw_reverse_counts: BTreeMap<String, u32>,
}

impl BamSnpPileupCounts {
    pub(super) fn evidence_lines(&self, locus: &str, target_pos: i64) -> Vec<String> {
        vec![
            format!(
                "observed BAM SNP pileup at {locus} target_pos={target_pos} filtered_depth={} ref_count={} alt_count={}",
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
                "filters applied: min_base_quality=13 min_mapping_quality=0 filtered_low_base_quality={} filtered_low_mapping_quality={} filtered_non_acgt={} filtered_unmapped={} filtered_secondary={} filtered_qc_fail={} filtered_duplicate={} filtered_improper_pair={}",
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

pub(super) fn infer_snp_genotype(
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
        let mut alleles = [
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase(),
        ];
        alleles.sort_by_key(|allele| match allele {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => 99,
        });
        Some(alleles.iter().collect())
    }
}

pub(super) fn describe_snp_decision_rule(
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

pub(super) fn infer_copy_number_genotype(
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
        let mut alleles = [
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase(),
        ];
        alleles.sort_by_key(|allele| {
            allele.chars().next().map_or(u8::MAX, |ch| match ch {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                'I' => 4,
                'D' => 5,
                _ => 99,
            })
        });
        Some(alleles.concat())
    }
}

pub(super) fn describe_copy_number_decision_rule(
    reference: &str,
    alternate: &str,
    _ref_count: u32,
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
        "copy-number genotype rule: alt_fraction={alt_fraction:.3} with thresholds ref<=0.200, het=(0.200,0.800), alt>=0.800; counts alt={alt_count} depth={depth} for {reference}->{alternate}"
    )
}
