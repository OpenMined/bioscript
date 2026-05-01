use bioscript_core::{GenomicLocus, RuntimeError};

use crate::alignment::{AlignmentOpKind, AlignmentRecord};

#[derive(Debug, Clone, Copy)]
pub(crate) struct IndelClassification {
    pub(crate) covering: bool,
    pub(crate) reference_like: bool,
    pub(crate) matches_alt: bool,
    pub(crate) observed_len: usize,
}

pub(crate) fn len_as_i64(len: usize) -> Option<i64> {
    i64::try_from(len).ok()
}

pub(crate) fn spans_position(record: &AlignmentRecord, pos: i64) -> bool {
    pos >= record.start.saturating_sub(1) && pos <= record.end
}

pub(crate) fn record_overlaps_locus(record: &AlignmentRecord, locus: &GenomicLocus) -> bool {
    record.end >= locus.start && record.start <= locus.end
}

pub(crate) fn indel_at_anchor(
    record: &AlignmentRecord,
    anchor_pos: i64,
) -> Option<(AlignmentOpKind, usize)> {
    let mut ref_pos = record.start;

    for op in &record.cigar {
        match op.kind {
            AlignmentOpKind::Match
            | AlignmentOpKind::SequenceMatch
            | AlignmentOpKind::SequenceMismatch
            | AlignmentOpKind::Skip => {
                ref_pos += len_as_i64(op.len)?;
            }
            AlignmentOpKind::Insertion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((AlignmentOpKind::Insertion, op.len));
                }
            }
            AlignmentOpKind::Deletion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((AlignmentOpKind::Deletion, op.len));
                }
                ref_pos += len_as_i64(op.len)?;
            }
            AlignmentOpKind::SoftClip | AlignmentOpKind::HardClip | AlignmentOpKind::Pad => {}
        }
    }

    None
}

pub(crate) fn classify_expected_indel(
    record: &AlignmentRecord,
    locus: &GenomicLocus,
    reference_len: usize,
    alternate: &str,
) -> Result<IndelClassification, RuntimeError> {
    let alt_len = alternate.len();
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
                AlignmentOpKind::Insertion => reference_len + len,
                AlignmentOpKind::Deletion => reference_len.saturating_sub(len),
                _ => reference_len,
            };

            return Ok(IndelClassification {
                covering: true,
                reference_like: false,
                matches_alt: observed_len == alt_len,
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
