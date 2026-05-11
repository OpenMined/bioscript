use crate::{LibError, LibResult};

use super::variant::NativeVariantCall;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentOp {
    Match(usize),
    Mismatch(usize),
    Insertion(usize),
    Deletion(usize),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeAlignment {
    pub reference: String,
    pub haplotype: String,
    pub ops: Vec<AlignmentOp>,
}

pub fn align_haplotype(reference: &str, haplotype: &str) -> LibResult<NativeAlignment> {
    validate_sequence(reference, "reference")?;
    validate_sequence(haplotype, "haplotype")?;
    let reference = reference.to_ascii_uppercase();
    let haplotype = haplotype.to_ascii_uppercase();
    let ref_bases = reference.as_bytes();
    let hap_bases = haplotype.as_bytes();
    let rows = ref_bases.len() + 1;
    let cols = hap_bases.len() + 1;
    let mut scores = vec![0u32; rows * cols];

    for row in 1..rows {
        scores[row * cols] = row as u32;
    }
    for col in 1..cols {
        scores[col] = col as u32;
    }
    for row in 1..rows {
        for col in 1..cols {
            let substitution = scores[(row - 1) * cols + col - 1]
                + u32::from(ref_bases[row - 1] != hap_bases[col - 1]);
            let deletion = scores[(row - 1) * cols + col] + 1;
            let insertion = scores[row * cols + col - 1] + 1;
            scores[row * cols + col] = substitution.min(deletion).min(insertion);
        }
    }

    let mut row = ref_bases.len();
    let mut col = hap_bases.len();
    let mut ops = Vec::new();
    while row > 0 || col > 0 {
        if row > 0 && col > 0 {
            let cost = u32::from(ref_bases[row - 1] != hap_bases[col - 1]);
            if scores[row * cols + col] == scores[(row - 1) * cols + col - 1] + cost {
                push_op(
                    &mut ops,
                    if cost == 0 {
                        AlignmentOp::Match(1)
                    } else {
                        AlignmentOp::Mismatch(1)
                    },
                );
                row -= 1;
                col -= 1;
                continue;
            }
        }
        if row > 0 && scores[row * cols + col] == scores[(row - 1) * cols + col] + 1 {
            push_op(&mut ops, AlignmentOp::Deletion(1));
            row -= 1;
        } else {
            push_op(&mut ops, AlignmentOp::Insertion(1));
            col -= 1;
        }
    }
    ops.reverse();
    Ok(NativeAlignment {
        reference,
        haplotype,
        ops: coalesce_ops(ops),
    })
}

pub fn call_alignment_variants(
    sample_name: impl Into<String>,
    alignment: &NativeAlignment,
    reference_start: u32,
    variant_depth: u32,
    locus_depth: u32,
) -> LibResult<Vec<NativeVariantCall>> {
    let sample_name = sample_name.into();
    let mut variants = Vec::new();
    let mut ref_pos = reference_start;
    let mut ref_index = 0usize;
    let mut hap_pos = 0usize;
    for op in &alignment.ops {
        match *op {
            AlignmentOp::Match(length) => {
                ref_pos += u32::try_from(length).unwrap_or(u32::MAX);
                ref_index += length;
                hap_pos += length;
            }
            AlignmentOp::Mismatch(length) => {
                for offset in 0..length {
                    variants.push(NativeVariantCall::snp(
                        sample_name.clone(),
                        ref_pos + u32::try_from(offset).unwrap_or(u32::MAX),
                        alignment.reference[ref_index + offset..ref_index + offset + 1].to_owned(),
                        alignment.haplotype[hap_pos + offset..hap_pos + offset + 1].to_owned(),
                        variant_depth,
                        locus_depth,
                    ));
                }
                ref_pos += u32::try_from(length).unwrap_or(u32::MAX);
                ref_index += length;
                hap_pos += length;
            }
            AlignmentOp::Insertion(length) => {
                variants.push(NativeVariantCall::insertion(
                    sample_name.clone(),
                    ref_pos,
                    alignment.haplotype[hap_pos..hap_pos + length].to_owned(),
                    variant_depth,
                    locus_depth,
                ));
                hap_pos += length;
            }
            AlignmentOp::Deletion(length) => {
                variants.push(NativeVariantCall::deletion(
                    sample_name.clone(),
                    ref_pos,
                    alignment.reference[ref_index..ref_index + length].to_owned(),
                    variant_depth,
                    locus_depth,
                ));
                ref_pos += u32::try_from(length).unwrap_or(u32::MAX);
                ref_index += length;
            }
        }
    }
    Ok(variants)
}

fn push_op(ops: &mut Vec<AlignmentOp>, op: AlignmentOp) {
    ops.push(op);
}

fn coalesce_ops(ops: Vec<AlignmentOp>) -> Vec<AlignmentOp> {
    let mut coalesced = Vec::new();
    for op in ops {
        match (coalesced.last_mut(), op) {
            (Some(AlignmentOp::Match(length)), AlignmentOp::Match(next))
            | (Some(AlignmentOp::Mismatch(length)), AlignmentOp::Mismatch(next))
            | (Some(AlignmentOp::Insertion(length)), AlignmentOp::Insertion(next))
            | (Some(AlignmentOp::Deletion(length)), AlignmentOp::Deletion(next)) => *length += next,
            _ => coalesced.push(op),
        }
    }
    coalesced
}

fn validate_sequence(sequence: &str, name: &str) -> LibResult<()> {
    if sequence.is_empty() {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel alignment {name} sequence cannot be empty"
        )));
    }
    if sequence.bytes().any(|base| {
        !matches!(
            base,
            b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'N' | b'n'
        )
    }) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel alignment {name} sequence contains unsupported bases"
        )));
    }
    Ok(())
}
