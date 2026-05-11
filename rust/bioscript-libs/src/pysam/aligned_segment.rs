use bioscript_formats::alignment::{AlignmentOp, AlignmentOpKind, AlignmentRecord};

use crate::{LibError, LibResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignedSegment {
    pub query_name: Option<String>,
    pub reference_name: Option<String>,
    pub reference_start: Option<u64>,
    pub reference_end: Option<u64>,
    pub query_sequence: Option<String>,
    pub mapping_quality: Option<u8>,
    pub cigarstring: Option<String>,
    pub is_unmapped: bool,
    pub is_reverse: bool,
}

impl AlignedSegment {
    pub fn from_alignment_record(contig: &str, record: &AlignmentRecord) -> Self {
        Self {
            query_name: None,
            reference_name: Some(contig.to_owned()),
            reference_start: u64::try_from(record.start.saturating_sub(1)).ok(),
            reference_end: u64::try_from(record.end).ok(),
            query_sequence: None,
            mapping_quality: None,
            cigarstring: cigar_string(&record.cigar),
            is_unmapped: record.is_unmapped,
            is_reverse: false,
        }
    }

    pub fn unmapped(query_name: Option<String>) -> Self {
        Self {
            query_name,
            reference_name: None,
            reference_start: None,
            reference_end: None,
            query_sequence: None,
            mapping_quality: None,
            cigarstring: None,
            is_unmapped: true,
            is_reverse: false,
        }
    }

    pub fn get_tag(&self, _tag: &str) -> LibResult<()> {
        Err(LibError::unsupported_feature(super::MODULE, "read tags"))
    }

    pub fn set_tag(&mut self, _tag: &str, _value: &str) -> LibResult<()> {
        Err(LibError::unsupported_feature(
            super::MODULE,
            "read mutation",
        ))
    }
}

fn cigar_string(ops: &[AlignmentOp]) -> Option<String> {
    if ops.is_empty() {
        return None;
    }
    let mut out = String::new();
    for op in ops {
        out.push_str(&op.len.to_string());
        out.push(cigar_op_char(op.kind));
    }
    Some(out)
}

fn cigar_op_char(kind: AlignmentOpKind) -> char {
    match kind {
        AlignmentOpKind::Match => 'M',
        AlignmentOpKind::Insertion => 'I',
        AlignmentOpKind::Deletion => 'D',
        AlignmentOpKind::Skip => 'N',
        AlignmentOpKind::SoftClip => 'S',
        AlignmentOpKind::HardClip => 'H',
        AlignmentOpKind::Pad => 'P',
        AlignmentOpKind::SequenceMatch => '=',
        AlignmentOpKind::SequenceMismatch => 'X',
    }
}
