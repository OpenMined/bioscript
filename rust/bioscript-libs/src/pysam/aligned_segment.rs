use htslib_rs::sam;

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
    pub fn from_hts_record<R>(contig: &str, record: &R) -> LibResult<Self>
    where
        R: sam::alignment::Record + ?Sized,
    {
        let flags = record
            .flags()
            .map_err(|err| LibError::InvalidArguments(err.to_string()))?;
        let alignment_start = record
            .alignment_start()
            .transpose()
            .map_err(|err| LibError::InvalidArguments(err.to_string()))?
            .map(usize::from);
        let cigar_ops = cigar_ops(record)?;
        let reference_span = reference_span(&cigar_ops);
        let query_sequence = record.sequence().iter().collect::<Vec<_>>();

        Ok(Self {
            query_name: record
                .name()
                .map(|name| String::from_utf8_lossy(name).into_owned()),
            reference_name: (!flags.is_unmapped()).then(|| contig.to_owned()),
            reference_start: alignment_start
                .and_then(|start| u64::try_from(start.saturating_sub(1)).ok()),
            reference_end: alignment_start.and_then(|start| {
                reference_span.and_then(|span| u64::try_from(start + span - 1).ok())
            }),
            query_sequence: (!query_sequence.is_empty())
                .then(|| String::from_utf8_lossy(&query_sequence).into_owned()),
            mapping_quality: record
                .mapping_quality()
                .transpose()
                .map_err(|err| LibError::InvalidArguments(err.to_string()))?
                .map(|mapping_quality| mapping_quality.get()),
            cigarstring: cigar_string(&cigar_ops),
            is_unmapped: flags.is_unmapped(),
            is_reverse: flags.is_reverse_complemented(),
        })
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

fn cigar_ops<R>(record: &R) -> LibResult<Vec<sam::alignment::record::cigar::Op>>
where
    R: sam::alignment::Record + ?Sized,
{
    record
        .cigar()
        .iter()
        .collect::<Result<Vec<_>, _>>()
        .map_err(|err| LibError::InvalidArguments(err.to_string()))
}

fn cigar_string(ops: &[sam::alignment::record::cigar::Op]) -> Option<String> {
    if ops.is_empty() {
        return None;
    }
    let mut out = String::new();
    for op in ops {
        out.push_str(&op.len().to_string());
        out.push(cigar_op_char(op.kind()));
    }
    Some(out)
}

fn cigar_op_char(kind: sam::alignment::record::cigar::op::Kind) -> char {
    match kind {
        sam::alignment::record::cigar::op::Kind::Match => 'M',
        sam::alignment::record::cigar::op::Kind::Insertion => 'I',
        sam::alignment::record::cigar::op::Kind::Deletion => 'D',
        sam::alignment::record::cigar::op::Kind::Skip => 'N',
        sam::alignment::record::cigar::op::Kind::SoftClip => 'S',
        sam::alignment::record::cigar::op::Kind::HardClip => 'H',
        sam::alignment::record::cigar::op::Kind::Pad => 'P',
        sam::alignment::record::cigar::op::Kind::SequenceMatch => '=',
        sam::alignment::record::cigar::op::Kind::SequenceMismatch => 'X',
    }
}

fn reference_span(ops: &[sam::alignment::record::cigar::Op]) -> Option<usize> {
    let span = ops
        .iter()
        .filter(|op| {
            matches!(
                op.kind(),
                sam::alignment::record::cigar::op::Kind::Match
                    | sam::alignment::record::cigar::op::Kind::Deletion
                    | sam::alignment::record::cigar::op::Kind::Skip
                    | sam::alignment::record::cigar::op::Kind::SequenceMatch
                    | sam::alignment::record::cigar::op::Kind::SequenceMismatch
            )
        })
        .map(|op| op.len())
        .sum::<usize>();
    (span > 0).then_some(span)
}
