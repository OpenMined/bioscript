use std::path::Path;

use noodles::{
    bam,
    core::{Position, Region},
    csi, sam,
};

use bioscript_core::{GenomicLocus, RuntimeError};

use crate::genotype::GenotypeLoadOptions;

use super::{AlignmentOp, AlignmentOpKind, AlignmentRecord};

pub fn query_bam_records(
    path: &Path,
    options: &GenotypeLoadOptions,
    locus: &GenomicLocus,
) -> Result<Vec<AlignmentRecord>, RuntimeError> {
    let mut reader = build_indexed_reader(path, options)?;
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header: {err}")))?;
    let region = build_region(locus)?;
    let query = reader
        .query(&header, &region)
        .map_err(|err| RuntimeError::Io(format!("failed to query BAM region {region}: {err}")))?;

    let mut records = Vec::new();
    for result in query.records() {
        let record =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?;
        records.push(convert_record(&record)?);
    }
    Ok(records)
}

fn build_indexed_reader(
    path: &Path,
    options: &GenotypeLoadOptions,
) -> Result<bam::io::IndexedReader<noodles::bgzf::io::Reader<std::fs::File>>, RuntimeError> {
    let builder = if let Some(index) = options.input_index.as_deref() {
        match index.extension().and_then(|ext| ext.to_str()) {
            Some("bai") => bam::io::indexed_reader::Builder::default().set_index(
                bam::bai::fs::read(index)
                    .map_err(|err| RuntimeError::Io(format!("failed to read BAM index: {err}")))?,
            ),
            Some("csi") => bam::io::indexed_reader::Builder::default().set_index(
                csi::fs::read(index)
                    .map_err(|err| RuntimeError::Io(format!("failed to read CSI index: {err}")))?,
            ),
            _ => {
                return Err(RuntimeError::InvalidArguments(format!(
                    "unsupported BAM index extension: {}",
                    index.display()
                )));
            }
        }
    } else {
        bam::io::indexed_reader::Builder::default()
    };

    builder
        .build_from_path(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open indexed BAM: {err}")))
}

fn build_region(locus: &GenomicLocus) -> Result<Region, RuntimeError> {
    let start = usize::try_from(locus.start)
        .ok()
        .and_then(Position::new)
        .ok_or_else(|| RuntimeError::InvalidArguments("BAM query start must be >= 1".to_owned()))?;
    let end = usize::try_from(locus.end)
        .ok()
        .and_then(Position::new)
        .ok_or_else(|| RuntimeError::InvalidArguments("BAM query end must be >= 1".to_owned()))?;
    Ok(Region::new(locus.chrom.clone(), start..=end))
}

fn convert_record(record: &bam::Record) -> Result<AlignmentRecord, RuntimeError> {
    let start = match record.alignment_start().transpose() {
        Ok(Some(position)) => i64::try_from(usize::from(position)).map_err(|_| {
            RuntimeError::Unsupported("BAM alignment start exceeds i64 range".to_owned())
        })?,
        Ok(None) => -1,
        Err(err) => {
            return Err(RuntimeError::Io(format!(
                "failed to read BAM alignment_start: {err}"
            )));
        }
    };
    let end = match sam::alignment::Record::alignment_end(record).transpose() {
        Ok(Some(position)) => i64::try_from(usize::from(position)).map_err(|_| {
            RuntimeError::Unsupported("BAM alignment end exceeds i64 range".to_owned())
        })?,
        Ok(None) => start,
        Err(err) => {
            return Err(RuntimeError::Io(format!(
                "failed to read BAM alignment_end: {err}"
            )));
        }
    };
    let cigar = record
        .cigar()
        .iter()
        .map(|result| {
            result
                .map(map_op)
                .map_err(|err| RuntimeError::Io(format!("failed to read BAM CIGAR: {err}")))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let is_unmapped = record.flags().is_unmapped();

    Ok(AlignmentRecord {
        start,
        end,
        is_unmapped,
        cigar,
    })
}

fn map_op(op: sam::alignment::record::cigar::Op) -> AlignmentOp {
    use sam::alignment::record::cigar::op::Kind;

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

#[cfg(test)]
mod tests {
    use std::{fs, num::NonZero};

    use noodles::{
        bam,
        core::Position,
        sam::{
            self,
            alignment::{
                RecordBuf,
                io::Write,
                record::{
                    Flags,
                    cigar::{Op, op::Kind},
                },
                record_buf::{Cigar, Sequence},
            },
            header::record::{
                value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
                value::{
                    Map,
                    map::{Header, ReferenceSequence},
                },
            },
        },
    };

    use super::*;

    #[test]
    fn query_bam_records_streams_indexed_region() -> Result<(), Box<dyn std::error::Error>> {
        let dir = std::env::temp_dir().join(format!("bioscript-bam-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;
        let bam_path = dir.join("mini.bam");
        let bai_path = dir.join("mini.bam.bai");
        write_fixture_bam(&bam_path)?;
        let index = bam::fs::index(&bam_path)?;
        bam::bai::fs::write(&bai_path, &index)?;

        let records = query_bam_records(
            &bam_path,
            &GenotypeLoadOptions {
                input_index: Some(bai_path),
                ..GenotypeLoadOptions::default()
            },
            &GenomicLocus {
                chrom: "chr_test".to_owned(),
                start: 1000,
                end: 1002,
            },
        )?;

        fs::remove_dir_all(&dir)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].start, 1000);
        assert_eq!(records[0].end, 1003);
        assert_eq!(records[0].cigar[0].kind, AlignmentOpKind::Match);
        assert_eq!(records[0].cigar[0].len, 4);
        Ok(())
    }

    fn write_fixture_bam(path: &Path) -> Result<(), Box<dyn std::error::Error>> {
        let header = sam::Header::builder()
            .set_header(
                Map::<Header>::builder()
                    .insert(SORT_ORDER, COORDINATE)
                    .build()?,
            )
            .add_reference_sequence(
                "chr_test",
                Map::<ReferenceSequence>::new(NonZero::new(2000).unwrap()),
            )
            .build();
        let mut writer = fs::File::create(path).map(bam::io::Writer::new)?;
        writer.write_header(&header)?;
        writer.write_alignment_record(&header, &record("hit", 1000)?)?;
        writer.write_alignment_record(&header, &record("miss", 1500)?)?;
        writer.try_finish()?;
        Ok(())
    }

    fn record(name: &str, start: usize) -> Result<RecordBuf, Box<dyn std::error::Error>> {
        Ok(RecordBuf::builder()
            .set_name(name)
            .set_flags(Flags::empty())
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start)?)
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, 4)]))
            .set_sequence(Sequence::from(b"ACGT".as_slice()))
            .build())
    }
}
