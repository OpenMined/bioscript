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

pub fn query_bam_depth_summary(
    path: &Path,
    options: &GenotypeLoadOptions,
    locus: &GenomicLocus,
) -> Result<DepthSummary, RuntimeError> {
    let records = query_bam_records(path, options, locus)?;
    let span = depth_span(locus)?;
    let mut depths = vec![0_u32; span];
    for record in &records {
        add_record_depth(record, locus.start, &mut depths);
    }
    Ok(DepthSummary::from_depths(depths))
}

pub fn write_bam_region(
    input_path: &Path,
    output_path: &Path,
    options: &GenotypeLoadOptions,
    locus: &GenomicLocus,
) -> Result<usize, RuntimeError> {
    let mut reader = build_indexed_reader(input_path, options)?;
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header: {err}")))?;
    let region = build_region(locus)?;
    let query = reader
        .query(&header, &region)
        .map_err(|err| RuntimeError::Io(format!("failed to query BAM region {region}: {err}")))?;

    let output = std::fs::File::create(output_path)
        .map_err(|err| RuntimeError::Io(format!("failed to create BAM slice: {err}")))?;
    let mut writer = bam::io::Writer::new(output);
    writer
        .write_header(&header)
        .map_err(|err| RuntimeError::Io(format!("failed to write BAM header: {err}")))?;

    let mut count = 0;
    for result in query.records() {
        let record =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?;
        writer
            .write_record(&header, &record)
            .map_err(|err| RuntimeError::Io(format!("failed to write BAM record: {err}")))?;
        count += 1;
    }
    writer
        .try_finish()
        .map_err(|err| RuntimeError::Io(format!("failed to finish BAM slice: {err}")))?;
    Ok(count)
}

pub(crate) fn build_indexed_reader(
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

#[derive(Debug, Clone, PartialEq)]
pub struct DepthSummary {
    pub mean: f64,
    pub median: f64,
    pub stdev: f64,
    pub min: u32,
    pub max: u32,
    pub region_length: usize,
    pub uncovered_bases: usize,
    pub percent_uncovered: f64,
}

impl DepthSummary {
    fn from_depths(mut depths: Vec<u32>) -> Self {
        if depths.is_empty() {
            return Self {
                mean: 0.0,
                median: 0.0,
                stdev: 0.0,
                min: 0,
                max: 0,
                region_length: 0,
                uncovered_bases: 0,
                percent_uncovered: 0.0,
            };
        }
        let region_length = depths.len();
        let region_length_f64 = usize_to_f64(region_length);
        let uncovered_bases = depths.iter().filter(|depth| **depth == 0).count();
        let sum = depths.iter().map(|depth| f64::from(*depth)).sum::<f64>();
        let mean = sum / region_length_f64;
        let stdev = (depths
            .iter()
            .map(|depth| {
                let delta = f64::from(*depth) - mean;
                delta * delta
            })
            .sum::<f64>()
            / region_length_f64)
            .sqrt();
        let min = depths.iter().copied().min().unwrap_or(0);
        let max = depths.iter().copied().max().unwrap_or(0);
        depths.sort_unstable();
        let median = if region_length.is_multiple_of(2) {
            let upper = region_length / 2;
            f64::midpoint(f64::from(depths[upper - 1]), f64::from(depths[upper]))
        } else {
            f64::from(depths[region_length / 2])
        };
        Self {
            mean,
            median,
            stdev,
            min,
            max,
            region_length,
            uncovered_bases,
            percent_uncovered: usize_to_f64(uncovered_bases) / region_length_f64 * 100.0,
        }
    }
}

fn usize_to_f64(value: usize) -> f64 {
    f64::from(u32::try_from(value).expect("BAM depth region length must fit in u32"))
}

pub(crate) fn build_region(locus: &GenomicLocus) -> Result<Region, RuntimeError> {
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

fn depth_span(locus: &GenomicLocus) -> Result<usize, RuntimeError> {
    if locus.end < locus.start {
        return Err(RuntimeError::InvalidArguments(
            "BAM depth end must be >= start".to_owned(),
        ));
    }
    usize::try_from(locus.end - locus.start + 1).map_err(|_| {
        RuntimeError::InvalidArguments("BAM depth region length is too large".to_owned())
    })
}

fn add_record_depth(record: &AlignmentRecord, locus_start: i64, depths: &mut [u32]) {
    if record.is_unmapped || record.start < 1 {
        return;
    }
    let mut reference_position = record.start;
    for op in &record.cigar {
        match op.kind {
            AlignmentOpKind::Match
            | AlignmentOpKind::SequenceMatch
            | AlignmentOpKind::SequenceMismatch => {
                for offset in 0..op.len {
                    let pos = reference_position + i64::try_from(offset).unwrap_or(i64::MAX);
                    if let Ok(index) = usize::try_from(pos - locus_start)
                        && let Some(depth) = depths.get_mut(index)
                    {
                        *depth = depth.saturating_add(1);
                    }
                }
                reference_position += i64::try_from(op.len).unwrap_or(i64::MAX);
            }
            AlignmentOpKind::Deletion | AlignmentOpKind::Skip => {
                reference_position += i64::try_from(op.len).unwrap_or(i64::MAX);
            }
            AlignmentOpKind::Insertion
            | AlignmentOpKind::SoftClip
            | AlignmentOpKind::HardClip
            | AlignmentOpKind::Pad => {}
        }
    }
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

    #[test]
    fn query_bam_depth_summary_counts_zero_coverage_positions()
    -> Result<(), Box<dyn std::error::Error>> {
        let dir =
            std::env::temp_dir().join(format!("bioscript-bam-depth-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;
        let bam_path = dir.join("mini.bam");
        let bai_path = dir.join("mini.bam.bai");
        write_fixture_bam(&bam_path)?;
        let index = bam::fs::index(&bam_path)?;
        bam::bai::fs::write(&bai_path, &index)?;

        let summary = query_bam_depth_summary(
            &bam_path,
            &GenotypeLoadOptions {
                input_index: Some(bai_path),
                ..GenotypeLoadOptions::default()
            },
            &GenomicLocus {
                chrom: "chr_test".to_owned(),
                start: 999,
                end: 1004,
            },
        )?;

        fs::remove_dir_all(&dir)?;
        assert_eq!(summary.region_length, 6);
        assert_eq!(summary.uncovered_bases, 2);
        assert_eq!(summary.min, 0);
        assert_eq!(summary.max, 1);
        assert!((summary.mean - (4.0 / 6.0)).abs() < f64::EPSILON);
        Ok(())
    }

    #[test]
    fn write_bam_region_creates_slice_with_matching_records()
    -> Result<(), Box<dyn std::error::Error>> {
        let dir =
            std::env::temp_dir().join(format!("bioscript-bam-slice-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;
        let bam_path = dir.join("mini.bam");
        let bai_path = dir.join("mini.bam.bai");
        let slice_path = dir.join("slice.bam");
        write_fixture_bam(&bam_path)?;
        let index = bam::fs::index(&bam_path)?;
        bam::bai::fs::write(&bai_path, &index)?;

        let count = write_bam_region(
            &bam_path,
            &slice_path,
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

        assert_eq!(count, 1);
        assert_eq!(count_bam_records(&slice_path)?, 1);
        fs::remove_dir_all(&dir)?;
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

    fn count_bam_records(path: &Path) -> Result<usize, Box<dyn std::error::Error>> {
        let mut reader = fs::File::open(path).map(bam::io::Reader::new)?;
        reader.read_header()?;
        let mut count = 0;
        for result in reader.records() {
            let _ = result?;
            count += 1;
        }
        Ok(count)
    }
}
