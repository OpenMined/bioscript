use std::{
    collections::HashSet,
    fs::File,
    io::{self, BufWriter, Write},
    path::Path,
};

use flate2::{Compression, write::GzEncoder};
use noodles::bam;

use bioscript_core::{GenomicLocus, RuntimeError};

use crate::genotype::GenotypeLoadOptions;

use super::bam_stream::{build_indexed_reader, build_region};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FastqPairSummary {
    pub read1_records: usize,
    pub read2_records: usize,
    pub skipped_records: usize,
}

pub fn write_bam_region_fastq_pair(
    input_path: &Path,
    read1_path: &Path,
    read2_path: &Path,
    options: &GenotypeLoadOptions,
    locus: &GenomicLocus,
) -> Result<FastqPairSummary, RuntimeError> {
    let target_names = collect_region_template_names(input_path, options, locus)?;
    let mut reader = File::open(input_path)
        .map(bam::io::Reader::new)
        .map_err(|err| RuntimeError::Io(format!("failed to open BAM: {err}")))?;
    reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header: {err}")))?;
    let mut read1 = FastqWriter::create(read1_path)?;
    let mut read2 = FastqWriter::create(read2_path)?;
    let mut summary = FastqPairSummary {
        read1_records: 0,
        read2_records: 0,
        skipped_records: 0,
    };

    for result in reader.records() {
        let record =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?;
        if !record_in_templates(&record, &target_names) {
            continue;
        }
        emit_fastq_record(&record, &mut read1, &mut read2, &mut summary)?;
    }

    read1.finish()?;
    read2.finish()?;
    Ok(summary)
}

fn collect_region_template_names(
    input_path: &Path,
    options: &GenotypeLoadOptions,
    locus: &GenomicLocus,
) -> Result<HashSet<Vec<u8>>, RuntimeError> {
    let mut reader = build_indexed_reader(input_path, options)?;
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header: {err}")))?;
    let region = build_region(locus)?;
    let query = reader
        .query(&header, &region)
        .map_err(|err| RuntimeError::Io(format!("failed to query BAM region {region}: {err}")))?;

    let mut names = HashSet::new();
    for result in query.records() {
        let record =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?;
        if let Some(name) = record.name() {
            let bytes: &[u8] = name.as_ref();
            names.insert(bytes.to_vec());
        }
    }
    Ok(names)
}

fn record_in_templates(record: &bam::Record, target_names: &HashSet<Vec<u8>>) -> bool {
    record
        .name()
        .is_some_and(|name| target_names.contains::<[u8]>(name.as_ref()))
}

fn emit_fastq_record(
    record: &bam::Record,
    read1: &mut FastqWriter,
    read2: &mut FastqWriter,
    summary: &mut FastqPairSummary,
) -> Result<(), RuntimeError> {
    let flags = record.flags();
    if flags.is_secondary() || flags.is_supplementary() {
        summary.skipped_records += 1;
    } else if flags.is_first_segment() {
        write_fastq_record(read1, record)?;
        summary.read1_records += 1;
    } else if flags.is_last_segment() {
        write_fastq_record(read2, record)?;
        summary.read2_records += 1;
    } else {
        summary.skipped_records += 1;
    }
    Ok(())
}

enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(GzEncoder<BufWriter<File>>),
}

impl FastqWriter {
    fn create(path: &Path) -> Result<Self, RuntimeError> {
        let file = File::create(path)
            .map_err(|err| RuntimeError::Io(format!("failed to create FASTQ: {err}")))?;
        let writer = BufWriter::new(file);
        if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
            Ok(Self::Gzip(GzEncoder::new(writer, Compression::default())))
        } else {
            Ok(Self::Plain(writer))
        }
    }

    fn finish(self) -> Result<(), RuntimeError> {
        match self {
            Self::Plain(mut writer) => writer
                .flush()
                .map_err(|err| RuntimeError::Io(format!("failed to flush FASTQ: {err}"))),
            Self::Gzip(writer) => writer
                .finish()
                .and_then(|mut writer| writer.flush())
                .map_err(|err| RuntimeError::Io(format!("failed to finish FASTQ gzip: {err}"))),
        }
    }
}

impl Write for FastqWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(writer) => writer.write(buf),
            Self::Gzip(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(writer) => writer.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

fn write_fastq_record(mut writer: impl Write, record: &bam::Record) -> Result<(), RuntimeError> {
    let name = record.name().map_or(b"*".as_slice(), |name| name.as_ref());
    let sequence = record.sequence().iter().collect::<Vec<_>>();
    let qualities = fastq_qualities(record, sequence.len())?;
    writer
        .write_all(b"@")
        .and_then(|()| writer.write_all(name))
        .and_then(|()| writer.write_all(b"\n"))
        .and_then(|()| writer.write_all(&sequence))
        .and_then(|()| writer.write_all(b"\n+\n"))
        .and_then(|()| writer.write_all(&qualities))
        .and_then(|()| writer.write_all(b"\n"))
        .map_err(|err| RuntimeError::Io(format!("failed to write FASTQ record: {err}")))
}

fn fastq_qualities(record: &bam::Record, sequence_len: usize) -> Result<Vec<u8>, RuntimeError> {
    let scores = record.quality_scores();
    if scores.is_empty() {
        return Ok(vec![b'I'; sequence_len]);
    }
    if scores.len() != sequence_len {
        return Err(RuntimeError::InvalidArguments(format!(
            "BAM record quality length {} does not match sequence length {sequence_len}",
            scores.len()
        )));
    }
    Ok(scores
        .iter()
        .map(|score| score.saturating_add(b'!'))
        .collect())
}

#[cfg(test)]
mod tests {
    use std::{fs, num::NonZero};

    use flate2::read::GzDecoder;
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
                record_buf::{Cigar, QualityScores, Sequence},
            },
            header::record::value::{Map, map::ReferenceSequence},
            header::record::{
                value::map::Header,
                value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER},
            },
        },
    };

    use super::*;

    #[test]
    fn write_bam_region_fastq_pair_rescues_mates() -> Result<(), Box<dyn std::error::Error>> {
        let dir =
            std::env::temp_dir().join(format!("bioscript-bam-fastq-test-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir)?;
        let bam_path = dir.join("mini.bam");
        let bai_path = dir.join("mini.bam.bai");
        let read1_path = dir.join("r1.fastq");
        let read2_path = dir.join("r2.fastq.gz");
        write_fixture_bam(&bam_path)?;
        let index = bam::fs::index(&bam_path)?;
        bam::bai::fs::write(&bai_path, &index)?;

        let summary = write_bam_region_fastq_pair(
            &bam_path,
            &read1_path,
            &read2_path,
            &GenotypeLoadOptions {
                input_index: Some(bai_path),
                ..GenotypeLoadOptions::default()
            },
            &GenomicLocus {
                chrom: "chr_test".to_owned(),
                start: 1000,
                end: 1004,
            },
        )?;

        assert_eq!(
            summary,
            FastqPairSummary {
                read1_records: 1,
                read2_records: 1,
                skipped_records: 1,
            }
        );
        assert_eq!(fs::read_to_string(read1_path)?, "@pair\nACGT\n+\nBCDE\n");
        let read2 = fs::File::open(read2_path).map(GzDecoder::new)?;
        assert_eq!(std::io::read_to_string(read2)?, "@pair\nTGCA\n+\nBCDE\n");
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
        writer.write_alignment_record(
            &header,
            &record(
                "pair",
                Flags::SEGMENTED | Flags::FIRST_SEGMENT,
                b"ACGT",
                1000,
            )?,
        )?;
        writer.write_alignment_record(
            &header,
            &record(
                "pair",
                Flags::SEGMENTED | Flags::LAST_SEGMENT,
                b"TGCA",
                1500,
            )?,
        )?;
        writer.write_alignment_record(&header, &record("skip", Flags::empty(), b"AAAA", 1002)?)?;
        writer.try_finish()?;
        Ok(())
    }

    fn record(
        name: &str,
        flags: Flags,
        sequence: &[u8],
        start: usize,
    ) -> Result<RecordBuf, Box<dyn std::error::Error>> {
        Ok(RecordBuf::builder()
            .set_name(name)
            .set_flags(flags)
            .set_reference_sequence_id(0)
            .set_alignment_start(Position::try_from(start)?)
            .set_cigar(Cigar::from(vec![Op::new(Kind::Match, sequence.len())]))
            .set_sequence(Sequence::from(sequence))
            .set_quality_scores(
                sequence
                    .iter()
                    .enumerate()
                    .map(|(i, _)| u8::try_from(i + 33).unwrap())
                    .collect::<QualityScores>(),
            )
            .build())
    }
}
