use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{self, BufWriter, Cursor, Write},
    path::Path,
};

use flate2::{Compression, write::GzEncoder};
use noodles::bam;

use bioscript_core::{GenomicLocus, RuntimeError};

use crate::genotype::GenotypeLoadOptions;

use super::{
    bam_stream::{build_indexed_reader, build_region},
    build_bam_indexed_reader_from_reader, parse_bai_bytes,
};

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
    let mut templates = TemplateFastqRecords::default();

    for result in reader.records() {
        let record =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?;
        if !record_in_templates(&record, &target_names) {
            continue;
        }
        templates.push(&record)?;
    }

    let mut read1 = FastqWriter::create(read1_path)?;
    let mut read2 = FastqWriter::create(read2_path)?;
    let summary = templates.write_paired(&mut read1, &mut read2)?;
    read1.finish()?;
    read2.finish()?;
    Ok(summary)
}

pub fn write_bam_region_fastq_pair_bytes(
    input_bytes: &[u8],
    index_bytes: &[u8],
    locus: &GenomicLocus,
    gzip: bool,
) -> Result<(Vec<u8>, Vec<u8>, FastqPairSummary), RuntimeError> {
    let target_names = collect_region_template_names_bytes(input_bytes, index_bytes, locus)?;
    let mut reader = bam::io::Reader::new(Cursor::new(input_bytes.to_vec()));
    reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header: {err}")))?;
    let mut templates = TemplateFastqRecords::default();

    for result in reader.records() {
        let record =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?;
        if !record_in_templates(&record, &target_names) {
            continue;
        }
        templates.push(&record)?;
    }

    let mut read1 = MemoryFastqWriter::new(gzip);
    let mut read2 = MemoryFastqWriter::new(gzip);
    let summary = templates.write_paired_memory(&mut read1, &mut read2)?;
    Ok((read1.finish()?, read2.finish()?, summary))
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

fn collect_region_template_names_bytes(
    input_bytes: &[u8],
    index_bytes: &[u8],
    locus: &GenomicLocus,
) -> Result<HashSet<Vec<u8>>, RuntimeError> {
    let bai = parse_bai_bytes(index_bytes)?;
    let mut reader = build_bam_indexed_reader_from_reader(Cursor::new(input_bytes.to_vec()), bai)?;
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
    if record.flags().is_unmapped() {
        return true;
    }
    record
        .name()
        .is_some_and(|name| target_names.contains::<[u8]>(name.as_ref()))
}

#[derive(Debug, Default)]
pub(crate) struct TemplateFastqRecords {
    order: Vec<Vec<u8>>,
    records: HashMap<Vec<u8>, TemplateFastqRecordPair>,
    skipped_records: usize,
}

impl TemplateFastqRecords {
    pub(crate) fn push<R>(&mut self, record: &R) -> Result<(), RuntimeError>
    where
        R: noodles::sam::alignment::Record,
    {
        let flags = record
            .flags()
            .map_err(|err| RuntimeError::Io(format!("failed to read alignment record flags: {err}")))?;
        if flags.is_secondary() || flags.is_supplementary() {
            self.skipped_records += 1;
            return Ok(());
        }
        let Some(name) = record.name() else {
            self.skipped_records += 1;
            return Ok(());
        };
        let bytes: &[u8] = name.as_ref();
        let key: Vec<u8> = bytes.to_vec();
        let fastq_record = FastqRecord::try_from_alignment_record(record)?;
        if let Some(pair) = self.records.get_mut(&key) {
            pair.push(fastq_record, &mut self.skipped_records);
        } else {
            let mut pair = TemplateFastqRecordPair::default();
            pair.push(fastq_record, &mut self.skipped_records);
            self.order.push(key.clone());
            self.records.insert(key, pair);
        }
        Ok(())
    }

    fn write_paired(
        self,
        read1: &mut FastqWriter,
        read2: &mut FastqWriter,
    ) -> Result<FastqPairSummary, RuntimeError> {
        let mut summary = FastqPairSummary {
            read1_records: 0,
            read2_records: 0,
            skipped_records: self.skipped_records,
        };
        for key in self.order {
            let pair = self.records.get(&key).expect("template order key exists");
            if let (Some(first), Some(last)) = (&pair.first, &pair.last) {
                first.write(&mut *read1)?;
                last.write(&mut *read2)?;
                summary.read1_records += 1;
                summary.read2_records += 1;
            } else {
                summary.skipped_records += pair.present_count();
            }
        }
        Ok(summary)
    }

    pub(crate) fn write_paired_memory(
        self,
        read1: &mut MemoryFastqWriter,
        read2: &mut MemoryFastqWriter,
    ) -> Result<FastqPairSummary, RuntimeError> {
        let mut summary = FastqPairSummary {
            read1_records: 0,
            read2_records: 0,
            skipped_records: self.skipped_records,
        };
        for key in self.order {
            let pair = self.records.get(&key).expect("template order key exists");
            if let (Some(first), Some(last)) = (&pair.first, &pair.last) {
                first.write(&mut *read1)?;
                last.write(&mut *read2)?;
                summary.read1_records += 1;
                summary.read2_records += 1;
            } else {
                summary.skipped_records += pair.present_count();
            }
        }
        Ok(summary)
    }
}

#[derive(Debug, Default)]
struct TemplateFastqRecordPair {
    first: Option<FastqRecord>,
    last: Option<FastqRecord>,
}

impl TemplateFastqRecordPair {
    fn push(&mut self, record: FastqRecord, skipped_records: &mut usize) {
        match record.segment {
            FastqSegment::First if self.first.is_none() => self.first = Some(record),
            FastqSegment::Last if self.last.is_none() => self.last = Some(record),
            _ => *skipped_records += 1,
        }
    }

    fn present_count(&self) -> usize {
        usize::from(self.first.is_some()) + usize::from(self.last.is_some())
    }
}

#[derive(Debug)]
struct FastqRecord {
    name: Vec<u8>,
    sequence: Vec<u8>,
    qualities: Vec<u8>,
    segment: FastqSegment,
}

impl FastqRecord {
    fn try_from_alignment_record<R>(record: &R) -> Result<Self, RuntimeError>
    where
        R: noodles::sam::alignment::Record,
    {
        let flags = record
            .flags()
            .map_err(|err| RuntimeError::Io(format!("failed to read alignment record flags: {err}")))?;
        let segment = if flags.is_first_segment() {
            FastqSegment::First
        } else if flags.is_last_segment() {
            FastqSegment::Last
        } else {
            FastqSegment::Other
        };
        let sequence = record.sequence().iter().collect::<Vec<_>>();
        Ok(Self {
            name: record.name().map_or_else(
                || b"*".to_vec(),
                |name| {
                    let bytes: &[u8] = name.as_ref();
                    bytes.to_vec()
                },
            ),
            qualities: fastq_qualities(record, sequence.len())?,
            sequence,
            segment,
        })
    }

    fn write(&self, mut writer: impl Write) -> Result<(), RuntimeError> {
        writer
            .write_all(b"@")
            .and_then(|()| writer.write_all(&self.name))
            .and_then(|()| writer.write_all(b"\n"))
            .and_then(|()| writer.write_all(&self.sequence))
            .and_then(|()| writer.write_all(b"\n+\n"))
            .and_then(|()| writer.write_all(&self.qualities))
            .and_then(|()| writer.write_all(b"\n"))
            .map_err(|err| RuntimeError::Io(format!("failed to write FASTQ record: {err}")))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FastqSegment {
    First,
    Last,
    Other,
}

enum FastqWriter {
    Plain(BufWriter<File>),
    Gzip(Box<GzEncoder<BufWriter<File>>>),
}

impl FastqWriter {
    fn create(path: &Path) -> Result<Self, RuntimeError> {
        let file = File::create(path)
            .map_err(|err| RuntimeError::Io(format!("failed to create FASTQ: {err}")))?;
        let writer = BufWriter::new(file);
        if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
            Ok(Self::Gzip(Box::new(GzEncoder::new(
                writer,
                Compression::default(),
            ))))
        } else {
            Ok(Self::Plain(writer))
        }
    }

    fn finish(self) -> Result<(), RuntimeError> {
        match self {
            Self::Plain(mut writer) => writer
                .flush()
                .map_err(|err| RuntimeError::Io(format!("failed to flush FASTQ: {err}"))),
            Self::Gzip(writer) => (*writer)
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

pub(crate) enum MemoryFastqWriter {
    Plain(Vec<u8>),
    Gzip(Box<GzEncoder<Vec<u8>>>),
}

impl MemoryFastqWriter {
    pub(crate) fn new(gzip: bool) -> Self {
        if gzip {
            Self::Gzip(Box::new(GzEncoder::new(Vec::new(), Compression::default())))
        } else {
            Self::Plain(Vec::new())
        }
    }

    pub(crate) fn finish(self) -> Result<Vec<u8>, RuntimeError> {
        match self {
            Self::Plain(bytes) => Ok(bytes),
            Self::Gzip(writer) => (*writer)
                .finish()
                .map_err(|err| RuntimeError::Io(format!("failed to finish FASTQ gzip: {err}"))),
        }
    }
}

impl Write for MemoryFastqWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(bytes) => bytes.write(buf),
            Self::Gzip(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(bytes) => bytes.flush(),
            Self::Gzip(writer) => writer.flush(),
        }
    }
}

fn fastq_qualities<R>(record: &R, sequence_len: usize) -> Result<Vec<u8>, RuntimeError>
where
    R: noodles::sam::alignment::Record,
{
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
    scores
        .iter()
        .map(|result| {
            result
                .map(|score| score.saturating_add(b'!'))
                .map_err(|err| RuntimeError::Io(format!("failed to read alignment record quality: {err}")))
        })
        .collect()
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
                read1_records: 2,
                read2_records: 2,
                skipped_records: 1,
            }
        );
        let read1 = fs::read_to_string(read1_path)?;
        assert!(read1.contains("@pair\nACGT\n+\nBCDE\n"));
        assert!(read1.contains("@unmapped\nTTTT\n+\nBCDE\n"));
        let read2 = fs::File::open(read2_path).map(GzDecoder::new)?;
        let read2 = std::io::read_to_string(read2)?;
        assert!(read2.contains("@pair\nTGCA\n+\nBCDE\n"));
        assert!(read2.contains("@unmapped\nCCCC\n+\nBCDE\n"));
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
        writer.write_alignment_record(
            &header,
            &unmapped_record("unmapped", Flags::SEGMENTED | Flags::FIRST_SEGMENT, b"TTTT")?,
        )?;
        writer.write_alignment_record(
            &header,
            &unmapped_record("unmapped", Flags::SEGMENTED | Flags::LAST_SEGMENT, b"CCCC")?,
        )?;
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

    fn unmapped_record(
        name: &str,
        flags: Flags,
        sequence: &[u8],
    ) -> Result<RecordBuf, Box<dyn std::error::Error>> {
        Ok(RecordBuf::builder()
            .set_name(name)
            .set_flags(flags | Flags::UNMAPPED)
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
