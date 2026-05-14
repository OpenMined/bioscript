use std::{
    cmp,
    collections::HashMap,
    io::{BufRead, BufReader, Cursor, Read, Seek},
    path::Path,
};

use noodles::core::Position;
use noodles::sam::alignment::Record as _;
use noodles::sam::header::record::value::map::header::{sort_order::COORDINATE, tag::SORT_ORDER};
use noodles::{
    bam, bgzf,
    cram::{self, container::ReferenceSequenceContext, crai, io::reader::Container},
    csi::{
        self,
        binning_index::{Indexer, index::reference_sequence::bin::Chunk},
    },
    fasta::{self, repository::adapters::IndexedReader as FastaIndexedReader},
    tabix, vcf,
    vcf::variant::Record as _,
};

use bioscript_core::RuntimeError;

use crate::genotype::GenotypeLoadOptions;

/// Build a CRAM `IndexedReader` over any `Read + Seek` source given a parsed
/// CRAI index and a reference repository. Mirrors `build_from_path` but with
/// an externally-provided reader — the wasm path uses this with a JS-backed
/// reader; native paths still go through the path-based helper below.
pub fn build_cram_indexed_reader_from_reader<R>(
    reader: R,
    crai_index: crai::Index,
    repository: fasta::Repository,
) -> Result<cram::io::indexed_reader::IndexedReader<R>, RuntimeError>
where
    R: Read,
{
    cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .set_index(crai_index)
        .build_from_reader(reader)
        .map_err(|err| RuntimeError::Io(format!("failed to build indexed CRAM reader: {err}")))
}

/// Build a FASTA `Repository` over any `BufRead + Seek + Send + Sync` source
/// given a parsed FAI index. The `Send + Sync + 'static` bounds come from
/// `fasta::Repository`'s internal `Arc<RwLock<dyn Adapter + Send + Sync>>`
/// cache — on single-threaded wasm32 these can be met via `unsafe impl`.
pub fn build_reference_repository_from_readers<R>(
    reader: R,
    fai_index: fasta::fai::Index,
) -> fasta::Repository
where
    R: BufRead + Seek + Send + Sync + 'static,
{
    let indexed = fasta::io::IndexedReader::new(reader, fai_index);
    fasta::Repository::new(FastaIndexedReader::new(indexed))
}

/// Parse a CRAM index (`.crai`) from an in-memory byte buffer. Used by wasm
/// callers that receive the small index inline while the big CRAM stays on a
/// JS-backed reader.
pub fn parse_crai_bytes(bytes: &[u8]) -> Result<crai::Index, RuntimeError> {
    crai::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse CRAM index bytes: {err}")))
}

/// Parse a BAM index (`.bai`) from an in-memory byte buffer.
pub fn parse_bai_bytes(bytes: &[u8]) -> Result<bam::bai::Index, RuntimeError> {
    bam::bai::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse BAM index bytes: {err}")))
}

/// Build a BAM `IndexedReader` over any `Read` source given a parsed BAI index.
pub fn build_bam_indexed_reader_from_reader<R>(
    reader: R,
    bai_index: bam::bai::Index,
) -> Result<bam::io::indexed_reader::IndexedReader<bgzf::io::Reader<R>>, RuntimeError>
where
    R: Read,
{
    bam::io::indexed_reader::Builder::default()
        .set_index(bai_index)
        .build_from_reader(reader)
        .map_err(|err| RuntimeError::Io(format!("failed to build indexed BAM reader: {err}")))
}

/// Parse a FASTA index (`.fai`) from an in-memory byte buffer.
pub fn parse_fai_bytes(bytes: &[u8]) -> Result<fasta::fai::Index, RuntimeError> {
    fasta::fai::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse FASTA index bytes: {err}")))
}

/// Parse a tabix index (`.tbi`) from an in-memory byte buffer. Used by wasm
/// callers that pass the small index inline while the bgzipped VCF stays on
/// a JS-backed `Read + Seek` reader.
pub fn parse_tbi_bytes(bytes: &[u8]) -> Result<tabix::Index, RuntimeError> {
    tabix::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse tabix index bytes: {err}")))
}

/// Generate a tabix index (`.tbi`) for an in-memory bgzipped VCF.
pub fn generate_vcf_tbi_bytes(bytes: &[u8]) -> Result<Vec<u8>, RuntimeError> {
    let cursor = Cursor::new(bytes);
    let mut reader = vcf::io::Reader::new(bgzf::io::Reader::new(cursor));
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read VCF header: {err}")))?;
    let mut indexer = tabix::index::Indexer::default();
    indexer.set_header(csi::binning_index::index::header::Builder::vcf().build());

    let mut record = vcf::Record::default();
    let mut start_position = reader.get_ref().virtual_position();
    while reader
        .read_record(&mut record)
        .map_err(|err| RuntimeError::Io(format!("failed to read VCF record: {err}")))?
        != 0
    {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);
        let reference_sequence_name = record.reference_sequence_name();
        let start = record
            .variant_start()
            .transpose()
            .map_err(|err| RuntimeError::Io(format!("invalid VCF record start: {err}")))?
            .ok_or_else(|| RuntimeError::Io("VCF record is missing a position".to_owned()))?;
        let end = record
            .variant_end(&header)
            .map_err(|err| RuntimeError::Io(format!("invalid VCF record end: {err}")))?;

        indexer
            .add_record(reference_sequence_name, start, end, chunk)
            .map_err(|err| {
                RuntimeError::Io(format!("failed to add VCF record to tabix index: {err}"))
            })?;
        start_position = end_position;
    }

    let index = indexer.build();
    let mut writer = tabix::io::Writer::new(Vec::new());
    writer
        .write_index(&index)
        .and_then(|()| writer.try_finish())
        .map_err(|err| RuntimeError::Io(format!("failed to write tabix index bytes: {err}")))?;
    Ok(writer.into_inner().into_inner())
}

/// Generate a CRAM CRAI index (`.crai`) for an in-memory CRAM.
pub fn generate_cram_crai_bytes(bytes: &[u8]) -> Result<Vec<u8>, RuntimeError> {
    generate_cram_crai_reader(Cursor::new(bytes))
}

/// Generate a CRAM CRAI index (`.crai`) from any sequential CRAM reader.
pub fn generate_cram_crai_reader<R>(reader: R) -> Result<Vec<u8>, RuntimeError>
where
    R: Read + Seek,
{
    let mut reader = cram::io::Reader::new(reader);
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read CRAM header: {err}")))?;
    let mut index = Vec::new();
    let mut container = Container::default();
    let mut container_position = reader
        .position()
        .map_err(|err| RuntimeError::Io(format!("failed to read CRAM position: {err}")))?;

    loop {
        let container_len = reader
            .read_container(&mut container)
            .map_err(|err| RuntimeError::Io(format!("failed to read CRAM container: {err}")))?;
        if container_len == 0 {
            break;
        }

        let compression_header = container.compression_header().map_err(|err| {
            RuntimeError::Io(format!("failed to read CRAM compression header: {err}"))
        })?;
        let landmarks = container.header().landmarks();
        let slice_count = landmarks.len();

        for (i, result) in container.slices().enumerate() {
            let slice = result
                .map_err(|err| RuntimeError::Io(format!("failed to read CRAM slice: {err}")))?;
            let landmark = landmarks[i];
            let slice_length = if i < slice_count - 1 {
                landmarks[i + 1] - landmark
            } else {
                container_len - landmark
            };
            if slice.header().reference_sequence_context().is_many() {
                let mut reference_sequence_ids: HashMap<
                    Option<usize>,
                    SliceReferenceSequenceAlignmentRangeInclusive,
                > = HashMap::new();

                let (core_data_src, external_data_srcs) = slice.decode_blocks().map_err(|err| {
                    RuntimeError::Io(format!("failed to decode CRAM slice blocks: {err}"))
                })?;

                for record in slice
                    .records(
                        fasta::Repository::default(),
                        &header,
                        &compression_header,
                        &core_data_src,
                        &external_data_srcs,
                    )
                    .map_err(|err| {
                        RuntimeError::Io(format!("failed to decode CRAM slice records: {err}"))
                    })?
                {
                    let range = reference_sequence_ids
                        .entry(record.reference_sequence_id(&header).transpose().map_err(
                            |err| {
                                RuntimeError::Io(format!(
                                    "failed to read CRAM record reference id: {err}"
                                ))
                            },
                        )?)
                        .or_default();

                    let alignment_start = record.alignment_start().transpose().map_err(|err| {
                        RuntimeError::Io(format!(
                            "failed to read CRAM record alignment start: {err}"
                        ))
                    })?;
                    let alignment_end = record.alignment_end().transpose().map_err(|err| {
                        RuntimeError::Io(format!("failed to read CRAM record alignment end: {err}"))
                    })?;
                    range.start = cmp::min(range.start, alignment_start);
                    range.end = cmp::max(range.end, alignment_end);
                }

                let mut sorted_reference_sequence_ids: Vec<_> =
                    reference_sequence_ids.keys().copied().collect();
                sorted_reference_sequence_ids.sort_unstable();

                for reference_sequence_id in sorted_reference_sequence_ids {
                    let (alignment_start, alignment_span) = if reference_sequence_id.is_some() {
                        let range = &reference_sequence_ids[&reference_sequence_id];
                        if let (Some(start), Some(end)) = (range.start, range.end) {
                            let span = usize::from(end) - usize::from(start) + 1;
                            (Some(start), span)
                        } else {
                            (None, 0)
                        }
                    } else {
                        (None, 0)
                    };

                    index.push(crai::Record::new(
                        reference_sequence_id,
                        alignment_start,
                        alignment_span,
                        container_position,
                        landmark as u64,
                        slice_length as u64,
                    ));
                }
            } else {
                let (reference_sequence_id, alignment_start, alignment_span) =
                    match slice.header().reference_sequence_context() {
                        ReferenceSequenceContext::Some(context) => {
                            let reference_sequence_id = Some(context.reference_sequence_id());
                            let alignment_start = Some(context.alignment_start());
                            let alignment_span = context.alignment_span();
                            (reference_sequence_id, alignment_start, alignment_span)
                        }
                        ReferenceSequenceContext::None => (None, None, 0),
                        ReferenceSequenceContext::Many => unreachable!(),
                    };

                index.push(crai::Record::new(
                    reference_sequence_id,
                    alignment_start,
                    alignment_span,
                    container_position,
                    landmark as u64,
                    slice_length as u64,
                ));
            }
        }

        container_position = reader
            .position()
            .map_err(|err| RuntimeError::Io(format!("failed to read CRAM position: {err}")))?;
    }

    let mut writer = crai::io::Writer::new(Vec::new());
    writer
        .write_index(&index)
        .map_err(|err| RuntimeError::Io(format!("failed to write CRAI index bytes: {err}")))?;
    writer
        .finish()
        .map_err(|err| RuntimeError::Io(format!("failed to finish CRAI index bytes: {err}")))
}

/// Generate a BAM BAI index (`.bai`) for an in-memory coordinate-sorted BAM.
pub fn generate_bam_bai_bytes(bytes: &[u8]) -> Result<Vec<u8>, RuntimeError> {
    generate_bam_bai_reader(Cursor::new(bytes))
}

/// Generate a BAM BAI index (`.bai`) from any coordinate-sorted BAM reader.
pub fn generate_bam_bai_reader<R>(reader: R) -> Result<Vec<u8>, RuntimeError>
where
    R: Read,
{
    let mut reader = bam::io::Reader::new(reader);
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header: {err}")))?;
    if header
        .header()
        .and_then(|hdr| hdr.other_fields().get(&SORT_ORDER))
        .is_none_or(|sort_order| sort_order != COORDINATE)
    {
        return Err(RuntimeError::Io(
            "BAM must be coordinate-sorted (SO:coordinate) before indexing".to_owned(),
        ));
    }

    let mut record = bam::Record::default();
    let mut builder = Indexer::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader
        .read_record(&mut record)
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM record: {err}")))?
        != 0
    {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);
        let alignment_context = match (
            record.reference_sequence_id().transpose().map_err(|err| {
                RuntimeError::Io(format!("failed to read BAM reference id: {err}"))
            })?,
            record.alignment_start().transpose().map_err(|err| {
                RuntimeError::Io(format!("failed to read BAM alignment start: {err}"))
            })?,
            record.alignment_end().transpose().map_err(|err| {
                RuntimeError::Io(format!("failed to read BAM alignment end: {err}"))
            })?,
        ) {
            (Some(id), Some(start), Some(end)) => {
                let flags = record.flags();
                Some((id, start, end, !flags.is_unmapped()))
            }
            _ => None,
        };

        builder
            .add_record(alignment_context, chunk)
            .map_err(|err| {
                RuntimeError::Io(format!("failed to add BAM record to BAI index: {err}"))
            })?;
        start_position = end_position;
    }

    let index = builder.build(header.reference_sequences().len());
    let mut writer = bam::bai::io::Writer::new(Vec::new());
    writer
        .write_index(&index)
        .map_err(|err| RuntimeError::Io(format!("failed to write BAI index bytes: {err}")))?;
    Ok(writer.into_inner())
}

/// Generate a FASTA FAI index (`.fai`) for an in-memory reference FASTA.
pub fn generate_fasta_fai_bytes(bytes: &[u8]) -> Result<Vec<u8>, RuntimeError> {
    generate_fasta_fai_reader(Cursor::new(bytes))
}

/// Generate a FASTA FAI index (`.fai`) from any FASTA reader.
pub fn generate_fasta_fai_reader<R>(reader: R) -> Result<Vec<u8>, RuntimeError>
where
    R: Read,
{
    let mut indexer = fasta::io::Indexer::new(BufReader::new(reader));
    let mut records = Vec::new();
    while let Some(record) = indexer
        .index_record()
        .map_err(|err| RuntimeError::Io(format!("failed to index FASTA record: {err}")))?
    {
        records.push(record);
    }

    let index = fasta::fai::Index::from(records);
    let mut writer = fasta::fai::io::Writer::new(Vec::new());
    writer
        .write_index(&index)
        .map_err(|err| RuntimeError::Io(format!("failed to write FAI index bytes: {err}")))?;
    Ok(writer.into_inner())
}

#[derive(Debug)]
struct SliceReferenceSequenceAlignmentRangeInclusive {
    start: Option<Position>,
    end: Option<Position>,
}

impl Default for SliceReferenceSequenceAlignmentRangeInclusive {
    fn default() -> Self {
        Self {
            start: Position::new(usize::MAX),
            end: None,
        }
    }
}

pub(crate) fn build_cram_indexed_reader_from_path(
    path: &Path,
    options: &GenotypeLoadOptions,
    repository: fasta::Repository,
) -> Result<cram::io::indexed_reader::IndexedReader<std::fs::File>, RuntimeError> {
    let mut builder =
        cram::io::indexed_reader::Builder::default().set_reference_sequence_repository(repository);

    if let Some(index_path) = options.input_index.as_ref() {
        let index = crai::fs::read(index_path).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read CRAM index {} for {}: {err}",
                index_path.display(),
                path.display()
            ))
        })?;
        builder = builder.set_index(index);
    }

    builder.build_from_path(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open indexed CRAM {}: {err}",
            path.display()
        ))
    })
}

pub(crate) fn build_reference_repository(
    reference_file: &Path,
) -> Result<fasta::Repository, RuntimeError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_file)
        .map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open indexed FASTA {}: {err}",
                reference_file.display()
            ))
        })?;

    Ok(fasta::Repository::new(FastaIndexedReader::new(reader)))
}
