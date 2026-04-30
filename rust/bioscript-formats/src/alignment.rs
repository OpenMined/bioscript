use std::{
    collections::{BTreeMap, HashSet},
    io::{BufRead, Read, Seek},
    path::Path,
};

use noodles::{
    core::{Position, Region},
    cram::{self, crai, io::reader::Container},
    fasta::{self, repository::adapters::IndexedReader as FastaIndexedReader},
    sam::{
        self,
        alignment::{Record as _, record::Cigar as _},
    },
    tabix,
};

use bioscript_core::{GenomicLocus, RuntimeError};

use crate::genotype::GenotypeLoadOptions;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentOpKind {
    Match,
    Insertion,
    Deletion,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    SequenceMatch,
    SequenceMismatch,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentOp {
    pub kind: AlignmentOpKind,
    pub len: usize,
}

#[derive(Debug, Clone)]
pub struct AlignmentRecord {
    pub start: i64,
    pub end: i64,
    pub is_unmapped: bool,
    pub cigar: Vec<AlignmentOp>,
}

#[derive(Debug, Clone)]
struct SelectedContainer {
    offset: u64,
    landmarks: HashSet<u64>,
}

pub(crate) fn for_each_cram_record<F>(
    path: &Path,
    options: &GenotypeLoadOptions,
    reference_file: &Path,
    locus: &GenomicLocus,
    on_record: F,
) -> Result<(), RuntimeError>
where
    F: FnMut(AlignmentRecord) -> Result<bool, RuntimeError>,
{
    let repository = build_reference_repository(reference_file)?;
    let mut reader = build_cram_indexed_reader_from_path(path, options, repository)?;
    let label = path.display().to_string();
    for_each_cram_record_with_reader_inner(
        &mut reader,
        &label,
        locus,
        options.allow_reference_md5_mismatch,
        on_record,
    )
}

pub(crate) fn query_cram_records(
    path: &Path,
    options: &GenotypeLoadOptions,
    reference_file: &Path,
    locus: &GenomicLocus,
) -> Result<Vec<AlignmentRecord>, RuntimeError> {
    let mut records = Vec::new();
    for_each_cram_record(path, options, reference_file, locus, |record| {
        records.push(record);
        Ok(true)
    })?;
    Ok(records)
}

/// Iterate decoded alignment records intersecting `locus`, streaming from an
/// already-built CRAM `IndexedReader`. This is the reader-based entry point
/// used by non-filesystem callers (e.g. wasm with a JS `ReadAt` shim).
pub fn for_each_cram_record_with_reader<R, F>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    on_record: F,
) -> Result<(), RuntimeError>
where
    R: Read + Seek,
    F: FnMut(AlignmentRecord) -> Result<bool, RuntimeError>,
{
    for_each_cram_record_with_reader_inner(reader, label, locus, false, on_record)
}

fn for_each_cram_record_with_reader_inner<R, F>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    allow_reference_md5_mismatch: bool,
    mut on_record: F,
) -> Result<(), RuntimeError>
where
    R: Read + Seek,
    F: FnMut(AlignmentRecord) -> Result<bool, RuntimeError>,
{
    // Same idempotent-rewind rationale as `for_each_raw_cram_record_with_reader`
    // — the CRAM header lives at offset 0; we must rewind before each call so
    // callers can iterate over multiple loci with the same `IndexedReader`.
    reader
        .get_mut()
        .seek(std::io::SeekFrom::Start(0))
        .map_err(|err| RuntimeError::Io(format!("failed to rewind CRAM {label}: {err}")))?;
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read CRAM header {label}: {err}")))?;

    let region = build_region(&header, locus).ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "indexed CRAM does not contain contig {} for {}:{}-{}",
            locus.chrom, locus.chrom, locus.start, locus.end
        ))
    })?;

    let selected_containers = select_query_containers(reader.index(), &header, &region)?;

    stream_selected_alignment_records(
        label,
        reader,
        &header,
        &region,
        locus.end,
        &selected_containers,
        allow_reference_md5_mismatch,
        &mut on_record,
    )
}

/// Iterate raw CRAM records intersecting `locus`, streaming from an
/// already-built CRAM `IndexedReader`. The raw variant preserves the
/// `cram::Record` handle so callers can pull base+quality at a specific
/// reference position (needed for SNP pileups).
pub fn for_each_raw_cram_record_with_reader<R, F>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    on_record: F,
) -> Result<(), RuntimeError>
where
    R: Read + Seek,
    F: FnMut(cram::Record<'_>) -> Result<bool, RuntimeError>,
{
    for_each_raw_cram_record_with_reader_inner(reader, label, locus, false, on_record)
}

pub(crate) fn for_each_raw_cram_record_with_reader_inner<R, F>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    allow_reference_md5_mismatch: bool,
    mut on_record: F,
) -> Result<(), RuntimeError>
where
    R: Read + Seek,
    F: FnMut(cram::Record<'_>) -> Result<bool, RuntimeError>,
{
    // Re-seeks to position 0 before reading the header so this helper is
    // idempotent across repeated calls on the same indexed reader (e.g. a
    // wasm caller looking up N variants in a loop). Otherwise the second
    // call reads garbage because the stream position is wherever the
    // previous container iteration left it.
    reader
        .get_mut()
        .seek(std::io::SeekFrom::Start(0))
        .map_err(|err| RuntimeError::Io(format!("failed to rewind CRAM {label}: {err}")))?;
    let header = reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read CRAM header {label}: {err}")))?;

    let region = build_region(&header, locus).ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "indexed CRAM does not contain contig {} for {}:{}-{}",
            locus.chrom, locus.chrom, locus.start, locus.end
        ))
    })?;

    let selected_containers = select_query_containers(reader.index(), &header, &region)?;

    stream_selected_cram_records(
        label,
        reader,
        &header,
        &region,
        locus.end,
        &selected_containers,
        allow_reference_md5_mismatch,
        &mut on_record,
    )
}

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

fn build_region(header: &sam::Header, locus: &GenomicLocus) -> Option<Region> {
    let chrom = resolve_reference_name(header, &locus.chrom)?;
    let start = Position::try_from(usize::try_from(locus.start).ok()?).ok()?;
    let end = Position::try_from(usize::try_from(locus.end).ok()?).ok()?;
    let raw = format!("{chrom}:{start}-{end}");
    raw.parse().ok()
}

fn select_query_containers(
    index: &crai::Index,
    header: &sam::Header,
    region: &Region,
) -> Result<Vec<SelectedContainer>, RuntimeError> {
    let reference_sequence_id =
        resolve_reference_sequence_id(header, region.name()).ok_or_else(|| {
            RuntimeError::Unsupported(format!(
                "indexed CRAM does not contain contig {}",
                String::from_utf8_lossy(region.name())
            ))
        })?;

    let interval = region.interval();
    let mut containers = BTreeMap::<u64, HashSet<u64>>::new();

    for record in index {
        if record.reference_sequence_id() != Some(reference_sequence_id) {
            continue;
        }

        if !record_intersects_interval(record, interval) {
            continue;
        }

        containers
            .entry(record.offset())
            .or_default()
            .insert(record.landmark());
    }

    Ok(containers
        .into_iter()
        .map(|(offset, landmarks)| SelectedContainer { offset, landmarks })
        .collect())
}

fn stream_selected_alignment_records<R, F>(
    label: &str,
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    header: &sam::Header,
    region: &Region,
    locus_end: i64,
    selected_containers: &[SelectedContainer],
    allow_reference_md5_mismatch: bool,
    on_record: &mut F,
) -> Result<(), RuntimeError>
where
    R: Read + Seek,
    F: FnMut(AlignmentRecord) -> Result<bool, RuntimeError>,
{
    stream_selected_cram_records(
        label,
        reader,
        header,
        region,
        locus_end,
        selected_containers,
        allow_reference_md5_mismatch,
        &mut |record| {
            let alignment_record = build_alignment_record_from_cram(label, &record)?;
            on_record(alignment_record)
        },
    )
}

fn stream_selected_cram_records<R, F>(
    label: &str,
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    header: &sam::Header,
    region: &Region,
    locus_end: i64,
    selected_containers: &[SelectedContainer],
    allow_reference_md5_mismatch: bool,
    on_record: &mut F,
) -> Result<(), RuntimeError>
where
    R: Read + Seek,
    F: FnMut(cram::Record<'_>) -> Result<bool, RuntimeError>,
{
    let interval = region.interval();

    for selected_container in selected_containers {
        let offset = selected_container.offset;
        reader
            .get_mut()
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to seek CRAM container at offset {offset} in {label}: {err}"
                ))
            })?;

        let mut container = Container::default();
        let len = reader.read_container(&mut container).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read CRAM container at offset {offset} in {label}: {err}"
            ))
        })?;

        if len == 0 {
            break;
        }

        let compression_header = container.compression_header().map_err(|err| {
            RuntimeError::Io(format!(
                "failed to decode CRAM compression header from {label}: {err}"
            ))
        })?;

        let landmarks = container.header().landmarks().to_vec();
        let reference_sequence_repository = reader.reference_sequence_repository().clone();

        let mut stop = false;

        for (index, slice_result) in container.slices().enumerate() {
            let slice = slice_result.map_err(|err| {
                RuntimeError::Io(format!("failed to read CRAM slice from {label}: {err}"))
            })?;

            let Some(&landmark_i32) = landmarks.get(index) else {
                return Err(RuntimeError::Io(format!(
                    "missing CRAM slice landmark {index} in {label}"
                )));
            };
            let Ok(landmark) = u64::try_from(landmark_i32) else {
                continue;
            };
            if !selected_container.landmarks.contains(&landmark) {
                continue;
            }

            let (core_data_src, external_data_srcs) = slice.decode_blocks().map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to decode CRAM slice blocks from {label}: {err}"
                ))
            })?;

            let mut callback_err: Option<RuntimeError> = None;
            let decode_result = slice.records_while(
                reference_sequence_repository.clone(),
                header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
                true,
                |record| {
                    let alignment_record = match build_alignment_record_from_cram(label, record) {
                        Ok(r) => r,
                        Err(e) => {
                            callback_err = Some(e);
                            return Ok(false);
                        }
                    };

                    if alignment_record.start > locus_end {
                        stop = true;
                        return Ok(false);
                    }

                    if !alignment_record_intersects_interval(&alignment_record, interval) {
                        return Ok(true);
                    }

                    match on_record(record.clone()) {
                        Ok(true) => Ok(true),
                        Ok(false) => {
                            stop = true;
                            Ok(false)
                        }
                        Err(e) => {
                            callback_err = Some(e);
                            Ok(false)
                        }
                    }
                },
            );

            match decode_result {
                Ok(()) => {}
                Err(err) if allow_reference_md5_mismatch && is_reference_md5_mismatch(&err) => {
                    eprintln!(
                        "[bioscript] warning: CRAM reference MD5 mismatch for {label} slice landmark {landmark} — \
                         retrying without checksum validation. Results may be incorrect if the \
                         supplied reference differs from the one used to encode this CRAM. \
                         Details: {err}"
                    );
                    callback_err = None;
                    stop = false;
                    slice
                        .records_while(
                            reference_sequence_repository.clone(),
                            header,
                            &compression_header,
                            &core_data_src,
                            &external_data_srcs,
                            false,
                            |record| {
                                let alignment_record =
                                    match build_alignment_record_from_cram(label, record) {
                                        Ok(r) => r,
                                        Err(e) => {
                                            callback_err = Some(e);
                                            return Ok(false);
                                        }
                                    };

                                if alignment_record.start > locus_end {
                                    stop = true;
                                    return Ok(false);
                                }

                                if !alignment_record_intersects_interval(
                                    &alignment_record,
                                    interval,
                                ) {
                                    return Ok(true);
                                }

                                match on_record(record.clone()) {
                                    Ok(true) => Ok(true),
                                    Ok(false) => {
                                        stop = true;
                                        Ok(false)
                                    }
                                    Err(e) => {
                                        callback_err = Some(e);
                                        Ok(false)
                                    }
                                }
                            },
                        )
                        .map_err(|err| {
                            RuntimeError::Io(format!(
                                "failed to decode CRAM slice records from {label} (unchecked): {err}"
                            ))
                        })?;
                }
                Err(err) if is_reference_md5_mismatch(&err) => {
                    return Err(RuntimeError::Io(format!(
                        "CRAM reference MD5 mismatch for {label} slice landmark {landmark}; rerun with --allow-md5-mismatch only if this lenient decode is intentional. Details: {err}"
                    )));
                }
                Err(err) => {
                    return Err(RuntimeError::Io(format!(
                        "failed to decode CRAM slice records from {label}: {err}"
                    )));
                }
            }

            if let Some(err) = callback_err {
                return Err(err);
            }

            if stop {
                break;
            }
        }

        if stop {
            break;
        }
    }

    Ok(())
}

fn is_reference_md5_mismatch(err: &std::io::Error) -> bool {
    err.to_string()
        .contains("reference sequence checksum mismatch")
}

fn build_alignment_record_from_cram(
    label: &str,
    record: &cram::Record<'_>,
) -> Result<AlignmentRecord, RuntimeError> {
    let flags = record.flags().map_err(|err| {
        RuntimeError::Io(format!(
            "failed to read CRAM record flags from {label}: {err}"
        ))
    })?;
    let is_unmapped = flags.is_unmapped();

    let start = match record.alignment_start() {
        Some(Ok(pos)) => i64::try_from(usize::from(pos)).map_err(|_| {
            RuntimeError::Unsupported(format!(
                "record alignment start exceeds i64 range in {label}"
            ))
        })?,
        Some(Err(err)) => {
            return Err(RuntimeError::Io(format!(
                "failed to read CRAM alignment_start from {label}: {err}"
            )));
        }
        None => 0,
    };

    let end = match record.alignment_end() {
        Some(Ok(pos)) => i64::try_from(usize::from(pos)).map_err(|_| {
            RuntimeError::Unsupported(format!("record alignment end exceeds i64 range in {label}"))
        })?,
        Some(Err(err)) => {
            return Err(RuntimeError::Io(format!(
                "failed to read CRAM alignment_end from {label}: {err}"
            )));
        }
        None => start,
    };

    let cigar = record
        .cigar()
        .iter()
        .map(|result| {
            result.map(map_op).map_err(|err| {
                RuntimeError::Io(format!("failed to read record CIGAR from {label}: {err}"))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(AlignmentRecord {
        start,
        end,
        is_unmapped,
        cigar,
    })
}

fn resolve_reference_sequence_id(header: &sam::Header, name: &[u8]) -> Option<usize> {
    header
        .reference_sequences()
        .iter()
        .position(|(candidate, _)| {
            let candidate_name: &[u8] = candidate.as_ref();
            candidate_name == name
        })
}

fn record_intersects_interval(
    record: &crai::Record,
    interval: noodles::core::region::Interval,
) -> bool {
    let Some(start) = record.alignment_start() else {
        return false;
    };

    if record.alignment_span() == 0 {
        return false;
    }

    let Some(end) = start.checked_add(record.alignment_span() - 1) else {
        return false;
    };

    interval.intersects((start..=end).into())
}

fn alignment_record_intersects_interval(
    record: &AlignmentRecord,
    interval: noodles::core::region::Interval,
) -> bool {
    let Ok(start) = usize::try_from(record.start).and_then(Position::try_from) else {
        return false;
    };
    let Ok(end) = usize::try_from(record.end).and_then(Position::try_from) else {
        return false;
    };

    interval.intersects((start..=end).into())
}

fn resolve_reference_name(header: &sam::Header, chrom: &str) -> Option<String> {
    let candidates = [
        chrom.to_owned(),
        format!("chr{chrom}"),
        chrom.trim_start_matches("chr").to_owned(),
    ];

    candidates.into_iter().find(|candidate| {
        header.reference_sequences().iter().any(|(name, _)| {
            let name_bytes: &[u8] = name.as_ref();
            name_bytes == candidate.as_bytes()
        })
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
    use super::*;
    use std::{fs::File, num::NonZero, path::PathBuf};

    use noodles::sam::{
        self,
        alignment::record::cigar::{Op, op::Kind},
        header::record::value::{Map, map::ReferenceSequence},
    };

    fn locus(chrom: &str, start: i64, end: i64) -> GenomicLocus {
        GenomicLocus {
            chrom: chrom.to_owned(),
            start,
            end,
        }
    }

    fn header() -> sam::Header {
        sam::Header::builder()
            .add_reference_sequence(
                "chr1",
                Map::<ReferenceSequence>::new(NonZero::new(100).unwrap()),
            )
            .add_reference_sequence(
                "2",
                Map::<ReferenceSequence>::new(NonZero::new(200).unwrap()),
            )
            .build()
    }

    fn mini_fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
    }

    #[test]
    fn alignment_helpers_cover_header_region_and_interval_logic() {
        let header = header();
        assert_eq!(
            resolve_reference_name(&header, "1").as_deref(),
            Some("chr1")
        );
        assert_eq!(
            resolve_reference_name(&header, "chr2").as_deref(),
            Some("2")
        );
        assert_eq!(resolve_reference_name(&header, "3"), None);
        assert_eq!(resolve_reference_sequence_id(&header, b"chr1"), Some(0));
        assert_eq!(resolve_reference_sequence_id(&header, b"2"), Some(1));
        assert_eq!(resolve_reference_sequence_id(&header, b"missing"), None);

        let region = build_region(&header, &locus("1", 10, 20)).unwrap();
        assert_eq!(region.name(), b"chr1");
        assert!(build_region(&header, &locus("missing", 10, 20)).is_none());
        assert!(build_region(&header, &locus("1", -1, 20)).is_none());

        let interval = region.interval();
        let hit = crai::Record::new(Some(0), Position::new(12), 3, 100, 1, 20);
        let miss_ref = crai::Record::new(Some(1), Position::new(12), 3, 100, 1, 20);
        let no_start = crai::Record::new(Some(0), None, 3, 100, 1, 20);
        let zero_span = crai::Record::new(Some(0), Position::new(12), 0, 100, 1, 20);
        assert!(record_intersects_interval(&hit, interval));
        assert!(record_intersects_interval(&miss_ref, interval));
        assert!(!record_intersects_interval(&no_start, interval));
        assert!(!record_intersects_interval(&zero_span, interval));

        let alignment_hit = AlignmentRecord {
            start: 11,
            end: 13,
            is_unmapped: false,
            cigar: Vec::new(),
        };
        let alignment_miss = AlignmentRecord {
            start: 30,
            end: 40,
            is_unmapped: false,
            cigar: Vec::new(),
        };
        let bad_start = AlignmentRecord {
            start: -1,
            end: 10,
            is_unmapped: false,
            cigar: Vec::new(),
        };
        assert!(alignment_record_intersects_interval(
            &alignment_hit,
            interval
        ));
        assert!(!alignment_record_intersects_interval(
            &alignment_miss,
            interval
        ));
        assert!(!alignment_record_intersects_interval(&bad_start, interval));
    }

    #[test]
    fn alignment_helpers_cover_index_selection_and_operation_mapping() {
        let header = header();
        let region = build_region(&header, &locus("1", 10, 20)).unwrap();
        let index = vec![
            crai::Record::new(Some(0), Position::new(8), 5, 100, 1, 20),
            crai::Record::new(Some(0), Position::new(19), 3, 100, 2, 20),
            crai::Record::new(Some(1), Position::new(12), 3, 200, 3, 20),
            crai::Record::new(Some(0), Position::new(30), 3, 300, 4, 20),
        ];
        let selected = select_query_containers(&index, &header, &region).unwrap();
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].offset, 100);
        assert!(selected[0].landmarks.contains(&1));
        assert!(selected[0].landmarks.contains(&2));

        let missing_region: Region = "missing:1-2".parse().unwrap();
        let err = select_query_containers(&index, &header, &missing_region).unwrap_err();
        assert!(err.to_string().contains("does not contain contig"));

        let cases = [
            (Kind::Match, AlignmentOpKind::Match),
            (Kind::Insertion, AlignmentOpKind::Insertion),
            (Kind::Deletion, AlignmentOpKind::Deletion),
            (Kind::Skip, AlignmentOpKind::Skip),
            (Kind::SoftClip, AlignmentOpKind::SoftClip),
            (Kind::HardClip, AlignmentOpKind::HardClip),
            (Kind::Pad, AlignmentOpKind::Pad),
            (Kind::SequenceMatch, AlignmentOpKind::SequenceMatch),
            (Kind::SequenceMismatch, AlignmentOpKind::SequenceMismatch),
        ];
        for (kind, expected) in cases {
            assert_eq!(
                map_op(Op::new(kind, 7)),
                AlignmentOp {
                    kind: expected,
                    len: 7
                }
            );
        }
    }

    #[test]
    fn alignment_helpers_cover_parser_and_builder_errors() {
        assert!(
            parse_crai_bytes(b"not a crai")
                .unwrap_err()
                .to_string()
                .contains("CRAM index")
        );
        assert!(
            parse_fai_bytes(b"not a fai")
                .unwrap_err()
                .to_string()
                .contains("FASTA index")
        );
        assert!(
            parse_tbi_bytes(b"not a tbi")
                .unwrap_err()
                .to_string()
                .contains("tabix index")
        );
        assert!(
            build_reference_repository(Path::new("/definitely/missing/reference.fa"))
                .unwrap_err()
                .to_string()
                .contains("failed to open indexed FASTA")
        );

        let repository = build_reference_repository_from_readers(
            std::io::Cursor::new(Vec::<u8>::new()),
            fasta::fai::Index::default(),
        );
        let options = GenotypeLoadOptions {
            input_index: Some(Path::new("/definitely/missing/input.crai").to_path_buf()),
            ..GenotypeLoadOptions::default()
        };
        let Err(err) =
            build_cram_indexed_reader_from_path(Path::new("sample.cram"), &options, repository)
        else {
            panic!("expected missing CRAM index to fail");
        };
        assert!(err.to_string().contains("failed to read CRAM index"));

        let repository = build_reference_repository_from_readers(
            std::io::Cursor::new(Vec::<u8>::new()),
            fasta::fai::Index::default(),
        );
        let reader = build_cram_indexed_reader_from_reader(
            std::io::Cursor::new(Vec::<u8>::new()),
            crai::Index::default(),
            repository,
        );
        assert!(reader.is_ok());

        let err = std::io::Error::other("reference sequence checksum mismatch: expected");
        assert!(is_reference_md5_mismatch(&err));
        let err = std::io::Error::other("other decode error");
        assert!(!is_reference_md5_mismatch(&err));
    }

    #[test]
    fn alignment_path_and_reader_wrappers_stream_mini_cram_records() {
        let dir = mini_fixtures_dir();
        let cram = dir.join("mini.cram");
        let cram_index = dir.join("mini.cram.crai");
        let reference = dir.join("mini.fa");
        let target = locus("chr_test", 1000, 1000);
        let options = GenotypeLoadOptions {
            input_index: Some(cram_index.clone()),
            ..GenotypeLoadOptions::default()
        };

        let mut path_seen = 0;
        for_each_cram_record(&cram, &options, &reference, &target, |record| {
            path_seen += 1;
            assert!(!record.is_unmapped);
            assert!(record.start <= 1000);
            assert!(record.end >= 1000);
            Ok(path_seen < 3)
        })
        .unwrap();
        assert_eq!(path_seen, 3);

        let records = query_cram_records(&cram, &options, &reference, &target).unwrap();
        assert_eq!(records.len(), 50);

        let missing = locus("missing", 1, 1);
        let err = query_cram_records(&cram, &options, &reference, &missing).unwrap_err();
        assert!(err.to_string().contains("does not contain contig missing"));

        let repository = build_reference_repository(&reference).unwrap();
        let index = parse_crai_bytes(&std::fs::read(cram_index).unwrap()).unwrap();
        let mut reader =
            build_cram_indexed_reader_from_reader(File::open(cram).unwrap(), index, repository)
                .unwrap();

        let mut reader_seen = 0;
        for_each_cram_record_with_reader(&mut reader, "mini.cram", &target, |_| {
            reader_seen += 1;
            Ok(reader_seen < 2)
        })
        .unwrap();
        assert_eq!(reader_seen, 2);

        let err = for_each_cram_record_with_reader(&mut reader, "mini.cram", &target, |_| {
            Err(RuntimeError::Unsupported("callback failed".to_owned()))
        })
        .unwrap_err();
        assert!(err.to_string().contains("callback failed"));

        let mut raw_seen = 0;
        for_each_raw_cram_record_with_reader(&mut reader, "mini.cram", &target, |_| {
            raw_seen += 1;
            Ok(raw_seen < 2)
        })
        .unwrap();
        assert_eq!(raw_seen, 2);
    }
}
