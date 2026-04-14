use std::{
    collections::{BTreeMap, HashSet},
    io::Seek as _,
    path::Path,
};

use noodles::{
    core::{Position, Region},
    cram::{self, crai, io::reader::Container},
    fasta::{self, repository::adapters::IndexedReader as FastaIndexedReader},
    sam::{
        self,
        alignment::{
            Record as _,
            record::{Cigar as _, Sequence as _},
        },
    },
};

use bioscript_core::{GenomicLocus, RuntimeError};

use crate::genotype::GenotypeLoadOptions;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AlignmentOpKind {
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
pub(crate) struct AlignmentOp {
    pub kind: AlignmentOpKind,
    pub len: usize,
}

#[derive(Debug, Clone)]
pub(crate) struct AlignmentRecord {
    pub start: i64,
    pub end: i64,
    pub is_unmapped: bool,
    pub sequence: Vec<u8>,
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
    mut on_record: F,
) -> Result<(), RuntimeError>
where
    F: FnMut(AlignmentRecord) -> Result<bool, RuntimeError>,
{
    let repository = build_reference_repository(reference_file)?;
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

    let mut reader = builder.build_from_path(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open indexed CRAM {}: {err}",
            path.display()
        ))
    })?;

    let header = reader.read_header().map_err(|err| {
        RuntimeError::Io(format!(
            "failed to read CRAM header {}: {err}",
            path.display()
        ))
    })?;

    let region = build_region(&header, locus).ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "indexed CRAM does not contain contig {} for {}:{}-{}",
            locus.chrom, locus.chrom, locus.start, locus.end
        ))
    })?;

    let selected_containers = select_query_containers(reader.index(), &header, &region)?;

    stream_selected_containers(
        path,
        &mut reader,
        &header,
        &region,
        locus.end,
        &selected_containers,
        &mut on_record,
    )?;

    Ok(())
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

fn build_reference_repository(reference_file: &Path) -> Result<fasta::Repository, RuntimeError> {
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

fn stream_selected_containers<F>(
    path: &Path,
    reader: &mut cram::io::indexed_reader::IndexedReader<std::fs::File>,
    header: &sam::Header,
    region: &Region,
    locus_end: i64,
    selected_containers: &[SelectedContainer],
    on_record: &mut F,
) -> Result<(), RuntimeError>
where
    F: FnMut(AlignmentRecord) -> Result<bool, RuntimeError>,
{
    let interval = region.interval();

    for selected_container in selected_containers {
        let offset = selected_container.offset;
        reader
            .get_mut()
            .seek(std::io::SeekFrom::Start(offset))
            .map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to seek CRAM container at offset {offset} in {}: {err}",
                    path.display()
                ))
            })?;

        let mut container = Container::default();
        let len = reader.read_container(&mut container).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read CRAM container at offset {offset} in {}: {err}",
                path.display()
            ))
        })?;

        if len == 0 {
            break;
        }

        let compression_header = container.compression_header().map_err(|err| {
            RuntimeError::Io(format!(
                "failed to decode CRAM compression header from {}: {err}",
                path.display()
            ))
        })?;

        let landmarks = container.header().landmarks().to_vec();
        let reference_sequence_repository = reader.reference_sequence_repository().clone();

        let mut stop = false;

        for (index, slice_result) in container.slices().enumerate() {
            let slice = slice_result.map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read CRAM slice from {}: {err}",
                    path.display()
                ))
            })?;

            let Some(&landmark_i32) = landmarks.get(index) else {
                return Err(RuntimeError::Io(format!(
                    "missing CRAM slice landmark {} in {}",
                    index,
                    path.display()
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
                    "failed to decode CRAM slice blocks from {}: {err}",
                    path.display()
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
                    let alignment_record = match build_alignment_record_from_cram(path, record) {
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

                    match on_record(alignment_record) {
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
                Err(err) if is_reference_md5_mismatch(&err) => {
                    eprintln!(
                        "[bioscript] warning: CRAM reference MD5 mismatch for {} slice landmark {} — \
                         retrying without checksum validation. Results may be incorrect if the \
                         supplied reference differs from the one used to encode this CRAM. \
                         Details: {}",
                        path.display(),
                        landmark,
                        err
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
                                    match build_alignment_record_from_cram(path, record) {
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

                                match on_record(alignment_record) {
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
                                "failed to decode CRAM slice records from {} (unchecked): {err}",
                                path.display()
                            ))
                        })?;
                }
                Err(err) => {
                    return Err(RuntimeError::Io(format!(
                        "failed to decode CRAM slice records from {}: {err}",
                        path.display()
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
    path: &Path,
    record: &cram::Record<'_>,
) -> Result<AlignmentRecord, RuntimeError> {
    let flags = record.flags().map_err(|err| {
        RuntimeError::Io(format!(
            "failed to read CRAM record flags from {}: {err}",
            path.display()
        ))
    })?;
    let is_unmapped = flags.is_unmapped();

    let start = match record.alignment_start() {
        Some(Ok(pos)) => i64::try_from(usize::from(pos)).map_err(|_| {
            RuntimeError::Unsupported(format!(
                "record alignment start exceeds i64 range in {}",
                path.display()
            ))
        })?,
        Some(Err(err)) => {
            return Err(RuntimeError::Io(format!(
                "failed to read CRAM alignment_start from {}: {err}",
                path.display()
            )));
        }
        None => 0,
    };

    let end = match record.alignment_end() {
        Some(Ok(pos)) => i64::try_from(usize::from(pos)).map_err(|_| {
            RuntimeError::Unsupported(format!(
                "record alignment end exceeds i64 range in {}",
                path.display()
            ))
        })?,
        Some(Err(err)) => {
            return Err(RuntimeError::Io(format!(
                "failed to read CRAM alignment_end from {}: {err}",
                path.display()
            )));
        }
        None => start,
    };

    let sequence: Vec<u8> = record.sequence().iter().collect();
    let cigar = record
        .cigar()
        .iter()
        .map(|result| {
            result.map(map_op).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read record CIGAR from {}: {err}",
                    path.display()
                ))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(AlignmentRecord {
        start,
        end,
        is_unmapped,
        sequence,
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
