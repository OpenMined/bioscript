use std::{
    collections::{BTreeMap, HashSet},
    io::{Read, Seek},
};

use noodles::{
    core::{Position, Region, region::Interval},
    cram::{self, crai, io::reader::Container},
    sam::{
        self,
        alignment::{Record as _, record::Cigar as _},
    },
};

use bioscript_core::{GenomicLocus, RuntimeError};

use super::{AlignmentOp, AlignmentOpKind, AlignmentRecord};

#[derive(Debug, Clone)]
pub(crate) struct SelectedContainer {
    pub(crate) offset: u64,
    pub(crate) landmarks: HashSet<u64>,
}

pub(crate) fn for_each_cram_record_with_reader_inner<R, F>(
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

pub(crate) fn build_region(header: &sam::Header, locus: &GenomicLocus) -> Option<Region> {
    let chrom = resolve_reference_name(header, &locus.chrom)?;
    let start = Position::try_from(usize::try_from(locus.start).ok()?).ok()?;
    let end = Position::try_from(usize::try_from(locus.end).ok()?).ok()?;
    let raw = format!("{chrom}:{start}-{end}");
    raw.parse().ok()
}

pub(crate) fn select_query_containers(
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

            let records = slice.records(
                reference_sequence_repository.clone(),
                header,
                &compression_header,
                &core_data_src,
                &external_data_srcs,
            );

            match records {
                Ok(records) => {
                    let mut callback_err: Option<RuntimeError> = None;
                    for record in &records {
                        if !handle_decoded_cram_record(
                            label,
                            record,
                            interval,
                            locus_end,
                            &mut stop,
                            &mut callback_err,
                            on_record,
                        ) {
                            break;
                        }
                    }
                    if let Some(err) = callback_err {
                        return Err(err);
                    }
                }
                Err(err) if allow_reference_md5_mismatch && is_reference_md5_mismatch(&err) => {
                    eprintln!(
                        "[bioscript] warning: CRAM reference MD5 mismatch for {label} slice landmark {landmark} — \
                         this noodles version cannot retry without checksum validation. \
                         Details: {err}"
                    );
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

fn handle_decoded_cram_record<F>(
    label: &str,
    record: &cram::Record<'_>,
    interval: Interval,
    locus_end: i64,
    stop: &mut bool,
    callback_err: &mut Option<RuntimeError>,
    on_record: &mut F,
) -> bool
where
    F: FnMut(cram::Record<'_>) -> Result<bool, RuntimeError>,
{
    let alignment_record = match build_alignment_record_from_cram(label, record) {
        Ok(record) => record,
        Err(err) => {
            *callback_err = Some(err);
            return false;
        }
    };

    if alignment_record.start > locus_end {
        *stop = true;
        return false;
    }

    if !alignment_record_intersects_interval(&alignment_record, interval) {
        return true;
    }

    match on_record(record.clone()) {
        Ok(true) => true,
        Ok(false) => {
            *stop = true;
            false
        }
        Err(err) => {
            *callback_err = Some(err);
            false
        }
    }
}

pub(crate) fn is_reference_md5_mismatch(err: &std::io::Error) -> bool {
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

pub(crate) fn resolve_reference_sequence_id(header: &sam::Header, name: &[u8]) -> Option<usize> {
    header
        .reference_sequences()
        .iter()
        .position(|(candidate, _)| {
            let candidate_name: &[u8] = candidate.as_ref();
            candidate_name == name
        })
}

pub(crate) fn record_intersects_interval(
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

pub(crate) fn alignment_record_intersects_interval(
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

pub(crate) fn resolve_reference_name(header: &sam::Header, chrom: &str) -> Option<String> {
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

pub(crate) fn map_op(op: sam::alignment::record::cigar::Op) -> AlignmentOp {
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
