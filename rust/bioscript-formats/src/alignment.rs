use std::{
    io::{Read, Seek},
    path::Path,
};

use noodles::cram;

use bioscript_core::{GenomicLocus, RuntimeError};

use crate::genotype::GenotypeLoadOptions;

mod cram_stream;
mod readers;

pub use readers::{
    build_cram_indexed_reader_from_reader, build_reference_repository_from_readers,
    parse_crai_bytes, parse_fai_bytes, parse_tbi_bytes,
};

pub(crate) use cram_stream::for_each_raw_cram_record_with_reader_inner;
pub(crate) use readers::{build_cram_indexed_reader_from_path, build_reference_repository};

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
    cram_stream::for_each_cram_record_with_reader_inner(
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
    cram_stream::for_each_cram_record_with_reader_inner(reader, label, locus, false, on_record)
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

#[cfg(test)]
mod tests {
    use super::cram_stream::*;
    use super::*;
    use std::{fs::File, num::NonZero, path::PathBuf};

    use noodles::sam::{
        self,
        alignment::record::cigar::{Op, op::Kind},
        header::record::value::{Map, map::ReferenceSequence},
    };
    use noodles::{
        core::{Position, Region},
        cram::crai,
        fasta,
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
