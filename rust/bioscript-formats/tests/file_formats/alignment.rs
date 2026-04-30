use super::*;

fn mini_fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

#[test]
fn alignment_index_parsers_handle_in_memory_bytes() {
    let fai = alignment::parse_fai_bytes(b"chr1\t4\t6\t4\t5\n").unwrap();
    let _repository = alignment::build_reference_repository_from_readers(
        std::io::BufReader::new(std::io::Cursor::new(b">chr1\nACGT\n".to_vec())),
        fai,
    );

    let err = alignment::parse_fai_bytes(b"not a fai").unwrap_err();
    assert!(format!("{err:?}").contains("failed to parse FASTA index bytes"));

    let err = alignment::parse_crai_bytes(b"not a crai").unwrap_err();
    assert!(format!("{err:?}").contains("failed to parse CRAM index bytes"));

    let err = alignment::parse_tbi_bytes(b"not a tbi").unwrap_err();
    assert!(format!("{err:?}").contains("failed to parse tabix index bytes"));
}

#[test]
fn alignment_reader_api_reports_invalid_cram_headers_without_real_fixtures() {
    let fai = alignment::parse_fai_bytes(b"chr1\t4\t6\t4\t5\n").unwrap();
    let repository = alignment::build_reference_repository_from_readers(
        std::io::BufReader::new(std::io::Cursor::new(b">chr1\nACGT\n".to_vec())),
        fai,
    );
    let locus = bioscript_core::GenomicLocus {
        chrom: "chr1".to_owned(),
        start: 1,
        end: 1,
    };
    let crai_bytes = fs::read(mini_fixtures_dir().join("mini.cram.crai")).unwrap();
    let mut reader = alignment::build_cram_indexed_reader_from_reader(
        std::io::Cursor::new(b"not a cram".to_vec()),
        alignment::parse_crai_bytes(&crai_bytes).unwrap(),
        repository,
    )
    .unwrap();

    let err =
        alignment::for_each_cram_record_with_reader(&mut reader, "bad.cram", &locus, |_| Ok(true))
            .unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to read CRAM header bad.cram"),
        "{err:?}"
    );

    let fai = alignment::parse_fai_bytes(b"chr1\t4\t6\t4\t5\n").unwrap();
    let repository = alignment::build_reference_repository_from_readers(
        std::io::BufReader::new(std::io::Cursor::new(b">chr1\nACGT\n".to_vec())),
        fai,
    );
    let mut raw_reader = alignment::build_cram_indexed_reader_from_reader(
        std::io::Cursor::new(b"still not a cram".to_vec()),
        alignment::parse_crai_bytes(&crai_bytes).unwrap(),
        repository,
    )
    .unwrap();

    let err = alignment::for_each_raw_cram_record_with_reader(
        &mut raw_reader,
        "raw-bad.cram",
        &locus,
        |_| Ok(true),
    )
    .unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to read CRAM header raw-bad.cram"),
        "{err:?}"
    );
}
