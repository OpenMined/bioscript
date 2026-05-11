use std::path::PathBuf;

use bioscript_libs::{
    LibError, ModuleName,
    pyfaidx::Fasta,
    pysam::{AlignedSegment, AlignmentFile},
    supported_modules,
    vcf::{VcfDirection, chosen_initial_surface},
};

#[test]
fn registry_lists_initial_bioscript_import_modules() {
    let modules = supported_modules();
    assert!(
        modules
            .iter()
            .any(|module| module.name == ModuleName::Pysam)
    );
    assert!(
        modules
            .iter()
            .any(|module| module.import_path == "from bioscript import pyfaidx")
    );
    assert_eq!(ModuleName::parse("pysam").unwrap(), ModuleName::Pysam);
    assert!(matches!(
        ModuleName::parse("numpy"),
        Err(LibError::UnknownModule(name)) if name == "numpy"
    ));
}

#[test]
fn pysam_alignment_file_accepts_read_modes_and_rejects_write_modes() {
    let file = AlignmentFile::open(
        "sample.cram",
        "rc",
        Some(PathBuf::from("ref.fa")),
        Some(PathBuf::from("sample.cram.crai")),
    )
    .unwrap();
    assert_eq!(file.path(), PathBuf::from("sample.cram").as_path());
    assert_eq!(
        file.reference_filename(),
        Some(PathBuf::from("ref.fa").as_path())
    );

    let err = AlignmentFile::open("out.bam", "wb", None, None).unwrap_err();
    assert!(matches!(
        err,
        LibError::UnsupportedMode {
            object: "AlignmentFile",
            ..
        }
    ));

    let err = AlignmentFile::open("https://example.org/sample.cram", "rc", None, None).unwrap_err();
    assert!(err.to_string().contains("remote alignment files"));
}

#[test]
fn pysam_fetch_validates_region_before_backend_exists() {
    let file = AlignmentFile::open("sample.cram", "rc", None, None).unwrap();
    let err = file.fetch("", Some(1), Some(2)).unwrap_err();
    assert!(err.to_string().contains("requires a contig"));

    let err = file.fetch("22", Some(10), Some(9)).unwrap_err();
    assert!(err.to_string().contains("stop must be >= start"));

    let err = file.fetch("22", Some(9), Some(10)).unwrap_err();
    assert!(err.to_string().contains("requires reference_filename"));
}

#[test]
fn pysam_fetch_streams_tiny_cram_fixture() {
    let fixtures =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../bioscript-formats/tests/fixtures");
    let cram = fixtures.join("mini.cram");
    let reference = fixtures.join("mini.fa");
    let index = fixtures.join("mini.cram.crai");
    let file = AlignmentFile::open(cram, "rc", Some(reference), Some(index)).unwrap();
    let fetched = file.fetch("chr_test", Some(999), Some(1001)).unwrap();
    assert_eq!(fetched.contig, "chr_test");
    assert!(fetched.records.iter().any(|record| {
        record.reference_name.as_deref() == Some("chr_test")
            && record.reference_start.is_some()
            && record.reference_end.is_some()
    }));
}

#[test]
fn pysam_read_tags_and_mutation_are_explicitly_unsupported() {
    let mut read = AlignedSegment::unmapped(Some("read1".to_owned()));
    assert!(
        read.get_tag("NM")
            .unwrap_err()
            .to_string()
            .contains("read tags")
    );
    assert!(
        read.set_tag("NM", "1")
            .unwrap_err()
            .to_string()
            .contains("read mutation")
    );
}

#[test]
fn pyfaidx_fasta_records_support_python_style_slicing() {
    let fasta = Fasta::open("ref.fa");
    assert_eq!(fasta.path(), PathBuf::from("ref.fa").as_path());

    let record = bioscript_libs::pyfaidx::FastaRecord {
        name: "22".to_owned(),
        sequence: "ACGT".to_owned(),
    };
    assert_eq!(record.slice(1, 3).unwrap(), "CG");
    assert!(record.slice(3, 1).is_err());
}

#[test]
fn pyfaidx_fasta_loads_fixture_and_fetches_contig_sequence() {
    let fixture = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../bioscript-formats/tests/fixtures/mini.fa");
    let fasta = Fasta::from_path(&fixture).unwrap();
    let record = fasta.get("chr_test").unwrap();
    assert_eq!(record.name, "chr_test");
    // Ported from the pyfaidx test_feature_bounds_check.py edge case:
    // seq[0:0] should return a blank string.
    assert_eq!(record.slice(0, 0).unwrap(), "");
    assert_eq!(record.slice(0, 6).unwrap(), "TGTACC");
    assert!(fasta.get("missing").is_err());
}

#[test]
fn vcf_direction_is_pysam_variant_file_first() {
    assert_eq!(chosen_initial_surface(), VcfDirection::PysamVariantFile);
    assert!(bioscript_libs::vcf::open_variant_file().is_err());
}
