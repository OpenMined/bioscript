#![allow(clippy::float_cmp)]

use std::{
    io::{Read, Write},
    path::PathBuf,
};

use bioscript_libs::{
    LibError, ModuleName, bcftools,
    kestrel::{
        KestrelRunConfig,
        native::{
            NativeKestrelRunOptions, NativeReferenceRegion, call_fastq_paths_to_vcf_references,
            call_sequences_to_vcf,
        },
    },
    pyfaidx::Fasta,
    pysam::{AlignedSegment, AlignmentFile},
    samtools, supported_modules,
    vcf::{VcfDirection, chosen_initial_surface, parse_kestrel_vcf},
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
    assert_eq!(ModuleName::parse("kestrel").unwrap(), ModuleName::Kestrel);
    assert_eq!(ModuleName::parse("samtools").unwrap(), ModuleName::Samtools);
    assert_eq!(ModuleName::parse("bcftools").unwrap(), ModuleName::Bcftools);
    assert!(matches!(
        ModuleName::parse("numpy"),
        Err(LibError::UnknownModule(name)) if name == "numpy"
    ));
}

#[test]
fn bcftools_vntyper_subset_builds_allowed_commands() {
    let sorted = bcftools::sort(
        PathBuf::from("calls.vcf").as_path(),
        PathBuf::from("calls.vcf.gz").as_path(),
    )
    .unwrap();
    assert_eq!(
        sorted.argv(),
        vec!["bcftools", "sort", "-Oz", "-o", "calls.vcf.gz", "calls.vcf"]
    );

    let filtered = bcftools::view_filter(
        PathBuf::from("calls.vcf").as_path(),
        PathBuf::from("pass.vcf.gz").as_path(),
        "FILTER=\"PASS\"",
    )
    .unwrap();
    assert_eq!(filtered.program(), "bcftools");
    assert_eq!(filtered.args()[0], "view");
    assert!(filtered.args().contains(&"FILTER=\"PASS\"".to_owned()));

    let viewed = bcftools::view(
        PathBuf::from("calls.vcf").as_path(),
        PathBuf::from("calls.bcf").as_path(),
        "b",
    )
    .unwrap();
    assert_eq!(
        viewed.argv(),
        vec![
            "bcftools",
            "view",
            "-O",
            "b",
            "-o",
            "calls.bcf",
            "calls.vcf"
        ]
    );

    let normalized = bcftools::norm(
        PathBuf::from("calls.vcf").as_path(),
        PathBuf::from("ref.fa").as_path(),
        PathBuf::from("norm.vcf.gz").as_path(),
    )
    .unwrap();
    assert_eq!(
        normalized.argv(),
        vec![
            "bcftools",
            "norm",
            "-f",
            "ref.fa",
            "-Oz",
            "-o",
            "norm.vcf.gz",
            "calls.vcf"
        ]
    );
}

#[test]
fn bcftools_native_view_header_uses_vendored_bcftools_rs() {
    let temp = tempfile::tempdir().unwrap();
    let input = temp.path().join("input.vcf");
    let output = temp.path().join("header.vcf");
    std::fs::write(
        &input,
        concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=16>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr1\t5\t.\tC\tT\t.\tPASS\t.\n",
        ),
    )
    .unwrap();

    bcftools::view_header_native(&input, &output).unwrap();
    let header = std::fs::read_to_string(output).unwrap();

    assert!(header.contains("##fileformat=VCFv4.2\n"));
    assert!(header.contains("##contig=<ID=chr1,length=16>\n"));
    assert!(header.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert!(!header.contains("chr1\t5\t.\tC\tT"));
    assert!(!header.contains("##bcftools_viewVersion="));
}

#[test]
fn bcftools_native_view_writes_bgzf_vcf_and_index_writes_tbi() {
    let temp = tempfile::tempdir().unwrap();
    let input = temp.path().join("input.vcf");
    let compressed = temp.path().join("output.vcf.gz");
    let index = temp.path().join("output.vcf.gz.tbi");
    std::fs::write(
        &input,
        concat!(
            "##fileformat=VCFv4.2\n",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
            "##contig=<ID=chr1,length=1000>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "chr1\t5\t.\tC\tT\t.\tPASS\t.\n",
            "chr1\t8\t.\tG\tA\t.\tPASS\t.\n",
        ),
    )
    .unwrap();

    bcftools::view_native(&input, &compressed, "z").unwrap();
    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(&compressed).unwrap());
    let mut vcf = String::new();
    decoder.read_to_string(&mut vcf).unwrap();
    assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"));
    assert!(vcf.contains("chr1\t5\t.\tC\tT"));
    assert!(!vcf.contains("##bcftools_viewVersion="));

    bcftools::index_native(&compressed, Some(&index), true, true).unwrap();
    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(index).unwrap());
    let mut magic = [0u8; 4];
    decoder.read_exact(&mut magic).unwrap();
    assert_eq!(&magic, b"TBI\x01");
}

#[test]
fn bcftools_native_sort_writes_bgzf_vcf_and_csi() {
    let temp = tempfile::tempdir().unwrap();
    let input = temp.path().join("unsorted.vcf");
    let output = temp.path().join("output_indel.vcf.gz");
    let index = temp.path().join("output_indel.vcf.gz.csi");
    std::fs::write(
        &input,
        concat!(
            "##fileformat=VCFv4.2\n",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">\n",
            "##contig=<ID=1,length=1000>\n",
            "##contig=<ID=2,length=1000>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "2\t25\t.\tA\tT\t100\tPASS\t.\n",
            "1\t20\t.\tC\tT\t100\tPASS\t.\n",
            "1\t10\t.\tA\tG\t100\tPASS\t.\n",
        ),
    )
    .unwrap();

    bcftools::sort_native(&input, &output, "z", true).unwrap();

    let mut decoder = flate2::read::MultiGzDecoder::new(std::fs::File::open(&output).unwrap());
    let mut vcf = String::new();
    decoder.read_to_string(&mut vcf).unwrap();
    let records = vcf
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(
        records,
        vec![
            "1\t10\t.\tA\tG\t100\tPASS\t.",
            "1\t20\t.\tC\tT\t100\tPASS\t.",
            "2\t25\t.\tA\tT\t100\tPASS\t.",
        ]
    );
    assert!(std::fs::metadata(index).unwrap().len() > 0);
}

#[test]
fn bcftools_native_sort_reports_invalid_input_errors() {
    let temp = tempfile::tempdir().unwrap();
    let input = temp.path().join("malformed.vcf");
    let output = temp.path().join("out.vcf.gz");
    std::fs::write(&input, "not a vcf\n").unwrap();

    let err = bcftools::sort_native(&input, &output, "z", true).unwrap_err();

    assert!(
        matches!(err, LibError::InvalidArguments(message) if message.contains("bcftools.sort failed"))
    );
    assert!(!output.exists());
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

    // Focused port of pysam AlignmentFile fetch coordinate behavior:
    // reversed coordinates are rejected before backend I/O.
    let err = file.fetch("22", Some(10), Some(9)).unwrap_err();
    assert!(err.to_string().contains("stop must be >= start"));

    // Focused port of pysam AlignmentFile fetch mode behavior:
    // CRAM fetches need an explicit reference source.
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

    // Focused port of pysam's invalid-contig fetch behavior: unknown
    // references surface as errors rather than empty successful iterators.
    let err = file
        .fetch("missing_chr", Some(999), Some(1001))
        .unwrap_err();
    assert!(
        err.to_string().contains("invalid reference sequence"),
        "{err}"
    );
}

#[test]
fn pysam_fetch_routes_bam_to_native_indexed_backend() {
    let file = AlignmentFile::open(
        "missing.bam",
        "rb",
        None,
        Some(PathBuf::from("missing.bam.bai")),
    )
    .unwrap();
    let err = file.fetch("chr_test", Some(999), Some(1001)).unwrap_err();
    assert!(
        err.to_string().contains("missing associated index"),
        "{err}"
    );
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
    // Ported from pyfaidx test_Fasta_integer_index.py's invalid-key behavior:
    // a missing contig should fail explicitly.
    assert!(fasta.get("missing").is_err());
}

#[test]
fn vcf_direction_is_pysam_variant_file_first() {
    assert_eq!(chosen_initial_surface(), VcfDirection::PysamVariantFile);
    assert!(bioscript_libs::vcf::open_variant_file().is_err());
}

#[test]
fn vcf_reads_kestrel_records_without_metadata() {
    let records = parse_kestrel_vcf(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nMUC1\t100\t.\tC\tCGGCA\t.\tPASS\t.\tGT\tDel:120:10000\n",
    )
    .unwrap();
    assert_eq!(records.len(), 1);
    assert_eq!(records[0].get("CHROM").map(String::as_str), Some("MUC1"));
    assert_eq!(
        records[0].get("Sample").map(String::as_str),
        Some("Del:120:10000")
    );
}

#[test]
fn kestrel_vntyper_command_uses_structured_argv() {
    let config = KestrelRunConfig::vntyper(
        "kestrel.jar",
        "muc1.fa",
        "out.vcf",
        "out.sam",
        "tmp",
        "sample1",
        "r1.fastq.gz",
        "r2.fastq.gz",
    );
    let command = config.command().unwrap();
    assert_eq!(command.program(), "java");
    assert_eq!(
        command.argv(),
        vec![
            "java",
            "-Xmx12g",
            "-jar",
            "kestrel.jar",
            "-k",
            "20",
            "--maxalignstates",
            "40",
            "--maxhapstates",
            "40",
            "-r",
            "muc1.fa",
            "-o",
            "out.vcf",
            "-ssample1",
            "r1.fastq.gz",
            "r2.fastq.gz",
            "--hapfmt",
            "sam",
            "-p",
            "out.sam",
            "--logstderr",
            "--logstdout",
            "--loglevel",
            "INFO",
            "--temploc",
            "tmp",
        ]
    );
}

#[test]
fn kestrel_native_adapter_calls_vendored_kestrel_rs_for_sequences() {
    let mut options = NativeKestrelRunOptions::new("sample1");
    options.minimum_difference = 1;
    options.max_haplotypes = 4;
    options.max_saved_states = 4;

    let vcf = call_sequences_to_vcf(
        "chr1",
        "AAAACCCCGGGGTTTT",
        ["AAAATCCCGGGGTTTT"; 5],
        4,
        &options,
    )
    .unwrap();

    assert!(vcf.contains("##fileformat=VCFv4.2\n"));
    assert!(vcf.contains("##contig=<ID=chr1,length=16"));
    assert!(vcf.contains("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1"));
    assert!(vcf.contains("chr1\t5\t.\tC\tT"), "{vcf}");
}

#[test]
fn kestrel_native_adapter_accepts_gzipped_fastq_for_kestrel_rs() {
    let temp = tempfile::tempdir().unwrap();
    let fastq = temp.path().join("reads.fastq.gz");
    let file = std::fs::File::create(&fastq).unwrap();
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    for index in 0..5 {
        writeln!(encoder, "@r{index}\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII").unwrap();
    }
    encoder.finish().unwrap();

    let mut options = NativeKestrelRunOptions::new("sample1");
    options.minimum_difference = 1;
    options.max_haplotypes = 4;
    options.max_saved_states = 4;

    let vcf = call_fastq_paths_to_vcf_references(
        &[NativeReferenceRegion::new(
            "chr1",
            "AAAACCCCGGGGTTTT",
            "2a9fd43653a81f9ec44e34c7ec038636",
        )],
        [fastq.as_path()],
        4,
        &options,
    )
    .unwrap();

    assert!(vcf.contains("##fileformat=VCFv4.2\n"));
    assert!(vcf.contains("##contig=<ID=chr1,length=16"));
    assert!(vcf.contains("chr1\t5\t.\tC\tT"), "{vcf}");
}

#[test]
fn samtools_vntyper_subset_builds_allowed_commands() {
    let view = samtools::view(
        PathBuf::from("sample.bam").as_path(),
        "chr1:1-10",
        PathBuf::from("slice.bam").as_path(),
    )
    .unwrap();
    assert_eq!(
        view.argv(),
        vec![
            "samtools",
            "view",
            "-b",
            "sample.bam",
            "chr1:1-10",
            "-o",
            "slice.bam"
        ]
    );

    let view = samtools::view_region(
        PathBuf::from("sample.bam").as_path(),
        "chr1:1-10",
        PathBuf::from("slice.bam").as_path(),
        false,
    )
    .unwrap();
    assert_eq!(
        view.argv(),
        vec![
            "samtools",
            "view",
            "-b",
            "sample.bam",
            "chr1:1-10",
            "-o",
            "slice.bam"
        ]
    );

    let fastq = samtools::fastq(
        PathBuf::from("slice.bam").as_path(),
        PathBuf::from("r1.fastq.gz").as_path(),
        PathBuf::from("r2.fastq.gz").as_path(),
    )
    .unwrap();
    assert_eq!(fastq.program(), "samtools");
    assert_eq!(fastq.args()[0], "fastq");

    let sorted = samtools::sort(
        PathBuf::from("slice.bam").as_path(),
        PathBuf::from("slice.name.bam").as_path(),
        true,
    )
    .unwrap();
    assert_eq!(
        sorted.argv(),
        vec![
            "samtools",
            "sort",
            "-n",
            "-o",
            "slice.name.bam",
            "slice.bam"
        ]
    );

    let faidx = samtools::faidx(PathBuf::from("ref.fa").as_path()).unwrap();
    assert_eq!(faidx.argv(), vec!["samtools", "faidx", "ref.fa"]);
}

#[test]
fn samtools_native_adapter_handles_tiny_indexed_bam() {
    let temp = tempfile::tempdir().unwrap();
    let sam = temp.path().join("tiny.sam");
    let bam = temp.path().join("tiny.bam");
    let slice = temp.path().join("slice.bam");
    let r1 = temp.path().join("r1.fastq.gz");
    let r2 = temp.path().join("r2.fastq.gz");
    std::fs::write(
        &sam,
        concat!(
            "@HD\tVN:1.6\tSO:coordinate\n",
            "@SQ\tSN:chr1\tLN:8\n",
            "pair\t65\tchr1\t1\t60\t4M\t=\t5\t8\tACGT\t!!!!\n",
            "pair\t129\tchr1\t5\t60\t4M\t=\t1\t-8\tTGCA\t####\n",
        ),
    )
    .unwrap();
    htslib_rs::alignment_compat::write_bam_from_sam_path(
        &sam,
        std::fs::File::create(&bam).unwrap(),
    )
    .unwrap();
    samtools_rs::native::index(&bam, Option::<&PathBuf>::None, Some(1)).unwrap();

    let records_written =
        samtools::view_region_native(&bam, None, None, None, "chr1:1-4", &slice).unwrap();
    assert_eq!(records_written, 0);
    assert!(std::fs::metadata(&slice).unwrap().len() > 0);

    let depth = samtools::depth_native(&bam, None, "chr1:1-8").unwrap();
    assert_eq!(depth.region_length, 8);
    assert_eq!(depth.uncovered_bases, 0);
    assert_eq!(depth.min, 1);
    assert_eq!(depth.max, 1);
    assert_eq!(depth.mean, 1.0);
    assert_eq!(depth.median, 1.0);

    let fastq = samtools::fastq_native(&bam, None, None, None, "chr1:1-4", &r1, &r2).unwrap();
    assert_eq!(fastq.read1_records, 1);
    assert_eq!(fastq.read2_records, 1);
    assert_eq!(fastq.skipped_records, 0);

    let err = samtools::depth_native(&bam, None, "chr1:8-1").unwrap_err();
    assert!(err.to_string().contains("region"), "{err}");
}
