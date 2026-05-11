use std::{fs, io::Write, path::PathBuf};

use bioscript_libs::{
    LibError, ModuleName, bcftools,
    kestrel::{
        KestrelRunConfig,
        native::{
            ActiveRegion, ActiveRegionDetectorConfig, AlignmentOp, AlignmentWeight,
            HaplotypeAssemblyConfig, HaplotypeEvidence, KestrelVcfWriter, KmerCountMap,
            NativeKestrelCallConfig, NativeReferenceRegion, NativeVariantCall, ReferenceRegion,
            ReferenceSequence, RegionStats, VariantCall, align_haplotype, assemble_haplotypes,
            call_alignment_variants, call_assembled_haplotypes_to_vcf,
            call_counted_kmers_to_vcf_references, call_explicit_haplotypes_to_vcf,
            call_fastq_paths_to_vcf, call_fastq_paths_to_vcf_references, call_sequences_to_vcf,
            count_fastq_kmers, count_sequence_kmers, detect_active_regions, difference_threshold,
            read_reference_records, recovery_threshold, reference_kmers, scan_limit_length,
            score_haplotype_alignment,
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
        err.to_string().contains("failed to read BAM index"),
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
fn kestrel_native_vcf_writer_matches_java_writer_surface() {
    let mut writer = KestrelVcfWriter::new(
        "1.0.2",
        vec![ReferenceSequence {
            name: "MUC1".to_owned(),
            length: 120,
            md5: "abc123".to_owned(),
        }],
    );
    writer.add_sample("sample1").unwrap();
    writer.add_sample("sample2").unwrap();
    writer
        .add_variant(VariantCall {
            sample_name: "sample2".to_owned(),
            chrom: "MUC1".to_owned(),
            pos: 21,
            ref_allele: "T".to_owned(),
            alt_allele: "G".to_owned(),
            variant_depth: 7,
            locus_depth: 100,
        })
        .unwrap();

    assert_eq!(
        writer.to_vcf_string(),
        concat!(
            "##fileformat=VCF4.2\n",
            "##source=Kestrel1.0.2\n",
            "##contig=<ID=MUC1,length=120,md5=abc123>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "##FORMAT=<ID=GDP,Number=A,Type=Integer,Description=\"Estimated depth of all haplotypes supporting the alternate variant\">\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Estimated depth of all haplotypes in the variant active region\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\tsample2\n",
            "MUC1\t21\t.\tT\tG\t.\t.\t.\tGT:GDP:DP\t0:.:.\t1:7:100\n",
        )
    );
    assert!(writer.add_sample("bad sample").is_err());
}

#[test]
fn kestrel_native_variants_use_java_vcf_normalization_rules() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTACGT".to_owned(),
    };
    let snp = NativeVariantCall::snp("sample1", 3, "G", "T", 4, 10)
        .to_vcf_call(&region)
        .unwrap();
    assert_eq!(
        (snp.pos, snp.ref_allele.as_str(), snp.alt_allele.as_str()),
        (3, "G", "T")
    );

    let insertion = NativeVariantCall::insertion("sample1", 4, "AA", 5, 10)
        .to_vcf_call(&region)
        .unwrap();
    assert_eq!(
        (
            insertion.pos,
            insertion.ref_allele.as_str(),
            insertion.alt_allele.as_str()
        ),
        (3, "G", "GAA")
    );

    let start_insertion = NativeVariantCall::insertion("sample1", 1, "TT", 5, 10)
        .to_vcf_call(&region)
        .unwrap();
    assert_eq!(
        (
            start_insertion.pos,
            start_insertion.ref_allele.as_str(),
            start_insertion.alt_allele.as_str()
        ),
        (1, "A", "TTA")
    );

    let deletion = NativeVariantCall::deletion("sample1", 4, "TA", 6, 10)
        .to_vcf_call(&region)
        .unwrap();
    assert_eq!(
        (
            deletion.pos,
            deletion.ref_allele.as_str(),
            deletion.alt_allele.as_str()
        ),
        (3, "GTA", "G")
    );
}

#[test]
fn kestrel_native_kmer_count_map_counts_canonical_bases() {
    let counts = count_sequence_kmers("ACGTACGTA", 3).unwrap();
    assert_eq!(counts.get("ACG"), Some(&2));
    assert_eq!(counts.get("CGT"), Some(&2));
    assert_eq!(counts.get("GTA"), Some(&2));
    assert_eq!(counts.get("TAC"), Some(&1));

    let map = KmerCountMap::from_sequences(["acgtnacgt", "ACGT"], 4).unwrap();
    assert_eq!(map.kmer_size(), 4);
    assert_eq!(map.get("ACGT").unwrap(), 3);
    assert_eq!(map.get("CGTA").unwrap(), 0);
    assert_eq!(map.transition_count("ACGT", "CGTN").is_err(), true);
    assert_eq!(map.transition_count("ACGT", "CGTA").unwrap(), 0);
    assert!(map.get("ACGN").is_err());

    let transitions = KmerCountMap::from_sequences(["AACCG"], 3).unwrap();
    assert_eq!(transitions.transition_count("AAC", "ACC").unwrap(), 1);
    assert_eq!(transitions.transition_count("AAC", "CCG").unwrap(), 0);
}

#[test]
fn kestrel_native_kmer_count_map_validates_inputs() {
    assert!(count_sequence_kmers("ACGT", 0).is_err());
    assert!(count_sequence_kmers("ACGX", 3).is_err());

    let map = KmerCountMap::from_sequences(["ACGT"], 3).unwrap();
    assert!(map.get("AC").is_err());
}

#[test]
fn kestrel_native_kmer_count_map_reads_fastq_inputs() {
    let dir = std::env::temp_dir().join(format!(
        "bioscript-kestrel-kmer-test-{}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let plain_path = dir.join("reads.fastq");
    fs::write(
        &plain_path,
        b"@r1\nACGTAC\n+\nIIIIII\n@r2\nTTNNAC\n+\nIIIIII\n",
    )
    .unwrap();
    let gz_path = dir.join("reads.fastq.gz");
    {
        let file = fs::File::create(&gz_path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        encoder.write_all(b"@r3\nACGT\n+\nIIII\n").unwrap();
        encoder.finish().unwrap();
    }

    let map = KmerCountMap::from_fastq_paths([plain_path.as_path(), gz_path.as_path()], 3).unwrap();
    assert_eq!(map.get("ACG").unwrap(), 2);
    assert_eq!(map.get("CGT").unwrap(), 2);
    assert_eq!(map.get("GTA").unwrap(), 1);
    assert_eq!(map.get("TAC").unwrap(), 1);
    assert_eq!(count_fastq_kmers(&plain_path, 3).unwrap().get("TTA"), None);

    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn kestrel_native_ports_upstream_reference_reader_resources() {
    let cases = [
        ("general.us-ascii.fasta", 10, 3000),
        ("general.us-ascii.fastq", 10, 3000),
        ("allchars.us-ascii.fasta", 20, 2000),
        ("allchars.us-ascii.fastq", 20, 2000),
    ];

    for (file_name, expected_records, expected_len) in cases {
        let records = read_reference_records(&kestrel_refreader_fixture(file_name)).unwrap();
        assert_eq!(records.len(), expected_records, "{file_name}");
        assert_eq!(records[0].name, "Seq-1", "{file_name}");
        assert_eq!(records[0].sequence.len(), expected_len, "{file_name}");
        assert_eq!(
            records.last().unwrap().sequence.len(),
            expected_len,
            "{file_name}"
        );

        for kmer_size in [1, 2, 21, 32, 64] {
            let kmers = reference_kmers(&records[0].sequence, kmer_size).unwrap();
            assert_eq!(kmers.len(), expected_len - kmer_size + 1, "{file_name}");
            assert!(kmers.iter().all(|kmer| kmer.len() == kmer_size));
            assert!(kmers.iter().all(|kmer| {
                kmer.bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            }));
        }
    }
}

#[test]
fn kestrel_native_reference_kmers_match_upstream_ambiguous_base_shape() {
    assert_eq!(
        reference_kmers("AUn.-r", 2).unwrap(),
        vec!["AT", "TA", "AC", "CG", "GT"]
    );
}

#[test]
fn kestrel_native_region_stats_match_java_percentiles() {
    let stats = RegionStats::from_counts(&[10, 4, 8, 2, 6], 0, 5).unwrap();
    assert_eq!(stats.min, 2);
    assert_eq!(stats.pct25, 4.0);
    assert_eq!(stats.pct50, 6.0);
    assert_eq!(stats.pct75, 8.0);
    assert_eq!(stats.max, 10);
    assert_eq!(stats.n, 5);

    let interpolated = RegionStats::from_counts(&[10, 20, 30, 40], 0, 4).unwrap();
    assert_eq!(interpolated.pct25, 17.5);
    assert_eq!(interpolated.pct50, 25.0);
    assert_eq!(interpolated.pct75, 32.5);
    assert!(RegionStats::from_counts(&[1], 1, 1).is_err());
}

#[test]
fn kestrel_native_active_region_tracks_anchors_and_stats() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTACGT".to_owned(),
    };
    let active = ActiveRegion::new(&region, Some(1), Some(4), &[5, 10, 20, 30, 40, 50], 3).unwrap();
    assert_eq!(active.reference_name, "MUC1");
    assert_eq!(active.start_index, 1);
    assert_eq!(active.end_index, 6);
    assert_eq!(active.left_end_kmer.as_deref(), Some("CGT"));
    assert_eq!(active.right_end_kmer.as_deref(), Some("ACG"));
    assert!(active.matches_left_end("CGT"));
    assert!(active.matches_right_end("ACG"));
    assert_eq!(active.stats.n, 3);
    assert_eq!(active.stats.min, 10);
    assert_eq!(active.stats.max, 30);

    let left_open = ActiveRegion::new(&region, None, Some(3), &[5, 10, 20, 30, 40, 50], 3).unwrap();
    assert!(left_open.left_end);
    assert_eq!(left_open.left_end_kmer, None);

    assert!(ActiveRegion::new(&region, Some(2), Some(2), &[5, 10, 20], 3).is_err());
}

#[test]
fn kestrel_native_reference_counts_support_detector_inputs() {
    let map = KmerCountMap::from_sequences(["AAAACCCCGGGGTTTT"], 4).unwrap();
    assert_eq!(
        map.reference_counts("AAAANCCCC", false).unwrap(),
        vec![1, 0, 0, 0, 0, 1]
    );

    let reverse = KmerCountMap::from_sequences(["AAAA"], 4).unwrap();
    assert_eq!(reverse.reference_counts("TTTT", true).unwrap(), vec![1]);
}

#[test]
fn kestrel_native_active_region_detector_finds_depth_drop_candidates() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(
        [
            "AAAA", "AAAC", "AACC", "ACCC", "GGGT", "GGTT", "GTTT", "TTTT",
        ],
        4,
    )
    .unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 1,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: true,
        decay_min: 1.0,
        decay_alpha: 0.80,
        peak_scan_length: 7,
        scan_limit_factor: 7.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    let detection = detect_active_regions(&region, &counts, &config).unwrap();
    assert_eq!(detection.difference_threshold, 1);
    assert_eq!(
        detection.reference_counts,
        vec![1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1]
    );
    assert_eq!(detection.regions.len(), 1);
    let active = &detection.regions[0];
    assert_eq!(active.start_kmer_index, 3);
    assert_eq!(active.end_kmer_index, 9);
    assert_eq!(active.left_end_kmer.as_deref(), Some("ACCC"));
    assert_eq!(active.right_end_kmer.as_deref(), Some("GGGT"));
}

#[test]
fn kestrel_native_active_region_detector_emits_right_open_candidates() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(["AAAA", "AAAC", "AACC", "ACCC"], 4).unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 1,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: false,
        decay_min: 1.0,
        decay_alpha: 0.80,
        peak_scan_length: 7,
        scan_limit_factor: 7.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    let detection = detect_active_regions(&region, &counts, &config).unwrap();
    assert_eq!(detection.regions.len(), 1);
    let active = &detection.regions[0];
    assert_eq!(active.start_kmer_index, 3);
    assert_eq!(active.end_kmer_index, 12);
    assert_eq!(active.left_end_kmer.as_deref(), Some("ACCC"));
    assert_eq!(active.right_end_kmer, None);
    assert_eq!(active.end_index, 15);
}

#[test]
fn kestrel_native_active_region_detector_respects_anchor_both_ends() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(["AAAA", "AAAC", "AACC", "ACCC"], 4).unwrap();

    let detection = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: true,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
    )
    .unwrap();
    assert!(detection.regions.is_empty());
}

#[test]
fn kestrel_native_active_region_detector_emits_left_open_candidates() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(["GGGT", "GGTT", "GTTT", "TTTT"], 4).unwrap();
    let detection = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: false,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
    )
    .unwrap();

    assert_eq!(
        detection.reference_counts,
        vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1]
    );
    assert_eq!(detection.regions.len(), 1);
    let active = &detection.regions[0];
    assert!(active.left_end);
    assert_eq!(active.left_end_kmer, None);
    assert_eq!(active.right_end_kmer.as_deref(), Some("GGGT"));
    assert_eq!(active.end_kmer_index, 9);
}

#[test]
fn kestrel_native_active_region_detector_scans_past_short_peaks() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let mut read_kmers = Vec::new();
    for kmer in [
        "AAAA", "ACCC", "CGGG", "GGGG", "GGGT", "GGTT", "GTTT", "TTTT",
    ] {
        for _ in 0..5 {
            read_kmers.push(kmer);
        }
    }
    let counts = KmerCountMap::from_sequences(read_kmers, 4).unwrap();

    let without_peak_scan = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: true,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 0,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
    )
    .unwrap();
    assert_eq!(without_peak_scan.regions.len(), 1);
    assert_eq!(without_peak_scan.regions[0].start_kmer_index, 3);

    let with_peak_scan = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: true,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
    )
    .unwrap();
    assert_eq!(with_peak_scan.regions.len(), 1);
    let active = &with_peak_scan.regions[0];
    assert_eq!(active.start_kmer_index, 0);
    assert_eq!(active.end_kmer_index, 7);
    assert_eq!(active.left_end_kmer.as_deref(), Some("AAAA"));
    assert_eq!(active.right_end_kmer.as_deref(), Some("CGGG"));
}

#[test]
fn kestrel_native_active_region_detector_discards_over_limit_scans() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(["AAAA"], 4).unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 1,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: false,
        decay_min: 1.0,
        decay_alpha: 0.80,
        peak_scan_length: 0,
        scan_limit_factor: 1.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    assert_eq!(scan_limit_length(4, &config).unwrap(), 4);
    assert_eq!(
        scan_limit_length(
            4,
            &ActiveRegionDetectorConfig {
                max_gap_size: 3,
                ..config.clone()
            }
        )
        .unwrap(),
        7
    );
    let detection = detect_active_regions(&region, &counts, &config).unwrap();
    assert!(detection.regions.is_empty());

    assert!(
        scan_limit_length(
            4,
            &ActiveRegionDetectorConfig {
                scan_limit_factor: f32::INFINITY,
                max_gap_size: 0,
                recover_right_anchor: true,
                call_ambiguous_regions: true,
                ..config
            }
        )
        .is_err()
    );
}

#[test]
fn kestrel_native_alignment_weight_matches_java_gap_limit_shape() {
    let default_weight = AlignmentWeight::default();
    assert_eq!(default_weight.initial_score(4).unwrap(), 40.0);
    assert_eq!(default_weight.max_exclusive_gap_size(4).unwrap(), 0);
    assert_eq!(default_weight.max_exclusive_gap_size(20).unwrap(), 40);
    assert_eq!(
        scan_limit_length(
            20,
            &ActiveRegionDetectorConfig {
                scan_limit_factor: 7.0,
                max_gap_size: default_weight.max_exclusive_gap_size(20).unwrap(),
                ..ActiveRegionDetectorConfig::default()
            }
        )
        .unwrap(),
        180
    );

    let custom_weight = AlignmentWeight::new(-8.0, 2.0, 12.0, 3.0, 0.0).unwrap();
    assert_eq!(custom_weight.match_weight, 8.0);
    assert_eq!(custom_weight.mismatch, -2.0);
    assert_eq!(custom_weight.gap_open, -12.0);
    assert_eq!(custom_weight.gap_extend, -3.0);
    assert_eq!(custom_weight.max_exclusive_gap_size(4).unwrap(), 6);
    assert!(AlignmentWeight::new(0.0, -1.0, -1.0, -1.0, 0.0).is_err());
}

#[test]
fn kestrel_native_alignment_weight_parses_java_weight_vectors() {
    assert_eq!(
        AlignmentWeight::parse(None).unwrap(),
        AlignmentWeight::default()
    );
    assert_eq!(
        AlignmentWeight::parse(Some("")).unwrap(),
        AlignmentWeight::default()
    );

    let parsed = AlignmentWeight::parse(Some("( -8, 2, 12, 3, -5 )")).unwrap();
    assert_eq!(
        parsed,
        AlignmentWeight {
            match_weight: 8.0,
            mismatch: -2.0,
            gap_open: -12.0,
            gap_extend: -3.0,
            init_score: 5.0,
        }
    );

    let partial = AlignmentWeight::parse(Some("[, -6, , -2]")).unwrap();
    assert_eq!(partial.match_weight, AlignmentWeight::DEFAULT_MATCH);
    assert_eq!(partial.mismatch, -6.0);
    assert_eq!(partial.gap_open, AlignmentWeight::DEFAULT_GAP_OPEN);
    assert_eq!(partial.gap_extend, -2.0);

    let integer_formats = AlignmentWeight::parse(Some("<0xA, 012, #28, 04, 0>")).unwrap();
    assert_eq!(integer_formats.match_weight, 10.0);
    assert_eq!(integer_formats.mismatch, -12.0);
    assert_eq!(integer_formats.gap_open, -40.0);
    assert_eq!(integer_formats.gap_extend, -4.0);
    assert!(AlignmentWeight::parse(Some("(1,2]")).is_err());
    assert!(AlignmentWeight::parse(Some("1,2,3,4,5,6")).is_err());
    assert!(AlignmentWeight::parse(Some("1,bad")).is_err());
}

#[test]
fn kestrel_native_active_region_detector_recovers_right_anchor() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let mut read_kmers = Vec::new();
    for _ in 0..20 {
        read_kmers.push("AAAA");
    }
    for _ in 0..8 {
        read_kmers.push("CCCG");
    }
    let counts = KmerCountMap::from_sequences(read_kmers, 4).unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 5,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: true,
        decay_min: 0.80,
        decay_alpha: 0.80,
        peak_scan_length: 0,
        scan_limit_factor: 7.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    let detection = detect_active_regions(&region, &counts, &config).unwrap();
    assert_eq!(detection.regions.len(), 1);
    let active = &detection.regions[0];
    assert_eq!(active.start_kmer_index, 0);
    assert_eq!(active.end_kmer_index, 5);
    assert_eq!(active.right_end_kmer.as_deref(), Some("CCCG"));

    let disabled = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            recover_right_anchor: false,
            call_ambiguous_regions: true,
            ..config
        },
    )
    .unwrap();
    assert!(disabled.regions.is_empty());
}

#[test]
fn kestrel_native_active_region_detector_skips_left_peak() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let mut read_kmers = Vec::new();
    for _ in 0..5 {
        read_kmers.push("CCCC");
    }
    for _ in 0..2 {
        read_kmers.push("CCCG");
    }
    let counts = KmerCountMap::from_sequences(read_kmers, 4).unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 5,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: false,
        decay_min: 1.0,
        decay_alpha: 0.80,
        peak_scan_length: 7,
        scan_limit_factor: 7.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    let detection = detect_active_regions(&region, &counts, &config).unwrap();
    assert!(detection.regions.is_empty());

    let without_peak_scan = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            peak_scan_length: 0,
            ..config
        },
    )
    .unwrap();
    assert_eq!(without_peak_scan.regions.len(), 1);
    assert!(without_peak_scan.regions[0].left_end);
    assert_eq!(without_peak_scan.regions[0].end_kmer_index, 4);
}

#[test]
fn kestrel_native_active_region_detector_limits_left_open_scans() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(["CCCG"], 4).unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 1,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: false,
        decay_min: 1.0,
        decay_alpha: 0.80,
        peak_scan_length: 0,
        scan_limit_factor: 1.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    let detection = detect_active_regions(&region, &counts, &config).unwrap();
    assert!(detection.regions.is_empty());

    let relaxed = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            ..config
        },
    )
    .unwrap();
    assert!(
        relaxed
            .regions
            .iter()
            .any(|region| region.left_end && region.end_kmer_index == 5)
    );
}

#[test]
fn kestrel_native_active_region_detector_discards_left_scan_recovery_before_left_end() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let mut read_kmers = Vec::new();
    for _ in 0..5 {
        read_kmers.push("AAAA");
        read_kmers.push("AAAC");
        read_kmers.push("ACCC");
    }
    let counts = KmerCountMap::from_sequences(read_kmers, 4).unwrap();
    let detection = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: false,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 0,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
    )
    .unwrap();

    assert_eq!(detection.reference_counts[..4], [5, 5, 0, 5]);
    assert!(
        detection
            .regions
            .iter()
            .all(|region| !(region.left_end && region.end_kmer_index == 3))
    );
}

#[test]
fn kestrel_native_active_region_detector_honors_ambiguous_region_flag() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCNGGGGTTTT".to_owned(),
    };
    let counts = KmerCountMap::from_sequences(
        [
            "AAAA", "AAAC", "AACC", "ACCC", "GGGG", "GGGT", "GGTT", "GTTT", "TTTT",
        ],
        4,
    )
    .unwrap();
    let config = ActiveRegionDetectorConfig {
        minimum_difference: 1,
        difference_quantile: 0.0,
        count_reverse_kmers: false,
        anchor_both_ends: true,
        decay_min: 1.0,
        decay_alpha: 0.80,
        peak_scan_length: 7,
        scan_limit_factor: 7.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
    };

    let allowed = detect_active_regions(&region, &counts, &config).unwrap();
    assert_eq!(allowed.regions.len(), 1);
    assert_eq!(allowed.regions[0].left_end_kmer.as_deref(), Some("ACCC"));
    assert_eq!(allowed.regions[0].right_end_kmer.as_deref(), Some("GGGG"));

    let rejected = detect_active_regions(
        &region,
        &counts,
        &ActiveRegionDetectorConfig {
            call_ambiguous_regions: false,
            ..config
        },
    )
    .unwrap();
    assert!(rejected.regions.is_empty());
}

#[test]
fn kestrel_native_difference_threshold_matches_java_quantile_shape() {
    assert_eq!(
        difference_threshold(&[10, 10, 1, 1, 10], 5, 0.90).unwrap(),
        6
    );
    assert_eq!(difference_threshold(&[10, 10], 5, 0.90).unwrap(), 5);
    assert!(difference_threshold(&[10, 10, 1], 0, 0.90).is_ok());
    assert!(difference_threshold(&[10, 10, 1], 1, 1.0).is_err());
}

#[test]
fn kestrel_native_recovery_threshold_matches_java_decay_shape() {
    let constant = ActiveRegionDetectorConfig {
        decay_min: 1.0,
        ..ActiveRegionDetectorConfig::default()
    };
    assert_eq!(
        recovery_threshold(200, 5, 48, 48, &constant).unwrap(),
        195.0
    );

    let decayed = ActiveRegionDetectorConfig {
        decay_min: 0.50,
        decay_alpha: 0.80,
        peak_scan_length: 7,
        scan_limit_factor: 7.0,
        max_gap_size: 0,
        recover_right_anchor: true,
        call_ambiguous_regions: true,
        ..ActiveRegionDetectorConfig::default()
    };
    assert_eq!(
        recovery_threshold(200, 5, 48, 48, &decayed).unwrap() as u32,
        180
    );
    assert_eq!(
        recovery_threshold(200, 5, 96, 48, &decayed).unwrap() as u32,
        164
    );
    assert!(
        recovery_threshold(
            200,
            5,
            48,
            48,
            &ActiveRegionDetectorConfig {
                decay_alpha: 1.0,
                ..decayed
            }
        )
        .is_err()
    );
}

#[test]
fn kestrel_native_alignment_emits_edit_operations() {
    let alignment = align_haplotype("ACGTAC", "ACGTTAC").unwrap();
    assert_eq!(
        alignment.ops,
        vec![
            AlignmentOp::Match(3),
            AlignmentOp::Insertion(1),
            AlignmentOp::Match(3)
        ]
    );

    let deletion = align_haplotype("ACGTAC", "ACAC").unwrap();
    assert_eq!(
        deletion.ops,
        vec![
            AlignmentOp::Match(2),
            AlignmentOp::Deletion(2),
            AlignmentOp::Match(2)
        ]
    );
    assert!(align_haplotype("ACGT", "ACGX").is_err());
}

#[test]
fn kestrel_native_alignment_scores_with_java_weight_shape() {
    let weight = AlignmentWeight::default();

    assert_eq!(
        score_haplotype_alignment("ACGTAC", "ACGTAC", &weight).unwrap(),
        60.0
    );
    assert_eq!(
        score_haplotype_alignment("ACGTAC", "ACGTTC", &weight).unwrap(),
        40.0
    );
    assert_eq!(
        score_haplotype_alignment("ACGTAC", "ACGTTAC", &weight).unwrap(),
        20.0
    );
    assert_eq!(
        score_haplotype_alignment("ACGTAC", "ACGTACAA", &weight).unwrap(),
        16.0
    );
}

#[test]
fn kestrel_native_alignment_calls_native_variants() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTACGT".to_owned(),
    };
    let alignment = align_haplotype("ACGTAC", "ATGTTAC").unwrap();
    let variants = call_alignment_variants("sample1", &alignment, 1, 6, 10).unwrap();
    assert_eq!(variants.len(), 2);

    let snp = variants[0].to_vcf_call(&region).unwrap();
    assert_eq!(
        (snp.pos, snp.ref_allele.as_str(), snp.alt_allele.as_str()),
        (2, "C", "T")
    );
    let insertion = variants[1].to_vcf_call(&region).unwrap();
    assert_eq!(
        (
            insertion.pos,
            insertion.ref_allele.as_str(),
            insertion.alt_allele.as_str()
        ),
        (3, "G", "GT")
    );
}

#[test]
fn kestrel_native_explicit_haplotype_engine_writes_vcf() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTAC".to_owned(),
    };
    let vcf = call_explicit_haplotypes_to_vcf(
        &region,
        &[HaplotypeEvidence {
            sequence: "ATGTTAC".to_owned(),
            variant_depth: 6,
            locus_depth: 10,
        }],
        &NativeKestrelCallConfig::new("native", "sample1", "md5"),
    )
    .unwrap();

    assert!(vcf.contains("##source=Kestrelnative\n"));
    assert!(vcf.contains("##contig=<ID=MUC1,length=6,md5=md5>\n"));
    assert!(vcf.contains("MUC1\t2\t.\tC\tT\t.\t.\t.\tGT:GDP:DP\t1:6:10\n"));
    assert!(vcf.contains("MUC1\t3\t.\tG\tGT\t.\t.\t.\tGT:GDP:DP\t1:6:10\n"));
}

#[test]
fn kestrel_native_haplotype_assembler_follows_counted_kmer_paths() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTAC".to_owned(),
    };
    let active = ActiveRegion::new(&region, Some(0), Some(3), &[10, 1, 1, 10], 3).unwrap();
    let counts = KmerCountMap::from_sequences(["ACGTTAC"], 3).unwrap();
    let haplotypes = assemble_haplotypes(
        &active,
        &counts,
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 20,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 10,
        },
    )
    .unwrap();

    assert_eq!(haplotypes.len(), 1);
    assert_eq!(haplotypes[0].sequence, "ACGTTAC");
    assert_eq!(haplotypes[0].variant_depth, 1);
    assert_eq!(haplotypes[0].locus_depth, 10);
}

#[test]
fn kestrel_native_haplotype_assembler_uses_total_active_region_depth() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTAC".to_owned(),
    };
    let active = ActiveRegion::new(&region, Some(0), Some(3), &[2, 2, 1, 2], 3).unwrap();
    let counts = KmerCountMap::from_sequences(["ACGTAC", "ACGTTAC"], 3).unwrap();
    let haplotypes = assemble_haplotypes(
        &active,
        &counts,
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 20,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 1,
        },
    )
    .unwrap();

    assert_eq!(haplotypes.len(), 2);
    assert!(
        haplotypes
            .iter()
            .all(|haplotype| haplotype.locus_depth == 2)
    );
    assert!(
        haplotypes
            .iter()
            .any(|haplotype| haplotype.sequence == "ACGTTAC" && haplotype.variant_depth == 1)
    );
}

#[test]
fn kestrel_native_haplotype_assembler_limits_repeated_kmers() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAAAA".to_owned(),
    };
    let active = ActiveRegion::new(&region, Some(0), Some(1), &[10, 10], 3).unwrap();
    let counts = KmerCountMap::from_sequences(["AAAAAA"], 3).unwrap();
    let no_repeats = assemble_haplotypes(
        &active,
        &counts,
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 8,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 10,
        },
    )
    .unwrap();
    assert!(no_repeats.is_empty());

    let one_repeat = assemble_haplotypes(
        &active,
        &counts,
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 8,
            max_repeat_count: 1,
            max_saved_states: 4,
            locus_depth: 10,
        },
    )
    .unwrap();
    assert_eq!(one_repeat.len(), 1);
    assert_eq!(one_repeat[0].sequence, "AAAA");
}

#[test]
fn kestrel_native_assembled_haplotype_engine_writes_vcf() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTAC".to_owned(),
    };
    let active = ActiveRegion::new(&region, Some(0), Some(3), &[10, 1, 1, 10], 3).unwrap();
    let counts = KmerCountMap::from_sequences(["ACGTTAC"], 3).unwrap();
    let vcf = call_assembled_haplotypes_to_vcf(
        &region,
        &active,
        &counts,
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 20,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 10,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "md5"),
    )
    .unwrap();

    assert!(vcf.contains("MUC1\t3\t.\tG\tGT\t.\t.\t.\tGT:GDP:DP\t1:1:10\n"));
}

#[test]
fn kestrel_native_assembled_haplotype_engine_prefers_alternate_over_reference_haplotype() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "ACGTAC".to_owned(),
    };
    let active = ActiveRegion::new(&region, Some(0), Some(3), &[2, 2, 1, 2], 3).unwrap();
    let counts = KmerCountMap::from_sequences(["ACGTAC", "ACGTTAC"], 3).unwrap();
    let vcf = call_assembled_haplotypes_to_vcf(
        &region,
        &active,
        &counts,
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 20,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 1,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "md5"),
    )
    .unwrap();

    assert!(vcf.contains("MUC1\t3\t.\tG\tGT\t.\t.\t.\tGT:GDP:DP\t1:1:2\n"));
}

#[test]
fn kestrel_native_sequences_engine_counts_detects_assembles_and_writes_vcf() {
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let vcf = call_sequences_to_vcf(
        &region,
        [
            "AAAA", "AAAC", "AACC", "ACCC", "CCCT", "CCTG", "CTGG", "TGGG", "GGGT", "GGTT", "GTTT",
            "TTTT",
        ],
        4,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: true,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 30,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 10,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "md5"),
    )
    .unwrap();

    assert!(vcf.contains("##contig=<ID=MUC1,length=16,md5=md5>\n"));
    assert!(vcf.contains("GT:GDP:DP\t1:1:10\n"));
}

#[test]
fn kestrel_native_fastq_engine_does_not_bridge_split_reads() {
    let dir = std::env::temp_dir().join(format!(
        "bioscript-kestrel-fastq-engine-test-{}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let fastq = dir.join("reads.fastq");
    fs::write(
        &fastq,
        b"@r1\nAAAACCC\n+\nIIIIIII\n@r2\nCCCTGGG\n+\nIIIIIII\n@r3\nGGGTTTT\n+\nIIIIIII\n",
    )
    .unwrap();
    let region = ReferenceRegion {
        reference_name: "MUC1".to_owned(),
        sequence: "AAAACCCCGGGGTTTT".to_owned(),
    };
    let vcf = call_fastq_paths_to_vcf(
        &region,
        [fastq.as_path()],
        4,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: true,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 4,
            max_bases: 30,
            max_repeat_count: 0,
            max_saved_states: 4,
            locus_depth: 10,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "md5"),
    )
    .unwrap();

    assert!(vcf.contains("##fileformat=VCF4.2\n"));
    assert!(
        !vcf.lines()
            .any(|line| !line.is_empty() && !line.starts_with('#'))
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn kestrel_native_multi_reference_engine_writes_all_contigs_and_calls_matching_region() {
    let references = vec![
        NativeReferenceRegion::new("REF1", "AAAACCCCGGGGTTTT", "md5-ref1"),
        NativeReferenceRegion::new("REF2", "ACAGTCCGTAAG", "md5-ref2"),
    ];
    let counts = KmerCountMap::from_sequences(["ACAGTTCGTAAG"; 5], 4).unwrap();
    let vcf = call_counted_kmers_to_vcf_references(
        &references,
        &counts,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: false,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 40,
            max_bases: 100,
            max_repeat_count: 0,
            max_saved_states: 40,
            locus_depth: 1,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "."),
    )
    .unwrap();

    assert!(vcf.contains("##contig=<ID=REF1,length=16,md5=md5-ref1>\n"));
    assert!(vcf.contains("##contig=<ID=REF2,length=12,md5=md5-ref2>\n"));
    assert!(vcf.contains("REF2\t6\t.\tC\tT\t.\t.\t.\tGT:GDP:DP\t1:5:5\n"));
    assert!(!vcf.contains("REF1\t"));
}

#[test]
fn kestrel_native_multi_reference_fastq_engine_reuses_counted_reads() {
    let dir = std::env::temp_dir().join(format!(
        "bioscript-kestrel-multiref-fastq-test-{}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    let fastq = dir.join("reads.fastq");
    fs::write(
        &fastq,
        b"@r1\nACAGTTCGTAAG\n+\nIIIIIIIIIIII\n@r2\nACAGTTCGTAAG\n+\nIIIIIIIIIIII\n",
    )
    .unwrap();
    let references = vec![NativeReferenceRegion::new("REF", "ACAGTCCGTAAG", "md5-ref")];
    let vcf = call_fastq_paths_to_vcf_references(
        &references,
        [fastq.as_path()],
        4,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: false,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: 40,
            max_bases: 100,
            max_repeat_count: 0,
            max_saved_states: 40,
            locus_depth: 1,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "."),
    )
    .unwrap();

    assert!(vcf.contains("##contig=<ID=REF,length=12,md5=md5-ref>\n"));
    assert!(vcf.contains("REF\t6\t.\tC\tT\t.\t.\t.\tGT:GDP:DP\t1:2:2\n"));
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn samtools_vntyper_subset_builds_allowed_commands() {
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
}

fn kestrel_refreader_fixture(file_name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../ports/vntyper/kestrel/bin/edu/gatech/kestrel/test/files/refreader")
        .join(file_name)
}
