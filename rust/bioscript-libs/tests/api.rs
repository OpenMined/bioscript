use std::{fs, io::Write, path::PathBuf};

use bioscript_libs::{
    LibError, ModuleName, bcftools,
    kestrel::{
        KestrelRunConfig,
        native::{
            ActiveRegion, ActiveRegionDetectorConfig, AlignmentOp, HaplotypeAssemblyConfig,
            HaplotypeEvidence, KestrelVcfWriter, KmerCountMap, NativeKestrelCallConfig,
            NativeVariantCall, ReferenceRegion, ReferenceSequence, RegionStats, VariantCall,
            align_haplotype, assemble_haplotypes, call_alignment_variants,
            call_assembled_haplotypes_to_vcf, call_explicit_haplotypes_to_vcf, count_fastq_kmers,
            count_sequence_kmers, detect_active_regions, difference_threshold,
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
    assert!(map.get("ACGN").is_err());
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
            locus_depth: 10,
        },
        &NativeKestrelCallConfig::new("native", "sample1", "md5"),
    )
    .unwrap();

    assert!(vcf.contains("MUC1\t3\t.\tG\tGT\t.\t.\t.\tGT:GDP:DP\t1:1:10\n"));
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
