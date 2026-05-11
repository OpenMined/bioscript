use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

use bioscript_libs::kestrel::native::{
    ActiveRegionDetectorConfig, AlignmentWeight, HaplotypeAssemblyConfig, NativeKestrelCallConfig,
    ReferenceRegion, call_fastq_paths_to_vcf,
};

const RUN_ENV: &str = "BIOSCRIPT_RUN_KESTREL_JAVA_PARITY";
const LONG_NONREPETITIVE_REFERENCE: &str =
    "ACGTTGCAACGAGTCCATGCTAGGCTAACCGTATCGGATCCGTAAGCTTGCAAGTCGATGCTAACGTTAGC";

#[test]
fn native_kestrel_fastq_output_matches_java_for_tiny_no_variant_fixture() {
    let dir = parity_temp_dir("tiny-no-variant");
    let fixture = KestrelParityFixture::new(
        "MUC1",
        "AAAACCCCGGGGTTTT",
        "2a9fd43653a81f9ec44e34c7ec038636",
        b"@r1\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n",
    );
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_tiny_snp_fixture() {
    let dir = parity_temp_dir("tiny-snp");
    let mut fastq = Vec::new();
    for read_index in 1..=5 {
        fastq.extend_from_slice(
            format!("@r{read_index}\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n").as_bytes(),
        );
    }
    let fixture = KestrelParityFixture::new(
        "MUC1",
        "AAAACCCCGGGGTTTT",
        "2a9fd43653a81f9ec44e34c7ec038636",
        &fastq,
    );
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_nonrepetitive_snp_fixture() {
    let dir = parity_temp_dir("nonrepetitive-snp");
    let mut fastq = Vec::new();
    for read_index in 1..=5 {
        fastq.extend_from_slice(
            format!("@r{read_index}\nACAGTTCGTAAG\n+\nIIIIIIIIIIII\n").as_bytes(),
        );
    }
    let fixture = KestrelParityFixture::new(
        "REF",
        "ACAGTCCGTAAG",
        "f17cc056a4c30b8661b5585d2641a37a",
        &fastq,
    );
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_adjacent_nonrepetitive_snps() {
    let dir = parity_temp_dir("adjacent-nonrepetitive-snps");
    let mut fastq = Vec::new();
    for read_index in 1..=5 {
        fastq.extend_from_slice(
            format!("@r{read_index}\nACAGTTTGTAAG\n+\nIIIIIIIIIIII\n").as_bytes(),
        );
    }
    let fixture = KestrelParityFixture::new(
        "REF",
        "ACAGTCCGTAAG",
        "f17cc056a4c30b8661b5585d2641a37a",
        &fastq,
    );
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_k20_nonrepetitive_snp() {
    let dir = parity_temp_dir("k20-nonrepetitive-snp");
    let reference = "ACGTTGCAACGAGTCCATGCTAGGCTAACCGTATCGGATCCGTAAGCT";
    let read = "ACGTTGCAACGAGTCCATGCTAGGTTAACCGTATCGGATCCGTAAGCT";
    let mut fastq = Vec::new();
    for read_index in 1..=5 {
        fastq.extend_from_slice(format!("@r{read_index}\n{read}\n+\n").as_bytes());
        fastq.extend_from_slice(format!("{}\n", "I".repeat(read.len())).as_bytes());
    }
    let fixture =
        KestrelParityFixture::new("REF", reference, "0f6b419f89dfba198188d4160b1c8329", &fastq)
            .with_kmer_size(20);
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_k20_nonrepetitive_deletion() {
    let dir = parity_temp_dir("k20-nonrepetitive-deletion");
    let read = "ACGTTGCAACGAGTCCATGCTAGGCTAACCGTACGGATCCGTAAGCTTGCAAGTCGATGCTAACGTTAGC";
    let fastq = repeated_fastq(read, 10);
    let fixture = KestrelParityFixture::new(
        "REF",
        LONG_NONREPETITIVE_REFERENCE,
        "e50386beaaf4c2113705c82a71502260",
        &fastq,
    )
    .with_kmer_size(20)
    .with_max_states(80);
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_k20_mixed_depth_deletion() {
    let dir = parity_temp_dir("k20-mixed-depth-deletion");
    let deletion_read = "ACGTTGCAACGAGTCCATGCTAGGCTAACCGTACGGATCCGTAAGCTTGCAAGTCGATGCTAACGTTAGC";
    let fastq = mixed_fastq(LONG_NONREPETITIVE_REFERENCE, 5, deletion_read, 5);
    let fixture = KestrelParityFixture::new(
        "REF",
        LONG_NONREPETITIVE_REFERENCE,
        "e50386beaaf4c2113705c82a71502260",
        &fastq,
    )
    .with_kmer_size(20)
    .with_max_states(80);
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_k20_nonrepetitive_insertion() {
    let dir = parity_temp_dir("k20-nonrepetitive-insertion");
    let read = "ACGTTGCAACGAGTCCATGCTAGGCTAACCGTTGATATCGGATCCGTAAGCTTGCAAGTCGATGCTAACGTTAGC";
    let fastq = repeated_fastq(read, 10);
    let fixture = KestrelParityFixture::new(
        "REF",
        LONG_NONREPETITIVE_REFERENCE,
        "e50386beaaf4c2113705c82a71502260",
        &fastq,
    )
    .with_kmer_size(20)
    .with_max_states(80);
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

#[test]
fn native_kestrel_fastq_output_matches_java_for_sparse_split_reads() {
    let dir = parity_temp_dir("sparse-split-reads");
    let fixture = KestrelParityFixture::new(
        "MUC1",
        "AAAACCCCGGGGTTTT",
        "2a9fd43653a81f9ec44e34c7ec038636",
        b"@r1\nAAAACCC\n+\nIIIIIII\n@r2\nCCCTGGG\n+\nIIIIIII\n@r3\nGGGTTTT\n+\nIIIIIII\n",
    );
    let (java_vcf, native_vcf) = run_java_and_native(&dir, &fixture);

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
}

struct KestrelParityFixture<'a> {
    reference_name: &'a str,
    reference_sequence: &'a str,
    reference_md5: &'a str,
    fastq_contents: &'a [u8],
    kmer_size: usize,
    max_states: usize,
}

impl<'a> KestrelParityFixture<'a> {
    fn new(
        reference_name: &'a str,
        reference_sequence: &'a str,
        reference_md5: &'a str,
        fastq_contents: &'a [u8],
    ) -> Self {
        Self {
            reference_name,
            reference_sequence,
            reference_md5,
            fastq_contents,
            kmer_size: 4,
            max_states: 40,
        }
    }

    fn with_kmer_size(mut self, kmer_size: usize) -> Self {
        self.kmer_size = kmer_size;
        self
    }

    fn with_max_states(mut self, max_states: usize) -> Self {
        self.max_states = max_states;
        self
    }
}

fn run_java_and_native(dir: &Path, fixture: &KestrelParityFixture<'_>) -> (String, String) {
    if std::env::var_os(RUN_ENV).is_none() {
        return (String::new(), String::new());
    }

    let jar = kestrel_jar();
    assert!(
        jar.exists(),
        "Kestrel Java parity gate requires {} or {} to exist: {}",
        RUN_ENV,
        "BIOSCRIPT_KESTREL_JAR",
        jar.display()
    );

    fs::create_dir_all(dir).unwrap();
    let reference_path = dir.join("ref.fa");
    let fastq_path = dir.join("reads.fq");
    let java_vcf_path = dir.join("java.vcf");
    let java_sam_path = dir.join("java.sam");

    fs::write(
        &reference_path,
        format!(
            ">{}\n{}\n",
            fixture.reference_name, fixture.reference_sequence
        ),
    )
    .unwrap();
    fs::write(&fastq_path, fixture.fastq_contents).unwrap();

    let status = Command::new("java")
        .arg("-Xmx512m")
        .arg("-jar")
        .arg(&jar)
        .arg("-k")
        .arg(fixture.kmer_size.to_string())
        .args([
            "--minsize",
            "4",
            "--mincount",
            "1",
            "--mindiff",
            "1",
            "--diffq",
            "0",
            "--decaymin",
            "1.0",
        ])
        .arg("--maxalignstates")
        .arg(fixture.max_states.to_string())
        .arg("--maxhapstates")
        .arg(fixture.max_states.to_string())
        .args(["--noanchorboth", "--nocountrev", "-r"])
        .arg(&reference_path)
        .arg("-o")
        .arg(&java_vcf_path)
        .arg("-ssample1")
        .arg(&fastq_path)
        .args(["--hapfmt", "sam", "-p"])
        .arg(&java_sam_path)
        .args(["--logstderr", "--loglevel", "ERROR", "--temploc"])
        .arg(&dir)
        .status()
        .unwrap();
    assert!(status.success(), "Java Kestrel exited with {status}");

    let java_vcf = fs::read_to_string(&java_vcf_path).unwrap();
    let native_vcf = call_fastq_paths_to_vcf(
        &ReferenceRegion {
            reference_name: fixture.reference_name.to_owned(),
            sequence: fixture.reference_sequence.to_owned(),
        },
        [fastq_path.as_path()],
        fixture.kmer_size,
        &ActiveRegionDetectorConfig {
            minimum_difference: 1,
            difference_quantile: 0.0,
            count_reverse_kmers: false,
            anchor_both_ends: false,
            decay_min: 1.0,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: AlignmentWeight::default()
                .max_exclusive_gap_size(fixture.kmer_size)
                .unwrap(),
            recover_right_anchor: true,
            call_ambiguous_regions: true,
        },
        &HaplotypeAssemblyConfig {
            min_kmer_count: 1,
            max_haplotypes: fixture.max_states,
            max_bases: 500,
            max_repeat_count: 0,
            max_saved_states: fixture.max_states,
            locus_depth: 1,
        },
        &NativeKestrelCallConfig::new("1.0.2", "sample1", fixture.reference_md5),
    )
    .unwrap();

    (java_vcf, native_vcf)
}

fn kestrel_jar() -> PathBuf {
    std::env::var_os("BIOSCRIPT_KESTREL_JAR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../..")
                .join("ports/vntyper/test-data/tools/kestrel/kestrel.jar")
        })
}

fn parity_temp_dir(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "bioscript-kestrel-java-parity-{name}-{}",
        std::process::id()
    ))
}

fn repeated_fastq(read: &str, copies: usize) -> Vec<u8> {
    let mut fastq = Vec::new();
    for read_index in 1..=copies {
        fastq.extend_from_slice(format!("@r{read_index}\n{read}\n+\n").as_bytes());
        fastq.extend_from_slice(format!("{}\n", "I".repeat(read.len())).as_bytes());
    }
    fastq
}

fn mixed_fastq(
    first_read: &str,
    first_copies: usize,
    second_read: &str,
    second_copies: usize,
) -> Vec<u8> {
    let mut fastq = Vec::new();
    for read_index in 1..=first_copies {
        fastq.extend_from_slice(format!("@ref{read_index}\n{first_read}\n+\n").as_bytes());
        fastq.extend_from_slice(format!("{}\n", "I".repeat(first_read.len())).as_bytes());
    }
    for read_index in 1..=second_copies {
        fastq.extend_from_slice(format!("@alt{read_index}\n{second_read}\n+\n").as_bytes());
        fastq.extend_from_slice(format!("{}\n", "I".repeat(second_read.len())).as_bytes());
    }
    fastq
}

fn variant_rows(vcf: &str) -> Vec<&str> {
    vcf.lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect()
}

fn header_without_source(vcf: &str) -> Vec<&str> {
    vcf.lines()
        .filter(|line| line.starts_with('#') && !line.starts_with("##source="))
        .collect()
}
