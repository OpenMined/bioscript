use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
};

use bioscript_libs::kestrel::native::{
    ActiveRegionDetectorConfig, HaplotypeAssemblyConfig, NativeKestrelCallConfig, ReferenceRegion,
    call_fastq_paths_to_vcf,
};

const RUN_ENV: &str = "BIOSCRIPT_RUN_KESTREL_JAVA_PARITY";

#[test]
fn native_kestrel_fastq_output_matches_java_for_tiny_no_variant_fixture() {
    if std::env::var_os(RUN_ENV).is_none() {
        return;
    }

    let jar = kestrel_jar();
    assert!(
        jar.exists(),
        "Kestrel Java parity gate requires {} or {} to exist: {}",
        RUN_ENV,
        "BIOSCRIPT_KESTREL_JAR",
        jar.display()
    );

    let dir = parity_temp_dir("tiny-no-variant");
    fs::create_dir_all(&dir).unwrap();
    let reference_path = dir.join("ref.fa");
    let fastq_path = dir.join("reads.fq");
    let java_vcf_path = dir.join("java.vcf");
    let java_sam_path = dir.join("java.sam");

    fs::write(&reference_path, b">MUC1\nAAAACCCCGGGGTTTT\n").unwrap();
    fs::write(&fastq_path, b"@r1\nAAAACCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n").unwrap();

    let status = Command::new("java")
        .arg("-Xmx512m")
        .arg("-jar")
        .arg(&jar)
        .args([
            "-k",
            "4",
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
            "--maxalignstates",
            "40",
            "--maxhapstates",
            "40",
            "--noanchorboth",
            "--nocountrev",
            "-r",
        ])
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
            reference_name: "MUC1".to_owned(),
            sequence: "AAAACCCCGGGGTTTT".to_owned(),
        },
        [fastq_path.as_path()],
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
            max_bases: 500,
            max_repeat_count: 0,
            max_saved_states: 40,
            locus_depth: 1,
        },
        &NativeKestrelCallConfig::new("1.0.2", "sample1", "2a9fd43653a81f9ec44e34c7ec038636"),
    )
    .unwrap();

    assert_eq!(variant_rows(&native_vcf), variant_rows(&java_vcf));
    assert_eq!(
        header_without_source(&native_vcf),
        header_without_source(&java_vcf)
    );
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
