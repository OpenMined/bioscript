use std::{
    env,
    io::Write as _,
    path::PathBuf,
    time::{Instant, SystemTime, UNIX_EPOCH},
};

use bioscript_core::Assembly;
use bioscript_formats::{DetectedKind, FileContainer, InspectOptions, inspect_bytes, inspect_file};
use zip::write::SimpleFileOptions;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("bioscript repo root")
        .to_path_buf()
}

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn shared_test_data_root() -> Option<PathBuf> {
    if let Some(path) = env::var_os("BIOSCRIPT_TEST_DATA_DIR") {
        let candidate = PathBuf::from(path);
        if candidate.exists() {
            return Some(candidate);
        }
    }

    let sibling = repo_root().join("test-data");
    if sibling.exists() {
        return Some(sibling);
    }

    env::var_os("HOME")
        .map(PathBuf::from)
        .map(|home| home.join(".bioscript/cache/test-data"))
        .filter(|path| path.exists())
}

fn shared_fixture_or_skip(test_name: &str, relative: &str) -> Option<PathBuf> {
    let root = shared_test_data_root()?;
    let path = root.join(relative);
    if !path.exists() {
        eprintln!("skipping {test_name}: missing {}", path.display());
        return None;
    }
    Some(path)
}

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-inspect-{label}-{}-{nanos}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn zip_bytes(entry_name: &str, contents: &[u8]) -> Vec<u8> {
    let cursor = std::io::Cursor::new(Vec::new());
    let mut writer = zip::ZipWriter::new(cursor);
    writer
        .start_file(entry_name, SimpleFileOptions::default())
        .unwrap();
    writer.write_all(contents).unwrap();
    writer.finish().unwrap().into_inner()
}

#[test]
fn inspect_bytes_handles_genotype_text() {
    let inspection = inspect_bytes(
        "sample.txt",
        b"rsid\tchromosome\tposition\tgenotype\n\
          rs73885319\t22\t36265860\tAG\n\
          rs60910145\t22\t36265900\tTG\n\
          rs71785313\t22\t36266005\tII\n",
        &InspectOptions::default(),
    )
    .unwrap();

    assert_eq!(inspection.container, FileContainer::Plain);
    assert_eq!(inspection.detected_kind, DetectedKind::GenotypeText);
}

#[test]
fn inspect_bytes_handles_vcf() {
    let inspection = inspect_bytes(
        "sample.vcf",
        b"##fileformat=VCFv4.2\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
          22\t36265860\trs73885319\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        &InspectOptions::default(),
    )
    .unwrap();

    assert_eq!(inspection.container, FileContainer::Plain);
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.phased, Some(false));
}

#[test]
fn inspect_bytes_handles_zip() {
    let bytes = zip_bytes(
        "nested/sample.txt",
        b"rsid\tchromosome\tposition\tgenotype\n\
          rs73885319\t22\t36265860\tAG\n\
          rs60910145\t22\t36265900\tTG\n\
          rs71785313\t22\t36266005\tII\n",
    );

    let inspection = inspect_bytes("sample.zip", &bytes, &InspectOptions::default()).unwrap();

    assert_eq!(inspection.container, FileContainer::Zip);
    assert_eq!(inspection.detected_kind, DetectedKind::GenotypeText);
    assert_eq!(
        inspection.selected_entry.as_deref(),
        Some("nested/sample.txt")
    );
}

#[test]
fn inspect_bytes_handles_unknown_bytes_conservatively() {
    let inspection = inspect_bytes(
        "sample.bin",
        b"this is not recognizable genotype or alignment data\n",
        &InspectOptions::default(),
    )
    .unwrap();

    assert_eq!(inspection.container, FileContainer::Plain);
    assert_eq!(inspection.detected_kind, DetectedKind::Unknown);
    assert!(
        inspection
            .warnings
            .contains(&"file did not match known textual heuristics".to_owned())
    );
}

#[test]
fn inspect_bytes_detects_alignment_reference_and_explicit_indexes() {
    let dir = temp_dir("bytes-explicit-indexes");
    let index_path = dir.join("sample.cram.crai");
    std::fs::write(&index_path, b"not crai").unwrap();
    let cram = inspect_bytes(
        "sample.cram",
        b"not actually decoded",
        &InspectOptions {
            input_index: Some(index_path.clone()),
            ..InspectOptions::default()
        },
    )
    .unwrap();
    assert_eq!(cram.detected_kind, DetectedKind::AlignmentCram);
    assert_eq!(cram.has_index, Some(index_path.exists()));
    assert_eq!(cram.index_path.as_deref(), Some(index_path.as_path()));
    assert!(cram.evidence.contains(&"extension .cram".to_owned()));

    let bam = inspect_bytes("sample.bam", b"bam bytes", &InspectOptions::default()).unwrap();
    assert_eq!(bam.detected_kind, DetectedKind::AlignmentBam);
    assert_eq!(
        bam.confidence,
        bioscript_formats::DetectionConfidence::Authoritative
    );

    let reference_index = PathBuf::from("ref.fa.fai");
    let reference = inspect_bytes(
        "ref.fa",
        b">chr1\nACGT\n",
        &InspectOptions {
            reference_index: Some(reference_index.clone()),
            ..InspectOptions::default()
        },
    )
    .unwrap();
    assert_eq!(reference.detected_kind, DetectedKind::ReferenceFasta);
    assert_eq!(
        reference.index_path.as_deref(),
        Some(reference_index.as_path())
    );
}

#[test]
fn inspect_bytes_zip_skips_macosx_entries_and_can_fallback_to_unknown_file() {
    let cursor = std::io::Cursor::new(Vec::new());
    let mut writer = zip::ZipWriter::new(cursor);
    writer
        .start_file("__MACOSX/._sample.txt", SimpleFileOptions::default())
        .unwrap();
    writer.write_all(b"ignored").unwrap();
    writer
        .start_file("nested/readme.bin", SimpleFileOptions::default())
        .unwrap();
    writer.write_all(b"not genotype text").unwrap();
    let bytes = writer.finish().unwrap().into_inner();

    let inspection = inspect_bytes("archive.zip", &bytes, &InspectOptions::default()).unwrap();

    assert_eq!(inspection.container, FileContainer::Zip);
    assert_eq!(inspection.detected_kind, DetectedKind::Unknown);
    assert_eq!(
        inspection.selected_entry.as_deref(),
        Some("nested/readme.bin")
    );
}

#[test]
fn inspect_bytes_rejects_empty_or_malformed_zip_archives() {
    let err = inspect_bytes("bad.zip", b"not a zip", &InspectOptions::default()).unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to read zip bytes"),
        "{err:?}"
    );

    let cursor = std::io::Cursor::new(Vec::new());
    let writer = zip::ZipWriter::new(cursor);
    let bytes = writer.finish().unwrap().into_inner();
    let err = inspect_bytes("empty.zip", &bytes, &InspectOptions::default()).unwrap_err();
    assert!(
        format!("{err:?}").contains("zip archive does not contain a supported file"),
        "{err:?}"
    );
}

#[test]
fn inspection_render_text_includes_empty_source_fields() {
    let inspection =
        inspect_bytes("sample.fa", b">chr1\nACGT\n", &InspectOptions::default()).unwrap();
    let text = inspection.render_text();

    assert!(text.contains("kind\treference_fasta"), "{text}");
    assert!(text.contains("vendor\t"), "{text}");
    assert!(text.contains("source_confidence\t"), "{text}");
    assert!(text.contains("reference_matches\t"), "{text}");
}

#[test]
fn ancestrydna_text_fixture_reports_vendor_platform_and_build() {
    let path = fixtures_dir().join("ancestrydna_v2_sample.txt");
    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();

    assert_eq!(inspection.container, FileContainer::Plain);
    assert_eq!(inspection.detected_kind, DetectedKind::GenotypeText);
    assert_eq!(inspection.assembly, Some(Assembly::Grch37));
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.vendor.as_deref()),
        Some("AncestryDNA")
    );
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.platform_version.as_deref()),
        Some("V2.0")
    );
    assert!(inspection.duration_ms < 1000);
}

#[test]
fn familytreedna_fixture_reports_vendor() {
    let path = fixtures_dir().join("familytreedna_sample.csv");
    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();

    assert_eq!(inspection.detected_kind, DetectedKind::GenotypeText);
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.vendor.as_deref()),
        Some("FamilyTreeDNA")
    );
}

#[test]
fn genesforgood_fixture_reports_vendor_platform_and_build() {
    let path = fixtures_dir().join("genesforgood_sample.txt");
    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();

    assert_eq!(inspection.detected_kind, DetectedKind::GenotypeText);
    assert_eq!(inspection.assembly, Some(Assembly::Grch37));
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.vendor.as_deref()),
        Some("Genes for Good")
    );
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.platform_version.as_deref()),
        Some("v1.2.1")
    );
}

#[test]
fn zipped_dynamicdna_fixture_reports_vendor_platform_and_build_quickly() {
    let Some(path) = shared_fixture_or_skip(
        "zipped_dynamicdna_fixture_reports_vendor_platform_and_build_quickly",
        "dynamicdna/100001-synthetic/100001_X_X_GSAv3-DTC_GRCh38-07-12-2025.txt.zip",
    ) else {
        return;
    };

    let started = Instant::now();
    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    let elapsed = started.elapsed().as_millis();

    assert_eq!(inspection.container, FileContainer::Zip);
    assert_eq!(inspection.detected_kind, DetectedKind::GenotypeText);
    assert_eq!(inspection.assembly, Some(Assembly::Grch38));
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.vendor.as_deref()),
        Some("Dynamic DNA")
    );
    assert_eq!(
        inspection
            .source
            .as_ref()
            .and_then(|source| source.platform_version.as_deref()),
        Some("GSAv3-DTC")
    );
    assert!(inspection.duration_ms < 1000);
    assert!(elapsed < 1000);
}

#[test]
#[allow(clippy::too_many_lines)]
fn shared_zipped_vendor_fixtures_report_expected_metadata() {
    struct Expectation {
        relative: &'static str,
        vendor: &'static str,
        platform_version: Option<&'static str>,
        assembly: Option<Assembly>,
    }

    let fixtures = [
        Expectation {
            relative: "23andme/v2/hu0199C8/23data20100526.txt.zip",
            vendor: "23andMe",
            platform_version: Some("v2"),
            assembly: None,
        },
        Expectation {
            relative: "23andme/v3/huE4DAE4/huE4DAE4_20120522224129.txt.zip",
            vendor: "23andMe",
            platform_version: Some("v3"),
            assembly: None,
        },
        Expectation {
            relative: "23andme/v4/huE18D82/genome__v4_Full_2016.txt.zip",
            vendor: "23andMe",
            platform_version: Some("v4"),
            assembly: Some(Assembly::Grch37),
        },
        Expectation {
            relative: "23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip",
            vendor: "23andMe",
            platform_version: Some("v5"),
            assembly: Some(Assembly::Grch37),
        },
        Expectation {
            relative: "dynamicdna/100001-synthetic/100001_X_X_GSAv3-DTC_GRCh38-07-12-2025.txt.zip",
            vendor: "Dynamic DNA",
            platform_version: Some("GSAv3-DTC"),
            assembly: Some(Assembly::Grch38),
        },
        Expectation {
            relative: "ancestrydna/huE922FC/AncestryDNA.txt.zip",
            vendor: "AncestryDNA",
            platform_version: Some("V2.0"),
            assembly: Some(Assembly::Grch37),
        },
        Expectation {
            relative: "familytreedna/hu17B792/2017-04-29_Family_Tree_DNA_Data.csv.zip",
            vendor: "FamilyTreeDNA",
            platform_version: None,
            assembly: None,
        },
        Expectation {
            relative: "genesforgood/hu80B047/GFG0_filtered_imputed_genotypes_noY_noMT_23andMe.txt.zip",
            vendor: "Genes for Good",
            platform_version: Some("v1.2.1"),
            assembly: Some(Assembly::Grch37),
        },
        Expectation {
            relative: "myheritage/hu33515F/MyHeritage_raw_dna_data.zip",
            vendor: "MyHeritage",
            platform_version: None,
            assembly: Some(Assembly::Grch37),
        },
    ];

    for fixture in fixtures {
        let Some(path) = shared_fixture_or_skip(
            "shared_zipped_vendor_fixtures_report_expected_metadata",
            fixture.relative,
        ) else {
            return;
        };

        let started = Instant::now();
        let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
        let elapsed = started.elapsed().as_millis();

        assert_eq!(
            inspection.container,
            FileContainer::Zip,
            "fixture {}",
            fixture.relative
        );
        assert_eq!(
            inspection.detected_kind,
            DetectedKind::GenotypeText,
            "fixture {}",
            fixture.relative
        );
        assert_eq!(
            inspection
                .source
                .as_ref()
                .and_then(|source| source.vendor.as_deref()),
            Some(fixture.vendor),
            "fixture {}",
            fixture.relative
        );
        assert_eq!(
            inspection
                .source
                .as_ref()
                .and_then(|source| source.platform_version.as_deref()),
            fixture.platform_version,
            "fixture {}",
            fixture.relative
        );
        assert_eq!(
            inspection.assembly, fixture.assembly,
            "fixture {}",
            fixture.relative
        );
        assert!(
            inspection.duration_ms < 1000,
            "fixture {}",
            fixture.relative
        );
        assert!(elapsed < 1000, "fixture {}", fixture.relative);
    }
}

#[test]
fn chr_y_cram_fixture_reports_index_without_decoding_entire_file() {
    let Some(path) = shared_fixture_or_skip(
        "chr_y_cram_fixture_reports_index_without_decoding_entire_file",
        "NA06985-chrY/aligned/NA06985.final.chrY.cram",
    ) else {
        return;
    };

    let started = Instant::now();
    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    let elapsed = started.elapsed().as_millis();

    assert_eq!(inspection.detected_kind, DetectedKind::AlignmentCram);
    assert_eq!(inspection.has_index, Some(true));
    assert!(inspection.index_path.is_some());
    assert!(inspection.duration_ms < 1000);
    assert!(elapsed < 1000);
}

#[test]
fn inspect_file_reports_bam_cram_and_reference_indexes_without_decoding() {
    let dir = temp_dir("index-detection");
    let cram = dir.join("sample.cram");
    let cram_index = dir.join("sample.cram.crai");
    let bam = dir.join("sample.bam");
    let bai = dir.join("sample.bam.bai");
    let fasta = dir.join("reference.fa");
    let fai = dir.join("reference.fa.fai");
    std::fs::write(&cram, b"not cram").unwrap();
    std::fs::write(&cram_index, b"not crai").unwrap();
    std::fs::write(&bam, b"not bam").unwrap();
    std::fs::write(&bai, b"not bai").unwrap();
    std::fs::write(&fasta, b">chr1\nACGT\n").unwrap();
    std::fs::write(&fai, b"chr1\t4\t6\t4\t5\n").unwrap();

    let cram_inspection = inspect_file(&cram, &InspectOptions::default()).unwrap();
    assert_eq!(cram_inspection.detected_kind, DetectedKind::AlignmentCram);
    assert_eq!(cram_inspection.has_index, Some(true));
    assert_eq!(
        cram_inspection.index_path.as_deref(),
        Some(cram_index.as_path())
    );

    let bam_inspection = inspect_file(&bam, &InspectOptions::default()).unwrap();
    assert_eq!(bam_inspection.detected_kind, DetectedKind::AlignmentBam);
    assert_eq!(bam_inspection.has_index, Some(true));
    assert_eq!(bam_inspection.index_path.as_deref(), Some(bai.as_path()));

    let fasta_inspection = inspect_file(&fasta, &InspectOptions::default()).unwrap();
    assert_eq!(fasta_inspection.detected_kind, DetectedKind::ReferenceFasta);
    assert_eq!(fasta_inspection.has_index, Some(true));
    assert_eq!(fasta_inspection.index_path.as_deref(), Some(fai.as_path()));
}

#[test]
fn inspect_file_reports_missing_short_cram_and_bam_indexes() {
    let dir = temp_dir("missing-index-detection");
    let cram = dir.join("sample.cram");
    let bam = dir.join("sample.bam");
    std::fs::write(&cram, b"not cram").unwrap();
    std::fs::write(&bam, b"not bam").unwrap();

    let cram_inspection = inspect_file(&cram, &InspectOptions::default()).unwrap();
    assert_eq!(cram_inspection.has_index, Some(false));
    assert_eq!(
        cram_inspection.index_path.as_deref(),
        Some(dir.join("sample.crai").as_path())
    );

    let bam_inspection = inspect_file(&bam, &InspectOptions::default()).unwrap();
    assert_eq!(bam_inspection.has_index, Some(false));
    assert_eq!(
        bam_inspection.index_path.as_deref(),
        Some(dir.join("sample.bai").as_path())
    );
}

#[test]
fn phased_vcf_reports_phasing() {
    let dir = temp_dir("phased-vcf");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0|1\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.phased, Some(true));
}

#[test]
fn unphased_vcf_reports_false_for_slash_genotypes() {
    let dir = temp_dir("unphased-vcf");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.phased, Some(false));
}

#[test]
fn dot_slash_genotype_is_treated_as_unphased() {
    let dir = temp_dir("dot-slash-vcf");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t100\trs1\tA\tG\t.\tPASS\t.\tGT:GQ\t./.:9\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.phased, Some(false));
}

#[test]
fn vcf_without_usable_gt_reports_unknown_phasing() {
    let dir = temp_dir("unknown-phasing-vcf");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t100\trs1\tA\tG\t.\tPASS\t.\tDP\t42\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.phased, None);
}

#[test]
fn mixed_vcf_sampled_rows_prefer_true_when_pipe_is_seen() {
    let dir = temp_dir("mixed-phasing-vcf");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
1\t101\trs2\tA\tG\t.\tPASS\t.\tGT\t0|1\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.phased, Some(true));
}

#[test]
fn vcf_b37_contig_headers_report_grch37_assembly() {
    let dir = temp_dir("b37-vcf");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=249250621,assembly=b37>\n\
##SentieonCommandLine.Haplotyper=<ID=Haplotyper,CommandLine=\"--reference human_g1k_v37_decoy.fasta\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.assembly, Some(Assembly::Grch37));
}

#[test]
fn vcf_grch38_reference_token_reports_grch38_assembly() {
    let dir = temp_dir("grch38-vcf-reference");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##reference=GCA_000001405.15_GRCh38_no_alt_analysis_set.fa\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.assembly, Some(Assembly::Grch38));
}

#[test]
fn vcf_grch38_contig_lengths_report_grch38_assembly() {
    let dir = temp_dir("grch38-vcf-contigs");
    let path = dir.join("sample.vcf");
    std::fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=248956422>\n\
##contig=<ID=chr22,length=50818468>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
chr1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let inspection = inspect_file(&path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.detected_kind, DetectedKind::Vcf);
    assert_eq!(inspection.assembly, Some(Assembly::Grch38));
}
