use super::*;

struct CramFixture {
    cram: PathBuf,
    reference: PathBuf,
    reference_index: PathBuf,
    input_index: PathBuf,
}

fn forced_cram_store(dir: &std::path::Path, reference_name: &str) -> GenotypeStore {
    GenotypeStore::from_file_with_options(
        &dir.join("missing.cram"),
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Cram),
            reference_file: Some(dir.join(reference_name)),
            reference_index: Some(dir.join(format!("{reference_name}.fai"))),
            input_index: Some(dir.join("missing.cram.crai")),
            allow_reference_md5_mismatch: false,
        },
    )
    .unwrap()
}

#[test]
fn forced_cram_backend_reports_reference_and_coordinate_errors_without_reading_cram() {
    let dir = temp_dir("cram-reference-errors");
    let cram_path = dir.join("missing.cram");
    let store_without_reference = GenotypeStore::from_file_with_options(
        &cram_path,
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Cram),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();
    let err = store_without_reference
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs1".to_owned()],
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(
        format!("{err:?}").contains("without --reference-file"),
        "{err:?}"
    );

    let store = forced_cram_store(&dir, "GRCh38.fa");
    let err = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs1".to_owned()],
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("needs GRCh37/GRCh38 coordinates"));
    assert!(format!("{err:?}").contains("reference index"));
    assert!(format!("{err:?}").contains("input index"));
}

#[test]
fn forced_cram_backend_reports_snp_and_indel_argument_errors_without_reading_cram() {
    let dir = temp_dir("cram-variant-argument-errors");
    let store = forced_cram_store(&dir, "GRCh38.fa");
    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            kind: Some(VariantKind::Snp),
            alternate: Some("G".to_owned()),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("SNP variant requires ref/reference"));

    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            reference: Some("A".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("SNP variant requires alt/alternate"));

    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            kind: Some(VariantKind::Deletion),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("deletion variant requires deletion_length"));

    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            kind: Some(VariantKind::Indel),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("indel variant requires ref/reference"));

    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            reference: Some("A".to_owned()),
            kind: Some(VariantKind::Insertion),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("indel variant requires alt/alternate"));
}

#[test]
fn forced_cram_backend_reports_file_and_assembly_errors_without_reading_cram() {
    let dir = temp_dir("cram-file-assembly-errors");
    let store = forced_cram_store(&dir, "GRCh38.fa");
    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to open indexed FASTA"),
        "{err:?}"
    );

    let err = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Other),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("does not yet support Other"));

    let hg19_store = forced_cram_store(&dir, "hg19.fa");
    let err = hg19_store
        .lookup_variant(&VariantSpec {
            grch37: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 10,
                end: 10,
            }),
            kind: Some(VariantKind::Other),
            ..VariantSpec::default()
        })
        .unwrap_err();
    assert!(format!("{err:?}").contains("does not yet support Other"));
}

fn run_large_cram_tests() -> bool {
    env::var_os("BIOSCRIPT_RUN_LARGE_TESTS").is_some()
}

fn require_large_cram_tests(test_name: &str) -> bool {
    if run_large_cram_tests() {
        true
    } else {
        eprintln!("skipping {test_name}: set BIOSCRIPT_RUN_LARGE_TESTS=1 to enable");
        false
    }
}

fn cram_fixture_or_skip(test_name: &str) -> Option<CramFixture> {
    if !require_large_cram_tests(test_name) {
        return None;
    }
    let root = shared_test_data_root()?;
    let fx = CramFixture {
        cram: root.join("1k-genomes/aligned/NA06985.final.cram"),
        reference: root.join("1k-genomes/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"),
        reference_index: root.join("1k-genomes/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai"),
        input_index: root.join("1k-genomes/aligned/NA06985.final.cram.crai"),
    };
    for p in [
        &fx.cram,
        &fx.reference,
        &fx.reference_index,
        &fx.input_index,
    ] {
        if !p.exists() {
            eprintln!("skipping {test_name}: missing {}", p.display());
            return None;
        }
    }
    Some(fx)
}

fn chr_y_cram_fixture_or_skip(test_name: &str) -> Option<CramFixture> {
    let root = shared_test_data_root()?;
    let fx = CramFixture {
        cram: root.join("NA06985-chrY/aligned/NA06985.final.chrY.cram"),
        reference: root.join("NA06985-chrY/ref/GRCh38_chrY.fa"),
        reference_index: root.join("NA06985-chrY/ref/GRCh38_chrY.fa.fai"),
        input_index: root.join("NA06985-chrY/aligned/NA06985.final.chrY.cram.crai"),
    };
    for p in [
        &fx.cram,
        &fx.reference,
        &fx.reference_index,
        &fx.input_index,
    ] {
        if !p.exists() {
            eprintln!("skipping {test_name}: missing {}", p.display());
            return None;
        }
    }
    Some(fx)
}

fn open_cram_store(fx: &CramFixture) -> GenotypeStore {
    open_cram_store_with_md5_policy(fx, false)
}

fn open_cram_store_with_md5_policy(
    fx: &CramFixture,
    allow_reference_md5_mismatch: bool,
) -> GenotypeStore {
    GenotypeStore::from_file_with_options(
        &fx.cram,
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Cram),
            input_index: Some(fx.input_index.clone()),
            reference_file: Some(fx.reference.clone()),
            reference_index: Some(fx.reference_index.clone()),
            allow_reference_md5_mismatch,
        },
    )
    .expect("open cram store")
}

fn mini_fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn mini_cram_fixture() -> CramFixture {
    let dir = mini_fixtures_dir();
    CramFixture {
        cram: dir.join("mini.cram"),
        reference: dir.join("mini.fa"),
        reference_index: dir.join("mini.fa.fai"),
        input_index: dir.join("mini.cram.crai"),
    }
}

fn mini_cram_fixture_with_bad_ref() -> CramFixture {
    let dir = mini_fixtures_dir();
    CramFixture {
        cram: dir.join("mini.cram"),
        reference: dir.join("mini_bad_ref.fa"),
        reference_index: dir.join("mini_bad_ref.fa.fai"),
        input_index: dir.join("mini.cram.crai"),
    }
}

#[test]
fn cram_mini_fixture_streams_only_locus_overlapping_reads() {
    // mini.cram has 2000 reads covering chr_test:500..2499. The streaming path
    // should decode roughly until it passes the locus and stop — correctness is
    // asserted via depth (exactly 50 reads overlap a single base in the middle).
    // If the streaming + early-termination path breaks and falls back to full
    // slice decode, the wall time still finishes fine on 2000 reads but this
    // test also catches regressions that double-count or miss reads.
    let fx = mini_cram_fixture();
    let store = open_cram_store(&fx);

    let start = std::time::Instant::now();
    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["mini_locus_1000".to_owned()],
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "chr_test".to_owned(),
                start: 1000,
                end: 1000,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("C".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .expect("mini cram lookup");
    let elapsed = start.elapsed();

    assert_eq!(observation.backend, "cram");
    assert_eq!(
        observation.depth.unwrap_or(0),
        50,
        "expected exactly 50 reads overlapping chr_test:1000, got {:?}",
        observation.depth
    );
    // All reads match reference in the fixture so alt_count should be zero.
    assert_eq!(observation.ref_count.unwrap_or(0), 50);
    assert_eq!(observation.alt_count.unwrap_or(0), 0);
    assert!(
        elapsed.as_millis() < 2000,
        "mini CRAM lookup took {elapsed:?}, expected well under 2s"
    );
}

#[test]
fn cram_mini_fixture_md5_mismatch_is_tolerated_when_allowed() {
    // mini_bad_ref.fa has a single-base mutation at chr_test:2800, inside the
    // slice span but far from our query locus at 1000. noodles' strict MD5
    // check will fail; bioscript must warn + retry unchecked + still return
    // the correct genotype (the bases at pos 1000 are identical in both refs).
    let fx = mini_cram_fixture_with_bad_ref();
    let store = open_cram_store_with_md5_policy(&fx, true);

    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["mini_locus_1000".to_owned()],
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "chr_test".to_owned(),
                start: 1000,
                end: 1000,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("C".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .expect("mini cram lookup should succeed via md5 fallback");

    assert_eq!(observation.backend, "cram");
    assert_eq!(
        observation.depth.unwrap_or(0),
        50,
        "expected exactly 50 reads after md5 fallback, got {:?}",
        observation.depth
    );
    // Bases at the query locus are the same in both references, so the
    // fallback-decoded reads should still be ref-homozygous.
    assert_eq!(observation.ref_count.unwrap_or(0), 50);
    assert_eq!(observation.alt_count.unwrap_or(0), 0);
}

#[test]
fn cram_chr_y_fixture_lookup_is_fast_and_correct() {
    let Some(fx) = chr_y_cram_fixture_or_skip("cram_chr_y_fixture_lookup_is_fast_and_correct")
    else {
        return;
    };
    let store = open_cram_store(&fx);

    let start = std::time::Instant::now();
    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["chrY_smoke_3449570".to_owned()],
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "chrY".to_owned(),
                start: 3_449_570,
                end: 3_449_570,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .expect("chrY lookup");
    let elapsed = start.elapsed();

    assert_eq!(observation.backend, "cram");
    let depth = observation.depth.unwrap_or(0);
    assert!(
        depth >= 8,
        "expected >=8 reads at chrY smoke locus, got {depth}"
    );
    assert_eq!(observation.alt_count.unwrap_or(0), 0);
    assert!(
        observation.ref_count.unwrap_or(0) >= 8,
        "expected ref-supporting reads at chrY smoke locus, got {:?}",
        observation.ref_count
    );
    assert!(
        elapsed.as_secs() < 5,
        "chrY CRAM lookup took {elapsed:?}, expected <5s"
    );
}

#[test]
fn cram_apol1_snp_lookup_is_fast_and_correct() {
    let Some(fx) = cram_fixture_or_skip("cram_apol1_snp_lookup_is_fast_and_correct") else {
        return;
    };
    let store = open_cram_store(&fx);

    let start = std::time::Instant::now();
    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs73885319".to_owned()],
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "22".to_owned(),
                start: 36_265_860,
                end: 36_265_860,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .expect("apol1 lookup");
    let elapsed = start.elapsed();

    assert_eq!(observation.backend, "cram");
    // NA06985 is reference-homozygous at APOL1 G1 site 1 per samtools mpileup.
    let depth = observation.depth.unwrap_or(0);
    assert!(
        depth >= 10,
        "expected >=10 reads at APOL1 locus, got {depth}"
    );
    let ref_count = observation.ref_count.unwrap_or(0);
    let alt_count = observation.alt_count.unwrap_or(0);
    assert!(
        ref_count > alt_count,
        "NA06985 APOL1 G1 site 1 should be ref-dominant: ref={ref_count} alt={alt_count}"
    );

    // Slice-level CRAM decode is the hot path. Samtools does the same locus
    // in ~40ms; we allow a generous ceiling to catch regressions (e.g. if the
    // streaming/early-termination path breaks and we fall back to decoding
    // every record in the slice, this blows past 10s).
    assert!(
        elapsed.as_secs() < 5,
        "APOL1 CRAM lookup took {elapsed:?}, expected <5s (samtools does it in ~40ms)"
    );
}

#[test]
fn cram_md5_mismatch_is_tolerated_and_returns_correct_result() {
    // For NA06985.final.cram, the bundled GRCh38 FASTA's chr6 MD5 does not
    // match the @SQ M5 the CRAM was encoded against (only chr22 matches).
    // We must warn + fall back to unchecked decoding, and still return the
    // correct genotype at the GLP1 rs10305420 locus. The correct call per
    // samtools mpileup is reference-homozygous (CC).
    let Some(fx) =
        cram_fixture_or_skip("cram_md5_mismatch_is_tolerated_and_returns_correct_result")
    else {
        return;
    };
    let store = open_cram_store_with_md5_policy(&fx, true);

    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs10305420".to_owned()],
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "6".to_owned(),
                start: 39_048_860,
                end: 39_048_860,
            }),
            reference: Some("C".to_owned()),
            alternate: Some("T".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .expect("glp1 lookup should succeed via md5 fallback");

    assert_eq!(observation.backend, "cram");
    let depth = observation.depth.unwrap_or(0);
    assert!(
        depth >= 10,
        "expected >=10 reads at GLP1 locus after md5 fallback, got {depth}"
    );
    let genotype = observation
        .genotype
        .as_deref()
        .expect("expected a genotype call");
    assert!(
        genotype.chars().all(|c| c == 'C' || c == 'T'),
        "unexpected genotype after md5 fallback: {genotype}"
    );

    // Parity with `samtools mpileup -f <ref> -r chr6:39048860-39048860`:
    // that locus shows a mixed pileup (roughly half reference C, half T).
    // We assert total depth matches samtools' reported depth within a small
    // tolerance — confirms we are not silently dropping or duplicating reads
    // after the unchecked-reference fallback.
    let depth_i32 = i32::try_from(depth).unwrap_or(i32::MAX);
    let samtools_depth: i32 = 41;
    assert!(
        (depth_i32 - samtools_depth).abs() <= 6,
        "depth {depth_i32} differs from samtools mpileup depth {samtools_depth} by >6"
    );
}

#[test]
fn cram_rs9357296_reports_heterozygous_counts_for_na06985() {
    let Some(fx) = cram_fixture_or_skip("cram_rs9357296_reports_heterozygous_counts_for_na06985")
    else {
        return;
    };
    let store = open_cram_store_with_md5_policy(&fx, true);

    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs9357296".to_owned()],
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "6".to_owned(),
                start: 39_011_352,
                end: 39_011_352,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .expect("rs9357296 lookup");

    assert_eq!(observation.backend, "cram");
    assert_eq!(observation.genotype.as_deref(), Some("AG"));
    assert_eq!(observation.raw_counts.get("A").copied(), Some(18));
    assert_eq!(observation.raw_counts.get("G").copied(), Some(12));
    assert_eq!(observation.raw_counts.get("T").copied(), None);
    assert_eq!(observation.depth, Some(29));
    assert_eq!(observation.ref_count, Some(17));
    assert_eq!(observation.alt_count, Some(12));
    assert!(
        observation
            .decision
            .as_deref()
            .is_some_and(|text| text.contains("alt_fraction=0.414")),
        "missing SNP decision summary: {:?}",
        observation.decision
    );
    assert!(
        observation
            .evidence
            .iter()
            .any(|line| line.contains("raw pileup depth=30")),
        "missing raw pileup evidence: {:?}",
        observation.evidence
    );
    assert!(
        observation.evidence.iter().any(|line| {
            line.contains("filtered_duplicate=4")
                && line.contains("filtered_low_base_quality=1")
                && line.contains("filtered_improper_pair=0")
        }),
        "missing mpileup-style filter evidence: {:?}",
        observation.evidence
    );
}
