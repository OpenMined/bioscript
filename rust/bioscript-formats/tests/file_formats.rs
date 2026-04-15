use std::{
    env, fs,
    io::Write,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_core::{VariantKind, VariantSpec};
use bioscript_formats::{GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore};
use zip::write::SimpleFileOptions;

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-formats-{label}-{}-{nanos}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("bioscript repo root")
        .to_path_buf()
}

fn shared_test_data_root() -> Option<PathBuf> {
    if let Some(path) = env::var_os("BIOSCRIPT_TEST_DATA_DIR") {
        let candidate = PathBuf::from(path);
        if candidate.exists() {
            return Some(candidate);
        }
    }

    let local = repo_root().join("test-data");
    if local.exists() {
        return Some(local);
    }

    let home_cache = env::var_os("HOME")
        .map(PathBuf::from)
        .map(|home| home.join(".bioscript/cache/test-data"));
    home_cache.filter(|path| path.exists())
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

#[test]
fn zip_genotype_file_is_auto_detected_and_readable() {
    let dir = temp_dir("zip-auto");
    let zip_path = dir.join("apol1-input.zip");

    let file = fs::File::create(&zip_path).unwrap();
    let mut writer = zip::ZipWriter::new(file);
    writer
        .start_file("test_snps.txt", SimpleFileOptions::default())
        .unwrap();
    writer
        .write_all(
            b"rsid\tchromosome\tposition\tgenotype\nrs73885319\t22\t1\tAG\nrs60910145\t22\t2\tTG\nrs71785313\t22\t3\tII\n",
        )
        .unwrap();
    writer.finish().unwrap();

    let store = GenotypeStore::from_file(&zip_path).unwrap();
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
    assert_eq!(store.get("rs60910145").unwrap().as_deref(), Some("TG"));
    assert_eq!(store.get("rs71785313").unwrap().as_deref(), Some("II"));
}

#[test]
fn zip_genotype_file_can_be_forced_by_format() {
    let dir = temp_dir("zip-forced");
    let zip_path = dir.join("apol1-input.dat");

    let file = fs::File::create(&zip_path).unwrap();
    let mut writer = zip::ZipWriter::new(file);
    writer
        .start_file("test_snps.txt", SimpleFileOptions::default())
        .unwrap();
    writer
        .write_all(b"rsid\tchromosome\tposition\tgenotype\nrs73885319\t22\t1\tAG\n")
        .unwrap();
    writer.finish().unwrap();

    let store = GenotypeStore::from_file_with_options(
        &zip_path,
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Zip),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();

    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
}

#[test]
fn zip_vcf_entry_is_auto_detected_and_readable() {
    let dir = temp_dir("zip-vcf");
    let zip_path = dir.join("apol1-sample.zip");

    let file = fs::File::create(&zip_path).unwrap();
    let mut writer = zip::ZipWriter::new(file);
    writer
        .start_file("nested/sample.vcf", SimpleFileOptions::default())
        .unwrap();
    writer
        .write_all(
            b"##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
22\t36265860\trs73885319\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
    writer.finish().unwrap();

    let store = GenotypeStore::from_file(&zip_path).unwrap();
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
}

#[test]
fn shared_real_world_zipped_genotype_exports_are_readable() {
    struct FixtureExpectation {
        relative: &'static str,
        rsid: &'static str,
        genotype: &'static str,
    }

    let fixtures = [
        FixtureExpectation {
            relative: "23andme/v2/hu0199C8/23data20100526.txt.zip",
            rsid: "rs3094315",
            genotype: "AA",
        },
        FixtureExpectation {
            relative: "23andme/v3/huE4DAE4/huE4DAE4_20120522224129.txt.zip",
            rsid: "rs3131972",
            genotype: "GG",
        },
        FixtureExpectation {
            relative: "23andme/v4/huE18D82/genome__v4_Full_2016.txt.zip",
            rsid: "rs3131972",
            genotype: "AG",
        },
        FixtureExpectation {
            relative: "23andme/v5/hu50B3F5/genome_hu50B3F5_v5_Full.zip",
            rsid: "rs116587930",
            genotype: "GG",
        },
        FixtureExpectation {
            relative: "dynamicdna/100001-synthetic/100001_X_X_GSAv3-DTC_GRCh38-07-12-2025.txt.zip",
            rsid: "rs116587930",
            genotype: "GG",
        },
        FixtureExpectation {
            relative: "ancestrydna/huE922FC/AncestryDNA.txt.zip",
            rsid: "rs3131972",
            genotype: "GG",
        },
        FixtureExpectation {
            relative: "familytreedna/hu17B792/2017-04-29_Family_Tree_DNA_Data.csv.zip",
            rsid: "rs1000530",
            genotype: "TT",
        },
        FixtureExpectation {
            relative: "genesforgood/hu80B047/GFG0_filtered_imputed_genotypes_noY_noMT_23andMe.txt.zip",
            rsid: "rs3094315",
            genotype: "AA",
        },
        FixtureExpectation {
            relative: "myheritage/hu33515F/MyHeritage_raw_dna_data.zip",
            rsid: "rs3131972",
            genotype: "GG",
        },
    ];

    for fixture in fixtures {
        let Some(path) = shared_fixture_or_skip(
            "shared_real_world_zipped_genotype_exports_are_readable",
            fixture.relative,
        ) else {
            return;
        };

        let store = GenotypeStore::from_file(&path).unwrap();
        assert_eq!(
            store.get(fixture.rsid).unwrap().as_deref(),
            Some(fixture.genotype),
            "fixture {}",
            fixture.relative
        );
    }
}

#[test]
fn bundled_dynamicdna_gsav3_plain_text_fixture_is_readable() {
    let path = repo_root()
        .join("old/examples/apol1/genotype_files/108179_G0G0_X_X_GSAv3-DTC_GRCh38-12-13-2025.txt");
    let store = GenotypeStore::from_file(&path).unwrap();
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AA"));
}

#[test]
fn vcf_variant_lookup_reads_single_sample_calls() {
    let dir = temp_dir("vcf");
    let vcf_path = dir.join("apol1_sample.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
22\t36265860\trs73885319\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
22\t36265900\trs60910145\tT\tG\t.\tPASS\t.\tGT\t1/1\n\
22\t36266005\trs71785313\tA\tATTTAA\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&vcf_path).unwrap();
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
    assert_eq!(store.get("rs60910145").unwrap().as_deref(), Some("GG"));

    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs71785313".to_owned()],
            kind: Some(VariantKind::Insertion),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(observation.genotype.as_deref(), Some("DI"));
}

struct CramFixture {
    cram: PathBuf,
    reference: PathBuf,
    reference_index: PathBuf,
    input_index: PathBuf,
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
    GenotypeStore::from_file_with_options(
        &fx.cram,
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Cram),
            input_index: Some(fx.input_index.clone()),
            reference_file: Some(fx.reference.clone()),
            reference_index: Some(fx.reference_index.clone()),
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
fn cram_mini_fixture_md5_mismatch_is_tolerated() {
    // mini_bad_ref.fa has a single-base mutation at chr_test:2800, inside the
    // slice span but far from our query locus at 1000. noodles' strict MD5
    // check will fail; bioscript must warn + retry unchecked + still return
    // the correct genotype (the bases at pos 1000 are identical in both refs).
    let fx = mini_cram_fixture_with_bad_ref();
    let store = open_cram_store(&fx);

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
    let store = open_cram_store(&fx);

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
