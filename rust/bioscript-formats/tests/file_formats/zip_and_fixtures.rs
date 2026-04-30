use super::*;

#[test]
fn batch_lookup_preserves_input_order_after_coordinate_sorting() {
    let dir = temp_dir("batch-order");
    let path = dir.join("sample.txt");
    fs::write(
        &path,
        "rsid\tchromosome\tposition\tgenotype\n\
         rs2\t1\t20\tCT\n\
         rs1\t1\t10\tAG\n",
    )
    .unwrap();
    let store = GenotypeStore::from_file(&path).unwrap();

    let results = store
        .lookup_variants(&[
            VariantSpec {
                grch38: Some(bioscript_core::GenomicLocus {
                    chrom: "1".to_owned(),
                    start: 20,
                    end: 20,
                }),
                ..VariantSpec::default()
            },
            VariantSpec {
                grch38: Some(bioscript_core::GenomicLocus {
                    chrom: "1".to_owned(),
                    start: 10,
                    end: 10,
                }),
                ..VariantSpec::default()
            },
        ])
        .unwrap();

    assert_eq!(results[0].genotype.as_deref(), Some("CT"));
    assert_eq!(
        results[0].evidence,
        vec!["resolved by locus 1:20".to_owned()]
    );
    assert_eq!(results[1].genotype.as_deref(), Some("AG"));
    assert_eq!(
        results[1].evidence,
        vec!["resolved by locus 1:10".to_owned()]
    );
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
fn zip_vcf_gz_entry_is_selected_and_read_as_vcf() {
    let dir = temp_dir("zip-vcf-gz-entry");
    let zip_path = dir.join("sample.zip");

    let file = fs::File::create(&zip_path).unwrap();
    let mut writer = zip::ZipWriter::new(file);
    writer
        .add_directory("nested/", SimpleFileOptions::default())
        .unwrap();
    writer
        .start_file("nested/sample.vcf.gz", SimpleFileOptions::default())
        .unwrap();
    writer
        .write_all(
            b"##fileformat=VCFv4.2\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
2\t22\trsZipVcfGz\tG\tA\t.\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
    writer.finish().unwrap();

    let store = GenotypeStore::from_file(&zip_path).unwrap();
    assert_eq!(store.backend_name(), "vcf");
    assert_eq!(store.get("rsZipVcfGz").unwrap().as_deref(), Some("GA"));
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
