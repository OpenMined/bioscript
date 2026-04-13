use std::{
    fs,
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
