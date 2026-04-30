use super::*;

#[test]
fn genotype_store_from_bytes_handles_genotype_text() {
    let store = GenotypeStore::from_bytes(
        "sample.txt",
        b"\xef\xbb\xbfrsid\tchromosome\tposition\tgenotype\n\
          # skipped comment\n\
          rs73885319\t22\t36265860\tag\n\
          rs60910145\t22\t36265900\tN/A\n",
    )
    .unwrap();

    assert_eq!(store.backend_name(), "text");
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
    assert_eq!(store.get("rs60910145").unwrap().as_deref(), Some("--"));
}

#[test]
fn genotype_source_format_parses_supported_values_and_rejects_unknowns() {
    assert_eq!(
        "txt".parse::<GenotypeSourceFormat>().unwrap(),
        GenotypeSourceFormat::Text
    );
    assert_eq!(
        "GENOTYPE".parse::<GenotypeSourceFormat>().unwrap(),
        GenotypeSourceFormat::Text
    );
    assert_eq!(
        "zip".parse::<GenotypeSourceFormat>().unwrap(),
        GenotypeSourceFormat::Zip
    );
    assert_eq!(
        "vcf".parse::<GenotypeSourceFormat>().unwrap(),
        GenotypeSourceFormat::Vcf
    );
    assert_eq!(
        "cram".parse::<GenotypeSourceFormat>().unwrap(),
        GenotypeSourceFormat::Cram
    );

    let err = "bam".parse::<GenotypeSourceFormat>().unwrap_err();
    assert_eq!(err, "unsupported input format: bam");
}

#[test]
fn backend_capabilities_match_query_backend_type() {
    let rsid_map = GenotypeStore::from_bytes(
        "sample.txt",
        b"rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\n",
    )
    .unwrap();
    assert_eq!(rsid_map.backend_name(), "text");
    assert!(rsid_map.supports(QueryKind::GenotypeByRsid));
    assert!(!rsid_map.supports(QueryKind::GenotypeByLocus));

    let dir = temp_dir("backend-capabilities");
    let text_path = dir.join("sample.txt");
    fs::write(
        &text_path,
        "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\n",
    )
    .unwrap();
    let delimited = GenotypeStore::from_file(&text_path).unwrap();
    assert_eq!(delimited.backend_name(), "text");
    assert!(delimited.supports(QueryKind::GenotypeByRsid));
    assert!(delimited.supports(QueryKind::GenotypeByLocus));

    let vcf_path = dir.join("sample.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();
    let vcf = GenotypeStore::from_file(&vcf_path).unwrap();
    assert_eq!(vcf.backend_name(), "vcf");
    assert!(vcf.supports(QueryKind::GenotypeByRsid));
    assert!(vcf.supports(QueryKind::GenotypeByLocus));

    let cram = GenotypeStore::from_file_with_options(
        &dir.join("sample.dat"),
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Cram),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();
    assert_eq!(cram.backend_name(), "cram");
    assert!(!cram.supports(QueryKind::GenotypeByRsid));
    assert!(cram.supports(QueryKind::GenotypeByLocus));
}

#[test]
fn genotype_store_from_bytes_handles_vcf() {
    let store = GenotypeStore::from_bytes(
        "sample.vcf",
        b"##fileformat=VCFv4.2\n\
          ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
          22\t36265860\trs73885319\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    assert_eq!(store.backend_name(), "vcf");
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
}

#[test]
fn vcf_bytes_skip_unusable_rows_and_decode_no_call_forms() {
    let store = GenotypeStore::from_bytes(
        "sample.vcf",
        b"##fileformat=VCFv4.2\n\
          #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
          1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
          1\t11\trsEmptyRef\t.\tG\t.\tPASS\t.\tGT\t0/1\n\
          1\t12\trsEmptyAlt\tA\t.\t.\tPASS\t.\tGT\t0/1\n\
          1\t13\trsShort\tA\tG\n\
          1\t14\trsNoCall\tA\tG\t.\tPASS\t.\tGT\t.\n\
          1\t15\trsPartialNoCall\tA\tG\t.\tPASS\t.\tGT\t./1\n\
          1\t16\trsOutOfRange\tA\tG\t.\tPASS\t.\tGT\t0/2\n\
          1\t17\trsValid\tC\tT\t.\tPASS\t.\tGT\t1|1\n",
    )
    .unwrap();

    assert_eq!(store.backend_name(), "vcf");
    assert_eq!(store.get("rsValid").unwrap().as_deref(), Some("TT"));
    assert_eq!(store.get("rsNoCall").unwrap().as_deref(), Some("--"));
    assert_eq!(store.get("rsPartialNoCall").unwrap().as_deref(), Some("--"));
    assert_eq!(store.get("rsOutOfRange").unwrap(), None);
    assert_eq!(store.get("rsEmptyRef").unwrap().as_deref(), Some(".G"));
    assert_eq!(store.get("rsEmptyAlt").unwrap(), None);
}

#[test]
fn extensionless_vcf_is_detected_by_content_and_can_be_forced() {
    let dir = temp_dir("extensionless-vcf");
    let path = dir.join("sample.data");
    fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
         ##reference=GRCh37\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         1\t10\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let detected = GenotypeStore::from_file(&path).unwrap();
    assert_eq!(detected.backend_name(), "vcf");
    assert_eq!(detected.get("rs1").unwrap().as_deref(), Some("AG"));

    let forced = GenotypeStore::from_file_with_options(
        &path,
        &GenotypeLoadOptions {
            format: Some(GenotypeSourceFormat::Vcf),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();
    assert_eq!(forced.backend_name(), "vcf");
    assert_eq!(forced.get("rs1").unwrap().as_deref(), Some("AG"));
}

#[test]
fn vcf_file_lookup_handles_gt_field_order_no_calls_and_bad_positions() {
    let dir = temp_dir("vcf-field-order");
    let path = dir.join("sample.vcf");
    fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
         ##reference=GRCh38\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         1\t10\trs1\tA\tG\t.\tPASS\t.\tDP:GT\t14:0|1\n\
         1\t11\trs2\tC\tT\t.\tPASS\t.\tGT:DP\t./.:9\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();
    assert_eq!(store.get("rs1").unwrap().as_deref(), Some("AG"));
    assert_eq!(store.get("rs2").unwrap().as_deref(), Some("--"));

    let bad_path = dir.join("bad.vcf");
    fs::write(
        &bad_path,
        "##fileformat=VCFv4.2\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         1\tnot-a-pos\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();
    let err = GenotypeStore::from_file(&bad_path)
        .unwrap()
        .get("rs1")
        .unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to parse VCF position 'not-a-pos'"),
        "{err:?}"
    );
}

#[test]
fn genotype_store_from_bytes_handles_zip() {
    let bytes = zip_bytes(
        "nested/sample.txt",
        b"rsid\tchromosome\tposition\tgenotype\nrs73885319\t22\t36265860\tAG\n",
    );

    let store = GenotypeStore::from_bytes("sample.zip", &bytes).unwrap();

    assert_eq!(store.backend_name(), "zip");
    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
}

#[test]
fn rsid_map_batch_lookup_preserves_order_and_reports_missing_rsids() {
    let store = GenotypeStore::from_bytes(
        "sample.txt",
        b"rsid\tchromosome\tposition\tgenotype\nrs2\t1\t20\tCT\nrs1\t1\t10\tAG\n",
    )
    .unwrap();

    let results = store
        .lookup_variants(&[
            VariantSpec {
                rsids: vec!["rs2".to_owned()],
                ..VariantSpec::default()
            },
            VariantSpec {
                rsids: vec!["rsMissing".to_owned()],
                ..VariantSpec::default()
            },
            VariantSpec {
                rsids: vec!["rs1".to_owned()],
                ..VariantSpec::default()
            },
        ])
        .unwrap();

    assert_eq!(results[0].genotype.as_deref(), Some("CT"));
    assert_eq!(results[1].genotype, None);
    assert_eq!(
        results[1].evidence,
        vec!["no matching rsid found".to_owned()]
    );
    assert_eq!(results[2].genotype.as_deref(), Some("AG"));
}

#[test]
fn genotype_store_from_bytes_rejects_malformed_zip() {
    let err = GenotypeStore::from_bytes("sample.zip", b"not a zip").unwrap_err();

    assert!(
        format!("{err:?}").contains("failed to read genotype zip sample.zip"),
        "{err:?}"
    );
}

#[test]
fn genotype_store_from_bytes_rejects_zip_without_supported_entry() {
    let bytes = zip_bytes("notes.bin", b"not genotype data");

    let err = GenotypeStore::from_bytes("sample.zip", &bytes).unwrap_err();

    assert!(
        format!("{err:?}")
            .contains("zip archive sample.zip does not contain a supported genotype file"),
        "{err:?}"
    );
}
