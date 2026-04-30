use std::{
    env, fs,
    io::Write,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_core::{VariantKind, VariantSpec};
use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, QueryKind, alignment,
};
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

#[test]
fn alignment_index_parsers_handle_in_memory_bytes() {
    let fai = alignment::parse_fai_bytes(b"chr1\t4\t6\t4\t5\n").unwrap();
    let _repository = alignment::build_reference_repository_from_readers(
        std::io::BufReader::new(std::io::Cursor::new(b">chr1\nACGT\n".to_vec())),
        fai,
    );

    let err = alignment::parse_fai_bytes(b"not a fai").unwrap_err();
    assert!(format!("{err:?}").contains("failed to parse FASTA index bytes"));

    let err = alignment::parse_crai_bytes(b"not a crai").unwrap_err();
    assert!(format!("{err:?}").contains("failed to parse CRAM index bytes"));

    let err = alignment::parse_tbi_bytes(b"not a tbi").unwrap_err();
    assert!(format!("{err:?}").contains("failed to parse tabix index bytes"));
}

#[test]
fn alignment_reader_api_reports_invalid_cram_headers_without_real_fixtures() {
    let fai = alignment::parse_fai_bytes(b"chr1\t4\t6\t4\t5\n").unwrap();
    let repository = alignment::build_reference_repository_from_readers(
        std::io::BufReader::new(std::io::Cursor::new(b">chr1\nACGT\n".to_vec())),
        fai,
    );
    let locus = bioscript_core::GenomicLocus {
        chrom: "chr1".to_owned(),
        start: 1,
        end: 1,
    };
    let crai_bytes = fs::read(mini_fixtures_dir().join("mini.cram.crai")).unwrap();
    let mut reader = alignment::build_cram_indexed_reader_from_reader(
        std::io::Cursor::new(b"not a cram".to_vec()),
        alignment::parse_crai_bytes(&crai_bytes).unwrap(),
        repository,
    )
    .unwrap();

    let err =
        alignment::for_each_cram_record_with_reader(&mut reader, "bad.cram", &locus, |_| Ok(true))
            .unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to read CRAM header bad.cram"),
        "{err:?}"
    );

    let fai = alignment::parse_fai_bytes(b"chr1\t4\t6\t4\t5\n").unwrap();
    let repository = alignment::build_reference_repository_from_readers(
        std::io::BufReader::new(std::io::Cursor::new(b">chr1\nACGT\n".to_vec())),
        fai,
    );
    let mut raw_reader = alignment::build_cram_indexed_reader_from_reader(
        std::io::Cursor::new(b"still not a cram".to_vec()),
        alignment::parse_crai_bytes(&crai_bytes).unwrap(),
        repository,
    )
    .unwrap();

    let err = alignment::for_each_raw_cram_record_with_reader(
        &mut raw_reader,
        "raw-bad.cram",
        &locus,
        |_| Ok(true),
    )
    .unwrap_err();
    assert!(
        format!("{err:?}").contains("failed to read CRAM header raw-bad.cram"),
        "{err:?}"
    );
}

#[test]
fn delimited_parser_handles_comments_blank_lines_csv_and_split_alleles() {
    let dir = temp_dir("csv-split-alleles");
    let path = dir.join("sample.csv");
    fs::write(
        &path,
        "\n\
         # rsid,chromosome,position,allele1,allele2\n\
         // ignored comment\n\
         rs73885319,chr22,36265860,a,g\n\
         rs60910145,22,36265900,n/a,\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();

    assert_eq!(store.get("rs73885319").unwrap().as_deref(), Some("AG"));
    assert_eq!(store.get("rs60910145").unwrap().as_deref(), Some("--"));

    let observation = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "22".to_owned(),
                start: 36_265_860,
                end: 36_265_861,
            }),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(observation.genotype.as_deref(), Some("AG"));
    assert_eq!(
        observation.evidence,
        vec!["resolved by locus chr22:36265860".to_owned()]
    );
}

#[test]
fn delimited_parser_uses_comment_headers_aliases_quotes_and_extra_columns() {
    let dir = temp_dir("comment-header-aliases");
    let path = dir.join("sample.csv");
    fs::write(
        &path,
        "# SNP ID, Chrom, Base Pair Position, Result, Ignored\n\
         \"rsQuoted\", \"chr3\", \"300\", \"a t\", \"unused, value\"\n\
         rsSlash,3,301,A/-,\n\
         rsNone,3,302,None,\n\
         no_position,3,,AG,\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();

    assert_eq!(store.get("rsQuoted").unwrap().as_deref(), Some("AT"));
    assert_eq!(store.get("rsSlash").unwrap().as_deref(), Some("ID"));
    assert_eq!(store.get("rsNone").unwrap().as_deref(), Some("--"));
    assert_eq!(store.get("no_position").unwrap().as_deref(), Some("AG"));
    let observation = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "3".to_owned(),
                start: 300,
                end: 300,
            }),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(observation.genotype.as_deref(), Some("AT"));
}

#[test]
fn delimited_parser_handles_space_delimited_rows_without_headers_and_inline_comments() {
    let dir = temp_dir("space-default-header");
    let path = dir.join("sample.txt");
    fs::write(
        &path,
        "\n\
         rsSpace chr2 200 tc # inline comment\n\
         chrOnly chr2 201 aa\n\
         badrow\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();

    assert_eq!(store.get("rsSpace").unwrap().as_deref(), Some("TC"));
    let observation = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "2".to_owned(),
                start: 201,
                end: 201,
            }),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(observation.genotype.as_deref(), Some("AA"));
    assert_eq!(
        observation.evidence,
        vec!["resolved by locus chr2:201".to_owned()]
    );
}

#[test]
fn vcf_coordinate_lookup_normalizes_chr_prefix_and_handles_multiallelic_gt() {
    let dir = temp_dir("vcf-chr-normalize");
    let path = dir.join("sample.vcf");
    fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
         ##reference=GRCh38\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         chr1\t1000\t.\tA\tC,G\t.\tPASS\t.\tGT\t2/1\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();
    let observation = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 1000,
                end: 1001,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .unwrap();

    assert_eq!(observation.genotype.as_deref(), Some("GC"));
    assert_eq!(observation.assembly, Some(bioscript_core::Assembly::Grch38));
    assert_eq!(
        observation.evidence,
        vec!["resolved by locus chr1:1000".to_owned()]
    );
}

#[test]
fn vcf_locus_lookup_handles_deletion_insertion_and_unresolved_evidence() {
    let dir = temp_dir("vcf-indel-locus");
    let path = dir.join("sample.hg19.vcf");
    fs::write(
        &path,
        "##fileformat=VCFv4.2\n\
         ##reference=hg19\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
         1\t99\t.\tAT\tA\t.\tPASS\t.\tGT\t0/1\n\
         chr1\t199\t.\tA\tATG\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();
    let deletion = store
        .lookup_variant(&VariantSpec {
            grch37: Some(bioscript_core::GenomicLocus {
                chrom: "chr1".to_owned(),
                start: 100,
                end: 100,
            }),
            reference: Some("AT".to_owned()),
            alternate: Some("A".to_owned()),
            kind: Some(VariantKind::Deletion),
            deletion_length: Some(1),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(deletion.genotype.as_deref(), Some("ID"));
    assert_eq!(deletion.assembly, Some(bioscript_core::Assembly::Grch37));
    assert_eq!(deletion.evidence, vec!["resolved by locus 1:99".to_owned()]);

    let insertion = store
        .lookup_variant(&VariantSpec {
            grch37: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 200,
                end: 200,
            }),
            reference: Some("A".to_owned()),
            alternate: Some("ATG".to_owned()),
            kind: Some(VariantKind::Insertion),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(insertion.genotype.as_deref(), Some("DI"));
    assert_eq!(
        insertion.evidence,
        vec!["resolved by locus chr1:199".to_owned()]
    );

    let unresolved = store
        .lookup_variant(&VariantSpec {
            grch37: Some(bioscript_core::GenomicLocus {
                chrom: "1".to_owned(),
                start: 300,
                end: 300,
            }),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(unresolved.genotype, None);
    assert_eq!(
        unresolved.evidence,
        vec!["no matching rsid or locus found for variant_by_locus".to_owned()]
    );
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

#[test]
fn vcf_variant_lookup_ignores_symbolic_non_ref_alt_when_decoding_gt() {
    let dir = temp_dir("vcf-non-ref");
    let vcf_path = dir.join("sample.g.vcf");
    fs::write(
        &vcf_path,
        "##fileformat=VCFv4.2\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
6\t39016636\trs10305420\tC\tT,<NON_REF>\t.\tPASS\t.\tGT\t0/1\n\
6\t38979128\trs9357296\tA\tG,<NON_REF>\t.\tPASS\t.\tGT\t0/1\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file(&vcf_path).unwrap();
    assert_eq!(store.get("rs10305420").unwrap().as_deref(), Some("CT"));
    assert_eq!(store.get("rs9357296").unwrap().as_deref(), Some("AG"));
}

#[test]
fn real_world_clean_vcf_supports_locus_lookup_without_rsids() {
    let test_name = "real_world_clean_vcf_supports_locus_lookup_without_rsids";
    let Some(clean_vcf) = shared_fixture_or_skip(test_name, "1k-genomes/vcf/NA06985.clean.vcf.gz")
    else {
        return;
    };
    let Some(original_vcf) = shared_fixture_or_skip(test_name, "1k-genomes/vcf/NA06985.vcf.gz")
    else {
        return;
    };

    let clean_store = GenotypeStore::from_file(&clean_vcf).expect("open cleaned VCF");
    let original_store = GenotypeStore::from_file(&original_vcf).expect("open original VCF");

    let queries = [
        (
            VariantSpec {
                grch37: Some(bioscript_core::GenomicLocus {
                    chrom: "1".to_owned(),
                    start: 12_783,
                    end: 12_783,
                }),
                reference: Some("G".to_owned()),
                alternate: Some("A".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            "GA",
        ),
        (
            VariantSpec {
                grch37: Some(bioscript_core::GenomicLocus {
                    chrom: "1".to_owned(),
                    start: 13_110,
                    end: 13_110,
                }),
                reference: Some("G".to_owned()),
                alternate: Some("A".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            "GA",
        ),
        (
            VariantSpec {
                rsids: vec!["rs78601809".to_owned()],
                grch37: Some(bioscript_core::GenomicLocus {
                    chrom: "1".to_owned(),
                    start: 15_211,
                    end: 15_211,
                }),
                reference: Some("T".to_owned()),
                alternate: Some("G".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            "GG",
        ),
    ];

    for (query, expected_genotype) in queries {
        let clean = clean_store.lookup_variant(&query).expect("clean lookup");
        let original = original_store
            .lookup_variant(&query)
            .expect("original lookup");

        assert_eq!(clean.backend, "vcf");
        assert_eq!(clean.genotype.as_deref(), Some(expected_genotype));
        assert_eq!(original.genotype, clean.genotype);
    }
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
