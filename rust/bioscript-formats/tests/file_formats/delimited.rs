use super::*;

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
