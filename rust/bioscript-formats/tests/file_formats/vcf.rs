use super::*;

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
