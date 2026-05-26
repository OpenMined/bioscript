use super::*;

use bioscript_core::GenomicLocus;
use bioscript_formats::{
    DetectedKind, FileContainer, InspectOptions, convert_23andme_grch37_to_grch38, inspect_file,
};
use noodles::core::Position;
use noodles::vcf::{
    self,
    header::record::value::{
        Map,
        map::{
            Contig, Filter, Format,
            format::{Number, Type},
        },
    },
    variant::{
        io::Write as _,
        record_buf::{
            Samples,
            samples::{Keys, sample::Value},
        },
    },
};
use std::io::{BufRead, BufReader};

fn locus(chrom: &str, start: i64, end: i64) -> GenomicLocus {
    GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    }
}

fn bcf_bytes(records: &[(&str, i64, &str, &str, &str, &[f32])]) -> Vec<u8> {
    let header = vcf::Header::builder()
        .add_filter("PASS", Map::<Filter>::pass())
        .add_contig("chr6", Map::<Contig>::new())
        .add_contig("chr19", Map::<Contig>::new())
        .add_format(
            "HDS",
            Map::<Format>::new(Number::Unknown, Type::Float, "haploid alternate dosage"),
        )
        .add_sample_name("1")
        .build();

    let mut data = Vec::new();
    let mut writer = noodles::bcf::io::Writer::new(&mut data);
    writer.write_header(&header).unwrap();

    let keys: Keys = ["HDS".to_owned()].into_iter().collect();
    for (chrom, pos, id, reference, alternate, hds) in records {
        let samples = Samples::new(
            keys.clone(),
            vec![vec![Some(Value::from(
                hds.iter().copied().map(Some).collect::<Vec<_>>(),
            ))]],
        );
        let mut builder = vcf::variant::RecordBuf::builder()
            .set_reference_sequence_name(*chrom)
            .set_variant_start(Position::try_from(usize::try_from(*pos).unwrap()).unwrap())
            .set_reference_bases(*reference)
            .set_alternate_bases(vec![(*alternate).to_owned()].into())
            .set_samples(samples);
        if *id != "." {
            builder = builder.set_ids([(*id).to_owned()].into_iter().collect());
        }
        let record = builder.build();
        writer.write_variant_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
    drop(writer);
    data
}

#[test]
fn individual_bcf_hds_records_are_readable_by_locus() {
    let dir = temp_dir("bcf-hds");
    let path = dir.join("chr6.bcf");
    fs::write(
        &path,
        bcf_bytes(&[
            ("chr6", 39_011_352, ".", "A", "G", &[0.0, 0.0]),
            ("chr6", 39_048_860, ".", "C", "T", &[1.0, 1.0]),
            ("chr6", 39_051_898, ".", "G", "T", &[0.0, 1.0]),
        ]),
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();
    assert_eq!(store.backend_name(), "bcf");

    let observations = store
        .lookup_variants(&[
            VariantSpec {
                grch38: Some(locus("6", 39_011_352, 39_011_352)),
                reference: Some("A".to_owned()),
                alternate: Some("G".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            VariantSpec {
                grch38: Some(locus("6", 39_048_860, 39_048_860)),
                reference: Some("C".to_owned()),
                alternate: Some("T".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            VariantSpec {
                grch38: Some(locus("6", 39_051_898, 39_051_898)),
                reference: Some("G".to_owned()),
                alternate: Some("T".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
        ])
        .unwrap();

    assert_eq!(observations[0].genotype.as_deref(), Some("AA"));
    assert_eq!(observations[1].genotype.as_deref(), Some("TT"));
    assert_eq!(observations[2].genotype.as_deref(), Some("GT"));
    assert_eq!(observations[0].assembly, Some(Assembly::Grch38));
    assert!(
        observations[2].evidence[1].contains("HDS=0,1"),
        "{:?}",
        observations[2].evidence
    );
}

#[test]
fn zip_of_bcf_shards_is_treated_as_one_logical_store() {
    let dir = temp_dir("bcf-zip");
    let zip_path = dir.join("imputed_genotype_data_r6_sample.zip");
    let chr6 = bcf_bytes(&[("chr6", 39_051_898, ".", "G", "T", &[0.0, 1.0])]);
    let chr19 = bcf_bytes(&[("chr19", 45_678_134, ".", "G", "C", &[1.0, 0.0])]);

    let file = fs::File::create(&zip_path).unwrap();
    let mut writer = zip::ZipWriter::new(file);
    writer
        .start_file("fbbdf42da6/chr6.bcf", SimpleFileOptions::default())
        .unwrap();
    writer.write_all(&chr6).unwrap();
    writer
        .start_file("fbbdf42da6/chr19.bcf", SimpleFileOptions::default())
        .unwrap();
    writer.write_all(&chr19).unwrap();
    writer.finish().unwrap();

    let store = GenotypeStore::from_file(&zip_path).unwrap();
    assert_eq!(store.backend_name(), "bcf");

    let observations = store
        .lookup_variants(&[
            VariantSpec {
                grch38: Some(locus("19", 45_678_134, 45_678_134)),
                reference: Some("G".to_owned()),
                alternate: Some("C".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            VariantSpec {
                grch38: Some(locus("6", 39_051_898, 39_051_898)),
                reference: Some("G".to_owned()),
                alternate: Some("T".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
        ])
        .unwrap();

    assert_eq!(observations[0].genotype.as_deref(), Some("CG"));
    assert_eq!(observations[1].genotype.as_deref(), Some("GT"));

    let inspection = inspect_file(&zip_path, &InspectOptions::default()).unwrap();
    assert_eq!(inspection.container, FileContainer::Zip);
    assert_eq!(inspection.detected_kind, DetectedKind::Bcf);
    assert_eq!(inspection.assembly, Some(Assembly::Grch38));
}

#[test]
fn bcf_left_anchored_indel_matches_23andme_style_catalog_variant() {
    let dir = temp_dir("bcf-indel");
    let path = dir.join("chr19.bcf");
    fs::write(
        &path,
        bcf_bytes(&[
            ("chr19", 45_679_773, ".", "A", "AT", &[1.0, 0.0]),
            ("chr19", 45_679_773, ".", "A", "ATT", &[0.0, 0.0]),
        ]),
    )
    .unwrap();

    let store = GenotypeStore::from_file(&path).unwrap();
    let observation = store
        .lookup_variant(&VariantSpec {
            rsids: vec!["rs71338792".to_owned()],
            grch38: Some(locus("19", 45_679_774, 45_679_786)),
            reference: Some("TTTTTTTTTTTTT".to_owned()),
            alternate: Some("TTTTTTTTTTTTTT".to_owned()),
            kind: Some(VariantKind::Indel),
            ..VariantSpec::default()
        })
        .unwrap();

    assert_eq!(observation.genotype.as_deref(), Some("DI"));
    assert!(
        observation.evidence[0].contains("resolved by locus chr19:45679773"),
        "{:?}",
        observation.evidence
    );
}

#[test]
fn local_23andme_lifted_text_matches_real_bcf_zip_for_known_variant_when_enabled() {
    if std::env::var_os("BIOSCRIPT_RUN_LOCAL_23ANDME_BCF").is_none() {
        eprintln!("skipping local 23andMe BCF concordance; set BIOSCRIPT_RUN_LOCAL_23ANDME_BCF=1");
        return;
    }

    let workspace = repo_root().parent().expect("workspace root").to_path_buf();
    let raw_txt = workspace.join("genome_Madhava_Jay_v4_Full_20250613052552.txt");
    let bcf_zip = workspace
        .join("exvitae")
        .join("imputed_genotype_data_r6_Madhava_Jay.zip");
    if !raw_txt.exists() || !bcf_zip.exists() {
        eprintln!(
            "skipping local 23andMe BCF concordance; missing {} or {}",
            raw_txt.display(),
            bcf_zip.display()
        );
        return;
    }

    let dir = temp_dir("real-23andme-bcf");
    let lifted = dir.join("genome.grch38.txt");
    let unmapped = dir.join("unmapped.tsv");
    convert_23andme_grch37_to_grch38(&raw_txt, &lifted, &unmapped).unwrap();

    let lifted_store = GenotypeStore::from_file(&lifted).unwrap();
    let bcf_store = GenotypeStore::from_file(&bcf_zip).unwrap();
    let variant = VariantSpec {
        rsids: vec!["rs1800437".to_owned()],
        grch38: Some(locus("19", 45_678_134, 45_678_134)),
        reference: Some("G".to_owned()),
        alternate: Some("C".to_owned()),
        kind: Some(VariantKind::Snp),
        ..VariantSpec::default()
    };

    let lifted_observation = lifted_store.lookup_variant(&variant).unwrap();
    let bcf_observation = bcf_store.lookup_variant(&variant).unwrap();
    assert_eq!(lifted_observation.genotype.as_deref(), Some("CG"));
    assert_eq!(bcf_observation.genotype.as_deref(), Some("CG"));
    assert_eq!(lifted_observation.genotype, bcf_observation.genotype);
}

#[test]
fn local_23andme_lifted_text_broadly_matches_real_bcf_zip_when_enabled() {
    if std::env::var_os("BIOSCRIPT_RUN_LOCAL_23ANDME_BCF").is_none() {
        eprintln!(
            "skipping local 23andMe BCF broad concordance; set BIOSCRIPT_RUN_LOCAL_23ANDME_BCF=1"
        );
        return;
    }

    let workspace = repo_root().parent().expect("workspace root").to_path_buf();
    let raw_txt = workspace.join("genome_Madhava_Jay_v4_Full_20250613052552.txt");
    let bcf_zip = workspace
        .join("exvitae")
        .join("imputed_genotype_data_r6_Madhava_Jay.zip");
    if !raw_txt.exists() || !bcf_zip.exists() {
        eprintln!(
            "skipping local 23andMe BCF broad concordance; missing {} or {}",
            raw_txt.display(),
            bcf_zip.display()
        );
        return;
    }

    let dir = temp_dir("real-23andme-bcf-broad");
    let lifted = dir.join("genome.grch38.txt");
    let unmapped = dir.join("unmapped.tsv");
    convert_23andme_grch37_to_grch38(&raw_txt, &lifted, &unmapped).unwrap();

    let mut txt_genotypes = Vec::new();
    let mut variants = Vec::new();
    let file = fs::File::open(&lifted).unwrap();
    for line in BufReader::new(file).lines() {
        let line = line.unwrap();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 4 {
            continue;
        }
        let chrom = fields[1];
        if matches!(chrom, "Y" | "MT" | "M") {
            continue;
        }
        let genotype = fields[3].trim().to_ascii_uppercase();
        if genotype.is_empty()
            || genotype == "--"
            || !genotype
                .chars()
                .all(|base| matches!(base, 'A' | 'C' | 'G' | 'T'))
        {
            continue;
        }
        let Ok(pos) = fields[2].parse::<i64>() else {
            continue;
        };
        txt_genotypes.push(normalized_acgt_genotype(&genotype));
        variants.push(VariantSpec {
            rsids: vec![fields[0].to_owned()],
            grch38: Some(locus(chrom, pos, pos)),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        });
    }

    let bcf_store = GenotypeStore::from_file(&bcf_zip).unwrap();
    let observations = bcf_store.lookup_variants(&variants).unwrap();

    let mut queried = 0usize;
    let mut found = 0usize;
    let mut matched = 0usize;
    let mut mismatched = 0usize;
    let mut examples = Vec::new();
    for (idx, observation) in observations.iter().enumerate() {
        queried += 1;
        let Some(bcf_genotype) = observation.genotype.as_deref() else {
            continue;
        };
        found += 1;
        if normalized_acgt_genotype(bcf_genotype) == txt_genotypes[idx] {
            matched += 1;
        } else {
            mismatched += 1;
            if examples.len() < 8 {
                let locus = variants[idx].grch38.as_ref().unwrap();
                examples.push(format!(
                    "{} {}:{} txt={} bcf={}",
                    variants[idx].rsids[0],
                    locus.chrom,
                    locus.start,
                    txt_genotypes[idx],
                    bcf_genotype
                ));
            }
        }
    }

    eprintln!(
        "local 23andMe BCF broad concordance: queried={queried} found={found} matched={matched} mismatched={mismatched} missing={}",
        queried - found
    );
    if !examples.is_empty() {
        eprintln!("mismatch examples: {}", examples.join("; "));
    }
    assert!(queried > 500_000, "expected broad 23andMe SNP comparison");
    assert!(
        found > 500_000,
        "expected most lifted SNPs to be present in BCF"
    );
}

fn normalized_acgt_genotype(genotype: &str) -> String {
    let mut chars: Vec<char> = genotype
        .chars()
        .filter(|base| matches!(base, 'A' | 'C' | 'G' | 'T'))
        .collect();
    chars.sort_unstable();
    chars.into_iter().collect()
}
