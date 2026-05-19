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

    let store = GenotypeStore::from_file_with_options(
        &path,
        &GenotypeLoadOptions {
            assembly: Some(Assembly::Grch38),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();

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
    assert_eq!(observation.evidence[0], "resolved by locus chr22:36265860");
    assert!(
        observation
            .evidence
            .get(1)
            .is_some_and(|line| line.contains("source line: rs73885319,chr22,36265860")),
        "{:?}",
        observation.evidence
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

    let store = GenotypeStore::from_file_with_options(
        &path,
        &GenotypeLoadOptions {
            assembly: Some(Assembly::Grch38),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();

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

    let store = GenotypeStore::from_file_with_options(
        &path,
        &GenotypeLoadOptions {
            assembly: Some(Assembly::Grch38),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();

    assert_eq!(store.get("rsSpace").unwrap().as_deref(), Some("CT"));
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
    assert_eq!(observation.evidence[0], "resolved by locus chr2:201");
    assert!(
        observation
            .evidence
            .get(1)
            .is_some_and(|line| line.contains("source line: chrOnly chr2 201")),
        "{:?}",
        observation.evidence
    );
}

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn gsgt_sample_bytes() -> Vec<u8> {
    fs::read(fixtures_dir().join("carigenetics_gsgt_sample.txt")).unwrap()
}

#[test]
fn gsgt_in_memory_path_parses_merges_and_skips() {
    let store =
        GenotypeStore::from_bytes("carigenetics_gsgt_sample.txt", &gsgt_sample_bytes()).unwrap();

    // Plain rs + Plus columns.
    assert_eq!(store.get("rs9000001").unwrap().as_deref(), Some("GG"));
    // Prefixed rs recovered from SNP Name.
    assert_eq!(store.get("rs9000002").unwrap().as_deref(), Some("AG"));
    // Replicate probes: a real call + a no-call replicate -> the call.
    assert_eq!(store.get("rs9000003").unwrap().as_deref(), Some("CC"));
    // Replicate disagreement -> no-call, never auto-picked.
    assert_eq!(store.get("rs9000004").unwrap().as_deref(), Some("--"));
    // Co-located distinct SNP stays addressable by its own rsid.
    assert_eq!(store.get("rs9000007").unwrap().as_deref(), Some("AT"));
    // No-call only.
    assert_eq!(store.get("rs9000011").unwrap().as_deref(), Some("--"));
    // Chr == 0 and Position == 0 rows are skipped entirely.
    assert_eq!(store.get("rs9000009").unwrap(), None);
    assert_eq!(store.get("rs9000010").unwrap(), None);

    // Non-rs probe resolves only by locus.
    let obs = store
        .lookup_variant(&VariantSpec {
            grch38: Some(bioscript_core::GenomicLocus {
                chrom: "2".to_owned(),
                start: 6000,
                end: 6000,
            }),
            ..VariantSpec::default()
        })
        .unwrap();
    assert_eq!(obs.genotype.as_deref(), Some("CT"));
}

#[test]
fn gsgt_streaming_scan_path_merges_replicates() {
    let dir = temp_dir("gsgt-scan");
    let path = dir.join("carigenetics_gsgt_sample.txt");
    fs::write(&path, gsgt_sample_bytes()).unwrap();

    let store =
        GenotypeStore::from_file_with_options(&path, &GenotypeLoadOptions::default()).unwrap();

    let by_rsid = |rsid: &str| {
        store
            .lookup_variant(&VariantSpec {
                rsids: vec![rsid.to_owned()],
                ..VariantSpec::default()
            })
            .unwrap()
            .genotype
    };
    assert_eq!(by_rsid("rs9000001").as_deref(), Some("GG"));
    // call + no-call replicate merged in a single streaming pass.
    assert_eq!(by_rsid("rs9000003").as_deref(), Some("CC"));
    // disagreeing replicates collapse to a no-call.
    assert_eq!(by_rsid("rs9000004").as_deref(), Some("--"));
}

#[test]
fn gsgt_matches_ddna_concordance_on_shared_fixture() {
    let gsgt =
        GenotypeStore::from_bytes("carigenetics_gsgt_sample.txt", &gsgt_sample_bytes()).unwrap();
    let ddna = GenotypeStore::from_bytes(
        "carigenetics_ddna_concordance_sample.txt",
        &fs::read(fixtures_dir().join("carigenetics_ddna_concordance_sample.txt")).unwrap(),
    )
    .unwrap();

    let mut compared = 0;
    for rsid in ["rs9000001", "rs9000002", "rs9000003", "rs9000007"] {
        let g = gsgt.get(rsid).unwrap().unwrap();
        let d = ddna.get(rsid).unwrap().unwrap();
        assert_eq!(g, d, "genotype mismatch for {rsid}: gsgt={g} ddna={d}");
        compared += 1;
    }
    assert_eq!(compared, 4);
}

#[test]
fn ddna_path_unchanged_placeholder_rsids_not_collapsed_by_rsid() {
    // Regression guard (illumina.md Phase 3): keying merges by rsid alone
    // silently collapses every placeholder rsid into one row. DDNA must not
    // go through the GSGT path and must keep distinct loci distinct.
    let dir = temp_dir("ddna-placeholder");
    let path = dir.join("ddna.txt");
    fs::write(
        &path,
        "# rsid\tchromosome\tposition\tgenotype\n\
         .\t1\t100\tAA\n\
         .\t1\t200\tGG\n",
    )
    .unwrap();

    let store = GenotypeStore::from_file_with_options(
        &path,
        &GenotypeLoadOptions {
            assembly: Some(Assembly::Grch38),
            ..GenotypeLoadOptions::default()
        },
    )
    .unwrap();

    let at = |start: i64| {
        store
            .lookup_variant(&VariantSpec {
                grch38: Some(bioscript_core::GenomicLocus {
                    chrom: "1".to_owned(),
                    start,
                    end: start,
                }),
                ..VariantSpec::default()
            })
            .unwrap()
            .genotype
    };
    assert_eq!(at(100).as_deref(), Some("AA"));
    assert_eq!(at(200).as_deref(), Some("GG"));
}

fn carika_dir() -> Option<PathBuf> {
    if let Some(p) = std::env::var_os("BIOVAULT_CARIKA_DIR") {
        let d = PathBuf::from(p);
        if d.exists() {
            return Some(d);
        }
    }
    let default = PathBuf::from("/Users/madhavajay/dev/my_private_data/carika");
    default.exists().then_some(default)
}

fn is_acgt_call(g: &str) -> bool {
    g.len() == 2 && g.chars().all(|c| matches!(c, 'A' | 'C' | 'G' | 'T'))
}

fn is_palindromic(g: &str) -> bool {
    let mut s: Vec<char> = g.chars().collect();
    s.sort_unstable();
    matches!(s.as_slice(), ['A', 'T'] | ['C', 'G'])
}

// Phase 4 (illumina.md): prove the real PC0001 GSGT export and the
// already-supported DDNA export are genotype-concordant. Skips when the
// private data dir is unavailable so CI stays green.
#[test]
#[allow(clippy::cast_precision_loss)]
fn pc0001_gsgt_matches_ddna_real_files() {
    let Some(dir) = carika_dir() else {
        eprintln!("skipping pc0001_gsgt_matches_ddna_real_files: BIOVAULT_CARIKA_DIR unset");
        return;
    };
    let gsgt_path = dir.join("PC0001_Raw Data_Carigenetics.txt");
    let ddna_path = dir.join("PC0001_X_X_GSAv3-DTC_GRCh38-07-29-2025.txt");
    if !gsgt_path.exists() || !ddna_path.exists() {
        eprintln!("skipping pc0001_gsgt_matches_ddna_real_files: sample files missing");
        return;
    }

    let gsgt = GenotypeStore::from_bytes(
        "PC0001_Raw Data_Carigenetics.txt",
        &fs::read(&gsgt_path).unwrap(),
    )
    .unwrap();
    let ddna = GenotypeStore::from_bytes(
        "PC0001_X_X_GSAv3-DTC_GRCh38-07-29-2025.txt",
        &fs::read(&ddna_path).unwrap(),
    )
    .unwrap();

    // rsids come from the DDNA ground-truth file (col 0 == rsNNN).
    let ddna_raw = fs::read_to_string(&ddna_path).unwrap();
    let mut compared = 0u64;
    let mut concordant = 0u64;
    let mut palindromic = 0u64;
    let mut mismatch = 0u64;
    for line in ddna_raw.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let Some(rsid) = line.split('\t').next() else {
            continue;
        };
        if !rsid.starts_with("rs") {
            continue;
        }
        let (Some(g), Some(d)) = (
            gsgt.get(rsid).unwrap(),
            ddna.get(rsid).unwrap(),
        ) else {
            continue;
        };
        if !is_acgt_call(&g) || !is_acgt_call(&d) {
            continue;
        }
        compared += 1;
        if g == d {
            concordant += 1;
        } else if is_palindromic(&g) && is_palindromic(&d) {
            palindromic += 1;
        } else {
            mismatch += 1;
        }
    }

    let pct = 100.0 * concordant as f64 / compared as f64;
    eprintln!(
        "PC0001 concordance: compared={compared} concordant={concordant} \
         palindromic={palindromic} mismatch={mismatch} ({pct:.4}%)",
    );
    assert!(compared > 100_000, "too few comparable SNPs: {compared}");
    let non_palindromic_mismatch_rate = mismatch as f64 / compared as f64;
    assert!(
        non_palindromic_mismatch_rate < 0.001,
        "non-palindromic mismatch rate {non_palindromic_mismatch_rate:.6} exceeds 0.1% \
         (compared={compared} mismatch={mismatch})",
    );
}
