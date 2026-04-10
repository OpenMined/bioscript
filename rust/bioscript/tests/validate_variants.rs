use std::{fs, path::PathBuf, process::Command};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

#[test]
fn validate_variants_reports_known_shape_issues() {
    let root = repo_root();
    let fixture_dir = root.join("rust/bioscript/tests/fixtures/variants");
    fs::create_dir_all(&fixture_dir).unwrap();
    let fixture = fixture_dir.join("variant.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:variant"
version: "1.0"
variant_id: "TEST_rs1"
identifiers:
  rsids:
    - "rs1"
coordinates:
  grch38:
    chrom: "1"
    pos: 100
alleles:
  kind: "snv"
  ref: "G"
  alts: ["A"]
research:
  tags:
    - "snp"
clinical:
  pgx:
    drug_labels:
      - source: "FDA"
        title: "Example"
        genes: ["GENE1"]
        drugs: ["drug1"]
        pgx_level: ""
        actionable: false
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("validate-variants")
        .arg(&fixture)
        .output()
        .unwrap();

    assert!(output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stderr.is_empty());
    assert!(stdout.contains("research.tags"));
    assert!(stdout.contains("empty pgx_level string"));
}
