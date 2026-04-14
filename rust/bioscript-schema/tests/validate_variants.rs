use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_schema::validate_variants_path;

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-schema-{label}-{}-{nanos}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn validate_variants_reports_known_shape_issues() {
    let dir = temp_dir("validate");
    let fixture = dir.join("variant.yaml");
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

    let report = validate_variants_path(&fixture).unwrap();
    let text = report.render_text();

    assert_eq!(report.total_errors(), 0);
    assert_eq!(report.total_warnings(), 2);
    assert!(text.contains("research.tags"));
    assert!(text.contains("empty pgx_level string"));
}
