use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_schema::{
    RemoteResourceKind, load_variant_manifest_text, load_variant_manifest_text_for_lookup,
    resolve_remote_resource_text, validate_panels_path, validate_variants_path,
};

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
fn validate_variants_reports_shape_issues() {
    let dir = temp_dir("validate");
    let fixture = dir.join("variant.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:variant"
version: "1.0"
variant_id: "TEST_rs1"
name: "test-rs1"
identifiers:
  rsids:
    - "bad-rsid"
coordinates:
  grch38:
    chrom: "HG7_PATCH"
    pos: 0
alleles:
  kind: "snv"
  ref: "I"
  alts: ["D"]
findings:
  - schema: ""
    alt: "A"
"#,
    )
    .unwrap();

    let report = validate_variants_path(&fixture).unwrap();
    let text = report.render_text();

    assert_eq!(report.total_errors(), 9);
    assert_eq!(report.total_warnings(), 3);
    assert!(text.contains("legacy schema value"));
    assert!(text.contains("invalid chromosome"));
    assert!(text.contains("symbolic I/D alleles are not allowed"));
}

#[test]
fn validate_variants_accepts_current_shape() {
    let dir = temp_dir("validate-tags");
    let fixture = dir.join("rs1.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "test-rs1"
tags:
  - "type:trait"
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
findings:
  - schema: "bioscript:trait:1.0"
    alt: "A"
    summary: "Example finding"
  - schema: "bioscript:pgx:1.0"
    alt: "*"
    summary: "Example multiallelic finding"
provenance:
  sources:
    - kind: "database"
      label: "dbSNP"
      url: "https://example.org/rs1"
"#,
    )
    .unwrap();

    let report = validate_variants_path(&fixture).unwrap();
    assert_eq!(report.total_errors(), 0);
    assert_eq!(report.total_warnings(), 0);
}

#[test]
fn load_variant_manifest_text_accepts_start_end_coordinates() {
    let manifest = load_variant_manifest_text(
        "rs71338792.yaml",
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "GLP1-nature-23andme-rs71338792-INDEL"
identifiers:
  rsids:
    - "rs71338792"
coordinates:
  grch38:
    chrom: "19"
    start: 45679774
    end: 45679786
alleles:
  kind: "indel"
  ref: "TTTTTTTTTTTTT"
  alts:
    - "TTTTTTTTTTTT"
"#,
    )
    .unwrap();

    let grch38 = manifest.spec.grch38.expect("grch38 locus");
    assert_eq!(grch38.chrom, "19");
    assert_eq!(grch38.start, 45_679_774);
    assert_eq!(grch38.end, 45_679_786);
}

#[test]
fn lookup_compile_allows_non_execution_metadata_issues() {
    let text = r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "GLP1-nature-23andme-rs10305420-C-T"
identifiers:
  rsids:
    - "rs10305420"
coordinates:
  grch38:
    chrom: "6"
    pos: 39048860
alleles:
  kind: "snv"
  ref: "C"
  alts:
    - "G"
    - "T"
findings:
  - alt: "T"
    notes: "finding metadata without schema should not block lookup compilation"
"#;

    assert!(load_variant_manifest_text("rs10305420.yaml", text).is_err());
    let manifest = load_variant_manifest_text_for_lookup("rs10305420.yaml", text).unwrap();
    assert_eq!(manifest.name, "GLP1-nature-23andme-rs10305420-C-T");
    assert_eq!(manifest.spec.grch38.unwrap().start, 39_048_860);
    assert_eq!(manifest.spec.alternate.as_deref(), Some("T"));
}

#[test]
fn validate_panels_checks_permissions_and_download_origins() {
    let dir = temp_dir("validate-panel");
    let fixture = dir.join("panel.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "traits-common"
permissions:
  domains:
    - "https://example.org"
downloads:
  - id: "remote-rs1"
    url: "https://cdn.example.org/variants/rs1.yaml"
    sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    version: "1.0"
members:
  - kind: "variant"
    download: "remote-rs1"
"#,
    )
    .unwrap();

    let report = validate_panels_path(&fixture).unwrap();
    let text = report.render_text();

    assert_eq!(report.total_errors(), 1);
    assert!(text.contains("not listed in permissions.domains"));
}

#[test]
fn validate_panels_accepts_local_variant_members() {
    let dir = temp_dir("validate-panel-ok");
    let fixture = dir.join("panel.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "traits-common"
tags:
  - "type:trait"
members:
  - kind: "variant"
    path: "variants/rs671.yaml"
    version: "1.0"
"#,
    )
    .unwrap();

    let report = validate_panels_path(&fixture).unwrap();
    assert_eq!(report.total_errors(), 0);
    assert_eq!(report.total_warnings(), 0);
}

#[test]
fn remote_resource_resolution_detects_panel_members() {
    let text = r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "snp-subset"
label: "SNPedia-backed Shared SNP Subset Panel"
members:
  - kind: "variant"
    path: "variants/rs2048683.yaml"
    version: "1.0"
  - kind: "variant"
    path: "variants/rs4240157.yaml"
    version: "1.0"
"#;

    let resolved = resolve_remote_resource_text(
        "https://github.com/madhavajay/exvitae/blob/main/assays/pgx/GLP1/panel.yaml",
        "panel.yaml",
        text,
    )
    .unwrap();

    assert_eq!(resolved.kind, RemoteResourceKind::Panel);
    assert_eq!(resolved.schema.as_deref(), Some("bioscript:panel:1.0"));
    assert_eq!(resolved.version.as_deref(), Some("1.0"));
    assert_eq!(resolved.dependencies.len(), 2);
    assert_eq!(resolved.dependencies[0].version.as_deref(), Some("1.0"));
    assert_eq!(
        resolved.dependencies[0].url,
        "https://github.com/madhavajay/exvitae/blob/main/assays/pgx/GLP1/variants/rs2048683.yaml"
    );
}
