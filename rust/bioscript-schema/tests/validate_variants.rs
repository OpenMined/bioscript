use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_schema::{
    RemoteResourceKind, load_panel_manifest, load_variant_manifest_text,
    load_variant_manifest_text_for_lookup, resolve_remote_resource_text, validate_panels_path,
    validate_variants_path,
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
fn validate_variants_scans_nested_yaml_files_and_ignores_other_files() {
    let dir = temp_dir("validate-dir");
    let nested = dir.join("nested");
    fs::create_dir_all(&nested).unwrap();
    fs::write(dir.join("notes.txt"), "not yaml").unwrap();
    fs::write(
        dir.join("valid.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "rs1"
identifiers:
  rsids: ["rs1"]
alleles:
  kind: "snv"
  ref: "A"
  alts: ["G"]
"#,
    )
    .unwrap();
    fs::write(
        nested.join("missing-schema.yml"),
        r#"
version: "1.0"
name: "rs2"
"#,
    )
    .unwrap();

    let report = validate_variants_path(&dir).unwrap();
    let text = report.render_text();

    assert_eq!(report.files_scanned, 2);
    assert!(report.has_errors());
    assert_eq!(report.total_errors(), 1);
    assert!(text.contains("missing-schema.yml"));
    assert!(text.contains("missing schema"));
    assert!(!text.contains("notes.txt"));
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
fn load_panel_manifest_parses_downloads_permissions_and_member_metadata() {
    let dir = temp_dir("load-panel");
    let fixture = dir.join("panel.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "traits-common"
tags: ["type:trait"]
permissions:
  domains:
    - "https://example.org"
downloads:
  - id: "remote-rs1"
    url: "https://example.org/variants/rs1.yaml"
    sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    version: "2026-01-01"
members:
  - kind: "variant"
    download: "remote-rs1"
    sha256: "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
    version: "2026-01-01"
"#,
    )
    .unwrap();

    let panel = load_panel_manifest(&fixture).unwrap();

    assert_eq!(panel.name, "traits-common");
    assert_eq!(panel.tags, vec!["type:trait"]);
    assert_eq!(panel.permissions.domains, vec!["https://example.org"]);
    assert_eq!(panel.downloads.len(), 1);
    assert_eq!(panel.downloads[0].origin, "https://example.org");
    assert_eq!(panel.members.len(), 1);
    assert_eq!(panel.members[0].download.as_deref(), Some("remote-rs1"));
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
fn validate_panels_reports_member_and_download_shape_issues() {
    let dir = temp_dir("validate-panel-shape");
    let fixture = dir.join("panel.yaml");
    fs::write(
        &fixture,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "traits-common"
permissions:
  domains:
    - "https://example.org/path"
    - "ftp://example.org"
    - "https://example.org"
    - "https://example.org"
downloads:
  - id: "remote-rs1"
    url: "https://example.org/variants/rs1.yaml"
    sha256: "not-a-sha"
    version: "1.0"
  - id: "remote-rs1"
    url: "file:///tmp/rs1.yaml"
    sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    version: ""
members:
  - kind: "script"
    path: "variants/rs1.yaml"
    download: "remote-rs1"
    sha256: "not-a-sha"
    version: ""
  - kind: "variant"
    download: "missing-download"
  - kind: "variant"
    path: ""
"#,
    )
    .unwrap();

    let report = validate_panels_path(&fixture).unwrap();
    let text = report.render_text();

    assert!(report.total_errors() >= 11, "{text}");
    assert!(report.total_warnings() >= 1, "{text}");
    assert!(text.contains("expected origin only"));
    assert!(text.contains("expected http or https origin"));
    assert!(text.contains("duplicate origin"));
    assert!(text.contains("duplicate download id"));
    assert!(text.contains("unknown download id"));
    assert!(text.contains("unsupported member kind"));
    assert!(text.contains("expected exactly one of path or download"));
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

#[test]
fn remote_resource_resolution_classifies_python_without_parsing() {
    let text = "print('hello from remote bioscript')\n";

    let resolved = resolve_remote_resource_text(
        "https://github.com/OpenMined/bioscript/blob/main/example.py",
        "example.py",
        text,
    )
    .unwrap();

    assert_eq!(resolved.kind, RemoteResourceKind::Python);
    assert_eq!(resolved.schema, None);
    assert_eq!(resolved.title, "example.py");
    assert_eq!(
        resolved.sha256,
        "b6d9c1ee20c7fb054ebd7defd271d7956b25d8d0c3ef451eaf6adcfda8a61b0f"
    );
}

#[test]
fn remote_resource_resolution_classifies_schema_kinds() {
    let cases = [
        (
            "variant.yaml",
            "bioscript:variant:1.0",
            RemoteResourceKind::Variant,
        ),
        (
            "panel.yaml",
            "bioscript:panel:1.0",
            RemoteResourceKind::Panel,
        ),
        (
            "catalogue.yaml",
            "bioscript:catalogue:1.0",
            RemoteResourceKind::Catalogue,
        ),
        (
            "assay.yaml",
            "bioscript:assay:1.0",
            RemoteResourceKind::Assay,
        ),
    ];

    for (name, schema, expected) in cases {
        let text = format!(
            r#"
schema: "{schema}"
version: "1.0"
name: "{name}"
"#
        );

        let resolved =
            resolve_remote_resource_text("https://example.com/resources/index.yaml", name, &text)
                .unwrap();

        assert_eq!(resolved.kind, expected, "{name}");
        assert_eq!(resolved.schema.as_deref(), Some(schema), "{name}");
    }
}

#[test]
fn remote_resource_resolution_infers_kind_from_fields_without_schema() {
    let cases = [
        ("members.yaml", "members: []\n", RemoteResourceKind::Panel),
        ("variants.yaml", "variants: []\n", RemoteResourceKind::Panel),
        (
            "catalogue.yaml",
            "assays: []\n",
            RemoteResourceKind::Catalogue,
        ),
        (
            "assay.yaml",
            "assay:\n  package_version: \"2026.1\"\n",
            RemoteResourceKind::Assay,
        ),
        (
            "variant.yaml",
            "variant_id: TEST_rs1\ncoordinates: {}\n",
            RemoteResourceKind::Variant,
        ),
        (
            "unknown.yaml",
            "name: just-a-file\n",
            RemoteResourceKind::Unknown,
        ),
    ];

    for (name, text, expected) in cases {
        let resolved =
            resolve_remote_resource_text("https://example.com/resources/index.yaml", name, text)
                .unwrap();

        assert_eq!(resolved.kind, expected, "{name}");
    }
}

#[test]
fn remote_resource_resolution_resolves_github_dependencies_and_dedupes_urls() {
    let text = r#"
schema: "bioscript:panel:1.0"
version: "1.0"
members:
  - path: "variants/rs1.yaml"
    version: "1.1"
  - path: "variants/rs1.yaml"
    version: "1.1"
  - path: "/shared/rs2.yaml"
downloads:
  - url: "https://example.com/reference.json"
"#;

    let resolved = resolve_remote_resource_text(
        "https://github.com/OpenMined/bioscript/blob/main/panels/panel.yaml",
        "panel.yaml",
        text,
    )
    .unwrap();

    let urls = resolved
        .dependencies
        .iter()
        .map(|dependency| dependency.url.as_str())
        .collect::<Vec<_>>();

    assert_eq!(
        urls,
        vec![
            "https://example.com/reference.json",
            "https://github.com/OpenMined/bioscript/blob/main/panels/variants/rs1.yaml",
            "https://github.com/OpenMined/bioscript/blob/main/shared/rs2.yaml",
        ]
    );
    assert_eq!(resolved.dependencies[0].kind, "download");
    assert_eq!(resolved.dependencies[1].kind, "member");
    assert_eq!(resolved.dependencies[1].version.as_deref(), Some("1.1"));
}

#[test]
fn remote_resource_resolution_reports_invalid_structured_text() {
    let err = resolve_remote_resource_text("https://example.com/bad.yaml", "bad.yaml", ":\n")
        .unwrap_err();

    assert!(
        err.contains("failed to parse YAML resource bad.yaml"),
        "{err}"
    );
}
