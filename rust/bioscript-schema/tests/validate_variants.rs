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
fn validate_variants_reports_type_and_metadata_issues() {
    let dir = temp_dir("validate-variant-edges");
    fs::write(
        dir.join("typed-shape.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "2.0"
name: ""
label: 42
gene: ""
summary: ""
tags: "type:trait"
identifiers:
  rsids: "rs1"
coordinates:
  grch37:
    chrom: "0"
    pos: "one"
alleles:
  kind: "other"
  canonical_alt: "G"
  ref: ""
  alts:
    - 1
    - ""
findings:
  - "not-a-map"
  - schema: "bioscript:trait:1.0"
    alt: "G"
provenance:
  sources:
    - "not-a-map"
    - kind: ""
      label: ""
      url: "mailto:example"
"#,
    )
    .unwrap();

    let report = validate_variants_path(&dir).unwrap();
    let text = report.render_text();

    assert_eq!(report.files_scanned, 1);
    assert!(report.total_errors() >= 17, "{text}");
    assert!(report.total_warnings() >= 4, "{text}");
    for expected in [
        "expected '1.0'",
        "expected string",
        "expected a sequence of strings",
        "expected integer",
        "canonical_alt is not part of the current schema",
        "expected one of snv, deletion, insertion, indel",
        "finding alt 'G' is not present",
        "expected http or https URL",
    ] {
        assert!(text.contains(expected), "{expected}\n{text}");
    }
}

#[test]
fn validate_variants_reports_coordinate_edge_cases() {
    let dir = temp_dir("validate-variant-coordinate-edges");
    fs::write(
        dir.join("coordinate-range.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "range"
identifiers:
  rsids:
    - "rs1"
    - "rs1"
coordinates:
  grch38:
    chrom: "MT"
    start: 20
    end: 10
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "N"
"#,
    )
    .unwrap();
    fs::write(
        dir.join("single-position-range.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "single-position-range"
coordinates:
  grch38:
    chrom: "X"
    start: 5
    end: 5
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();
    fs::write(
        dir.join("pos-and-range.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "both-coordinate-styles"
coordinates:
  grch38:
    chrom: "Y"
    pos: 5
    start: 5
    end: 6
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();

    let report = validate_variants_path(&dir).unwrap();
    let text = report.render_text();

    assert_eq!(report.files_scanned, 3);
    assert!(report.total_errors() >= 3, "{text}");
    assert!(report.total_warnings() >= 2, "{text}");
    for expected in [
        "duplicate identifier 'rs1'",
        "expected end >= start",
        "single-position coordinate uses start/end",
        "use either pos or start/end",
    ] {
        assert!(text.contains(expected), "{expected}\n{text}");
    }
}

#[test]
fn validate_panels_reports_missing_empty_and_type_issues() {
    let dir = temp_dir("validate-panel-edges");
    fs::write(
        dir.join("missing-members.yaml"),
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "missing-members"
label: 7
summary: ""
tags: "type:trait"
permissions:
  domains: "https://example.org"
downloads:
  - "not-a-map"
"#,
    )
    .unwrap();
    fs::write(
        dir.join("empty-members.yaml"),
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "empty-members"
permissions:
  domains:
    - 3
    - "https://"
    - "https://example.org:8443"
members: []
"#,
    )
    .unwrap();

    let report = validate_panels_path(&dir).unwrap();
    let text = report.render_text();

    assert_eq!(report.files_scanned, 2);
    assert!(report.total_errors() >= 7, "{text}");
    assert!(report.total_warnings() >= 1, "{text}");
    for expected in [
        "expected a sequence of strings",
        "downloads[0]: expected mapping",
        "members: missing required field",
        "members: expected at least one member",
        "permissions.domains[0]: expected string",
        "invalid URL",
        "expected string",
        "empty string",
    ] {
        assert!(text.contains(expected), "{expected}\n{text}");
    }
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

#[test]
fn remote_resource_resolution_handles_json_versions_and_plain_relative_urls() {
    let resolved = resolve_remote_resource_text(
        "https://example.com/catalogues/index.json",
        "assay.json",
        r#"
{
  "name": "json-assay",
  "assay": {
    "version": "2026.4",
    "panel": "panels/common.yaml"
  },
  "artifact_url": "../artifacts/compiled.json"
}
"#,
    )
    .unwrap();

    assert_eq!(resolved.kind, RemoteResourceKind::Assay);
    assert_eq!(resolved.version.as_deref(), Some("2026.4"));
    let urls = resolved
        .dependencies
        .iter()
        .map(|dependency| dependency.url.as_str())
        .collect::<Vec<_>>();
    assert_eq!(
        urls,
        vec![
            "https://example.com/artifacts/compiled.json",
            "https://example.com/catalogues/panels/common.yaml",
        ]
    );

    let err =
        resolve_remote_resource_text("https://example.com/bad.json", "bad.json", "{").unwrap_err();
    assert!(
        err.contains("failed to parse JSON resource bad.json"),
        "{err}"
    );
}

#[test]
fn validate_variants_covers_remaining_identity_coordinate_and_allele_edges() {
    let dir = temp_dir("validate-variant-more-edges");
    fs::write(
        dir.join("not-a-variant.yaml"),
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "panel-shape"
members: []
"#,
    )
    .unwrap();
    fs::write(
        dir.join("many-errors.yaml"),
        r#"
schema: "bioscript:variant:1.0"
tags:
  - 7
  - ""
identifiers:
  aliases:
    - 7
    - "bad-alias"
    - "rs22"
    - "rs22"
coordinates:
  grch37:
    pos: 12
  grch38:
    chrom: "2"
alleles:
  kind: "snv"
provenance:
  sources:
    - kind: "database"
"#,
    )
    .unwrap();
    fs::write(
        dir.join("range-and-alleles.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "range-and-alleles"
coordinates:
  grch37:
    chrom: "3"
    start: 0
    end: 0
  grch38:
    chrom: "4"
    start: "bad"
    end: 9
alleles:
  kind: "deletion"
  ref: "A"
  alts: "T"
"#,
    )
    .unwrap();
    fs::write(
        dir.join("empty-alts.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "empty-alts"
coordinates:
  grch38:
    chrom: "5"
    pos: 11
alleles:
  kind: "insertion"
  ref: "A"
  alts: []
"#,
    )
    .unwrap();

    let report = validate_variants_path(&dir).unwrap();
    let text = report.render_text();

    assert_eq!(report.files_scanned, 4);
    assert!(report.total_issues() >= 19, "{text}");
    for expected in [
        "missing required field",
        "tags[0]: expected string",
        "tags[1]: empty tag string",
        "identifiers.aliases[0]: expected string",
        "expected rsid like rs123, found 'bad-alias'",
        "duplicate identifier 'rs22'",
        "coordinates.grch37.chrom: missing chrom",
        "coordinates.grch38: expected either pos or start/end",
        "coordinates.grch37.start: expected integer >= 1",
        "coordinates.grch37.end: expected integer >= 1",
        "coordinates.grch38: expected integer start/end",
        "alleles.ref: missing required field",
        "alleles.alts: expected a non-empty sequence of strings",
        "alleles.alts: expected at least one alternate allele",
        "provenance.sources[0].label: missing required field",
        "provenance.sources[0].url: missing required field",
    ] {
        assert!(text.contains(expected), "{expected}\n{text}");
    }
    assert!(!text.contains("panel-shape"));
}

#[test]
fn validate_panels_and_loaders_cover_parse_error_edges() {
    let dir = temp_dir("validate-panel-more-edges");
    let non_panel = dir.join("variant.yaml");
    fs::write(
        &non_panel,
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "rs1"
"#,
    )
    .unwrap();
    let missing_schema = dir.join("missing-schema.yaml");
    fs::write(
        &missing_schema,
        r#"
version: "1.0"
name: "missing-schema"
"#,
    )
    .unwrap();

    let report = validate_panels_path(&dir).unwrap();
    let text = report.render_text();
    assert_eq!(report.files_scanned, 2);
    assert_eq!(report.total_errors(), 1, "{text}");
    assert!(text.contains("missing schema"));
    assert!(!text.contains("rs1"));

    let invalid_panel = dir.join("invalid-panel.yaml");
    fs::write(
        &invalid_panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
downloads:
  - id: ""
    url: "http://"
    sha256: ""
members:
  - download: ""
  - "not-a-map"
"#,
    )
    .unwrap();
    let err = load_panel_manifest(&invalid_panel).unwrap_err();
    assert!(err.contains("name: missing required field"), "{err}");
    assert!(err.contains("downloads[0].id: empty string"), "{err}");
    assert!(err.contains("downloads[0].version: missing required field"), "{err}");
    assert!(err.contains("members[0].kind: missing required field"), "{err}");
    assert!(err.contains("members[0].download: empty string"), "{err}");
    assert!(err.contains("members[1]: expected mapping"), "{err}");

    let downloads_not_mapping = dir.join("downloads-not-mapping.yaml");
    fs::write(
        &downloads_not_mapping,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "bad-download"
downloads:
  - "not-a-map"
members:
  - kind: "variant"
    path: "rs1.yaml"
"#,
    )
    .unwrap();
    let err = load_panel_manifest(&downloads_not_mapping).unwrap_err();
    assert!(err.contains("downloads[0]: expected mapping"), "{err}");

    let members_not_mapping = dir.join("members-not-mapping.yaml");
    fs::write(
        &members_not_mapping,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "bad-member"
members:
  - "not-a-map"
"#,
    )
    .unwrap();
    let err = load_panel_manifest(&members_not_mapping).unwrap_err();
    assert!(err.contains("members[0]: expected mapping"), "{err}");

    let invalid_lookup = load_variant_manifest_text_for_lookup(
        "bad-lookup.yaml",
        r#"
schema: "wrong"
version: "1.0"
name: "bad-lookup"
coordinates:
  grch38:
    chrom: "1"
    pos: 0
alleles:
  kind: "snv"
  ref: "A"
  alts: ["G"]
"#,
    )
    .unwrap_err();
    assert!(invalid_lookup.contains("schema: expected schema"), "{invalid_lookup}");
    assert!(
        invalid_lookup.contains("coordinates.grch38.pos: expected integer >= 1"),
        "{invalid_lookup}"
    );
}
