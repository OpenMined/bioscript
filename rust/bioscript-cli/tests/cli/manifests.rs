use super::*;

#[test]
fn variant_manifest_runs_directly_via_cli() {
    let root = repo_root();
    let dir = temp_dir("variant-manifest");
    let manifest = dir.join("rs1.yaml");
    fs::write(
        &manifest,
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(&manifest)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("kind\tname\tpath"));
    assert!(stdout.contains("example-rs73885319"));
    assert!(stdout.contains("AG"));
}

#[test]
fn variant_manifest_requires_input_file() {
    let root = repo_root();
    let dir = temp_dir("variant-manifest-missing-input");
    let manifest = dir.join("rs1.yaml");
    fs::write(
        &manifest,
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs1"
identifiers:
  rsids:
    - "rs1"
coordinates:
  grch38:
    chrom: "1"
    pos: 10
alleles:
  kind: "snv"
  ref: "A"
  alts: ["G"]
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg(&manifest)
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        stderr_text(&output).contains("manifest execution requires --input-file"),
        "{}",
        stderr_text(&output)
    );
}

#[test]
fn variant_manifest_writes_output_trace_and_participant_id() {
    let root = repo_root();
    let dir = temp_dir("variant-manifest-output");
    let manifest = dir.join("rs1.yaml");
    let output_path = dir.join("reports/variant.tsv");
    let trace_path = dir.join("reports/variant.trace.tsv");
    fs::write(
        &manifest,
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg("--output-file")
        .arg(&output_path)
        .arg("--participant-id")
        .arg("participant-1")
        .arg("--trace-report")
        .arg(&trace_path)
        .arg(&manifest)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).is_empty());
    let table = fs::read_to_string(output_path).unwrap();
    assert!(table.contains("participant-1"), "{table}");
    assert!(table.contains("example-rs73885319"), "{table}");
    let trace = fs::read_to_string(trace_path).unwrap();
    assert!(trace.contains("step\tline\tcode"), "{trace}");
    assert!(trace.contains("rs1.yaml"), "{trace}");
}

#[test]
fn panel_manifest_runs_directly_via_cli() {
    let root = repo_root();
    let dir = temp_dir("panel-manifest");
    let variants_dir = dir.join("variants");
    fs::create_dir_all(&variants_dir).unwrap();
    fs::write(
        variants_dir.join("rs73885319.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();
    fs::write(
        variants_dir.join("rs60910145.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs60910145"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs60910145"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265988
alleles:
  kind: "snv"
  ref: "T"
  alts:
    - "G"
"#,
    )
    .unwrap();
    let panel = dir.join("panel.yaml");
    fs::write(
        &panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "example-panel"
tags:
  - "type:trait"
members:
  - kind: "variant"
    path: "variants/rs73885319.yaml"
    version: "1.0"
  - kind: "variant"
    path: "variants/rs60910145.yaml"
    version: "1.0"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg("--filter")
        .arg("name=rs73885319")
        .arg(&panel)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("example-rs73885319"));
    assert!(!stdout.contains("example-rs60910145"));
}

#[test]
fn panel_manifest_filters_by_kind_tag_path_and_rejects_unknown_filter_keys() {
    let root = repo_root();
    let dir = temp_dir("panel-filters");
    let variants_dir = dir.join("variants");
    fs::create_dir_all(&variants_dir).unwrap();
    fs::write(
        variants_dir.join("rs73885319.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();
    let panel = dir.join("panel.yaml");
    fs::write(
        &panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "example-panel"
members:
  - kind: "variant"
    path: "variants/rs73885319.yaml"
    version: "1.0"
"#,
    )
    .unwrap();

    let matched = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg("--filter")
        .arg("kind=variant")
        .arg("--filter")
        .arg("tag=type:trait")
        .arg("--filter")
        .arg("path=rs73885319")
        .arg(&panel)
        .output()
        .unwrap();

    assert!(
        matched.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&matched.stderr)
    );
    let stdout = String::from_utf8_lossy(&matched.stdout);
    assert!(stdout.contains("example-rs73885319"), "{stdout}");

    let filtered_out = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg("--filter")
        .arg("unknown=value")
        .arg(&panel)
        .output()
        .unwrap();

    assert!(
        filtered_out.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&filtered_out.stderr)
    );
    let stdout = String::from_utf8_lossy(&filtered_out.stdout);
    assert!(stdout.starts_with("kind\tname\tpath"), "{stdout}");
    assert!(!stdout.contains("example-rs73885319"), "{stdout}");
}

#[test]
fn panel_manifest_reports_remote_members_as_not_executable_yet() {
    let root = repo_root();
    let dir = temp_dir("panel-remote-member");
    let panel = dir.join("panel.yaml");
    fs::write(
        &panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "remote-panel"
permissions:
  domains:
    - "https://example.com"
downloads:
  - id: "remote-rs73885319"
    url: "https://example.com/rs73885319.yaml"
    sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
    version: "1.0"
members:
  - kind: "variant"
    download: "remote-rs73885319"
    version: "1.0"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(&panel)
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("remote panel members are not executable yet"),
        "{stderr}"
    );
}

#[test]
fn panel_manifest_reports_non_variant_members_as_not_executable_yet() {
    let root = repo_root();
    let dir = temp_dir("panel-nonvariant-member");
    let panel = dir.join("panel.yaml");
    fs::write(
        &panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "mixed-panel"
members:
  - kind: "script"
    path: "script.py"
    version: "1.0"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(&panel)
        .output()
        .unwrap();

    assert!(!output.status.success());
    assert!(
        stderr_text(&output).contains("unsupported member kind 'script'"),
        "{}",
        stderr_text(&output)
    );
}
