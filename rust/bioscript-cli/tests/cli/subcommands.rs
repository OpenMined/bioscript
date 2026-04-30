use super::*;

#[test]
fn inspect_subcommand_reports_detected_vendor_and_platform() {
    let root = repo_root();
    let path = root.join("rust/bioscript-formats/tests/fixtures/ancestrydna_v2_sample.txt");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("inspect")
        .arg(path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("kind\tgenotype_text"));
    assert!(stdout.contains("vendor\tAncestryDNA"));
    assert!(stdout.contains("platform_version\tV2.0"));
    assert!(stdout.contains("assembly\tgrch37"));
    assert!(stdout.contains("duration_ms\t"));
}

#[test]
fn prepare_subcommand_reports_reference_index_flags() {
    let root = repo_root();
    let dir = temp_dir("prepare-cli");
    let reference = dir.join("ref.fa");
    fs::write(&reference, b">chr1\nACGT\n").unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("prepare")
        .arg("--root")
        .arg(&dir)
        .arg("--reference-file")
        .arg("ref.fa")
        .arg("--cache-dir")
        .arg("cache")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("--reference-file"));
    assert!(stdout.contains("--reference-index"));
    assert!(stdout.contains("cache"));
}

#[test]
fn prepare_subcommand_reports_nothing_to_index_for_noop_auto_request() {
    let root = repo_root();
    let dir = temp_dir("prepare-noop-cli");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("prepare")
        .arg("--root")
        .arg(&dir)
        .arg("--input-format")
        .arg("auto")
        .output()
        .unwrap();

    assert!(output.status.success(), "stderr: {}", stderr_text(&output));
    assert!(String::from_utf8_lossy(&output.stdout).is_empty());
    assert!(
        stderr_text(&output).contains("bioscript prepare: nothing to index"),
        "{}",
        stderr_text(&output)
    );
}

#[test]
fn validate_variants_cli_returns_nonzero_and_writes_report() {
    let root = repo_root();
    let dir = temp_dir("validate-variants-cli");
    let manifest = dir.join("bad-variant.yaml");
    let report = dir.join("reports/variants.txt");
    fs::write(
        &manifest,
        r#"
schema: "bioscript:variant"
version: "1.0"
variant_id: "TEST_bad"
name: "bad"
identifiers:
  rsids:
    - "bad-rsid"
coordinates:
  grch38:
    chrom: "chrUn"
    pos: 0
alleles:
  kind: "snv"
  ref: "AA"
  alts: []
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("validate-variants")
        .arg(&manifest)
        .arg("--report")
        .arg(&report)
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("validation found"), "{stderr}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("errors"), "{stdout}");
    let report_text = fs::read_to_string(report).unwrap();
    assert!(report_text.contains("bad-rsid"));
}

#[test]
fn validate_panels_cli_returns_nonzero_and_writes_report() {
    let root = repo_root();
    let dir = temp_dir("validate-panels-cli");
    let panel = dir.join("bad-panel.yaml");
    let report = dir.join("reports/panels.txt");
    fs::write(
        &panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "bad-panel"
members:
  - kind: "variant"
    path: "../outside.yaml"
    sha256: "not-a-sha"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("validate-panels")
        .arg(&panel)
        .arg("--report")
        .arg(&report)
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("validation found"), "{stderr}");
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("errors"), "{stdout}");
    let report_text = fs::read_to_string(report).unwrap();
    assert!(report_text.contains("members[0].sha256"));
}
