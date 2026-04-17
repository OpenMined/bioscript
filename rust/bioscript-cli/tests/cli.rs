use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-cli-tests-tmp-{label}-{}-{nanos}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn hello_world_script_runs_via_cli_and_writes_within_root() {
    let root = repo_root();
    let output_path = root.join("bioscripts/output/hello-world.txt");
    if output_path.exists() {
        fs::remove_file(&output_path).unwrap();
    }

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("bioscripts/hello-world.py")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("hello from bioscript"));
    assert!(stdout.contains("2 + 3 = 5"));

    let written = fs::read_to_string(output_path).unwrap();
    assert!(written.contains("hello from bioscript"));
    assert!(written.contains("loaded: sample input for bioscript"));
}

#[test]
fn path_escape_is_rejected() {
    let root = repo_root();
    let script = root.join("rust/bioscript-cli/tests/fixtures/path_escape.py");
    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg(script)
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("path escapes bioscript root"));
}

#[test]
fn trace_report_is_written_for_hello_world() {
    let root = repo_root();
    let trace_path = root.join("bioscripts/output/hello-world.trace.tsv");
    if trace_path.exists() {
        fs::remove_file(&trace_path).unwrap();
    }

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--trace-report")
        .arg("bioscripts/output/hello-world.trace.tsv")
        .arg("bioscripts/hello-world.py")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let trace = fs::read_to_string(trace_path).unwrap();
    assert!(trace.contains("step\tline\tcode"));
    assert!(trace.contains("hello from bioscript"));
}

#[test]
fn batch_lookup_query_plan_runs_and_preserves_requested_result_order() {
    let root = repo_root();
    let script = root.join("rust/bioscript-cli/tests/fixtures/batch_lookup.py");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(script)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("AG"));
    assert!(stdout.contains("TC"));
    assert!(stdout.contains("II"));
}

#[test]
fn lookup_variant_details_returns_counts_and_decision_fields() {
    let root = repo_root();
    let script = root.join("rust/bioscript-cli/tests/fixtures/lookup_details.py");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(script)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("VariantObservation"));
    assert!(stdout.contains("genotype='AG'"));
    assert!(stdout.contains("raw_counts={"));
    assert!(stdout.contains("decision="));
    assert!(stdout.contains("evidence=["));
}

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
