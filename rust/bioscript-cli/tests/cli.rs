use std::{
    ffi::OsStr,
    fs,
    path::PathBuf,
    process::{Command, Output},
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

fn run_bioscript<I, S>(root: &PathBuf, args: I) -> Output
where
    I: IntoIterator<Item = S>,
    S: AsRef<OsStr>,
{
    Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(root)
        .args(args)
        .output()
        .unwrap()
}

fn stderr_text(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn cli_reports_usage_when_no_script_or_subcommand_is_provided() {
    let root = repo_root();

    let output = run_bioscript(&root, std::iter::empty::<&str>());

    assert!(!output.status.success());
    let stderr = stderr_text(&output);
    assert!(stderr.contains("usage: bioscript"), "{stderr}");
    assert!(stderr.contains("validate-variants"), "{stderr}");
    assert!(stderr.contains("inspect <path>"), "{stderr}");
}

#[test]
fn cli_rejects_missing_values_and_unexpected_arguments() {
    let root = repo_root();

    for (args, expected) in [
        (vec!["--root"], "--root requires a directory"),
        (vec!["--input-file"], "--input-file requires a path"),
        (vec!["--output-file"], "--output-file requires a path"),
        (
            vec!["--participant-id"],
            "--participant-id requires a value",
        ),
        (vec!["--trace-report"], "--trace-report requires a path"),
        (vec!["--timing-report"], "--timing-report requires a path"),
        (vec!["--filter"], "--filter requires key=value"),
        (vec!["--input-index"], "--input-index requires a path"),
        (vec!["--reference-file"], "--reference-file requires a path"),
        (
            vec!["--reference-index"],
            "--reference-index requires a path",
        ),
        (vec!["--cache-dir"], "--cache-dir requires a path"),
        (
            vec!["bioscripts/hello-world.py", "extra"],
            "unexpected argument: extra",
        ),
        (vec!["inspect"], "usage: bioscript inspect"),
        (
            vec!["inspect", "bioscripts/hello-world.py", "extra"],
            "unexpected argument: extra",
        ),
        (vec!["prepare", "--root"], "--root requires a directory"),
        (
            vec!["prepare", "--cache-dir"],
            "--cache-dir requires a path",
        ),
        (
            vec!["validate-variants", "--report"],
            "--report requires a path",
        ),
        (
            vec!["validate-panels", "--report"],
            "--report requires a path",
        ),
    ] {
        let output = run_bioscript(&root, args);
        assert!(!output.status.success(), "expected failure for {expected}");
        let stderr = stderr_text(&output);
        assert!(stderr.contains(expected), "{stderr}");
    }
}

#[test]
fn cli_rejects_invalid_numeric_limits_and_input_formats() {
    let root = repo_root();

    for (args, expected) in [
        (
            vec!["--input-format", "bam", "bioscripts/hello-world.py"],
            "invalid --input-format value bam",
        ),
        (
            vec!["--max-duration-ms", "soon", "bioscripts/hello-world.py"],
            "invalid --max-duration-ms value soon",
        ),
        (
            vec!["--max-memory-bytes", "large", "bioscripts/hello-world.py"],
            "invalid --max-memory-bytes value large",
        ),
        (
            vec!["--max-allocations", "many", "bioscripts/hello-world.py"],
            "invalid --max-allocations value many",
        ),
        (
            vec!["--max-recursion-depth", "deep", "bioscripts/hello-world.py"],
            "invalid --max-recursion-depth value deep",
        ),
        (
            vec!["prepare", "--input-format", "bam"],
            "invalid --input-format: unsupported input format: bam",
        ),
    ] {
        let output = run_bioscript(&root, args);
        assert!(!output.status.success(), "expected failure for {expected}");
        let stderr = stderr_text(&output);
        assert!(stderr.contains(expected), "{stderr}");
    }
}

#[test]
fn cli_rejects_unsupported_manifest_schema() {
    let root = repo_root();
    let dir = temp_dir("unsupported-manifest");
    let manifest = dir.join("unsupported.yaml");
    fs::write(
        &manifest,
        r#"
schema: "bioscript:catalogue:1.0"
version: "1.0"
name: "catalogue"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg(&manifest)
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = stderr_text(&output);
    assert!(
        stderr.contains("unsupported manifest schema 'bioscript:catalogue:1.0'"),
        "{stderr}"
    );
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
fn timing_report_is_written_for_hello_world() {
    let root = repo_root();
    let timing_path = root.join("bioscripts/output/hello-world.timing.tsv");
    if timing_path.exists() {
        fs::remove_file(&timing_path).unwrap();
    }

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--timing-report")
        .arg("bioscripts/output/hello-world.timing.tsv")
        .arg("bioscripts/hello-world.py")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let timing = fs::read_to_string(timing_path).unwrap();
    assert!(timing.contains("stage\tduration_ms\tdetail"));
    assert!(timing.contains("run_file_total\t"));
    assert!(timing.contains("script=bioscripts/hello-world.py"));
}

#[test]
fn auto_index_adds_reference_index_timing_for_script_runs() {
    let root = repo_root();
    let dir = temp_dir("auto-index-script");
    let cache_dir = dir.join("cache");
    let timing_path = dir.join("reports/timing.tsv");
    fs::write(dir.join("ref.fa"), b">chr1\nACGT\n").unwrap();
    fs::write(
        dir.join("script.py"),
        r#"
def main():
    print("indexed")


if __name__ == "__main__":
    main()
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--root")
        .arg(&dir)
        .arg("--reference-file")
        .arg("ref.fa")
        .arg("--auto-index")
        .arg("--cache-dir")
        .arg(&cache_dir)
        .arg("--timing-report")
        .arg(&timing_path)
        .arg(dir.join("script.py"))
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("bioscript: auto-indexed reference ->"),
        "{stderr}"
    );
    assert!(fs::read_dir(&cache_dir).unwrap().any(|entry| {
        entry
            .unwrap()
            .path()
            .extension()
            .is_some_and(|ext| ext == "fai")
    }));
    let timing = fs::read_to_string(timing_path).unwrap();
    assert!(timing.contains("auto_index\t"), "{timing}");
    assert!(timing.contains("run_file_total\t"), "{timing}");
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
