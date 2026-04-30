use super::*;

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
