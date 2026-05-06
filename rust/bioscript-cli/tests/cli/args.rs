use super::*;

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
        (
            vec!["inspect", "--input-index"],
            "--input-index requires a path",
        ),
        (
            vec!["inspect", "--reference-file"],
            "--reference-file requires a path",
        ),
        (
            vec!["inspect", "--reference-index"],
            "--reference-index requires a path",
        ),
        (vec!["prepare", "--root"], "--root requires a directory"),
        (
            vec!["prepare", "--input-file"],
            "--input-file requires a path",
        ),
        (
            vec!["prepare", "--reference-file"],
            "--reference-file requires a path",
        ),
        (vec!["prepare", "extra"], "unexpected argument: extra"),
        (
            vec!["prepare", "--cache-dir"],
            "--cache-dir requires a path",
        ),
        (
            vec!["validate-variants", "one.yaml", "two.yaml"],
            "unexpected argument: two.yaml",
        ),
        (
            vec!["validate-variants", "--report"],
            "--report requires a path",
        ),
        (
            vec!["validate-panels", "one.yaml", "two.yaml"],
            "unexpected argument: two.yaml",
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
fn cli_accepts_auto_format_and_explicit_loader_paths_for_script_runs() {
    let root = repo_root();
    let dir = temp_dir("loader-args");
    fs::write(
        dir.join("script.py"),
        r#"
def main():
    print("loader args accepted")

if __name__ == "__main__":
    main()
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--root")
        .arg(&dir)
        .arg("--input-format")
        .arg("auto")
        .arg("--input-index")
        .arg("input.crai")
        .arg("--reference-file")
        .arg("ref.fa")
        .arg("--reference-index")
        .arg("ref.fa.fai")
        .arg("--allow-md5-mismatch")
        .arg(dir.join("script.py"))
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(
        String::from_utf8_lossy(&output.stdout).contains("loader args accepted"),
        "stdout: {}",
        String::from_utf8_lossy(&output.stdout)
    );
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
