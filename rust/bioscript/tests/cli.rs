use std::{fs, path::PathBuf, process::Command};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
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

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
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
    let script = root.join("rust/bioscript/tests/fixtures/path_escape.py");
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

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
    let trace = fs::read_to_string(trace_path).unwrap();
    assert!(trace.contains("step\tline\tcode"));
    assert!(trace.contains("hello from bioscript"));
}
