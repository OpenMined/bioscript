use std::{path::PathBuf, process::Command};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn fixture_path(name: &str) -> PathBuf {
    repo_root()
        .join("rust/bioscript/tests/fixtures/resources")
        .join(name)
}

fn run_fixture(name: &str, extra_args: &[&str]) -> std::process::Output {
    let root = repo_root();
    let script = fixture_path(name);
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_bioscript"));
    cmd.current_dir(&root);
    for arg in extra_args {
        cmd.arg(arg);
    }
    cmd.arg(script);
    cmd.output().unwrap()
}

#[test]
fn infinite_loop_times_out() {
    let output = run_fixture("infinite_loop.py", &["--max-duration-ms", "10"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("time limit exceeded"));
}

#[test]
fn recursion_limit_trips() {
    let output = run_fixture("deep_recursion.py", &["--max-recursion-depth", "32"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("maximum recursion depth exceeded"));
}

#[test]
fn large_allocation_fails() {
    let output = run_fixture("large_allocation.py", &["--max-memory-bytes", "65536"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("memory limit exceeded"));
}

#[test]
fn giant_string_amplification_fails() {
    let output = run_fixture("giant_string_amplification.py", &["--max-memory-bytes", "65536"]);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("memory limit exceeded"));
}

#[test]
fn giant_list_growth_fails() {
    let output = run_fixture(
        "giant_list_growth.py",
        &["--max-memory-bytes", "65536", "--max-duration-ms", "50"],
    );
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("memory limit exceeded") || stderr.contains("time limit exceeded"));
}
