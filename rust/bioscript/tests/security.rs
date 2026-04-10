use std::{path::PathBuf, process::Command};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn run_fixture(name: &str) -> std::process::Output {
    let root = repo_root();
    let script = root.join("rust/bioscript/tests/fixtures/security").join(name);
    Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg(script)
        .output()
        .unwrap()
}

#[test]
fn open_builtin_is_not_available() {
    let output = run_fixture("open_builtin.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("unknown bioscript host function: open"));
}

#[test]
fn eval_builtin_is_not_available() {
    let output = run_fixture("eval_builtin.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("unknown bioscript host function: eval"));
}

#[test]
fn exec_builtin_is_not_available() {
    let output = run_fixture("exec_builtin.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("unknown bioscript host function: exec"));
}

#[test]
fn dunder_import_is_not_available() {
    let output = run_fixture("dunder_import.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("unknown bioscript host function: __import__"));
}

#[test]
fn os_getenv_is_blocked() {
    let output = run_fixture("os_getenv.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("OS call os.getenv is blocked"));
}

#[test]
fn pathlib_read_text_is_blocked() {
    let output = run_fixture("pathlib_read_text.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("OS call Path.read_text is blocked"));
}

#[test]
fn unsupported_networkish_import_fails() {
    let output = run_fixture("import_urllib.py");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("No module named 'urllib'"));
}
