use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_formats::GenotypeLoadOptions;
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-runtime-{label}-{}-{nanos}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn run_script(code: &str) -> Result<(), String> {
    let dir = temp_dir("security");
    let script = dir.join("script.py");
    fs::write(&script, code).unwrap();

    let runtime = BioscriptRuntime::with_config(
        &dir,
        RuntimeConfig {
            loader: GenotypeLoadOptions::default(),
            ..RuntimeConfig::default()
        },
    )
    .unwrap();

    runtime
        .run_file(&script, None, Vec::new())
        .map(|_| ())
        .map_err(|err| err.to_string())
}

#[test]
fn open_builtin_is_not_available() {
    let err = run_script("open('secret.txt')\n").unwrap_err();
    assert!(err.contains("unknown bioscript host function: open"));
}

#[test]
fn eval_builtin_is_not_available() {
    let err = run_script("eval('1 + 1')\n").unwrap_err();
    assert!(err.contains("unknown bioscript host function: eval"));
}

#[test]
fn exec_builtin_is_not_available() {
    let err = run_script("exec('x = 1')\n").unwrap_err();
    assert!(err.contains("unknown bioscript host function: exec"));
}

#[test]
fn dunder_import_is_not_available() {
    let err = run_script("__import__('os')\n").unwrap_err();
    assert!(err.contains("unknown bioscript host function: __import__"));
}

#[test]
fn os_getenv_is_blocked() {
    let err = run_script("import os\nprint(os.getenv('HOME'))\n").unwrap_err();
    assert!(err.contains("OS call os.getenv is blocked"));
}

#[test]
fn pathlib_read_text_is_blocked() {
    let err = run_script("from pathlib import Path\nprint(Path('sample.txt').read_text())\n")
        .unwrap_err();
    assert!(err.contains("OS call Path.read_text is blocked"));
}

#[test]
fn unsupported_networkish_import_fails() {
    let err = run_script("import urllib\n").unwrap_err();
    assert!(err.contains("No module named 'urllib'"));
}
