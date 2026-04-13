use std::{
    fs,
    path::PathBuf,
    time::{Duration, SystemTime, UNIX_EPOCH},
};

use bioscript_formats::GenotypeLoadOptions;
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use monty::ResourceLimits;

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

fn run_script(code: &str, limits: ResourceLimits) -> Result<(), String> {
    let dir = temp_dir("resources");
    let script = dir.join("script.py");
    fs::write(&script, code).unwrap();

    let runtime = BioscriptRuntime::with_config(
        &dir,
        RuntimeConfig {
            limits,
            loader: GenotypeLoadOptions::default(),
        },
    )
    .unwrap();

    runtime
        .run_file(&script, None, Vec::new())
        .map(|_| ())
        .map_err(|err| err.to_string())
}

#[test]
fn infinite_loop_times_out() {
    let err = run_script(
        "while True:\n    pass\n",
        ResourceLimits::new().max_duration(Duration::from_millis(10)),
    )
    .unwrap_err();
    assert!(err.contains("time limit exceeded"));
}

#[test]
fn recursion_limit_trips() {
    let err = run_script(
        "def f():\n    return f()\n\nf()\n",
        ResourceLimits::new().max_recursion_depth(Some(32)),
    )
    .unwrap_err();
    assert!(err.contains("maximum recursion depth exceeded"));
}

#[test]
fn large_allocation_fails() {
    let err = run_script(
        "x = 'a' * 1_000_000\nprint(len(x))\n",
        ResourceLimits::new().max_memory(65_536),
    )
    .unwrap_err();
    assert!(err.contains("memory limit exceeded"));
}

#[test]
fn giant_string_amplification_fails() {
    let err = run_script(
        "text = 'a'\nwhile True:\n    text = text + text\n",
        ResourceLimits::new().max_memory(65_536),
    )
    .unwrap_err();
    assert!(err.contains("memory limit exceeded"));
}

#[test]
fn giant_list_growth_fails() {
    let err = run_script(
        "items = []\nwhile True:\n    items.append('x' * 1000)\n",
        ResourceLimits::new()
            .max_memory(65_536)
            .max_duration(Duration::from_millis(50)),
    )
    .unwrap_err();
    assert!(err.contains("memory limit exceeded") || err.contains("time limit exceeded"));
}
