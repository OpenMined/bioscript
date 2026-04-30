use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
    time::{Duration, SystemTime, UNIX_EPOCH},
};

use bioscript_formats::GenotypeLoadOptions;
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use monty::ResourceLimits;

static TEMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!(
        "bioscript-runtime-coverage-{label}-{}-{nanos}-{counter}",
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
fn coverage_infinite_loop_times_out() {
    let err = run_script(
        "while True:\n    pass\n",
        ResourceLimits::new().max_duration(Duration::from_millis(10)),
    )
    .unwrap_err();

    assert!(err.contains("time limit exceeded"), "{err}");
}

#[test]
fn coverage_large_allocation_fails() {
    let err = run_script(
        "x = 'a' * 1_000_000\nprint(len(x))\n",
        ResourceLimits::new().max_memory(65_536),
    )
    .unwrap_err();

    assert!(err.contains("memory limit exceeded"), "{err}");
}

#[test]
fn coverage_giant_string_amplification_fails() {
    let err = run_script(
        "text = 'a'\nwhile True:\n    text = text + text\n",
        ResourceLimits::new().max_memory(65_536),
    )
    .unwrap_err();

    assert!(err.contains("memory limit exceeded"), "{err}");
}

#[test]
fn coverage_giant_list_growth_fails() {
    let err = run_script(
        "items = []\nwhile True:\n    items.append('x' * 1000)\n",
        ResourceLimits::new()
            .max_memory(65_536)
            .max_duration(Duration::from_millis(50)),
    )
    .unwrap_err();

    assert!(
        err.contains("memory limit exceeded") || err.contains("time limit exceeded"),
        "{err}"
    );
}
