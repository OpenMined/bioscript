use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_formats::GenotypeLoadOptions;
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use monty::MontyObject;

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

fn run_script_with_inputs(
    root: &PathBuf,
    code: &str,
    inputs: Vec<(&str, MontyObject)>,
) -> Result<BioscriptRuntime, String> {
    let script = root.join("script.py");
    fs::write(&script, code).unwrap();

    let runtime = BioscriptRuntime::with_config(
        root,
        RuntimeConfig {
            loader: GenotypeLoadOptions::default(),
            ..RuntimeConfig::default()
        },
    )
    .unwrap();

    runtime
        .run_file(&script, None, inputs)
        .map(|_| runtime.clone())
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
    assert!(err.contains("eval"), "{err}");
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

#[test]
fn host_read_write_text_allows_nested_safe_paths() {
    let dir = temp_dir("nested-safe-paths");
    fs::create_dir_all(dir.join("inputs/nested")).unwrap();
    fs::write(dir.join("inputs/nested/source.txt"), "hello nested").unwrap();

    run_script_with_inputs(
        &dir,
        r#"
def main():
    text = bioscript.read_text("inputs/nested/source.txt")
    bioscript.write_text("outputs/nested/result.txt", text + " output")

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    )
    .unwrap();

    let written = fs::read_to_string(dir.join("outputs/nested/result.txt")).unwrap();
    assert_eq!(written, "hello nested output");
}

#[test]
fn runtime_lookup_details_reports_missing_variant_and_no_call() {
    let dir = temp_dir("lookup-missing-no-call");
    fs::write(
        dir.join("genotypes.txt"),
        "rsid\tchromosome\tposition\tgenotype\nrsNoCall\t1\t100\tN/A\n",
    )
    .unwrap();

    let runtime = run_script_with_inputs(
        &dir,
        r#"
MISSING = bioscript.variant(rsid="rsMissing")
NO_CALL = bioscript.variant(rsid="rsNoCall")

def main():
    genotypes = bioscript.load_genotypes(input_file)
    missing = genotypes.lookup_variant_details(MISSING)
    no_call = genotypes.lookup_variant_details(NO_CALL)
    bioscript.write_text("outputs/details.txt", str(missing) + "\n" + str(no_call))

if __name__ == "__main__":
    main()
"#,
        vec![(
            "input_file",
            MontyObject::String("genotypes.txt".to_owned()),
        )],
    )
    .unwrap();

    let timings = runtime.timing_snapshot();
    assert!(
        timings
            .iter()
            .any(|timing| timing.stage == "lookup_variant_details")
    );
    let details = fs::read_to_string(dir.join("outputs/details.txt")).unwrap();
    assert!(details.contains("genotype=None"), "{details}");
    assert!(details.contains("genotype='--'"), "{details}");
}
