use std::{
    fs,
    path::PathBuf,
    sync::atomic::{AtomicUsize, Ordering},
    time::{Duration, SystemTime, UNIX_EPOCH},
};

use bioscript_formats::GenotypeLoadOptions;
use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use monty::{MontyObject, ResourceLimits};

static TEMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let counter = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!(
        "bioscript-runtime-{label}-{}-{nanos}-{counter}",
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

#[cfg(unix)]
#[test]
fn host_read_text_rejects_symlink_escape() {
    let dir = temp_dir("read-symlink-escape");
    let outside = temp_dir("read-symlink-outside");
    fs::write(outside.join("secret.txt"), "secret").unwrap();
    std::os::unix::fs::symlink(outside.join("secret.txt"), dir.join("linked-secret.txt")).unwrap();

    let Err(err) = run_script_with_inputs(
        &dir,
        r#"
def main():
    bioscript.read_text("linked-secret.txt")

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    ) else {
        panic!("expected symlink read to fail");
    };
    assert!(err.contains("path escapes bioscript root"), "{err}");
}

#[cfg(unix)]
#[test]
fn host_write_text_rejects_symlink_escape() {
    let dir = temp_dir("write-symlink-escape");
    let outside = temp_dir("write-symlink-outside");
    fs::write(outside.join("target.txt"), "before").unwrap();
    std::os::unix::fs::symlink(outside.join("target.txt"), dir.join("linked-target.txt")).unwrap();

    let Err(err) = run_script_with_inputs(
        &dir,
        r#"
def main():
    bioscript.write_text("linked-target.txt", "after")

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    ) else {
        panic!("expected symlink write to fail");
    };
    assert!(err.contains("path escapes bioscript root"), "{err}");
    assert_eq!(
        fs::read_to_string(outside.join("target.txt")).unwrap(),
        "before"
    );
}

#[test]
fn host_read_text_rejects_oversized_file() {
    let dir = temp_dir("oversized-read");
    fs::write(dir.join("large.txt"), vec![b'a'; 16 * 1024 * 1024 + 1]).unwrap();

    let Err(err) = run_script_with_inputs(
        &dir,
        r#"
def main():
    bioscript.read_text("large.txt")

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    ) else {
        panic!("expected oversized read to fail");
    };
    assert!(err.contains("exceeds 16777216 bytes"), "{err}");
}

#[test]
fn host_read_text_rejects_invalid_utf8() {
    let dir = temp_dir("invalid-utf8-read");
    fs::write(dir.join("invalid.txt"), [0xff, 0xfe, b'a']).unwrap();

    let Err(err) = run_script_with_inputs(
        &dir,
        r#"
def main():
    bioscript.read_text("invalid.txt")

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    ) else {
        panic!("expected invalid UTF-8 read to fail");
    };
    assert!(err.contains("failed to decode"), "{err}");
}

#[test]
fn host_write_text_rejects_oversized_content() {
    let dir = temp_dir("oversized-write");
    let script = dir.join("script.py");
    fs::write(
        &script,
        r#"
def main():
    bioscript.write_text("large.txt", content)

if __name__ == "__main__":
    main()
"#,
    )
    .unwrap();
    let huge = "a".repeat(16 * 1024 * 1024 + 1);
    let runtime = BioscriptRuntime::with_config(
        &dir,
        RuntimeConfig {
            limits: ResourceLimits::new()
                .max_duration(Duration::from_millis(100))
                .max_memory(64 * 1024 * 1024)
                .max_allocations(200_000)
                .gc_interval(1000)
                .max_recursion_depth(Some(200)),
            ..RuntimeConfig::default()
        },
    )
    .unwrap();

    let Err(err) = runtime.run_file(&script, None, vec![("content", MontyObject::String(huge))])
    else {
        panic!("expected oversized write to fail");
    };
    let err = err.to_string();
    assert!(
        err.contains("write_text content exceeds 16777216 bytes"),
        "{err}"
    );
    assert!(!dir.join("large.txt").exists());
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

#[test]
fn runtime_trace_report_records_rsid_and_coordinate_lookup_metadata() {
    let dir = temp_dir("trace-metadata");
    let script = dir.join("script.py");
    let trace = dir.join("reports/trace.tsv");
    fs::write(
        &script,
        r#"
RSID = bioscript.variant(rsid="rs73885319")
COORD = bioscript.variant(grch38="chr22:36265860-36265860", ref="A", alt="G", kind="snp")

def main():
    text = str(RSID) + str(COORD)

if __name__ == "__main__":
    main()
"#,
    )
    .unwrap();

    let runtime = BioscriptRuntime::new(&dir).unwrap();
    runtime.run_file(&script, Some(&trace), Vec::new()).unwrap();

    let report = fs::read_to_string(trace).unwrap();
    assert!(report.contains("lookup_key\tlookup_url"), "{report}");
    assert!(report.contains("rs73885319"), "{report}");
    assert!(
        report.contains("https://www.ncbi.nlm.nih.gov/snp/rs73885319"),
        "{report}"
    );
    assert!(report.contains("22:36265860-36265860"), "{report}");
    assert!(
        report.contains("https://www.ensembl.org/Homo_sapiens/Location/View"),
        "{report}"
    );
}

#[test]
fn runtime_write_tsv_serializes_rows_and_records_timing() {
    let dir = temp_dir("write-tsv");
    let runtime = run_script_with_inputs(
        &dir,
        r#"
def main():
    bioscript.write_tsv("outputs/table.tsv", [
        {"name": "alpha", "count": 2, "ok": True, "empty": None},
        {"name": "beta", "count": 3, "ok": False},
    ])

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    )
    .unwrap();

    let table = fs::read_to_string(dir.join("outputs/table.tsv")).unwrap();
    assert!(table.contains("count\tempty\tname\tok"), "{table}");
    assert!(table.contains("2\t\talpha\ttrue"), "{table}");
    assert!(table.contains("3\t\tbeta\tfalse"), "{table}");
    assert!(
        runtime
            .timing_snapshot()
            .iter()
            .any(|timing| timing.stage == "write_tsv")
    );
}

#[test]
fn runtime_reports_host_method_argument_errors() {
    for (code, expected) in [
        (
            "bioscript.variant(rsid=123)\n",
            "expected string or list of strings",
        ),
        (
            "bioscript.variant(foo='bar')\n",
            "bioscript.variant does not accept keyword 'foo'",
        ),
        (
            "bioscript.variant(grch38='not-a-locus')\n",
            "invalid locus string",
        ),
        (
            "bioscript.variant(kind='structural')\n",
            "invalid variant kind",
        ),
        (
            "bioscript.query_plan('not a plan')\n",
            "expected a list of Variant objects or a VariantPlan",
        ),
        (
            "bioscript.query_plan()\n",
            "bioscript.query_plan expects self and a list of variants",
        ),
        (
            "bioscript.query_plan(variants=[])\n",
            "bioscript.query_plan does not accept keyword arguments",
        ),
        (
            "bioscript.variant('rs1')\n",
            "bioscript.variant expects only self as a positional argument",
        ),
        (
            "bioscript.load_genotypes()\n",
            "bioscript.load_genotypes expects self and path",
        ),
        (
            "bioscript.load_genotypes(path='genotypes.txt')\n",
            "bioscript.load_genotypes does not accept keyword arguments",
        ),
        (
            "bioscript.missing_method()\n",
            "'Bioscript' object has no attribute 'missing_method'",
        ),
        (
            "bioscript.write_tsv('outputs/table.tsv', 'not rows')\n",
            "write_tsv expects a list of dict rows",
        ),
        (
            "bioscript.read_text(path='inputs/source.txt')\n",
            "read_text does not accept keyword arguments",
        ),
    ] {
        let err = run_script(code).unwrap_err();
        assert!(err.contains(expected), "{err}");
    }
}

#[test]
fn runtime_reports_genotype_method_argument_errors() {
    let dir = temp_dir("genotype-method-errors");
    fs::write(
        dir.join("genotypes.txt"),
        "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\n",
    )
    .unwrap();

    for (call, expected) in [
        ("genotypes.get()", "GenotypeFile.get expects self and rsid"),
        (
            "genotypes.get(123)",
            "GenotypeFile.get expected str at position 1",
        ),
        (
            "genotypes.get(rsid='rs1')",
            "GenotypeFile.get does not accept keyword arguments",
        ),
        (
            "genotypes.lookup_variant()",
            "GenotypeFile.lookup_variant expects self and variant",
        ),
        (
            "genotypes.lookup_variant_details()",
            "GenotypeFile.lookup_variant_details expects self and variant",
        ),
        (
            "genotypes.lookup_variants()",
            "GenotypeFile.lookup_variants expects self and a variant plan",
        ),
        (
            "genotypes.lookup_variants_details()",
            "GenotypeFile.lookup_variants_details expects self and a variant plan",
        ),
        (
            "genotypes.lookup_variants_details(plan=[])",
            "GenotypeFile.lookup_variants_details does not accept keyword arguments",
        ),
        (
            "genotypes.missing_method()",
            "'GenotypeFile' object has no attribute 'missing_method'",
        ),
    ] {
        let code = format!(
            r#"
def main():
    genotypes = bioscript.load_genotypes(input_file)
    {call}

if __name__ == "__main__":
    main()
"#
        );
        let Err(err) = run_script_with_inputs(
            &dir,
            &code,
            vec![(
                "input_file",
                MontyObject::String("genotypes.txt".to_owned()),
            )],
        ) else {
            panic!("expected genotype method call to fail: {call}");
        };
        assert!(err.contains(expected), "{err}");
    }
}

#[test]
fn runtime_batch_lookup_methods_return_values_details_and_timings() {
    let dir = temp_dir("batch-lookup-methods");
    fs::write(
        dir.join("genotypes.txt"),
        "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\nrs2\t1\t20\tCT\n",
    )
    .unwrap();

    let runtime = run_script_with_inputs(
        &dir,
        r#"
RS1 = bioscript.variant(rsid="rs1")
RS2 = bioscript.variant(rsid="rs2")
MISSING = bioscript.variant(rsid="rsMissing")

def main():
    genotypes = bioscript.load_genotypes(input_file)
    plan = bioscript.query_plan([RS2, MISSING, RS1])
    values = genotypes.lookup_variants(plan)
    details = genotypes.lookup_variants_details(plan)
    bioscript.write_text("outputs/batch.txt", str(values) + "\n" + str(details))

if __name__ == "__main__":
    main()
"#,
        vec![(
            "input_file",
            MontyObject::String("genotypes.txt".to_owned()),
        )],
    )
    .unwrap();

    let output = fs::read_to_string(dir.join("outputs/batch.txt")).unwrap();
    assert!(output.contains("CT"), "{output}");
    assert!(output.contains("AG"), "{output}");
    assert!(output.contains("VariantObservation"), "{output}");
    assert!(
        output.contains("no matching rsid or locus found"),
        "{output}"
    );
    let timings = runtime.timing_snapshot();
    assert!(
        timings
            .iter()
            .any(|timing| timing.stage == "lookup_variants")
    );
    assert!(
        timings
            .iter()
            .any(|timing| timing.stage == "lookup_variants_details")
    );
}

#[test]
fn runtime_single_lookup_methods_return_values_and_none() {
    let dir = temp_dir("single-lookup-methods");
    fs::write(
        dir.join("genotypes.txt"),
        "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\n",
    )
    .unwrap();

    let runtime = run_script_with_inputs(
        &dir,
        r#"
RS1 = bioscript.variant(rsid="rs1")
MISSING = bioscript.variant(rsid="rsMissing")

def main():
    genotypes = bioscript.load_genotypes(input_file)
    values = [
        genotypes.get("rs1"),
        genotypes.get("rsMissing"),
        genotypes.lookup_variant(RS1),
        genotypes.lookup_variant(MISSING),
    ]
    bioscript.write_text("outputs/single.txt", str(values))

if __name__ == "__main__":
    main()
"#,
        vec![(
            "input_file",
            MontyObject::String("genotypes.txt".to_owned()),
        )],
    )
    .unwrap();

    let output = fs::read_to_string(dir.join("outputs/single.txt")).unwrap();
    assert!(output.contains("AG"), "{output}");
    assert!(output.contains("None"), "{output}");
    assert!(
        runtime
            .timing_snapshot()
            .iter()
            .any(|timing| timing.stage == "lookup_variant")
    );
}

#[test]
fn runtime_variant_objects_preserve_optional_fields() {
    let dir = temp_dir("variant-optional-fields");
    run_script_with_inputs(
        &dir,
        r#"
def main():
    insertion = bioscript.variant(
        rsids=["rs1", "rs2"],
        grch37="chr1:10-11",
        grch38="2:20",
        ref="A",
        alt="AT",
        kind="insertion",
        motifs=["AT", "TA"],
    )
    deletion = bioscript.variant(
        grch38="3:30-32",
        kind="deletion",
        deletion_length=3,
    )
    indel = bioscript.variant(kind="indel")
    other = bioscript.variant(kind="other")
    plan = bioscript.query_plan([insertion, deletion, indel, other])
    bioscript.write_text("outputs/variants.txt", str(plan))

if __name__ == "__main__":
    main()
"#,
        Vec::new(),
    )
    .unwrap();

    let output = fs::read_to_string(dir.join("outputs/variants.txt")).unwrap();
    assert!(output.contains("grch37='1:10-11'"), "{output}");
    assert!(output.contains("grch38='2:20-20'"), "{output}");
    assert!(output.contains("reference='A'"), "{output}");
    assert!(output.contains("alternate='AT'"), "{output}");
    assert!(output.contains("kind='insertion'"), "{output}");
    assert!(output.contains("deletion_length=3"), "{output}");
    assert!(output.contains("motifs=['AT', 'TA']"), "{output}");
    assert!(output.contains("kind='indel'"), "{output}");
    assert!(output.contains("kind='other'"), "{output}");
}

#[test]
fn runtime_direct_run_script_and_accessors_are_usable() {
    let dir = temp_dir("direct-run-script");
    let runtime = BioscriptRuntime::with_config(
        &dir,
        RuntimeConfig {
            loader: GenotypeLoadOptions::default(),
            ..RuntimeConfig::default()
        },
    )
    .unwrap();

    assert_eq!(runtime.root(), dir.canonicalize().unwrap().as_path());
    assert!(runtime.config().loader.input_index.is_none());
    let result = runtime
        .run_script("result = 2 + 3\nresult\n", "inline.py", Vec::new())
        .unwrap();
    assert!(matches!(result, MontyObject::Int(5)));
}

#[test]
fn runtime_reports_filesystem_setup_errors() {
    let dir = temp_dir("filesystem-errors");
    let missing_root = dir.join("missing-root");
    let Err(err) = BioscriptRuntime::new(&missing_root) else {
        panic!("expected missing root to fail");
    };
    assert!(
        err.to_string()
            .contains("failed to canonicalize bioscript root"),
        "{err}"
    );

    let runtime = BioscriptRuntime::new(&dir).unwrap();
    let err = runtime
        .run_file(dir.join("missing-script.py"), None, Vec::new())
        .unwrap_err();
    assert!(err.to_string().contains("failed to read script"), "{err}");
}

#[test]
fn runtime_loader_paths_are_resolved_and_escape_checks_apply() {
    let dir = temp_dir("loader-paths");
    fs::write(
        dir.join("genotypes.txt"),
        "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\n",
    )
    .unwrap();

    let runtime = BioscriptRuntime::with_config(
        &dir,
        RuntimeConfig {
            loader: GenotypeLoadOptions {
                format: Some(bioscript_formats::GenotypeSourceFormat::Cram),
                input_index: Some(PathBuf::from("indexes/input.crai")),
                reference_file: Some(PathBuf::from("refs/ref.fa")),
                reference_index: Some(PathBuf::from("refs/ref.fa.fai")),
                ..GenotypeLoadOptions::default()
            },
            ..RuntimeConfig::default()
        },
    )
    .unwrap();
    let err = runtime
        .run_file(
            {
                let script = dir.join("script.py");
                fs::write(
                    &script,
                    r#"
SNP = bioscript.variant(grch38="1:10-10", ref="A", alt="G", kind="snp")

def main():
    genotypes = bioscript.load_genotypes("genotypes.txt")
    genotypes.lookup_variant(SNP)

if __name__ == "__main__":
    main()
"#,
                )
                .unwrap();
                script
            },
            None,
            Vec::new(),
        )
        .unwrap_err();
    assert!(
        err.to_string().contains("failed to open indexed FASTA"),
        "{err}"
    );

    let err = run_script("bioscript.read_text('/tmp/outside.txt')\n").unwrap_err();
    assert!(err.contains("absolute paths are not allowed"), "{err}");
    let err = run_script("bioscript.write_text('../outside.txt', 'x')\n").unwrap_err();
    assert!(err.contains("path escapes bioscript root"), "{err}");
}

#[test]
fn runtime_write_tsv_reports_row_shape_errors() {
    for (code, expected) in [
        (
            "bioscript.write_tsv('outputs/table.tsv', ['not a dict'])\n",
            "write_tsv row must be a dict",
        ),
        (
            "bioscript.write_tsv('outputs/table.tsv', [{1: 'bad key'}])\n",
            "write_tsv dict keys must be strings",
        ),
        (
            "bioscript.write_tsv('outputs/table.tsv')\n",
            "bioscript.write_tsv expects self, path, rows",
        ),
    ] {
        let err = run_script(code).unwrap_err();
        assert!(err.contains(expected), "{err}");
    }
}

#[test]
fn runtime_variant_argument_conversions_cover_optional_and_list_errors() {
    for (code, expected) in [
        (
            "bioscript.variant(rsids=['rs1', 2])\n",
            "expected list of strings",
        ),
        (
            "bioscript.variant(grch38=123)\n",
            "expected optional string",
        ),
        (
            "bioscript.variant(deletion_length='long')\n",
            "expected optional int",
        ),
        (
            "bioscript.variant(motifs=123)\n",
            "expected string or list of strings",
        ),
        (
            "bioscript.variant(grch38='1:not-a-start-20')\n",
            "invalid locus start",
        ),
        (
            "bioscript.variant(grch38='1:10-not-an-end')\n",
            "invalid locus end",
        ),
    ] {
        let err = run_script(code).unwrap_err();
        assert!(err.contains(expected), "{err}");
    }
}
