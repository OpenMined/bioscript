use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_runtime::{BioscriptRuntime, RuntimeConfig};
use monty::MontyObject;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn unique_output_path(root: &std::path::Path) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    root.join(format!(
        "target/vntyper-runtime-plan-{}-{nanos}.tsv",
        std::process::id()
    ))
}

#[test]
fn vntyper_bioscript_program_runs_through_runtime() {
    let root = repo_root();
    let output_path = unique_output_path(&root);
    let output_arg = output_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let runtime = BioscriptRuntime::with_config(&root, RuntimeConfig::default()).unwrap();

    runtime
        .run_file(
            root.join("ports/vntyper/bioscript/vntyper.bs"),
            None,
            vec![
                (
                    "input_file",
                    MontyObject::String(
                        "ports/vntyper/test-data/example_6449_hg19_subset.bam".to_owned(),
                    ),
                ),
                ("output_file", MontyObject::String(output_arg)),
                ("participant_id", MontyObject::String("positive".to_owned())),
            ],
        )
        .unwrap();

    let plan = fs::read_to_string(&output_path).unwrap();
    assert!(plan.contains("samtools_view_command"));
    assert!(plan.contains("kestrel_command"));
    assert!(plan.contains("bcftools_sort_command"));
    fs::remove_file(output_path).unwrap();
}

#[test]
fn vntyper_fastq_bioscript_program_runs_through_runtime() {
    let root = repo_root();
    let output_path = unique_output_path(&root);
    let output_arg = output_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let runtime = BioscriptRuntime::with_config(&root, RuntimeConfig::default()).unwrap();

    runtime
        .run_file(
            root.join("ports/vntyper/bioscript/vntyper-fastq.bs"),
            None,
            vec![
                (
                    "fastq_1",
                    MontyObject::String(
                        "ports/vntyper/test-data/example_6449_hg19_subset_R1.fastq.gz".to_owned(),
                    ),
                ),
                (
                    "fastq_2",
                    MontyObject::String(
                        "ports/vntyper/test-data/example_6449_hg19_subset_R2.fastq.gz".to_owned(),
                    ),
                ),
                ("output_file", MontyObject::String(output_arg)),
                ("participant_id", MontyObject::String("positive".to_owned())),
            ],
        )
        .unwrap();

    let plan = fs::read_to_string(&output_path).unwrap();
    assert!(plan.contains("fastq_1"));
    assert!(plan.contains("kestrel_command"));
    assert!(plan.contains("bcftools_sort_command"));
    fs::remove_file(output_path).unwrap();
}
