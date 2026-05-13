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
    let fixture_dir = root.join(format!(
        "target/vntyper-runtime-native-{}",
        std::process::id()
    ));
    fs::create_dir_all(&fixture_dir).unwrap();
    let reference_path = fixture_dir.join("ref.fa");
    let fastq_1_path = fixture_dir.join("r1.fastq");
    let fastq_2_path = fixture_dir.join("r2.fastq");
    let output_dir = fixture_dir.join("out");
    fs::write(&reference_path, ">chr1\nAAAACCCCGGGGTTTT\n").unwrap();
    fs::write(
        &fastq_1_path,
        "@r1\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n@r2\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n@r3\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    fs::write(
        &fastq_2_path,
        "@r4\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n@r5\nAAAATCCCGGGGTTTT\n+\nIIIIIIIIIIIIIIII\n",
    )
    .unwrap();
    let output_arg = output_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let reference_arg = reference_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let fastq_1_arg = fastq_1_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let fastq_2_arg = fastq_2_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let output_dir_arg = output_dir
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
                ("fastq_1", MontyObject::String(fastq_1_arg)),
                ("fastq_2", MontyObject::String(fastq_2_arg)),
                ("reference_fasta", MontyObject::String(reference_arg)),
                ("output_dir", MontyObject::String(output_dir_arg)),
                ("output_file", MontyObject::String(output_arg)),
                ("participant_id", MontyObject::String("positive".to_owned())),
            ],
        )
        .unwrap();

    let plan = fs::read_to_string(&output_path).unwrap();
    assert!(plan.contains("fastq_1"));
    assert!(plan.contains("kestrel_vcf"));
    assert!(plan.contains("first_variant_alt"));
    assert!(plan.contains("first_variant_confidence"));
    assert!(plan.contains("Low_Precision"));
    assert!(plan.contains("\tT"));
    assert!(output_dir.join("positive/kestrel/output.vcf").exists());
    assert!(
        output_dir
            .join("positive/kestrel/output.sorted.vcf.gz")
            .exists()
    );
    let kestrel_tsv = output_dir.join("positive/kestrel_result.tsv");
    assert!(kestrel_tsv.exists());
    let tsv = fs::read_to_string(&kestrel_tsv).unwrap();
    assert!(tsv.contains("Confidence"));
    assert!(tsv.contains("passes_vntyper_filters"));
    let report_json = output_dir.join("positive/report.json");
    assert!(report_json.exists());
    let report = fs::read_to_string(&report_json).unwrap();
    assert!(report.contains("\"algorithm_results\""));
    assert!(report.contains("\"kestrel\""));
    assert!(report.contains("\"Low_Precision\""));
    fs::remove_file(output_path).unwrap();
    fs::remove_dir_all(fixture_dir).unwrap();
}
