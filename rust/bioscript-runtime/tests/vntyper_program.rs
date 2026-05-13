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
    let fixture_dir = root.join(format!(
        "target/vntyper-runtime-main-{}",
        std::process::id()
    ));
    fs::create_dir_all(&fixture_dir).unwrap();
    let bam_source = root.join("vendor/rust/samtools-rs/samtools/test/stat/11_target.bam");
    let bai_source = root.join("vendor/rust/samtools-rs/samtools/test/stat/11_target.bam.bai");
    let bam_path = fixture_dir.join("input.bam");
    let bai_path = fixture_dir.join("input.bam.bai");
    let reference_path = fixture_dir.join("ref.fa");
    let output_dir = fixture_dir.join("out");
    fs::create_dir_all(&output_dir).unwrap();
    fs::copy(bam_source, &bam_path).unwrap();
    fs::copy(bai_source, &bai_path).unwrap();
    fs::write(&reference_path, ">ref1\nAAAACCCCGGGGTTTT\n").unwrap();
    let output_arg = output_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let bam_arg = bam_path.strip_prefix(&root).unwrap().display().to_string();
    let bai_arg = bai_path.strip_prefix(&root).unwrap().display().to_string();
    let reference_arg = reference_path
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
            root.join("ports/vntyper/bioscript/vntyper.bs"),
            None,
            vec![
                ("input_file", MontyObject::String(bam_arg)),
                ("input_bai", MontyObject::String(bai_arg)),
                ("bam_region", MontyObject::String("ref1:1-10".to_owned())),
                ("vntr_region", MontyObject::String("ref1:1-10".to_owned())),
                ("reference_fasta", MontyObject::String(reference_arg)),
                ("kmer_size", MontyObject::Int(4)),
                ("minimum_difference", MontyObject::Int(1)),
                ("max_haplotypes", MontyObject::Int(4)),
                ("max_saved_states", MontyObject::Int(4)),
                ("output_dir", MontyObject::String(output_dir_arg)),
                ("output_file", MontyObject::String(output_arg)),
                ("participant_id", MontyObject::String("main-bam".to_owned())),
            ],
        )
        .unwrap();

    let summary = fs::read_to_string(&output_path).unwrap();
    assert!(summary.contains("sliced_bam"));
    assert!(summary.contains("fastq_read1_records"));
    assert!(summary.contains("report_json"));
    assert!(output_dir.join("main-bam_kestrel_result.tsv").exists());
    assert!(output_dir.join("main-bam_report.json").exists());
    fs::remove_file(output_path).unwrap();
    fs::remove_dir_all(fixture_dir).unwrap();
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
                ("kmer_size", MontyObject::Int(4)),
                ("minimum_difference", MontyObject::Int(1)),
                ("max_haplotypes", MontyObject::Int(4)),
                ("max_saved_states", MontyObject::Int(4)),
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
    assert!(report.contains("\"native bioscript kestrel from FASTQ\""));
    fs::remove_file(output_path).unwrap();
    fs::remove_dir_all(fixture_dir).unwrap();
}

#[test]
fn vntyper_bam_native_bioscript_program_runs_through_runtime() {
    let root = repo_root();
    let output_path = unique_output_path(&root);
    let fixture_dir = root.join(format!("target/vntyper-runtime-bam-{}", std::process::id()));
    fs::create_dir_all(&fixture_dir).unwrap();
    let bam_source = root.join("vendor/rust/samtools-rs/samtools/test/stat/11_target.bam");
    let bai_source = root.join("vendor/rust/samtools-rs/samtools/test/stat/11_target.bam.bai");
    let bam_path = fixture_dir.join("input.bam");
    let bai_path = fixture_dir.join("input.bam.bai");
    let reference_path = fixture_dir.join("ref.fa");
    let output_dir = fixture_dir.join("out");
    fs::create_dir_all(&output_dir).unwrap();
    fs::copy(bam_source, &bam_path).unwrap();
    fs::copy(bai_source, &bai_path).unwrap();
    fs::write(&reference_path, ">ref1\nAAAACCCCGGGGTTTT\n").unwrap();
    let output_arg = output_path
        .strip_prefix(&root)
        .unwrap()
        .display()
        .to_string();
    let bam_arg = bam_path.strip_prefix(&root).unwrap().display().to_string();
    let bai_arg = bai_path.strip_prefix(&root).unwrap().display().to_string();
    let reference_arg = reference_path
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
            root.join("ports/vntyper/bioscript/vntyper-bam-native.bs"),
            None,
            vec![
                ("input_file", MontyObject::String(bam_arg)),
                ("input_bai", MontyObject::String(bai_arg)),
                ("bam_region", MontyObject::String("ref1:1-10".to_owned())),
                ("vntr_region", MontyObject::String("ref1:1-10".to_owned())),
                ("reference_fasta", MontyObject::String(reference_arg)),
                ("kmer_size", MontyObject::Int(4)),
                ("minimum_difference", MontyObject::Int(1)),
                ("max_haplotypes", MontyObject::Int(4)),
                ("max_saved_states", MontyObject::Int(4)),
                ("output_dir", MontyObject::String(output_dir_arg)),
                ("output_file", MontyObject::String(output_arg)),
                ("participant_id", MontyObject::String("tiny-bam".to_owned())),
            ],
        )
        .unwrap();

    let summary = fs::read_to_string(&output_path).unwrap();
    assert!(summary.contains("sliced_bam"));
    assert!(summary.contains("fastq_read1_records"));
    assert!(summary.contains("depth_region_length"));
    assert!(output_dir.join("tiny-bam_sliced.bam").exists());
    assert!(output_dir.join("tiny-bam_R1.fastq.gz").exists());
    assert!(output_dir.join("tiny-bam_R2.fastq.gz").exists());
    assert!(output_dir.join("tiny-bam_kestrel.vcf").exists());
    assert!(output_dir.join("tiny-bam_kestrel.sorted.vcf.gz").exists());
    assert!(output_dir.join("tiny-bam_kestrel_result.tsv").exists());
    let report_json = output_dir.join("tiny-bam_report.json");
    assert!(report_json.exists());
    let report = fs::read_to_string(&report_json).unwrap();
    assert!(report.contains("\"native bioscript samtools/kestrel\""));
    assert!(report.contains("\"region_length\""));
    fs::remove_file(output_path).unwrap();
    fs::remove_dir_all(fixture_dir).unwrap();
}
