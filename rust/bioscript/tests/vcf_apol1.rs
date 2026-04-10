use std::{fs, path::PathBuf, process::Command};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

#[test]
fn apol1_script_accepts_single_sample_vcf() {
    let root = repo_root();
    let output_rel = "bioscripts/output/apol1-from-vcf.tsv";
    let trace_rel = "bioscripts/output/apol1-from-vcf.trace.tsv";
    let output_path = root.join(output_rel);
    let trace_path = root.join(trace_rel);
    let _ = fs::remove_file(&output_path);
    let _ = fs::remove_file(&trace_path);

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("rust/bioscript/tests/fixtures/vcf/apol1_sample.vcf")
        .arg("--output-file")
        .arg(output_rel)
        .arg("--participant-id")
        .arg("VCF001")
        .arg("--trace-report")
        .arg(trace_rel)
        .arg("bioscripts/apol1.py")
        .output()
        .unwrap();

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("G2/G1"));

    let report = fs::read_to_string(&output_path).unwrap();
    assert!(report.contains("VCF001"));
    assert!(report.contains("G2/G1"));

    let trace = fs::read_to_string(&trace_path).unwrap();
    assert!(trace.contains("classify_apol1(genotypes)"));
}
