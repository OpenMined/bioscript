use std::{fs, path::PathBuf, process::Command};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn run_apol1_case(input_file: &str, participant_id: &str, expected_status: &str) {
    let root = repo_root();
    let output_rel = format!("bioscripts/output/{participant_id}.apol1.tsv");
    let trace_rel = format!("bioscripts/output/{participant_id}.apol1.trace.tsv");
    let output_path = root.join(&output_rel);
    let trace_path = root.join(&trace_rel);
    let _ = fs::remove_file(&output_path);
    let _ = fs::remove_file(&trace_path);

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg(input_file)
        .arg("--output-file")
        .arg(&output_rel)
        .arg("--participant-id")
        .arg(participant_id)
        .arg("--trace-report")
        .arg(&trace_rel)
        .arg("bioscripts/apol1.py")
        .output()
        .unwrap();

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains(expected_status));

    let report = fs::read_to_string(&output_path).unwrap();
    assert!(report.contains("participant_id"));
    assert!(report.contains("apol1_status"));
    assert!(report.contains(participant_id));
    assert!(report.contains(expected_status));

    let trace = fs::read_to_string(&trace_path).unwrap();
    assert!(trace.contains("classify_apol1(genotypes)"));
    assert!(trace.contains("status = classify_apol1(genotypes)"));
    assert!(trace.contains("https://www.ncbi.nlm.nih.gov/snp/rs73885319"));
}

#[test]
fn apol1_script_matches_simple_example_file() {
    run_apol1_case("examples/apol1/test_snps.txt", "P001", "G1/G0");
}

#[test]
fn apol1_script_matches_real_ddna_files() {
    run_apol1_case(
        "examples/apol1/genotype_files/108179_G0G0_X_X_GSAv3-DTC_GRCh38-12-13-2025.txt",
        "108179_G0G0",
        "G0/G0",
    );
    run_apol1_case(
        "examples/apol1/genotype_files/108187_G2G1_X_X_GSAv3-DTC_GRCh38-12-13-2025.txt",
        "108187_G2G1",
        "G2/G1",
    );
    run_apol1_case(
        "examples/apol1/genotype_files/108189_G2G2_X_X_GSAv3-DTC_GRCh38-12-13-2025.txt",
        "108189_G2G2",
        "G2/G2",
    );
}
