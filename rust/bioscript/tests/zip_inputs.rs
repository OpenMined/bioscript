use std::{fs, io::Write, path::PathBuf, process::Command};

use zip::write::SimpleFileOptions;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn create_zip_file(zip_path: &PathBuf, entry_name: &str, source_path: &str) {
    let root = repo_root();
    if let Some(parent) = zip_path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    let source = fs::read(root.join(source_path)).unwrap();
    let file = fs::File::create(zip_path).unwrap();
    let mut writer = zip::ZipWriter::new(file);
    writer
        .start_file(entry_name, SimpleFileOptions::default())
        .unwrap();
    writer.write_all(&source).unwrap();
    writer.finish().unwrap();
}

#[test]
fn apol1_script_accepts_zip_genotype_file_via_autodetect() {
    let root = repo_root();
    let zip_rel = "bioscripts/output/apol1-input.zip";
    let zip_path = root.join(zip_rel);
    create_zip_file(&zip_path, "test_snps.txt", "examples/apol1/test_snps.txt");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg(zip_rel)
        .arg("--output-file")
        .arg("bioscripts/output/apol1-from-zip.tsv")
        .arg("--participant-id")
        .arg("ZIP001")
        .arg("bioscripts/apol1.py")
        .output()
        .unwrap();

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("G1/G0"));
}

#[test]
fn apol1_script_accepts_zip_genotype_file_via_forced_format() {
    let root = repo_root();
    let zip_rel = "bioscripts/output/apol1-input.dat";
    let zip_path = root.join(zip_rel);
    create_zip_file(&zip_path, "test_snps.txt", "examples/apol1/test_snps.txt");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-format")
        .arg("zip")
        .arg("--input-file")
        .arg(zip_rel)
        .arg("--output-file")
        .arg("bioscripts/output/apol1-from-zip-flag.tsv")
        .arg("--participant-id")
        .arg("ZIP002")
        .arg("bioscripts/apol1.py")
        .output()
        .unwrap();

    assert!(output.status.success(), "stderr: {}", String::from_utf8_lossy(&output.stderr));
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("G1/G0"));
}
