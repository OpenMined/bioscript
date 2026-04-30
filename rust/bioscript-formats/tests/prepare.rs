use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use bioscript_formats::{
    GenotypeSourceFormat, PrepareRequest, PreparedPaths, prepare_indexes, shell_flags,
};

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-prepare-{label}-{}-{nanos}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn request(root: PathBuf, cwd: PathBuf, cache_dir: PathBuf) -> PrepareRequest {
    PrepareRequest {
        root,
        cwd,
        cache_dir,
        input_file: None,
        input_format: None,
        reference_file: None,
    }
}

#[test]
fn relative_input_path_resolves_inside_root() {
    let root = temp_dir("relative-input-root");
    let cwd = temp_dir("relative-input-cwd");
    fs::write(
        root.join("sample.txt"),
        "rsid\tchromosome\tposition\tgenotype\n",
    )
    .unwrap();

    let mut req = request(root.clone(), cwd.clone(), PathBuf::from("cache"));
    req.input_file = Some("sample.txt".to_owned());

    let prepared = prepare_indexes(&req).unwrap();

    let expected_input = root.join("sample.txt").canonicalize().unwrap();
    assert_eq!(
        prepared.input_file.as_deref(),
        Some(expected_input.as_path())
    );
    assert_eq!(prepared.input_index, None);
    assert_eq!(prepared.cache_dir, cwd.join("cache"));
}

#[test]
fn relative_path_escape_is_rejected() {
    let root = temp_dir("relative-escape-root");
    let cwd = temp_dir("relative-escape-cwd");
    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("../outside.txt".to_owned());

    let err = prepare_indexes(&req).unwrap_err();

    assert!(err.contains("path escapes bioscript root"), "{err}");
}

#[test]
fn absolute_path_outside_root_is_rejected() {
    let root = temp_dir("absolute-escape-root");
    let cwd = temp_dir("absolute-escape-cwd");
    let outside = temp_dir("absolute-escape-outside").join("sample.txt");
    fs::write(&outside, "rsid\tchromosome\tposition\tgenotype\n").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some(outside.to_string_lossy().into_owned());

    let err = prepare_indexes(&req).unwrap_err();

    assert!(err.contains("path escapes bioscript root"), "{err}");
}

#[test]
fn absolute_path_inside_root_is_allowed() {
    let root = temp_dir("absolute-inside-root");
    let cwd = temp_dir("absolute-inside-cwd");
    let input = root.join("sample.txt");
    fs::write(&input, "rsid\tchromosome\tposition\tgenotype\n").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some(input.to_string_lossy().into_owned());

    let prepared = prepare_indexes(&req).unwrap();

    assert_eq!(
        prepared.input_file.as_deref(),
        Some(input.canonicalize().unwrap().as_path())
    );
}

#[test]
fn missing_input_file_returns_clear_error() {
    let root = temp_dir("missing-input-root");
    let cwd = temp_dir("missing-input-cwd");

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("missing.txt".to_owned());

    let err = prepare_indexes(&req).unwrap_err();

    assert!(err.contains("failed to resolve"), "{err}");
    assert!(err.contains("missing.txt"), "{err}");
}

#[test]
fn absolute_cache_dir_is_preserved() {
    let root = temp_dir("absolute-cache-root");
    let cwd = temp_dir("absolute-cache-cwd");
    let cache = temp_dir("absolute-cache-target");

    let req = request(root, cwd, cache.clone());
    let prepared = prepare_indexes(&req).unwrap();

    assert_eq!(prepared.cache_dir, cache);
}

#[test]
fn adjacent_cram_index_is_detected_without_rebuilding() {
    let root = temp_dir("adjacent-cram-root");
    let cwd = temp_dir("adjacent-cram-cwd");
    let input_path = root.join("sample.cram");
    let index_path = root.join("sample.cram.crai");
    fs::write(&input_path, b"not a real cram").unwrap();
    fs::write(&index_path, b"not a real crai").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("sample.cram".to_owned());

    let prepared = prepare_indexes(&req).unwrap();

    let expected_index = index_path.canonicalize().unwrap();
    assert_eq!(
        prepared.input_index.as_deref(),
        Some(expected_index.as_path())
    );
}

#[test]
fn adjacent_short_cram_index_is_detected_without_rebuilding() {
    let root = temp_dir("adjacent-short-cram-root");
    let cwd = temp_dir("adjacent-short-cram-cwd");
    let input_path = root.join("sample.cram");
    let index_path = root.join("sample.crai");
    fs::write(&input_path, b"not a real cram").unwrap();
    fs::write(&index_path, b"not a real crai").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("sample.cram".to_owned());

    let prepared = prepare_indexes(&req).unwrap();

    assert_eq!(
        prepared.input_index.as_deref(),
        Some(index_path.canonicalize().unwrap().as_path())
    );
}

#[test]
fn bam_without_adjacent_index_returns_clear_error() {
    let root = temp_dir("bam-root");
    let cwd = temp_dir("bam-cwd");
    fs::write(root.join("sample.bam"), b"not a real bam").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("sample.bam".to_owned());

    let err = prepare_indexes(&req).unwrap_err();

    assert!(
        err.contains("alignment indexing only supports CRAM"),
        "{err}"
    );
}

#[test]
fn cram_without_adjacent_index_reports_build_failure() {
    let root = temp_dir("invalid-cram-root");
    let cwd = temp_dir("invalid-cram-cwd");
    fs::write(root.join("sample.cram"), b"not a real cram").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("sample.cram".to_owned());

    let err = prepare_indexes(&req).unwrap_err();

    assert!(err.contains("failed to build alignment index"), "{err}");
    assert!(err.contains("sample.cram"), "{err}");
}

#[test]
fn adjacent_fasta_index_is_detected() {
    let root = temp_dir("adjacent-fasta-root");
    let cwd = temp_dir("adjacent-fasta-cwd");
    let fasta = root.join("ref.fa");
    let fai = root.join("ref.fa.fai");
    fs::write(&fasta, b">chr1\nACGT\n").unwrap();
    fs::write(&fai, b"chr1\t4\t6\t4\t5\n").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.reference_file = Some("ref.fa".to_owned());

    let prepared = prepare_indexes(&req).unwrap();

    let expected_fasta = fasta.canonicalize().unwrap();
    let expected_fai = fai.canonicalize().unwrap();
    assert_eq!(
        prepared.reference_file.as_deref(),
        Some(expected_fasta.as_path())
    );
    assert_eq!(
        prepared.reference_index.as_deref(),
        Some(expected_fai.as_path())
    );
}

#[test]
fn invalid_fasta_reference_reports_index_build_failure() {
    let root = temp_dir("invalid-fasta-root");
    let cwd = temp_dir("invalid-fasta-cwd");
    fs::write(root.join("ref.fa"), b"not fasta\n").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.reference_file = Some("ref.fa".to_owned());

    let err = prepare_indexes(&req).unwrap_err();

    assert!(err.contains("failed to build FASTA index"), "{err}");
    assert!(err.contains("ref.fa"), "{err}");
}

#[test]
fn fasta_index_is_generated_in_cache_when_missing() {
    let root = temp_dir("generated-fasta-root");
    let cwd = temp_dir("generated-fasta-cwd");
    let cache = cwd.join("cache");
    fs::write(root.join("ref.fa"), b">chr1\nACGT\n").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.reference_file = Some("ref.fa".to_owned());

    let prepared = prepare_indexes(&req).unwrap();
    let cached_reference = prepared.reference_file.expect("cached reference");
    let cached_index = prepared.reference_index.expect("cached reference index");

    assert!(cached_reference.starts_with(&cache));
    assert!(cached_reference.exists());
    assert!(cached_index.starts_with(&cache));
    assert!(cached_index.exists());
}

#[test]
fn fasta_without_extension_uses_fai_extension() {
    let root = temp_dir("fasta-no-extension-root");
    let cwd = temp_dir("fasta-no-extension-cwd");
    let cache = cwd.join("cache");
    fs::write(root.join("reference"), b">chr1\nACGT\n").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.reference_file = Some("reference".to_owned());

    let prepared = prepare_indexes(&req).unwrap();
    let cached_reference = prepared.reference_file.expect("cached reference");
    let cached_index = prepared.reference_index.expect("cached reference index");

    assert!(cached_reference.starts_with(&cache));
    assert_eq!(
        cached_index.extension().and_then(|ext| ext.to_str()),
        Some("fai")
    );
    assert!(cached_index.exists());
}

#[test]
fn shell_flags_quote_paths_with_spaces_and_single_quotes() {
    let prepared = PreparedPaths {
        input_file: Some(PathBuf::from("/tmp/input files/sample's.cram")),
        input_index: Some(PathBuf::from("/tmp/input files/sample.cram.crai")),
        reference_file: Some(PathBuf::from("/tmp/ref files/ref's.fa")),
        reference_index: Some(PathBuf::from("/tmp/ref files/ref.fa.fai")),
        cache_dir: PathBuf::from("/tmp/cache"),
    };

    let flags = shell_flags(&prepared);

    assert!(flags.contains("--input-file '/tmp/input files/sample'\"'\"'s.cram'"));
    assert!(flags.contains("--input-index '/tmp/input files/sample.cram.crai'"));
    assert!(flags.contains("--reference-file '/tmp/ref files/ref'\"'\"'s.fa'"));
    assert!(flags.contains("--reference-index '/tmp/ref files/ref.fa.fai'"));
}

#[test]
fn shell_flags_are_empty_when_nothing_was_prepared() {
    let prepared = PreparedPaths {
        cache_dir: PathBuf::from("/tmp/cache"),
        ..PreparedPaths::default()
    };

    assert_eq!(shell_flags(&prepared), "");
}

#[test]
fn explicit_cram_format_triggers_index_detection_for_non_cram_extension() {
    let root = temp_dir("forced-cram-root");
    let cwd = temp_dir("forced-cram-cwd");
    let input = root.join("sample.dat");
    let crai = root.join("sample.crai");
    fs::write(&input, b"not a real cram").unwrap();
    fs::write(&crai, b"not a real crai").unwrap();

    let mut req = request(root, cwd, PathBuf::from("cache"));
    req.input_file = Some("sample.dat".to_owned());
    req.input_format = Some(GenotypeSourceFormat::Cram);

    let err = prepare_indexes(&req).unwrap_err();

    assert!(
        err.contains("alignment indexing only supports CRAM")
            || err.contains("failed to build alignment index"),
        "{err}"
    );
}
