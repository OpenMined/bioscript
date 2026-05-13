use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace rust dir")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn temp_dir(label: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock drift")
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "bioscript-cli-tests-tmp-{label}-{}-{nanos}",
        std::process::id()
    ));
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn hello_world_script_runs_via_cli_and_writes_within_root() {
    let root = repo_root();
    let output_path = root.join("bioscripts/output/hello-world.txt");
    if output_path.exists() {
        fs::remove_file(&output_path).unwrap();
    }

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("bioscripts/hello-world.py")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("hello from bioscript"));
    assert!(stdout.contains("2 + 3 = 5"));

    let written = fs::read_to_string(output_path).unwrap();
    assert!(written.contains("hello from bioscript"));
    assert!(written.contains("loaded: sample input for bioscript"));
}

#[test]
fn path_escape_is_rejected() {
    let root = repo_root();
    let script = root.join("rust/bioscript-cli/tests/fixtures/path_escape.py");
    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg(script)
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("path escapes bioscript root"));
}

#[test]
fn trace_report_is_written_for_hello_world() {
    let root = repo_root();
    let trace_path = root.join("bioscripts/output/hello-world.trace.tsv");
    if trace_path.exists() {
        fs::remove_file(&trace_path).unwrap();
    }

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--trace-report")
        .arg("bioscripts/output/hello-world.trace.tsv")
        .arg("bioscripts/hello-world.py")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let trace = fs::read_to_string(trace_path).unwrap();
    assert!(trace.contains("step\tline\tcode"));
    assert!(trace.contains("hello from bioscript"));
}

#[test]
fn batch_lookup_query_plan_runs_and_preserves_requested_result_order() {
    let root = repo_root();
    let script = root.join("rust/bioscript-cli/tests/fixtures/batch_lookup.py");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(script)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("AG"));
    assert!(stdout.contains("CT"));
    assert!(stdout.contains("II"));
}

#[test]
fn lookup_variant_details_returns_counts_and_decision_fields() {
    let root = repo_root();
    let script = root.join("rust/bioscript-cli/tests/fixtures/lookup_details.py");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(script)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("VariantObservation"));
    assert!(stdout.contains("genotype='AG'"));
    assert!(stdout.contains("raw_counts={"));
    assert!(stdout.contains("decision="));
    assert!(stdout.contains("evidence=["));
}

#[test]
fn vntyper_bioscript_program_runs_via_cli_and_writes_command_plan() {
    let root = repo_root();
    let output_path = root.join("target/vntyper-bs-plan.tsv");
    if output_path.exists() {
        fs::remove_file(&output_path).unwrap();
    }

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("ports/vntyper/test-data/example_6449_hg19_subset.bam")
        .arg("--output-file")
        .arg("target/vntyper-bs-plan.tsv")
        .arg("--participant-id")
        .arg("positive")
        .arg("ports/vntyper/bioscript/vntyper.bs")
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let plan = fs::read_to_string(&output_path).unwrap();
    assert!(plan.contains("participant_id"));
    assert!(plan.contains("positive"));
    assert!(plan.contains("samtools_view_command"));
    assert!(plan.contains("chr1:155158000-155163000"));
    assert!(plan.contains("kestrel_command"));
    assert!(plan.contains("bcftools_sort_command"));
    fs::remove_file(output_path).unwrap();
}

#[test]
fn inspect_subcommand_reports_detected_vendor_and_platform() {
    let root = repo_root();
    let path = root.join("rust/bioscript-formats/tests/fixtures/ancestrydna_v2_sample.txt");

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("inspect")
        .arg(path)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("kind\tgenotype_text"));
    assert!(stdout.contains("vendor\tAncestryDNA"));
    assert!(stdout.contains("platform_version\tV2.0"));
    assert!(stdout.contains("assembly\tgrch37"));
    assert!(stdout.contains("duration_ms\t"));
}

#[test]
fn variant_manifest_runs_directly_via_cli() {
    let root = repo_root();
    let dir = temp_dir("variant-manifest");
    let manifest = dir.join("rs1.yaml");
    fs::write(
        &manifest,
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(&manifest)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("kind\tname\tpath"));
    assert!(stdout.contains("example-rs73885319"));
    assert!(stdout.contains("AG"));
}

#[test]
fn panel_manifest_runs_directly_via_cli() {
    let root = repo_root();
    let dir = temp_dir("panel-manifest");
    let variants_dir = dir.join("variants");
    fs::create_dir_all(&variants_dir).unwrap();
    fs::write(
        variants_dir.join("rs73885319.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();
    fs::write(
        variants_dir.join("rs60910145.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs60910145"
tags:
  - "type:trait"
identifiers:
  rsids:
    - "rs60910145"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265988
alleles:
  kind: "snv"
  ref: "T"
  alts:
    - "G"
"#,
    )
    .unwrap();
    let panel = dir.join("panel.yaml");
    fs::write(
        &panel,
        r#"
schema: "bioscript:panel:1.0"
version: "1.0"
name: "example-panel"
tags:
  - "type:trait"
members:
  - kind: "variant"
    path: "variants/rs73885319.yaml"
    version: "1.0"
  - kind: "variant"
    path: "variants/rs60910145.yaml"
    version: "1.0"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg("--filter")
        .arg("name=rs73885319")
        .arg(&panel)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("example-rs73885319"));
    assert!(!stdout.contains("example-rs60910145"));
}

#[test]
fn assay_manifest_runs_directly_via_cli() {
    let root = repo_root();
    let dir = temp_dir("assay-manifest");
    fs::write(
        dir.join("rs73885319.yaml"),
        r#"
schema: "bioscript:variant:1.0"
version: "1.0"
name: "example-rs73885319"
identifiers:
  rsids:
    - "rs73885319"
coordinates:
  grch38:
    chrom: "22"
    pos: 36265860
alleles:
  kind: "snv"
  ref: "A"
  alts:
    - "G"
"#,
    )
    .unwrap();
    let assay = dir.join("assay.yaml");
    fs::write(
        &assay,
        r#"
schema: "bioscript:assay:1.0"
version: "1.0"
name: "example-assay"
members:
  - kind: "variant"
    path: "rs73885319.yaml"
    version: "1.0"
interpretations:
  - id: "example_status"
    kind: "bioscript"
    path: "example.py"
    derived_from:
      - "rs73885319.yaml"
    emits:
      - key: "example_status"
        label: "Example status"
        value_type: "string"
        format: "badge"
"#,
    )
    .unwrap();

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg(&assay)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("example-rs73885319"));
    assert!(stdout.contains("AG"));
}

#[test]
fn report_runs_variant_catalogue_analysis_with_sandboxed_assets() {
    let root = repo_root();
    let dir = temp_dir("catalogue-report");
    let (assay, output_dir) = write_catalogue_report_fixture(&dir);

    let output = Command::new(env!("CARGO_BIN_EXE_bioscript"))
        .current_dir(&root)
        .arg("report")
        .arg(&assay)
        .arg("--input-file")
        .arg("old/examples/apol1/test_snps.txt")
        .arg("--output-dir")
        .arg(&output_dir)
        .output()
        .unwrap();

    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let participant = "test_snps";
    let analysis_tsv =
        fs::read_to_string(output_dir.join(format!("analysis/{participant}/catalogue_status.tsv")))
            .unwrap();
    assert!(analysis_tsv.contains("catalogue_assets_loaded"));
    let observations = fs::read_to_string(output_dir.join("observations.tsv")).unwrap();
    assert!(observations.contains("example-rs73885319"));
    assert!(observations.contains("AG"));
}

fn write_catalogue_report_fixture(dir: &Path) -> (PathBuf, PathBuf) {
    let catalogue_dir = dir.join("catalogue");
    fs::create_dir_all(&catalogue_dir).unwrap();
    write_catalogue_manifest(&catalogue_dir);
    write_catalogue_tables(&catalogue_dir);
    write_catalogue_analysis(dir);
    let assay = write_catalogue_assay(dir);
    (assay, dir.join("report"))
}

fn write_catalogue_manifest(catalogue_dir: &Path) {
    fs::write(
        catalogue_dir.join("variants.yaml"),
        r#"
schema: "bioscript:variant-catalogue:1.0"
version: "1.0"
name: "apol1-catalogue"
variants:
  source: "variants.tsv"
  format: "tsv"
  key: "variant_id"
  columns:
    variant_id:
      role: "id"
      required: true
    name:
      role: "name"
    rsid:
      role: "identifier.rsid"
    kind:
      role: "alleles.kind"
    ref:
      role: "alleles.ref"
    alts:
      role: "alleles.alts"
    grch38_chrom:
      role: "coordinates.grch38.chrom"
    grch38_pos:
      role: "coordinates.grch38.pos"
      type: "integer"
provenance:
  sources:
    - id: "dbsnp"
      kind: "database"
      label: "dbSNP"
      url_template: "https://www.ncbi.nlm.nih.gov/snp/{rsid}"
"#,
    )
    .unwrap();
}

fn write_catalogue_tables(catalogue_dir: &Path) {
    fs::write(
        catalogue_dir.join("variants.tsv"),
        "variant_id\tname\trsid\tkind\tref\talts\tgrch38_chrom\tgrch38_pos\napol1-g1-site-1\texample-rs73885319\trs73885319\tsnv\tA\tG\t22\t36265860\n",
    )
    .unwrap();
    fs::write(
        catalogue_dir.join("findings.tsv"),
        "variant_id\tlabel\napol1-g1-site-1\tcatalogue_assets_loaded\n",
    )
    .unwrap();
}

fn write_catalogue_analysis(dir: &Path) {
    fs::write(
        dir.join("analysis.py"),
        r#"
observations = bioscript.read_tsv(bioscript.context["observations_file"])
variants = bioscript.read_tsv(bioscript.context["assets"]["variants"])
findings = bioscript.read_tsv(bioscript.context["assets"]["findings"])

status = "missing"
if len(observations) > 0 and len(variants) > 0 and len(findings) > 0:
    status = findings[0]["label"]

bioscript.write_tsv(bioscript.context["output_file"], [
    {
        "participant_id": bioscript.context["participant_id"],
        "status": status,
    }
])
"#,
    )
    .unwrap();
}

fn write_catalogue_assay(dir: &Path) -> PathBuf {
    let assay = dir.join("assay.yaml");
    fs::write(
        &assay,
        r#"
schema: "bioscript:assay:1.0"
version: "1.0"
name: "catalogue-assay"
members:
  - kind: "variant-catalogue"
    path: "catalogue/variants.yaml"
    version: "1.0"
analyses:
  - id: "catalogue_status"
    kind: "bioscript"
    path: "analysis.py"
    output_format: "tsv"
    derived_from:
      - "catalogue/variants.yaml"
    assets:
      - id: "variants"
        path: "catalogue/variants.tsv"
      - id: "findings"
        path: "catalogue/findings.tsv"
    emits:
      - key: "status"
        label: "Status"
        value_type: "string"
        format: "badge"
"#,
    )
    .unwrap();
    assay
}
