use bioscript_libs::vcf::{VcfRecord, parse_kestrel_vcf, vntyper::vntyper_kestrel_rows};
use serde_json::Value;

#[test]
fn parses_kestrel_vcf_sample_depth_fields_for_vntyper() {
    let records = parse_kestrel_vcf(concat!(
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnegative\n",
        "MUC1\t59\t.\tG\tGT\t.\tPASS\t.\tGT\tIns:491:18434\n",
    ))
    .unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].get("CHROM").map(String::as_str), Some("MUC1"));
    assert_eq!(records[0].get("POS").map(String::as_str), Some("59"));
    assert_eq!(records[0].get("REF").map(String::as_str), Some("G"));
    assert_eq!(records[0].get("ALT").map(String::as_str), Some("GT"));
    assert_eq!(
        records[0].get("Sample").map(String::as_str),
        Some("Ins:491:18434")
    );
}

#[test]
fn builds_vntyper_kestrel_call_rows_for_fixture() {
    let records = parse_kestrel_vcf(include_str!(
        "../../../ports/vntyper/tests/fixtures/kestrel_minimal.vcf"
    ))
    .unwrap();
    let rows = vntyper_kestrel_rows(&records);

    let selected = rows
        .iter()
        .map(|row| {
            [
                "CHROM",
                "POS",
                "REF",
                "ALT",
                "Estimated_Depth_AlternateVariant",
                "Estimated_Depth_Variant_ActiveRegion",
                "Depth_Score",
                "Confidence",
                "is_valid_frameshift",
                "alt_filter_pass",
                "passes_vntyper_filters",
            ]
            .into_iter()
            .map(|key| row.get(key).cloned().unwrap_or_default())
            .collect::<Vec<_>>()
            .join("\t")
        })
        .collect::<Vec<_>>();

    assert_eq!(
        selected,
        vec![
            "C-Q\t100\tC\tCGGCA\t120.0\t10000.0\t0.012\tHigh_Precision*\tTrue\tTrue\tTrue",
            "C-Q\t160\tATG\tA\t50.0\t10000.0\t0.005\tLow_Precision\tTrue\tTrue\tTrue",
            "C-Q\t220\tC\tCGG\t5.0\t10000.0\t0.0005\tNegative\tFalse\tTrue\tFalse",
        ]
    );
}

#[test]
fn annotates_and_filters_vntyper_motif_fields_like_python_port() {
    let records = parse_kestrel_vcf(concat!(
        "##fileformat=VCFv4.2\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n",
        "6-M\t61\t.\tG\tGT\t.\tPASS\t.\tGT:GDP:DP\t1:80:1000\n",
        "5C-M\t61\t.\tG\tGT\t.\tPASS\t.\tGT:GDP:DP\t1:80:1000\n",
        "5C-M\t61\t.\tG\tGG\t.\tPASS\t.\tGT:GDP:DP\t1:80:1000\n",
        "5C-M\t61\t.\tG\tGCCGCC\t.\tPASS\t.\tGT:GDP:DP\t1:80:1000\n",
    ))
    .unwrap();

    let rows = vntyper_kestrel_rows(&records);

    assert_eq!(rows[0].get("Motif").map(String::as_str), Some("6"));
    assert_eq!(
        rows[0].get("motif_filter_pass").map(String::as_str),
        Some("False")
    );
    assert_eq!(
        rows[1].get("motif_filter_pass").map(String::as_str),
        Some("True")
    );
    // Upstream-faithful motif_correction keeps a `G>GG` right-motif
    // insertion in a non-excluded motif (MOTIFS_FOR_ALT_GG is empty, so the
    // legacy GG `.any()` guard does not restrict). This is the canonical
    // MUC1 dup; the old per-row approximation wrongly rejected it.
    assert_eq!(
        rows[2].get("motif_filter_pass").map(String::as_str),
        Some("True")
    );
    // Not a valid frameshift (delta = +5), so it fails regardless of motif.
    assert_eq!(
        rows[3].get("motif_filter_pass").map(String::as_str),
        Some("False")
    );
    assert_eq!(rows[1].get("Motifs").map(String::as_str), Some("5C-M"));
    assert_eq!(rows[1].get("Motif_fasta").map(String::as_str), Some("5C-M"));
    assert_eq!(rows[1].get("POS_fasta").map(String::as_str), Some("61"));
}

#[test]
fn builds_vntyper_report_summary_for_fixture() {
    let records = parse_kestrel_vcf(include_str!(
        "../../../ports/vntyper/tests/fixtures/kestrel_minimal.vcf"
    ))
    .unwrap();
    let rows = vntyper_kestrel_rows(&records);
    let mut input_files = VcfRecord::new();
    input_files.insert("vcf".to_owned(), "kestrel_minimal.vcf".to_owned());
    let report: Value = serde_json::from_str(
        &bioscript_libs::vcf::vntyper_report_json("fixture", &input_files, &rows).unwrap(),
    )
    .unwrap();
    let expected: Value = serde_json::from_str(include_str!(
        "../../../ports/vntyper/tests/fixtures/kestrel_minimal_expected_report.json"
    ))
    .unwrap();

    assert_eq!(report["sample_name"], "fixture");
    assert_eq!(
        report["algorithm_results"]["kestrel"],
        expected["algorithm_results"]["kestrel"]
    );
    assert_eq!(
        report["algorithm_results"]["advntr"],
        expected["algorithm_results"]["advntr"]
    );
    assert_eq!(
        report["algorithm_results"]["quality_metrics_pass"],
        expected["algorithm_results"]["quality_metrics_pass"]
    );
    assert_eq!(report["coverage"]["status"], expected["coverage"]["status"]);
    assert_eq!(
        report["coverage"]["quality_pass"],
        expected["coverage"]["quality_pass"]
    );
    assert_eq!(report["screening_summary"], expected["screening_summary"]);
    assert_eq!(
        report["kestrel_variant_count"],
        expected["kestrel_variant_count"]
    );
    assert_eq!(report["best_call"], expected["best_call"]);
}

#[test]
fn ignores_metadata_and_blank_lines_until_header() {
    let records = parse_kestrel_vcf(concat!(
        "\n",
        "##fileformat=VCFv4.2\n",
        "MUC1\t1\t.\tA\tT\t.\tPASS\t.\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        "MUC1\t21\t.\tT\tG\t.\tPASS\t.\n",
    ))
    .unwrap();

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].get("POS").map(String::as_str), Some("21"));
}
