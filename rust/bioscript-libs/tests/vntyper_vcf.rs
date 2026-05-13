use bioscript_libs::vcf::{parse_kestrel_vcf, vntyper::vntyper_kestrel_rows};

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
            "MUC1\t100\tC\tCGGCA\t120.0\t10000.0\t0.012\tHigh_Precision*\tTrue\tTrue\tTrue",
            "MUC1\t160\tATG\tA\t50.0\t10000.0\t0.005\tLow_Precision\tTrue\tTrue\tTrue",
            "MUC1\t220\tC\tCGG\t5.0\t10000.0\t0.0005\tNegative\tFalse\tTrue\tFalse",
        ]
    );
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
