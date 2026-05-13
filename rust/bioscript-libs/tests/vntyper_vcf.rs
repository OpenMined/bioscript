use bioscript_libs::vcf::parse_kestrel_vcf;

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
