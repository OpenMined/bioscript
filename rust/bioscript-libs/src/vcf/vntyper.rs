use crate::{LibError, LibResult};

use super::VcfRecord;

const NEGATIVE_LABEL: &str = "Negative";
const LOW_DEPTH_SCORE: f64 = 0.00469;
const HIGH_DEPTH_SCORE: f64 = 0.00515;
const ALT_DEPTH_LOW: f64 = 20.0;
const ALT_DEPTH_MID_LOW: f64 = 21.0;
const ALT_DEPTH_MID_HIGH: f64 = 100.0;
const VAR_ACTIVE_REGION_THRESHOLD: f64 = 200.0;
const MOTIF_POSITION_THRESHOLD: i64 = 60;
const EXCLUDE_MOTIFS_RIGHT: &[&str] = &["8", "9", "7", "6p", "6"];
const ALT_FOR_MOTIF_RIGHT_GG: &str = "GG";
const MOTIFS_FOR_ALT_GG: &[&str] = &[];
const EXCLUDE_ALTS_COMBINED: &[&str] = &["CCGCC", "CGGCG", "CGGCC"];
const EXCLUDE_MOTIFS_COMBINED: &[&str] = &["6", "6p", "7"];

pub fn vntyper_kestrel_rows(records: &[VcfRecord]) -> Vec<VcfRecord> {
    records.iter().map(vntyper_kestrel_row).collect()
}

pub fn vntyper_report_json(
    sample_name: &str,
    input_files: &VcfRecord,
    rows: &[VcfRecord],
) -> LibResult<String> {
    vntyper_report_json_with_context(
        sample_name,
        input_files,
        rows,
        &VcfRecord::new(),
        &VcfRecord::new(),
    )
}

pub fn vntyper_report_json_with_context(
    sample_name: &str,
    input_files: &VcfRecord,
    rows: &[VcfRecord],
    metadata: &VcfRecord,
    coverage: &VcfRecord,
) -> LibResult<String> {
    let coverage_qc = coverage_json(coverage);
    let quality_pass = coverage_quality_pass(coverage);
    let kestrel_result = compute_kestrel_result(rows);
    let screening_summary = screening_summary(&kestrel_result, quality_pass);
    let best_call = best_kestrel_call(rows).map(best_call_json);
    let report_date = metadata_value(metadata, "report_date", "runtime-generated");
    let alignment_pipeline = metadata_value(
        metadata,
        "alignment_pipeline",
        "native bioscript kestrel from FASTQ",
    );
    let value = serde_json::json!({
        "sample_name": sample_name,
        "version": "bioscript-vntyper-port",
        "report_date": report_date,
        "metadata": {
            "sample_name": sample_name,
            "vntyper_version": "bioscript-vntyper-port",
            "report_date": report_date,
            "input_files": input_files,
            "alignment_pipeline": alignment_pipeline,
            "detected_assembly": metadata_value(metadata, "detected_assembly", "unknown"),
            "detected_contig": metadata_value(metadata, "detected_contig", "unknown"),
            "bam_header_warnings": [],
        },
        "input_files": input_files,
        "coverage": coverage_qc,
        "fastp": {
            "available": false,
        },
        "algorithm_results": {
            "kestrel": kestrel_result,
            "advntr": "none",
            "quality_metrics_pass": quality_pass,
        },
        "screening_summary": screening_summary,
        "kestrel_variants": rows,
        "advntr_variants": [],
        "cross_match_summary": {
            "available": false,
            "status": "not_performed",
            "message": "adVNTR genotyping was not performed.",
        },
        "pipeline_log": [],
        "best_call": best_call,
        "kestrel_variant_count": rows.len(),
    });
    serde_json::to_string_pretty(&value)
        .map_err(|err| LibError::InvalidArguments(format!("failed to build VNtyper report: {err}")))
}

fn vntyper_kestrel_row(record: &VcfRecord) -> VcfRecord {
    let mut row = record.clone();
    let sample = row.get("Sample").cloned().unwrap_or_default();
    let parts = sample.split(':').collect::<Vec<_>>();
    let alt_depth = parts
        .get(1)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0);
    let region_depth = parts
        .get(2)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0);
    let ref_len = row.get("REF").map_or(0, String::len);
    let alt_len = row.get("ALT").map_or(0, String::len);
    let delta = alt_len as isize - ref_len as isize;
    let frame_score = delta as f64 / 3.0;
    let direction = delta.signum();
    let frameshift_amount = delta.unsigned_abs() % 3;
    let is_frameshift = delta % 3 != 0;
    let is_valid_frameshift =
        (direction > 0 && frameshift_amount == 1) || (direction < 0 && frameshift_amount == 2);
    let depth_score = if region_depth == 0.0 {
        None
    } else {
        Some(alt_depth / region_depth)
    };
    let confidence = confidence(alt_depth, region_depth, depth_score);
    let depth_confidence_pass = confidence != NEGATIVE_LABEL;
    let alt_filter_pass = alt_filter_pass(row.get("ALT").map(String::as_str), depth_score);
    let motif_filter = motif_filter(&row, is_valid_frameshift);
    for (key, value) in &motif_filter.annotations {
        row.insert((*key).to_owned(), value.clone());
    }
    let passes_vntyper_filters =
        is_valid_frameshift && depth_confidence_pass && alt_filter_pass && motif_filter.passes;

    row.insert(
        "Estimated_Depth_AlternateVariant".to_owned(),
        decimal(alt_depth),
    );
    row.insert(
        "Estimated_Depth_Variant_ActiveRegion".to_owned(),
        decimal(region_depth),
    );
    row.insert(
        "Depth_Score".to_owned(),
        depth_score.map_or_else(|| "None".to_owned(), compact_float),
    );
    row.insert("Frame_Score".to_owned(), compact_float(frame_score));
    row.insert("Confidence".to_owned(), confidence.to_owned());
    row.insert("Flag".to_owned(), flags(&row, depth_score));
    row.insert("is_frameshift".to_owned(), title_bool(is_frameshift));
    row.insert(
        "is_valid_frameshift".to_owned(),
        title_bool(is_valid_frameshift),
    );
    row.insert("alt_filter_pass".to_owned(), title_bool(alt_filter_pass));
    row.insert(
        "motif_filter_pass".to_owned(),
        title_bool(motif_filter.passes),
    );
    row.insert(
        "passes_vntyper_filters".to_owned(),
        title_bool(passes_vntyper_filters),
    );
    row
}

fn confidence(alt_depth: f64, region_depth: f64, depth_score: Option<f64>) -> &'static str {
    let Some(depth_score) = depth_score else {
        return NEGATIVE_LABEL;
    };
    let mut confidence = NEGATIVE_LABEL;
    if depth_score >= LOW_DEPTH_SCORE {
        if region_depth <= VAR_ACTIVE_REGION_THRESHOLD || depth_score == LOW_DEPTH_SCORE {
            confidence = "Low_Precision";
        }
        if alt_depth >= ALT_DEPTH_MID_HIGH && depth_score >= HIGH_DEPTH_SCORE {
            confidence = "High_Precision*";
        }
        if (ALT_DEPTH_MID_LOW..ALT_DEPTH_MID_HIGH).contains(&alt_depth)
            && (LOW_DEPTH_SCORE..=HIGH_DEPTH_SCORE).contains(&depth_score)
        {
            confidence = "Low_Precision";
        }
        if alt_depth <= ALT_DEPTH_LOW {
            confidence = "Low_Precision";
        }
        if (ALT_DEPTH_MID_LOW..ALT_DEPTH_MID_HIGH).contains(&alt_depth)
            && depth_score >= HIGH_DEPTH_SCORE
        {
            confidence = "High_Precision";
        }
        if depth_score > LOW_DEPTH_SCORE && depth_score < HIGH_DEPTH_SCORE {
            confidence = "Low_Precision";
        }
    }
    confidence
}

fn alt_filter_pass(alt: Option<&str>, depth_score: Option<f64>) -> bool {
    alt != Some("GG") || depth_score.is_some_and(|score| score >= LOW_DEPTH_SCORE)
}

struct MotifFilter {
    passes: bool,
    annotations: Vec<(&'static str, String)>,
}

fn motif_filter(row: &VcfRecord, is_valid_frameshift: bool) -> MotifFilter {
    let motifs = row
        .get("Motifs")
        .or_else(|| row.get("CHROM"))
        .cloned()
        .unwrap_or_default();
    let parts = motifs.split('-').collect::<Vec<_>>();
    if parts.len() != 2 {
        return MotifFilter {
            passes: is_valid_frameshift,
            annotations: Vec::new(),
        };
    }

    let pos = parse_row_float(row, "POS") as i64;
    let is_right_motif = pos >= MOTIF_POSITION_THRESHOLD;
    let motif = if is_right_motif { parts[0] } else { parts[1] }.to_owned();
    let alt = row.get("ALT").map(String::as_str).unwrap_or_default();

    let mut passes = is_valid_frameshift;
    if is_right_motif && EXCLUDE_MOTIFS_RIGHT.contains(&motif.as_str()) {
        passes = false;
    }
    if is_right_motif
        && alt == ALT_FOR_MOTIF_RIGHT_GG
        && !MOTIFS_FOR_ALT_GG.contains(&motif.as_str())
    {
        passes = false;
    }
    if EXCLUDE_ALTS_COMBINED.contains(&alt) {
        passes = false;
    }
    if EXCLUDE_MOTIFS_COMBINED.contains(&motif.as_str()) {
        passes = false;
    }

    MotifFilter {
        passes,
        annotations: vec![
            ("Motifs", motifs.clone()),
            ("Motif_fasta", motifs),
            ("POS_fasta", pos.to_string()),
            ("Motif", motif),
        ],
    }
}

fn flags(row: &VcfRecord, depth_score: Option<f64>) -> String {
    let mut flags = Vec::new();
    if row.get("REF").map(String::as_str) == Some("C")
        && row.get("ALT").map(String::as_str) == Some("CGGCA")
    {
        flags.push("False_Positive_4bp_Insertion");
    }
    if depth_score.is_some_and(|score| score < 0.4)
        && matches!(
            row.get("Motif").map(String::as_str),
            Some("1" | "2" | "3" | "4" | "6" | "7" | "8" | "9")
        )
    {
        flags.push("Low_Depth_Conserved_Motifs");
    }
    if flags.is_empty() {
        "Not flagged".to_owned()
    } else {
        flags.join(", ")
    }
}

fn compute_kestrel_result(rows: &[VcfRecord]) -> String {
    for row in rows {
        if row.get("passes_vntyper_filters").map(String::as_str) == Some("False") {
            continue;
        }
        let confidence = row.get("Confidence").map(String::as_str);
        let flagged = row.get("Flag").map(String::as_str) != Some("Not flagged");
        match (confidence, flagged) {
            (Some("High_Precision" | "High_Precision*"), false) => {
                return "High_Precision".to_owned();
            }
            (Some("Low_Precision"), false) => return "Low_Precision".to_owned(),
            (Some("High_Precision" | "High_Precision*"), true) => {
                return "High_Precision_flagged".to_owned();
            }
            (Some("Low_Precision"), true) => return "Low_Precision_flagged".to_owned(),
            _ => {}
        }
    }
    "negative".to_owned()
}

fn screening_summary(kestrel_result: &str, quality_pass: bool) -> &'static str {
    match (kestrel_result, quality_pass) {
        ("High_Precision", true) => {
            "Kestrel detected a high-precision pathogenic variant.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the result using orthogonal methods (e.g., SNaPshot, long-read sequencing)."
        }
        ("High_Precision", false) => {
            "Kestrel detected a high-precision pathogenic variant with quality metrics below threshold, and adVNTR genotyping was not performed.<br>Further validation using alternative methods (e.g., SNaPshot, long-read sequencing) is strongly recommended."
        }
        ("High_Precision_flagged", true) => {
            "Kestrel detected a high-precision pathogenic variant with a flagged result.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the finding using orthogonal methods (e.g., SNaPshot, long-read sequencing)."
        }
        ("Low_Precision", true) => {
            "Kestrel detected a pathogenic variant with low precision.<br>Note: adVNTR genotyping was not performed.<br>It is recommended to perform adVNTR and validate the result using alternative methods (e.g., SNaPshot, long-read sequencing)."
        }
        ("negative", true) => "No variant detected.<br>Note: adVNTR genotyping was not performed.",
        _ => "The screening was negative (no valid Kestrel or adVNTR data).",
    }
}

fn best_kestrel_call(rows: &[VcfRecord]) -> Option<&VcfRecord> {
    rows.iter().max_by(|left, right| {
        parse_row_float(left, "Depth_Score").total_cmp(&parse_row_float(right, "Depth_Score"))
    })
}

fn best_call_json(row: &VcfRecord) -> serde_json::Value {
    serde_json::json!({
        "CHROM": row.get("CHROM").cloned().unwrap_or_default(),
        "POS": row.get("POS").cloned().unwrap_or_default(),
        "REF": row.get("REF").cloned().unwrap_or_default(),
        "ALT": row.get("ALT").cloned().unwrap_or_default(),
        "Estimated_Depth_AlternateVariant": parse_row_float(row, "Estimated_Depth_AlternateVariant"),
        "Estimated_Depth_Variant_ActiveRegion": parse_row_float(row, "Estimated_Depth_Variant_ActiveRegion"),
        "Depth_Score": parse_row_float(row, "Depth_Score"),
        "Confidence": row.get("Confidence").cloned().unwrap_or_default(),
        "passes_vntyper_filters": row.get("passes_vntyper_filters").map(String::as_str) == Some("True"),
    })
}

fn parse_row_float(row: &VcfRecord, key: &str) -> f64 {
    row.get(key)
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0)
}

fn metadata_value<'a>(metadata: &'a VcfRecord, key: &str, default: &'a str) -> &'a str {
    metadata.get(key).map_or(default, String::as_str)
}

fn coverage_json(coverage: &VcfRecord) -> serde_json::Value {
    let quality_pass = coverage_quality_pass(coverage);
    serde_json::json!({
        "mean": numeric_or_null(coverage, "mean"),
        "median": numeric_or_null(coverage, "median"),
        "stdev": numeric_or_null(coverage, "stdev"),
        "min": numeric_or_null(coverage, "min"),
        "max": numeric_or_null(coverage, "max"),
        "region_length": numeric_or_null(coverage, "region_length"),
        "uncovered_bases": numeric_or_null(coverage, "uncovered_bases"),
        "percent_uncovered": numeric_or_null(coverage, "percent_uncovered"),
        "threshold": 100,
        "quality_pass": quality_pass,
        "status": if quality_pass { "pass" } else { "warning" },
    })
}

fn coverage_quality_pass(coverage: &VcfRecord) -> bool {
    coverage
        .get("mean")
        .and_then(|value| value.parse::<f64>().ok())
        .is_none_or(|mean| mean >= 100.0)
}

fn numeric_or_null(coverage: &VcfRecord, key: &str) -> serde_json::Value {
    coverage
        .get(key)
        .and_then(|value| value.parse::<f64>().ok())
        .map_or(serde_json::Value::Null, serde_json::Value::from)
}

fn title_bool(value: bool) -> String {
    if value { "True" } else { "False" }.to_owned()
}

fn decimal(value: f64) -> String {
    format!("{value:.1}")
}

fn compact_float(value: f64) -> String {
    let mut text = value.to_string();
    if text.contains('.') {
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.push('0');
        }
    }
    text
}
