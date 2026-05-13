use super::VcfRecord;

const NEGATIVE_LABEL: &str = "Negative";
const LOW_DEPTH_SCORE: f64 = 0.00469;
const HIGH_DEPTH_SCORE: f64 = 0.00515;
const ALT_DEPTH_LOW: f64 = 20.0;
const ALT_DEPTH_MID_LOW: f64 = 21.0;
const ALT_DEPTH_MID_HIGH: f64 = 100.0;
const VAR_ACTIVE_REGION_THRESHOLD: f64 = 200.0;

pub fn vntyper_kestrel_rows(records: &[VcfRecord]) -> Vec<VcfRecord> {
    records.iter().map(vntyper_kestrel_row).collect()
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
    let motif_filter_pass = motif_filter_pass(&row, is_valid_frameshift);
    let passes_vntyper_filters =
        is_valid_frameshift && depth_confidence_pass && alt_filter_pass && motif_filter_pass;

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
        title_bool(motif_filter_pass),
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

fn motif_filter_pass(row: &VcfRecord, is_valid_frameshift: bool) -> bool {
    let Some(chrom) = row.get("CHROM") else {
        return is_valid_frameshift;
    };
    let parts = chrom.split('-').collect::<Vec<_>>();
    if parts.len() != 2 {
        return true;
    }
    is_valid_frameshift
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
