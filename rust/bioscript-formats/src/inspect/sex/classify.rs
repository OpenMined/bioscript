use crate::inspect::DetectedKind;

use super::{InferredSex, SexDetectionConfidence, SexInference, SexStats};

pub(super) fn classify_stats(stats: &SexStats, kind: DetectedKind) -> SexInference {
    match kind {
        DetectedKind::Vcf => classify_vcf_stats(stats),
        DetectedKind::GenotypeText | DetectedKind::Unknown => classify_y_fingerprint_stats(stats),
        _ => unsupported_sex_inference(),
    }
}

pub(super) fn supports_sex_detection(kind: DetectedKind) -> bool {
    matches!(
        kind,
        DetectedKind::Vcf | DetectedKind::GenotypeText | DetectedKind::Unknown
    )
}

pub(super) fn unsupported_sex_inference() -> SexInference {
    SexInference {
        sex: InferredSex::Unknown,
        confidence: SexDetectionConfidence::Low,
        method: "unsupported_source_type".to_owned(),
        evidence: vec!["sex detection currently supports genotype text and VCF inputs".to_owned()],
    }
}

fn classify_y_fingerprint_stats(stats: &SexStats) -> SexInference {
    let called_y_pct = called_y_percent(stats);
    let (sex, confidence) = if stats.total_y_snps > 1000
        && stats.male_markers_found >= 3
        && stats.male_markers_called == 0
        && called_y_pct < 20.0
    {
        (InferredSex::Female, SexDetectionConfidence::High)
    } else if stats.called_y_snps > 500 && called_y_pct >= 50.0 {
        (InferredSex::Male, SexDetectionConfidence::High)
    } else if stats.called_y_snps > 100 && stats.male_markers_called > 10 {
        (InferredSex::Male, SexDetectionConfidence::Medium)
    } else if stats.total_y_snps > 1000 && stats.called_y_snps < 10 && stats.male_markers_called < 3
    {
        (
            InferredSex::Female,
            if stats.called_y_snps == 0 {
                SexDetectionConfidence::High
            } else {
                SexDetectionConfidence::Medium
            },
        )
    } else if stats.called_y_snps > 50 || stats.male_markers_called > 5 {
        (InferredSex::Male, SexDetectionConfidence::Medium)
    } else if stats.called_y_snps < 10 && stats.total_y_snps > 0 {
        (InferredSex::Female, SexDetectionConfidence::Medium)
    } else {
        (InferredSex::Unknown, SexDetectionConfidence::Low)
    };
    SexInference {
        sex,
        confidence,
        method: "y_fingerprint".to_owned(),
        evidence: y_fingerprint_evidence(stats),
    }
}

fn classify_vcf_stats(stats: &SexStats) -> SexInference {
    let x_het_pct = if stats.x_diploid_gt_sites == 0 {
        0.0
    } else {
        stats.x_het_gt_sites as f64 * 100.0 / stats.x_diploid_gt_sites as f64
    };
    let y_to_x_pct = if stats.x_non_par_sites == 0 {
        0.0
    } else {
        stats.called_y_snps as f64 * 100.0 / stats.x_non_par_sites as f64
    };
    let female_like_x =
        stats.x_non_par_sites >= 50 && stats.x_diploid_gt_sites > 0 && x_het_pct >= 2.0;
    let male_like_x = stats.x_haploid_gt_sites >= 20
        || (stats.x_non_par_sites >= 50 && stats.x_diploid_gt_sites > 0 && x_het_pct <= 0.2);
    let male_like_y = stats.called_y_snps > 500;
    let strong_y_signal =
        stats.called_y_snps >= 10_000 || (stats.x_non_par_sites >= 1000 && y_to_x_pct >= 10.0);
    let (sex, confidence) = if strong_y_signal {
        (InferredSex::Male, SexDetectionConfidence::High)
    } else if female_like_x {
        (
            InferredSex::Female,
            if stats.x_non_par_sites >= 1000 && x_het_pct >= 10.0 {
                SexDetectionConfidence::High
            } else {
                SexDetectionConfidence::Medium
            },
        )
    } else if male_like_x {
        (InferredSex::Male, SexDetectionConfidence::High)
    } else if male_like_y {
        (InferredSex::Unknown, SexDetectionConfidence::Medium)
    } else {
        (InferredSex::Unknown, SexDetectionConfidence::Low)
    };
    SexInference {
        sex,
        confidence,
        method: "vcf_non_par_x_gt_y_count".to_owned(),
        evidence: vec![
            format!("x_non_par_sites={}", stats.x_non_par_sites),
            format!("x_haploid_gt_sites={}", stats.x_haploid_gt_sites),
            format!("x_diploid_gt_sites={}", stats.x_diploid_gt_sites),
            format!("x_het_gt_sites={}", stats.x_het_gt_sites),
            format!("x_het_pct={x_het_pct:.2}"),
            format!("called_y_snps={}", stats.called_y_snps),
            format!("y_to_x_pct={y_to_x_pct:.2}"),
        ],
    }
}

fn y_fingerprint_evidence(stats: &SexStats) -> Vec<String> {
    let mut evidence = vec![
        format!("total_y_snps={}", stats.total_y_snps),
        format!("called_y_snps={}", stats.called_y_snps),
        format!("called_y_pct={:.1}", called_y_percent(stats)),
        format!("male_markers_found={}", stats.male_markers_found),
        format!("male_markers_called={}", stats.male_markers_called),
    ];
    if !stats.y_examples.is_empty() {
        evidence.push(format!("y_examples={}", stats.y_examples.join(",")));
    }
    evidence
}

fn called_y_percent(stats: &SexStats) -> f64 {
    if stats.total_y_snps == 0 {
        0.0
    } else {
        let called_y_snps = u32::try_from(stats.called_y_snps).unwrap_or(u32::MAX);
        let total_y_snps = u32::try_from(stats.total_y_snps).unwrap_or(u32::MAX);
        f64::from(called_y_snps) * 100.0 / f64::from(total_y_snps)
    }
}
