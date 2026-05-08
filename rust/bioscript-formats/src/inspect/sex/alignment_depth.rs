use std::{io::Read, path::Path};

use bioscript_core::RuntimeError;

use crate::{alignment, genotype::GenotypeLoadOptions};

use super::{
    DetectedKind, InferredSex, InspectOptions, SexDetectionConfidence, SexInference,
    unsupported_sex_inference,
};

const ALIGNMENT_SEX_WINDOW_LEN: i64 = 1000;
const ALIGNMENT_SEX_WINDOW_RECORD_CAP: usize = 400;
const ALIGNMENT_AUTOSOME_WINDOWS: &[(&str, i64)] = &[
    ("1", 50_000_000),
    ("2", 50_000_000),
    ("3", 50_000_000),
    ("4", 50_000_000),
    ("5", 50_000_000),
    ("6", 50_000_000),
    ("7", 50_000_000),
    ("8", 50_000_000),
];
const ALIGNMENT_X_NON_PAR_WINDOWS: &[i64] = &[
    10_000_000,
    20_000_000,
    40_000_000,
    60_000_000,
    80_000_000,
    100_000_000,
    120_000_000,
    140_000_000,
];
const ALIGNMENT_Y_WINDOWS: &[i64] = &[
    3_500_000, 8_000_000, 12_000_000, 16_000_000, 20_000_000, 24_000_000, 28_000_000, 40_000_000,
];

#[derive(Debug, Default)]
struct AlignmentSexStats {
    autosome_windows: usize,
    autosome_records: usize,
    x_windows: usize,
    x_records: usize,
    y_windows: usize,
    y_records: usize,
}

pub(crate) fn infer_sex_from_alignment_path(
    path: &Path,
    options: &InspectOptions,
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    if kind != DetectedKind::AlignmentCram {
        return Ok(unsupported_sex_inference());
    }
    let Some(reference_file) = options.reference_file.as_ref() else {
        return Ok(SexInference {
            sex: InferredSex::Unknown,
            confidence: SexDetectionConfidence::Low,
            method: "alignment_y_x_coverage".to_owned(),
            evidence: vec!["CRAM sex detection requires --reference-file".to_owned()],
        });
    };

    let load_options = GenotypeLoadOptions {
        input_index: options.input_index.clone(),
        reference_file: options.reference_file.clone(),
        reference_index: options.reference_index.clone(),
        allow_reference_md5_mismatch: true,
        ..GenotypeLoadOptions::default()
    };

    let stats = sample_alignment_sex_windows(path, &load_options, reference_file)?;
    let autosome_mean = mean_records(stats.autosome_records, stats.autosome_windows);
    let x_mean = mean_records(stats.x_records, stats.x_windows);
    let y_mean = mean_records(stats.y_records, stats.y_windows);
    let x_ratio = ratio_to_autosome(x_mean, autosome_mean);
    let y_ratio = ratio_to_autosome(y_mean, autosome_mean);

    let (sex, confidence) = if autosome_mean < 5.0 {
        (InferredSex::Unknown, SexDetectionConfidence::Low)
    } else if x_ratio >= 0.75 && y_ratio < 0.15 {
        (InferredSex::Female, SexDetectionConfidence::High)
    } else if x_ratio < 0.75 && y_ratio >= 0.08 {
        (InferredSex::Male, SexDetectionConfidence::High)
    } else if x_ratio >= 0.75 && y_ratio < 0.25 {
        (InferredSex::Female, SexDetectionConfidence::Medium)
    } else if x_ratio < 0.85 && y_ratio >= 0.03 {
        (InferredSex::Male, SexDetectionConfidence::Medium)
    } else {
        (InferredSex::Unknown, SexDetectionConfidence::Low)
    };

    let mut evidence = vec![
        format!("autosome_windows={}", stats.autosome_windows),
        format!("autosome_records={}", stats.autosome_records),
        format!("x_windows={}", stats.x_windows),
        format!("x_records={}", stats.x_records),
        format!("y_windows={}", stats.y_windows),
        format!("y_records={}", stats.y_records),
        format!("autosome_mean_records={autosome_mean:.2}"),
        format!("x_mean_records={x_mean:.2}"),
        format!("y_mean_records={y_mean:.2}"),
        format!("x_to_autosome_ratio={x_ratio:.3}"),
        format!("y_to_autosome_ratio={y_ratio:.3}"),
    ];
    if options.input_index.is_none() {
        evidence.push("CRAM sex detection ran without explicit --input-index".to_owned());
    }

    Ok(SexInference {
        sex,
        confidence,
        method: "alignment_autosome_x_y_depth_ratio".to_owned(),
        evidence,
    })
}

fn sample_alignment_sex_windows(
    path: &Path,
    options: &GenotypeLoadOptions,
    reference_file: &Path,
) -> Result<AlignmentSexStats, RuntimeError> {
    let repository = alignment::build_reference_repository(reference_file)?;
    let mut reader = alignment::build_cram_indexed_reader_from_path(path, options, repository)?;
    let label = path.display().to_string();
    let mut stats = AlignmentSexStats::default();

    for (chrom, center) in ALIGNMENT_AUTOSOME_WINDOWS {
        stats.autosome_records += count_alignment_records_in_window(
            &mut reader,
            &label,
            chrom,
            *center,
            options.allow_reference_md5_mismatch,
        )?;
        stats.autosome_windows += 1;
    }
    for center in ALIGNMENT_X_NON_PAR_WINDOWS {
        stats.x_records += count_alignment_records_in_window(
            &mut reader,
            &label,
            "X",
            *center,
            options.allow_reference_md5_mismatch,
        )?;
        stats.x_windows += 1;
    }
    for center in ALIGNMENT_Y_WINDOWS {
        stats.y_records += count_alignment_records_in_window(
            &mut reader,
            &label,
            "Y",
            *center,
            options.allow_reference_md5_mismatch,
        )?;
        stats.y_windows += 1;
    }

    Ok(stats)
}

fn count_alignment_records_in_window<R: Read + std::io::Seek>(
    reader: &mut noodles::cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    chrom: &str,
    center: i64,
    allow_reference_md5_mismatch: bool,
) -> Result<usize, RuntimeError> {
    let half_window = ALIGNMENT_SEX_WINDOW_LEN / 2;
    let locus = bioscript_core::GenomicLocus {
        chrom: chrom.to_owned(),
        start: center.saturating_sub(half_window).max(1),
        end: center.saturating_add(half_window),
    };
    let mut count = 0usize;
    alignment::for_each_cram_record_with_reader_allow_md5_mismatch(
        reader,
        label,
        &locus,
        allow_reference_md5_mismatch,
        |record| {
            if !record.is_unmapped {
                count += 1;
            }
            Ok(count < ALIGNMENT_SEX_WINDOW_RECORD_CAP)
        },
    )?;
    Ok(count)
}

fn mean_records(records: usize, windows: usize) -> f64 {
    if windows == 0 {
        0.0
    } else {
        f64::from(u32::try_from(records).unwrap_or(u32::MAX))
            / f64::from(u32::try_from(windows).unwrap_or(u32::MAX))
    }
}

fn ratio_to_autosome(value: f64, autosome_mean: f64) -> f64 {
    if autosome_mean <= f64::EPSILON {
        0.0
    } else {
        value / autosome_mean
    }
}
