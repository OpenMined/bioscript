use std::{
    io::{BufRead, BufReader, Cursor, Read},
    path::Path,
};

use bioscript_core::RuntimeError;
use flate2::read::MultiGzDecoder;
use zip::ZipArchive;

use crate::genotype::{DelimitedColumnIndexes, Delimiter, detect_delimiter, parse_streaming_row};

use super::{DetectedKind, InspectOptions};
use crate::{alignment, genotype::GenotypeLoadOptions};

mod classify;

use classify::{classify_stats, supports_sex_detection, unsupported_sex_inference};

const MAX_SEX_DETECTION_LINES: usize = 50_000_000;
const MAX_ZIP_ENTRY_BYTES: u64 = 256 * 1024 * 1024;
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
const MALE_SPECIFIC_Y_MARKERS: &[&str] = &[
    "rs11575897",
    "rs2534636",
    "i3000043",
    "i3000045",
    "i4000162",
    "rs13303871",
    "rs35284970",
    "rs3895",
    "i4000120",
    "i4000121",
    "i4000123",
    "rs13447361",
    "rs2267801",
    "rs2267802",
    "rs9786142",
    "i4000099",
    "i4000174",
    "i4000095",
    "i4000052",
    "i4000102",
];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InferredSex {
    Male,
    Female,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SexDetectionConfidence {
    High,
    Medium,
    Low,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SexInference {
    pub sex: InferredSex,
    pub confidence: SexDetectionConfidence,
    pub method: String,
    pub evidence: Vec<String>,
}

#[derive(Debug, Default)]
struct SexStats {
    total_y_snps: usize,
    called_y_snps: usize,
    male_markers_found: usize,
    male_markers_called: usize,
    y_examples: Vec<String>,
    x_non_par_sites: usize,
    x_haploid_gt_sites: usize,
    x_diploid_gt_sites: usize,
    x_het_gt_sites: usize,
}

#[derive(Debug, Default)]
struct AlignmentSexStats {
    autosome_windows: usize,
    autosome_records: usize,
    x_windows: usize,
    x_records: usize,
    y_windows: usize,
    y_records: usize,
}

pub(crate) fn infer_sex_from_path(
    path: &Path,
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    if !supports_sex_detection(kind) {
        return Ok(unsupported_sex_inference());
    }
    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".zip") {
        let file = std::fs::File::open(path)
            .map_err(|err| RuntimeError::Io(format!("failed to open {}: {err}", path.display())))?;
        let mut archive = ZipArchive::new(file).map_err(|err| {
            RuntimeError::Io(format!("failed to read zip {}: {err}", path.display()))
        })?;
        let entry_name = select_sex_detection_zip_entry(&mut archive)?;
        let mut entry = archive.by_name(&entry_name).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open zip entry {entry_name} in {}: {err}",
                path.display()
            ))
        })?;
        let mut bytes = Vec::new();
        std::io::Read::by_ref(&mut entry)
            .take(MAX_ZIP_ENTRY_BYTES.saturating_add(1))
            .read_to_end(&mut bytes)
            .map_err(|err| {
                RuntimeError::Io(format!("failed to read zip entry {entry_name}: {err}"))
            })?;
        if u64::try_from(bytes.len()).unwrap_or(u64::MAX) > MAX_ZIP_ENTRY_BYTES {
            return Err(RuntimeError::InvalidArguments(format!(
                "zip entry {entry_name} exceeds sex detection limit of {MAX_ZIP_ENTRY_BYTES} bytes"
            )));
        }
        return infer_sex_from_bytes(&entry_name, &bytes, kind);
    }

    let file = std::fs::File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open {}: {err}", path.display())))?;
    if lower.ends_with(".vcf.gz") {
        return infer_sex_from_reader(BufReader::new(MultiGzDecoder::new(file)), kind);
    }
    infer_sex_from_reader(BufReader::new(file), kind)
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

pub(crate) fn infer_sex_from_bytes(
    name: &str,
    bytes: &[u8],
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    if !supports_sex_detection(kind) {
        return Ok(unsupported_sex_inference());
    }
    let lower = name.to_ascii_lowercase();
    if lower.ends_with(".vcf.gz") {
        return infer_sex_from_reader(
            BufReader::new(MultiGzDecoder::new(Cursor::new(bytes))),
            kind,
        );
    }
    infer_sex_from_reader(BufReader::new(Cursor::new(bytes)), kind)
}

pub(crate) fn infer_sex_from_zip_bytes(
    bytes: &[u8],
    selected_entry: &str,
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    let mut archive = ZipArchive::new(Cursor::new(bytes))
        .map_err(|err| RuntimeError::Io(format!("failed to read zip bytes: {err}")))?;
    let mut entry = archive.by_name(selected_entry).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open zip entry {selected_entry} from bytes: {err}"
        ))
    })?;
    let mut entry_bytes = Vec::new();
    Read::by_ref(&mut entry)
        .take(MAX_ZIP_ENTRY_BYTES.saturating_add(1))
        .read_to_end(&mut entry_bytes)
        .map_err(|err| {
            RuntimeError::Io(format!("failed to read zip entry {selected_entry}: {err}"))
        })?;
    if u64::try_from(entry_bytes.len()).unwrap_or(u64::MAX) > MAX_ZIP_ENTRY_BYTES {
        return Err(RuntimeError::InvalidArguments(format!(
            "zip entry {selected_entry} exceeds sex detection limit of {MAX_ZIP_ENTRY_BYTES} bytes"
        )));
    }
    infer_sex_from_bytes(selected_entry, &entry_bytes, kind)
}

pub(crate) fn infer_sex_from_text_lines(
    lines: &[String],
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    let mut stats = SexStats::default();
    let delimiter = detect_delimiter(lines);
    let mut column_indexes = None;
    let mut comment_header = None;
    for line in lines {
        update_stats_from_line(
            &mut stats,
            line,
            kind,
            delimiter,
            &mut column_indexes,
            &mut comment_header,
        )?;
    }
    Ok(classify_stats(&stats, kind))
}

fn infer_sex_from_reader<R: BufRead>(
    mut reader: R,
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    let mut stats = SexStats::default();
    let mut probe_lines = Vec::new();
    let mut line = String::new();
    for _ in 0..64 {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|err| RuntimeError::Io(format!("failed to scan sex markers: {err}")))?;
        if bytes == 0 {
            let delimiter = detect_delimiter(&probe_lines);
            let mut column_indexes = None;
            let mut comment_header = None;
            for probe_line in &probe_lines {
                update_stats_from_line(
                    &mut stats,
                    probe_line,
                    kind,
                    delimiter,
                    &mut column_indexes,
                    &mut comment_header,
                )?;
            }
            return Ok(classify_stats(&stats, kind));
        }
        probe_lines.push(line.trim_end_matches(['\n', '\r']).to_owned());
    }
    let delimiter = detect_delimiter(&probe_lines);
    let mut column_indexes = None;
    let mut comment_header = None;
    for probe_line in &probe_lines {
        update_stats_from_line(
            &mut stats,
            probe_line,
            kind,
            delimiter,
            &mut column_indexes,
            &mut comment_header,
        )?;
    }
    for _ in probe_lines.len()..MAX_SEX_DETECTION_LINES {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|err| RuntimeError::Io(format!("failed to scan sex markers: {err}")))?;
        if bytes == 0 {
            break;
        }
        update_stats_from_line(
            &mut stats,
            line.trim_end_matches(['\n', '\r']),
            kind,
            delimiter,
            &mut column_indexes,
            &mut comment_header,
        )?;
    }
    Ok(classify_stats(&stats, kind))
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
            *chrom,
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
        records as f64 / windows as f64
    }
}

fn ratio_to_autosome(value: f64, autosome_mean: f64) -> f64 {
    if autosome_mean <= f64::EPSILON {
        0.0
    } else {
        value / autosome_mean
    }
}

fn update_stats_from_line(
    stats: &mut SexStats,
    line: &str,
    kind: DetectedKind,
    delimiter: Delimiter,
    column_indexes: &mut Option<DelimitedColumnIndexes>,
    comment_header: &mut Option<Vec<String>>,
) -> Result<(), RuntimeError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(());
    }
    match kind {
        DetectedKind::Vcf => update_vcf_stats(stats, trimmed),
        DetectedKind::GenotypeText | DetectedKind::Unknown => {
            update_genotype_text_stats(stats, trimmed, delimiter, column_indexes, comment_header)?;
        }
        _ => {}
    }
    Ok(())
}

fn update_genotype_text_stats(
    stats: &mut SexStats,
    line: &str,
    delimiter: Delimiter,
    column_indexes: &mut Option<DelimitedColumnIndexes>,
    comment_header: &mut Option<Vec<String>>,
) -> Result<(), RuntimeError> {
    let Some(row) = parse_streaming_row(line, delimiter, column_indexes, comment_header)? else {
        return Ok(());
    };
    let rsid = row.rsid.as_deref().unwrap_or_default();
    let chrom = normalize_chrom(row.chrom.as_deref().unwrap_or_default());
    let genotype = row.genotype.as_str();
    if chrom != "Y" {
        return Ok(());
    }
    stats.total_y_snps += 1;
    let called = is_called_genotype_text(genotype);
    if called {
        stats.called_y_snps += 1;
        if stats.y_examples.len() < 5 {
            stats.y_examples.push(format!("{rsid}:{genotype}"));
        }
    }
    if MALE_SPECIFIC_Y_MARKERS.contains(&rsid) {
        stats.male_markers_found += 1;
        if called {
            stats.male_markers_called += 1;
        }
    }
    Ok(())
}

fn update_vcf_stats(stats: &mut SexStats, line: &str) {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 10 {
        return;
    }
    let chrom = normalize_chrom(fields[0]);
    let Ok(pos) = fields[1].parse::<u32>() else {
        return;
    };
    let gt = fields[9].split(':').next().unwrap_or_default();
    if chrom == "Y" {
        stats.total_y_snps += 1;
        if is_called_vcf_gt(gt) {
            stats.called_y_snps += 1;
        }
        return;
    }
    if chrom != "X" || !is_non_par_x(pos) || !is_called_vcf_gt(gt) {
        return;
    }
    stats.x_non_par_sites += 1;
    let allele_count = vcf_gt_allele_count(gt);
    if allele_count == 1 {
        stats.x_haploid_gt_sites += 1;
    } else if allele_count == 2 {
        stats.x_diploid_gt_sites += 1;
        if is_vcf_gt_het(gt) {
            stats.x_het_gt_sites += 1;
        }
    }
}

fn normalize_chrom(value: &str) -> String {
    value
        .trim()
        .trim_start_matches("chr")
        .trim_start_matches("CHR")
        .to_ascii_uppercase()
}

fn is_called_genotype_text(value: &str) -> bool {
    let value = value.trim();
    if value.is_empty() || matches!(value, "--" | "00" | "." | "./." | ".|.") {
        return false;
    }
    value
        .chars()
        .all(|ch| matches!(ch.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'))
}

fn is_called_vcf_gt(value: &str) -> bool {
    let value = value.trim();
    !value.is_empty()
        && !value.contains('.')
        && (value != "0" || matches!(value, "0" | "1" | "2" | "3"))
}

fn vcf_gt_allele_count(gt: &str) -> usize {
    gt.split(['/', '|'])
        .filter(|part| !part.is_empty() && *part != ".")
        .count()
}

fn is_vcf_gt_het(gt: &str) -> bool {
    let alleles: Vec<&str> = gt
        .split(['/', '|'])
        .filter(|part| !part.is_empty() && *part != ".")
        .collect();
    alleles.len() == 2 && alleles[0] != alleles[1]
}

fn is_non_par_x(pos: u32) -> bool {
    // Human GRCh38 non-PAR X used by bcftools +guess-ploidy.
    // GRCh37 differs slightly, but these bounds cover the common non-PAR body
    // and avoid both pseudoautosomal ends for this QC heuristic.
    (2_781_480..=154_931_043).contains(&pos)
}

fn select_sex_detection_zip_entry<R: std::io::Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
) -> Result<String, RuntimeError> {
    for idx in 0..archive.len() {
        let entry = archive.by_index(idx).map_err(|err| {
            RuntimeError::Io(format!("failed to inspect zip for sex detection: {err}"))
        })?;
        if entry.is_dir() || entry.name().starts_with("__MACOSX/") {
            continue;
        }
        let lower = entry.name().to_ascii_lowercase();
        if lower.ends_with(".txt")
            || lower.ends_with(".tsv")
            || lower.ends_with(".csv")
            || lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
        {
            return Ok(entry.name().to_owned());
        }
    }
    Err(RuntimeError::Unsupported(
        "zip archive does not contain a supported sex detection input".to_owned(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn y_fingerprint_detects_male_and_female_text_exports() {
        let male = [
            "rs11575897\tY\t1\tG".to_owned(),
            "rs2534636\tY\t2\tC".to_owned(),
            "i3000043\tY\t3\tG".to_owned(),
            "i3000045\tY\t4\tG".to_owned(),
            "i4000162\tY\t5\tT".to_owned(),
            "rs13303871\tY\t6\tG".to_owned(),
        ];
        let result = infer_sex_from_text_lines(&male, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Male);
        assert_eq!(result.confidence, SexDetectionConfidence::Medium);

        let female: Vec<String> = (0..1001)
            .map(|idx| format!("rs{idx}\tY\t{idx}\t--"))
            .collect();
        let result = infer_sex_from_text_lines(&female, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.confidence, SexDetectionConfidence::High);
    }

    #[test]
    fn y_fingerprint_uses_genotype_column_not_array_metrics() {
        let lines: Vec<String> = (0..1001)
            .map(|idx| format!("rs{idx}\tY\t{idx}\t--\t0\t0.5\t-4.0"))
            .chain([
                "rs11575897\tY\t2001\t--\t0\t0.5\t-4.2".to_owned(),
                "rs2534636\tY\t2002\t--\t0\t0.9\t-2.6".to_owned(),
            ])
            .collect();
        let result = infer_sex_from_text_lines(&lines, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.confidence, SexDetectionConfidence::High);
        assert!(result.evidence.iter().any(|item| item == "called_y_snps=0"));
    }

    #[test]
    fn y_fingerprint_respects_header_column_order() {
        let mut lines = vec!["# rsid\tchromosome\tposition\tgs\tbaf\tlrr\tgenotype".to_owned()];
        lines.extend((0..1001).map(|idx| format!("rs{idx}\tY\t{idx}\t0\t0.5\t-4.0\t--")));
        lines.push("rs11575897\tY\t2001\t0\t0.5\t-4.2\t--".to_owned());
        let result = infer_sex_from_text_lines(&lines, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.confidence, SexDetectionConfidence::High);
        assert!(result.evidence.iter().any(|item| item == "called_y_snps=0"));
    }

    #[test]
    fn y_fingerprint_treats_sparse_y_calls_without_male_markers_as_female() {
        let mut lines: Vec<String> = (0..5500)
            .map(|idx| format!("rs{idx}\tY\t{idx}\t--\t0\t0.5\t-4.0"))
            .collect();
        for idx in 0..900 {
            lines.push(format!("noise{idx}\tY\t{}\tAA\t0.2\t0\t0.0", idx + 6000));
        }
        for marker in [
            "rs11575897",
            "rs2534636",
            "rs35284970",
            "rs13447361",
            "rs2267801",
            "rs2267802",
            "rs9786142",
        ] {
            lines.push(format!("{marker}\tY\t9000\t--\t0\t0.5\t-4.0"));
        }
        let result = infer_sex_from_text_lines(&lines, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.confidence, SexDetectionConfidence::High);
    }

    #[test]
    fn vcf_non_par_x_gt_detects_diploid_het_signal() {
        let lines: Vec<String> = (0..100)
            .map(|idx| {
                let gt = if idx % 3 == 0 { "0/1" } else { "0/0" };
                format!("chrX\t{}\t.\tC\tT\t.\tPASS\t.\tGT\t{gt}", 3_000_000 + idx)
            })
            .collect();
        let result = infer_sex_from_text_lines(&lines, DetectedKind::Vcf).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.method, "vcf_non_par_x_gt_y_count");
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "x_het_gt_sites=34")
        );
    }

    #[test]
    fn vcf_non_par_x_gt_beats_low_y_call_presence_when_signals_conflict() {
        let mut lines: Vec<String> = (0..20)
            .map(|idx| format!("chrY\t{}\t.\tC\tT\t.\tPASS\t.\tGT\t1", 3_000_000 + idx))
            .collect();
        lines.extend((0..1000).map(|idx| {
            let gt = if idx % 2 == 0 { "0/1" } else { "0/0" };
            format!("chrX\t{}\t.\tC\tT\t.\tPASS\t.\tGT\t{gt}", 3_000_000 + idx)
        }));

        let result = infer_sex_from_text_lines(&lines, DetectedKind::Vcf).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.confidence, SexDetectionConfidence::High);
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "called_y_snps=20")
        );
    }

    #[test]
    fn vcf_strong_y_signal_beats_diploid_x_het_from_wgs_callers() {
        let mut lines: Vec<String> = (0..600)
            .map(|idx| format!("chrY\t{}\t.\tC\tT\t.\tPASS\t.\tGT\t1", 3_000_000 + idx))
            .collect();
        lines.extend((0..2000).map(|idx| {
            let gt = if idx % 10 < 3 { "0/1" } else { "0/0" };
            format!("chrX\t{}\t.\tC\tT\t.\tPASS\t.\tGT\t{gt}", 3_000_000 + idx)
        }));

        let result = infer_sex_from_text_lines(&lines, DetectedKind::Vcf).unwrap();
        assert_eq!(result.sex, InferredSex::Male);
        assert_eq!(result.confidence, SexDetectionConfidence::High);
        assert!(result.evidence.iter().any(|item| item == "x_het_pct=30.00"));
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "y_to_x_pct=30.00")
        );
    }
}
