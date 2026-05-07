use std::{
    io::{BufRead, BufReader, Cursor, Read},
    path::Path,
};

use bioscript_core::RuntimeError;
use flate2::read::MultiGzDecoder;
use zip::ZipArchive;

use super::{DetectedKind, split_fields};

const MAX_SEX_DETECTION_LINES: usize = 50_000_000;
const MAX_ZIP_ENTRY_BYTES: u64 = 256 * 1024 * 1024;
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
    for line in lines {
        update_stats_from_line(&mut stats, line, kind);
    }
    Ok(classify_stats(&stats, kind))
}

fn infer_sex_from_reader<R: BufRead>(
    mut reader: R,
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    let mut stats = SexStats::default();
    let mut line = String::new();
    for _ in 0..MAX_SEX_DETECTION_LINES {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|err| RuntimeError::Io(format!("failed to scan sex markers: {err}")))?;
        if bytes == 0 {
            break;
        }
        update_stats_from_line(&mut stats, line.trim_end_matches(['\n', '\r']), kind);
    }
    Ok(classify_stats(&stats, kind))
}

fn update_stats_from_line(stats: &mut SexStats, line: &str, kind: DetectedKind) {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("//") {
        return;
    }
    match kind {
        DetectedKind::Vcf => update_vcf_stats(stats, trimmed),
        DetectedKind::GenotypeText | DetectedKind::Unknown => {
            update_genotype_text_stats(stats, trimmed);
        }
        _ => {}
    }
}

fn update_genotype_text_stats(stats: &mut SexStats, line: &str) {
    let fields = split_fields(line);
    if fields.len() < 4 {
        return;
    }
    let rsid = fields[0].trim();
    let chrom = normalize_chrom(fields.get(1).map(String::as_str).unwrap_or_default());
    let genotype = genotype_text_field(&fields);
    if chrom != "Y" {
        return;
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

fn classify_stats(stats: &SexStats, kind: DetectedKind) -> SexInference {
    match kind {
        DetectedKind::Vcf => classify_vcf_stats(stats),
        DetectedKind::GenotypeText | DetectedKind::Unknown => classify_y_fingerprint_stats(stats),
        _ => unsupported_sex_inference(),
    }
}

fn supports_sex_detection(kind: DetectedKind) -> bool {
    matches!(
        kind,
        DetectedKind::Vcf | DetectedKind::GenotypeText | DetectedKind::Unknown
    )
}

fn unsupported_sex_inference() -> SexInference {
    SexInference {
        sex: InferredSex::Unknown,
        confidence: SexDetectionConfidence::Low,
        method: "unsupported_source_type".to_owned(),
        evidence: vec!["sex detection currently supports genotype text and VCF inputs".to_owned()],
    }
}

fn classify_y_fingerprint_stats(stats: &SexStats) -> SexInference {
    let (sex, confidence) = if stats.called_y_snps > 500 {
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
    let (sex, confidence) = if stats.called_y_snps > 500
        || (stats.x_haploid_gt_sites >= 20 && stats.x_diploid_gt_sites == 0)
    {
        (InferredSex::Male, SexDetectionConfidence::High)
    } else if stats.x_non_par_sites >= 50
        && stats.x_diploid_gt_sites > 0
        && stats.x_het_gt_sites * 100 / stats.x_diploid_gt_sites.max(1) >= 2
    {
        (InferredSex::Female, SexDetectionConfidence::Medium)
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
            format!("called_y_snps={}", stats.called_y_snps),
        ],
    }
}

fn y_fingerprint_evidence(stats: &SexStats) -> Vec<String> {
    let mut evidence = vec![
        format!("total_y_snps={}", stats.total_y_snps),
        format!("called_y_snps={}", stats.called_y_snps),
        format!("male_markers_found={}", stats.male_markers_found),
        format!("male_markers_called={}", stats.male_markers_called),
    ];
    if !stats.y_examples.is_empty() {
        evidence.push(format!("y_examples={}", stats.y_examples.join(",")));
    }
    evidence
}

fn normalize_chrom(value: &str) -> String {
    value
        .trim()
        .trim_start_matches("chr")
        .trim_start_matches("CHR")
        .to_ascii_uppercase()
}

fn genotype_text_field(fields: &[String]) -> &str {
    if fields.len() >= 4 {
        fields[3].trim()
    } else {
        fields.last().map(String::as_str).unwrap_or_default().trim()
    }
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
    fn vcf_non_par_x_gt_detects_diploid_het_signal() {
        let lines = vec![
            "chrX\t3000000\t.\tC\tT\t.\tPASS\t.\tGT\t0/1".to_owned(),
            "chrX\t4000000\t.\tC\tT\t.\tPASS\t.\tGT\t0/0".to_owned(),
            "chrX\t5000000\t.\tC\tT\t.\tPASS\t.\tGT\t1/1".to_owned(),
        ];
        let result = infer_sex_from_text_lines(&lines, DetectedKind::Vcf).unwrap();
        assert_eq!(result.method, "vcf_non_par_x_gt_y_count");
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "x_het_gt_sites=1")
        );
    }
}
