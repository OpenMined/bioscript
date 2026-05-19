use std::{
    io::{BufRead, BufReader, Cursor, Read},
    path::Path,
};

use bioscript_core::RuntimeError;
use flate2::read::MultiGzDecoder;
use zip::ZipArchive;

use crate::genotype::{
    DelimitedColumnIndexes, Delimiter, GsgtParser, detect_delimiter, lines_look_like_gsgt,
    parse_streaming_row,
};

use super::{DetectedKind, InspectOptions};

mod alignment_depth;
mod calls;
mod classify;

pub use alignment_depth::infer_sex_from_alignment_reader;

pub(crate) use alignment_depth::infer_sex_from_alignment_path;
use calls::{
    genotype_allele_count, is_called_genotype_text, is_called_vcf_gt, is_genotype_text_het,
    is_non_par_x, is_vcf_gt_het, normalize_chrom, vcf_gt_allele_count,
};
use classify::{classify_stats, supports_sex_detection, unsupported_sex_inference};

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

pub fn infer_sex_from_named_reader<R: Read>(
    name: &str,
    reader: R,
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    if !supports_sex_detection(kind) {
        return Ok(unsupported_sex_inference());
    }
    if name.to_ascii_lowercase().ends_with(".vcf.gz") {
        return infer_sex_from_reader(BufReader::new(MultiGzDecoder::new(reader)), kind);
    }
    infer_sex_from_reader(BufReader::new(reader), kind)
}

pub(crate) fn infer_sex_from_bytes(
    name: &str,
    bytes: &[u8],
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    if !supports_sex_detection(kind) {
        return Ok(unsupported_sex_inference());
    }
    infer_sex_from_named_reader(name, Cursor::new(bytes), kind)
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

pub fn infer_sex_from_text_lines(
    lines: &[String],
    kind: DetectedKind,
) -> Result<SexInference, RuntimeError> {
    let mut stats = SexStats::default();
    let delimiter = detect_delimiter(lines);
    let is_gsgt = lines_look_like_gsgt(lines);
    let mut gsgt = None;
    let mut column_indexes = None;
    let mut comment_header = None;
    for line in lines {
        update_stats_from_line(
            &mut stats,
            line,
            kind,
            delimiter,
            is_gsgt,
            &mut gsgt,
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
    // Treat any I/O error mid-stream (e.g. truncated bgzf head when the
    // wasm caller only loaded the first N MiB) as end-of-data: classify
    // whatever we got rather than failing the whole inspection. The CLI
    // path streams the full file so this is effectively unchanged for it.
    for _ in 0..64 {
        line.clear();
        let bytes = reader.read_line(&mut line).unwrap_or_default();
        if bytes == 0 {
            let delimiter = detect_delimiter(&probe_lines);
            let is_gsgt = lines_look_like_gsgt(&probe_lines);
            let mut gsgt = None;
            let mut column_indexes = None;
            let mut comment_header = None;
            for probe_line in &probe_lines {
                update_stats_from_line(
                    &mut stats,
                    probe_line,
                    kind,
                    delimiter,
                    is_gsgt,
                    &mut gsgt,
                    &mut column_indexes,
                    &mut comment_header,
                )?;
            }
            return Ok(classify_stats(&stats, kind));
        }
        probe_lines.push(line.trim_end_matches(['\n', '\r']).to_owned());
    }
    let delimiter = detect_delimiter(&probe_lines);
    let is_gsgt = lines_look_like_gsgt(&probe_lines);
    let mut gsgt = None;
    let mut column_indexes = None;
    let mut comment_header = None;
    for probe_line in &probe_lines {
        update_stats_from_line(
            &mut stats,
            probe_line,
            kind,
            delimiter,
            is_gsgt,
            &mut gsgt,
            &mut column_indexes,
            &mut comment_header,
        )?;
    }
    for _ in probe_lines.len()..MAX_SEX_DETECTION_LINES {
        line.clear();
        let bytes = reader.read_line(&mut line).unwrap_or_default();
        if bytes == 0 {
            break;
        }
        update_stats_from_line(
            &mut stats,
            line.trim_end_matches(['\n', '\r']),
            kind,
            delimiter,
            is_gsgt,
            &mut gsgt,
            &mut column_indexes,
            &mut comment_header,
        )?;
    }
    Ok(classify_stats(&stats, kind))
}

#[allow(clippy::too_many_arguments)]
fn update_stats_from_line(
    stats: &mut SexStats,
    line: &str,
    kind: DetectedKind,
    delimiter: Delimiter,
    is_gsgt: bool,
    gsgt: &mut Option<GsgtParser>,
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
            update_genotype_text_stats(
                stats,
                trimmed,
                delimiter,
                is_gsgt,
                gsgt,
                column_indexes,
                comment_header,
            )?;
        }
        _ => {}
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn update_genotype_text_stats(
    stats: &mut SexStats,
    line: &str,
    delimiter: Delimiter,
    is_gsgt: bool,
    gsgt: &mut Option<GsgtParser>,
    column_indexes: &mut Option<DelimitedColumnIndexes>,
    comment_header: &mut Option<Vec<String>>,
) -> Result<(), RuntimeError> {
    // GSGT SNP-array exports must normalize into the same
    // rsid/chrom/pos/genotype row as every other SNP text export so the
    // identical X/Y fingerprint heuristic applies regardless of source.
    let row = if is_gsgt {
        match gsgt.get_or_insert_with(GsgtParser::new).consume(line)? {
            Some(row) => row,
            None => return Ok(()),
        }
    } else {
        let Some(row) = parse_streaming_row(line, delimiter, column_indexes, comment_header)?
        else {
            return Ok(());
        };
        row
    };
    let rsid = row.rsid.as_deref().unwrap_or_default();
    let chrom = normalize_chrom(row.chrom.as_deref().unwrap_or_default());
    let genotype = row.genotype.as_str();
    if chrom == "X" {
        let Some(position) = row.position.and_then(|pos| u32::try_from(pos).ok()) else {
            return Ok(());
        };
        if !is_non_par_x(position) || !is_called_genotype_text(genotype) {
            return Ok(());
        }
        stats.x_non_par_sites += 1;
        let allele_count = genotype_allele_count(genotype);
        if allele_count == 1 {
            stats.x_haploid_gt_sites += 1;
        } else if allele_count == 2 {
            stats.x_diploid_gt_sites += 1;
            if is_genotype_text_het(genotype) {
                stats.x_het_gt_sites += 1;
            }
        }
        return Ok(());
    }
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
    use std::fmt::Write as _;
    use std::io::Write as _;

    fn zip_bytes(entries: &[(&str, &str)]) -> Vec<u8> {
        let cursor = Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        let options = zip::write::SimpleFileOptions::default();
        for (name, body) in entries {
            if name.ends_with('/') {
                writer.add_directory(*name, options).unwrap();
            } else {
                writer.start_file(*name, options).unwrap();
                writer.write_all(body.as_bytes()).unwrap();
            }
        }
        writer.finish().unwrap().into_inner()
    }

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
    fn gsgt_snp_array_goes_through_the_same_y_fingerprint() {
        // A male SNP array wrapped in the Illumina GSGT [Header]/[Data]
        // layout must reach the identical X/Y fingerprint and classify the
        // same as any other text export — source format is irrelevant.
        let mut gsgt = vec![
            "[Header]".to_owned(),
            "GSGT Version\t2.0.5".to_owned(),
            "[Data]".to_owned(),
            "SNP Name\tSNP\tChr\tPosition\tAllele1 - Plus\tAllele2 - Plus".to_owned(),
        ];
        let mut flat = Vec::new();
        for i in 0..600 {
            let rsid = format!("rs{}", 700_000 + i);
            gsgt.push(format!("{rsid}\t[A/G]\tY\t{}\tG\tG", i + 1));
            flat.push(format!("{rsid}\tY\t{}\tG", i + 1));
        }

        let gsgt_result = infer_sex_from_text_lines(&gsgt, DetectedKind::GenotypeText).unwrap();
        let flat_result = infer_sex_from_text_lines(&flat, DetectedKind::GenotypeText).unwrap();

        assert_eq!(gsgt_result.sex, InferredSex::Male);
        assert_eq!(gsgt_result.method, "snp_array_x_y_fingerprint");
        // Identical heuristic outcome regardless of which export format.
        assert_eq!(gsgt_result.sex, flat_result.sex);
        assert_eq!(gsgt_result.confidence, flat_result.confidence);
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
    fn snp_array_numeric_x_and_y_chromosomes_feed_sex_inference() {
        let mut lines = vec!["rsid\tchromosome\tposition\tallele1\tallele2".to_owned()];
        lines.extend((0..2000).map(|idx| {
            let (a, b) = if idx % 4 == 0 { ("A", "G") } else { ("A", "A") };
            format!("rsX{idx}\t23\t{}\t{a}\t{b}", 3_000_000 + idx)
        }));
        lines.extend((0..1000).map(|idx| format!("rsY{idx}\t24\t{}\t0\t0", 3_000_000 + idx)));

        let result = infer_sex_from_text_lines(&lines, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert_eq!(result.confidence, SexDetectionConfidence::Medium);
        assert_eq!(result.method, "snp_array_x_y_fingerprint");
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "x_het_gt_sites=500")
        );
    }

    #[test]
    fn snp_array_haploid_non_par_x_detects_male_without_y_rows() {
        let mut lines = vec!["rsid\tchromosome\tposition\tgenotype".to_owned()];
        lines.extend((0..2000).map(|idx| format!("rsX{idx}\tX\t{}\tA", 3_000_000 + idx)));
        lines.extend((0..20).map(|idx| format!("rsDiploidTail{idx}\tX\t{}\tAG", 4_000_000 + idx)));

        let result = infer_sex_from_text_lines(&lines, DetectedKind::GenotypeText).unwrap();
        assert_eq!(result.sex, InferredSex::Male);
        assert_eq!(result.confidence, SexDetectionConfidence::Medium);
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "x_haploid_gt_sites=2000")
        );
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

    #[test]
    fn sex_inference_bytes_and_zip_paths_cover_entry_selection_and_unsupported_kinds() {
        let text = "rsid\tchromosome\tposition\tgenotype\nrs11575897\tY\t1\tG\n";
        let unsupported =
            infer_sex_from_bytes("sample.txt", text.as_bytes(), DetectedKind::ReferenceFasta)
                .unwrap();
        assert_eq!(unsupported.sex, InferredSex::Unknown);
        assert_eq!(unsupported.method, "unsupported_source_type");

        let result =
            infer_sex_from_bytes("sample.txt", text.as_bytes(), DetectedKind::GenotypeText)
                .unwrap();
        assert_eq!(result.method, "snp_array_x_y_fingerprint");

        let archive = zip_bytes(&[
            ("__MACOSX/._sample.txt", "ignored"),
            ("notes.md", "ignored"),
            ("nested/sample.txt", text),
        ]);
        let result =
            infer_sex_from_zip_bytes(&archive, "nested/sample.txt", DetectedKind::GenotypeText)
                .unwrap();
        assert_eq!(result.method, "snp_array_x_y_fingerprint");

        let err = infer_sex_from_zip_bytes(&archive, "missing.txt", DetectedKind::GenotypeText)
            .unwrap_err();
        assert!(
            err.to_string()
                .contains("failed to open zip entry missing.txt")
        );

        let bad_zip =
            infer_sex_from_zip_bytes(b"not a zip", "sample.txt", DetectedKind::GenotypeText)
                .unwrap_err();
        assert!(bad_zip.to_string().contains("failed to read zip bytes"));

        let mut zip = ZipArchive::new(Cursor::new(archive)).unwrap();
        assert_eq!(
            select_sex_detection_zip_entry(&mut zip).unwrap(),
            "nested/sample.txt"
        );

        let unsupported_zip = zip_bytes(&[("docs/readme.md", "ignored")]);
        let mut zip = ZipArchive::new(Cursor::new(unsupported_zip)).unwrap();
        let err = select_sex_detection_zip_entry(&mut zip).unwrap_err();
        assert!(
            err.to_string()
                .contains("does not contain a supported sex detection input")
        );
    }

    #[test]
    fn sex_inference_reader_probe_and_late_stream_paths_handle_vcf_edges() {
        let mut text = String::new();
        text.push_str("chrY\tbad\t.\tC\tT\t.\tPASS\t.\tGT\t1\n");
        text.push_str("chrM\t1\t.\tC\tT\t.\tPASS\t.\tGT\t1\n");
        for idx in 0..70 {
            let gt = if idx % 2 == 0 { "0|1" } else { "0|0" };
            let _ = writeln!(
                text,
                "23\t{}\t.\tC\tT\t.\tPASS\t.\tGT\t{gt}:99",
                3_000_000 + idx
            );
        }
        text.push_str("24\t1\t.\tC\tT\t.\tPASS\t.\tGT\t.\n");
        text.push_str("chrX\t60000\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n");
        text.push_str("chrX\t155000000\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n");

        let result =
            infer_sex_from_bytes("sample.vcf", text.as_bytes(), DetectedKind::Vcf).unwrap();
        assert_eq!(result.sex, InferredSex::Female);
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "x_non_par_sites=70")
        );
        assert!(
            result
                .evidence
                .iter()
                .any(|item| item == "x_het_gt_sites=35")
        );
    }
}
