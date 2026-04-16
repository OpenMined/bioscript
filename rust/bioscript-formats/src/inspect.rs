use std::{
    fs::File,
    io::{BufRead, BufReader, Cursor, Read},
    path::{Path, PathBuf},
};

#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

// std::time::Instant::now() panics on wasm32-unknown-unknown ("time not
// implemented on this platform"). duration_ms is diagnostic-only, so on wasm
// we stub the timer instead of pulling in a perf/date shim crate.
#[cfg(target_arch = "wasm32")]
struct Instant;

#[cfg(target_arch = "wasm32")]
impl Instant {
    fn now() -> Self {
        Self
    }
    fn elapsed(&self) -> StubDuration {
        StubDuration
    }
}

#[cfg(target_arch = "wasm32")]
struct StubDuration;

#[cfg(target_arch = "wasm32")]
impl StubDuration {
    fn as_millis(&self) -> u128 {
        0
    }
}

use bioscript_core::{Assembly, RuntimeError};
use noodles::bgzf;
use zip::ZipArchive;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FileContainer {
    Plain,
    Zip,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedKind {
    GenotypeText,
    Vcf,
    AlignmentCram,
    AlignmentBam,
    ReferenceFasta,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectionConfidence {
    Authoritative,
    StrongHeuristic,
    WeakHeuristic,
    Unknown,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SourceMetadata {
    pub vendor: Option<String>,
    pub platform_version: Option<String>,
    pub confidence: DetectionConfidence,
    pub evidence: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct InspectOptions {
    pub input_index: Option<PathBuf>,
    pub reference_file: Option<PathBuf>,
    pub reference_index: Option<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileInspection {
    pub path: PathBuf,
    pub container: FileContainer,
    pub detected_kind: DetectedKind,
    pub confidence: DetectionConfidence,
    pub source: Option<SourceMetadata>,
    pub assembly: Option<Assembly>,
    pub phased: Option<bool>,
    pub selected_entry: Option<String>,
    pub has_index: Option<bool>,
    pub index_path: Option<PathBuf>,
    pub reference_matches: Option<bool>,
    pub evidence: Vec<String>,
    pub warnings: Vec<String>,
    pub duration_ms: u128,
}

impl FileInspection {
    #[must_use]
    pub fn render_text(&self) -> String {
        let mut lines = Vec::new();
        lines.push(format!("path\t{}", self.path.display()));
        lines.push(format!("container\t{}", render_container(self.container)));
        lines.push(format!("kind\t{}", render_kind(self.detected_kind)));
        lines.push(format!(
            "confidence\t{}",
            render_confidence(self.confidence)
        ));
        lines.push(format!("assembly\t{}", render_assembly(self.assembly)));
        lines.push(format!("phased\t{}", render_bool(self.phased)));
        lines.push(format!(
            "selected_entry\t{}",
            self.selected_entry.as_deref().unwrap_or("")
        ));
        lines.push(format!("has_index\t{}", render_bool(self.has_index)));
        lines.push(format!(
            "index_path\t{}",
            self.index_path
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default()
        ));
        lines.push(format!(
            "reference_matches\t{}",
            render_bool(self.reference_matches)
        ));
        if let Some(source) = &self.source {
            lines.push(format!(
                "vendor\t{}",
                source.vendor.as_deref().unwrap_or_default()
            ));
            lines.push(format!(
                "platform_version\t{}",
                source.platform_version.as_deref().unwrap_or_default()
            ));
            lines.push(format!(
                "source_confidence\t{}",
                render_confidence(source.confidence)
            ));
            lines.push(format!("source_evidence\t{}", source.evidence.join(" | ")));
        } else {
            lines.push("vendor\t".to_owned());
            lines.push("platform_version\t".to_owned());
            lines.push("source_confidence\t".to_owned());
            lines.push("source_evidence\t".to_owned());
        }
        lines.push(format!("evidence\t{}", self.evidence.join(" | ")));
        lines.push(format!("warnings\t{}", self.warnings.join(" | ")));
        lines.push(format!("duration_ms\t{}", self.duration_ms));
        lines.join("\n")
    }
}

/// Classify a file from in-memory bytes. Mirrors `inspect_file` but sources
/// its sample lines / zip entries from a byte buffer instead of the
/// filesystem. Needed by wasm targets where `std::fs` isn't available.
///
/// `name` is used for extension-based detection (.cram / .bam / .fa / .zip /
/// .vcf.gz) and vendor sniffing from the filename. `bytes` is the file
/// contents; for zips we scan the central directory out of these bytes.
pub fn inspect_bytes(
    name: &str,
    bytes: &[u8],
    options: &InspectOptions,
) -> Result<FileInspection, RuntimeError> {
    let started = Instant::now();
    let lower = name.to_ascii_lowercase();
    let mut evidence = Vec::new();
    let mut warnings = Vec::new();
    let path = Path::new(name);

    if lower.ends_with(".zip") {
        let selected_entry = select_zip_entry_from_bytes(bytes)?;
        let sample_lines = read_zip_sample_lines_from_bytes(bytes, &selected_entry)?;
        let mut inspection = inspect_from_textual_sample(
            path,
            FileContainer::Zip,
            &selected_entry,
            &sample_lines,
            options,
        );
        inspection.duration_ms = started.elapsed().as_millis();
        return Ok(inspection);
    }

    let detected_kind = if lower.ends_with(".cram") {
        evidence.push("extension .cram".to_owned());
        DetectedKind::AlignmentCram
    } else if lower.ends_with(".bam") {
        evidence.push("extension .bam".to_owned());
        DetectedKind::AlignmentBam
    } else if is_reference_path(path) {
        evidence.push("reference fasta extension".to_owned());
        DetectedKind::ReferenceFasta
    } else {
        let sample_lines = read_plain_sample_lines_from_bytes(&lower, bytes)?;
        let sample_lower = sample_lines.join("\n").to_ascii_lowercase();
        if looks_like_vcf_lines(&sample_lines) {
            evidence.push("vcf header markers".to_owned());
            DetectedKind::Vcf
        } else if looks_like_genotype_text(&sample_lines) {
            if sample_lower.contains("rsid") || sample_lower.contains("allele1") {
                evidence.push("genotype-like sampled rows and headers".to_owned());
            } else {
                evidence.push("genotype-like sampled rows".to_owned());
            }
            DetectedKind::GenotypeText
        } else {
            warnings.push("file did not match known textual heuristics".to_owned());
            DetectedKind::Unknown
        }
    };

    let sample_lines = match detected_kind {
        DetectedKind::AlignmentCram | DetectedKind::AlignmentBam | DetectedKind::ReferenceFasta => {
            Vec::new()
        }
        _ => read_plain_sample_lines_from_bytes(&lower, bytes)?,
    };
    let source = detect_source(&lower, &sample_lines, detected_kind);
    let assembly = detect_assembly(&lower, &sample_lines);
    let phased = (detected_kind == DetectedKind::Vcf)
        .then(|| detect_vcf_phasing(&sample_lines))
        .flatten();
    // Index discovery is filesystem-only; wasm callers pass indexes separately
    // through `InspectOptions.input_index` / `reference_index` and we surface
    // whichever is provided.
    let has_index = options
        .input_index
        .as_ref()
        .or(options.reference_index.as_ref())
        .map(|_| true);
    let index_path = options
        .input_index
        .clone()
        .or_else(|| options.reference_index.clone());
    let confidence = classify_confidence(detected_kind, &sample_lines, source.as_ref());

    Ok(FileInspection {
        path: path.to_path_buf(),
        container: FileContainer::Plain,
        detected_kind,
        confidence,
        source,
        assembly,
        phased,
        selected_entry: None,
        has_index,
        index_path,
        reference_matches: None,
        evidence,
        warnings,
        duration_ms: started.elapsed().as_millis(),
    })
}

fn read_plain_sample_lines_from_bytes(
    lower_name: &str,
    bytes: &[u8],
) -> Result<Vec<String>, RuntimeError> {
    if lower_name.ends_with(".vcf.gz") {
        return read_sample_lines_from_reader(BufReader::new(bgzf::io::Reader::new(Cursor::new(
            bytes,
        ))));
    }
    read_sample_lines_from_reader(BufReader::new(Cursor::new(bytes)))
}

fn read_zip_sample_lines_from_bytes(
    bytes: &[u8],
    selected_entry: &str,
) -> Result<Vec<String>, RuntimeError> {
    let mut archive = ZipArchive::new(Cursor::new(bytes))
        .map_err(|err| RuntimeError::Io(format!("failed to read zip bytes: {err}")))?;
    let mut entry = archive.by_name(selected_entry).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open zip entry {selected_entry} from bytes: {err}"
        ))
    })?;
    if selected_entry.to_ascii_lowercase().ends_with(".vcf.gz") {
        let mut inner = Vec::new();
        entry.read_to_end(&mut inner).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read compressed zip entry {selected_entry}: {err}"
            ))
        })?;
        let reader = bgzf::io::Reader::new(Cursor::new(inner));
        return read_sample_lines_from_reader(BufReader::new(reader));
    }
    read_sample_lines_from_reader(BufReader::new(entry))
}

fn select_zip_entry_from_bytes(bytes: &[u8]) -> Result<String, RuntimeError> {
    let mut archive = ZipArchive::new(Cursor::new(bytes))
        .map_err(|err| RuntimeError::Io(format!("failed to read zip bytes: {err}")))?;
    let mut fallback = None;
    for idx in 0..archive.len() {
        let entry = archive
            .by_index(idx)
            .map_err(|err| RuntimeError::Io(format!("failed to inspect zip bytes: {err}")))?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        if name.starts_with("__MACOSX/") {
            continue;
        }
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
            || lower.ends_with(".txt")
            || lower.ends_with(".tsv")
            || lower.ends_with(".csv")
        {
            return Ok(name);
        }
        if fallback.is_none() {
            fallback = Some(name);
        }
    }
    fallback.ok_or_else(|| {
        RuntimeError::Unsupported("zip archive does not contain a supported file".to_owned())
    })
}

pub fn inspect_file(path: &Path, options: &InspectOptions) -> Result<FileInspection, RuntimeError> {
    let started = Instant::now();
    let lower = path.to_string_lossy().to_ascii_lowercase();
    let mut evidence = Vec::new();
    let mut warnings = Vec::new();

    if lower.ends_with(".zip") {
        let selected_entry = select_zip_entry(path)?;
        let sample_lines = read_zip_sample_lines(path, &selected_entry)?;
        let mut inspection = inspect_from_textual_sample(
            path,
            FileContainer::Zip,
            &selected_entry,
            &sample_lines,
            options,
        );
        inspection.duration_ms = started.elapsed().as_millis();
        return Ok(inspection);
    }

    let detected_kind = if lower.ends_with(".cram") {
        evidence.push("extension .cram".to_owned());
        DetectedKind::AlignmentCram
    } else if lower.ends_with(".bam") {
        evidence.push("extension .bam".to_owned());
        DetectedKind::AlignmentBam
    } else if is_reference_path(path) {
        evidence.push("reference fasta extension".to_owned());
        DetectedKind::ReferenceFasta
    } else {
        let sample_lines = read_plain_sample_lines(path)?;
        let sample_lower = sample_lines.join("\n").to_ascii_lowercase();
        if looks_like_vcf_lines(&sample_lines) {
            evidence.push("vcf header markers".to_owned());
            DetectedKind::Vcf
        } else if looks_like_genotype_text(&sample_lines) {
            if sample_lower.contains("rsid") || sample_lower.contains("allele1") {
                evidence.push("genotype-like sampled rows and headers".to_owned());
            } else {
                evidence.push("genotype-like sampled rows".to_owned());
            }
            DetectedKind::GenotypeText
        } else {
            warnings.push("file did not match known textual heuristics".to_owned());
            DetectedKind::Unknown
        }
    };

    let sample_lines = match detected_kind {
        DetectedKind::AlignmentCram | DetectedKind::AlignmentBam | DetectedKind::ReferenceFasta => {
            Vec::new()
        }
        _ => read_plain_sample_lines(path)?,
    };
    let source = detect_source(&lower, &sample_lines, detected_kind);
    let assembly = detect_assembly(&lower, &sample_lines);
    let phased = (detected_kind == DetectedKind::Vcf)
        .then(|| detect_vcf_phasing(&sample_lines))
        .flatten();
    let (has_index, index_path) = detect_index(path, detected_kind, options);
    let confidence = classify_confidence(detected_kind, &sample_lines, source.as_ref());

    Ok(FileInspection {
        path: path.to_path_buf(),
        container: FileContainer::Plain,
        detected_kind,
        confidence,
        source,
        assembly,
        phased,
        selected_entry: None,
        has_index,
        index_path,
        reference_matches: None,
        evidence,
        warnings,
        duration_ms: started.elapsed().as_millis(),
    })
}

fn inspect_from_textual_sample(
    path: &Path,
    container: FileContainer,
    selected_entry: &str,
    sample_lines: &[String],
    options: &InspectOptions,
) -> FileInspection {
    let lower = selected_entry.to_ascii_lowercase();
    let detected_kind = if lower.ends_with(".vcf")
        || lower.ends_with(".vcf.gz")
        || looks_like_vcf_lines(sample_lines)
    {
        DetectedKind::Vcf
    } else if looks_like_genotype_text(sample_lines) {
        DetectedKind::GenotypeText
    } else {
        DetectedKind::Unknown
    };
    let path_lower = path.to_string_lossy().to_ascii_lowercase();
    let combined_name = format!("{path_lower}\n{lower}");
    let source = detect_source(&combined_name, sample_lines, detected_kind);
    let assembly = detect_assembly(&combined_name, sample_lines);
    let phased = (detected_kind == DetectedKind::Vcf)
        .then(|| detect_vcf_phasing(sample_lines))
        .flatten();
    let confidence = classify_confidence(detected_kind, sample_lines, source.as_ref());
    let mut evidence = vec![format!("selected zip entry {selected_entry}")];
    if detected_kind == DetectedKind::Vcf {
        evidence.push("vcf entry markers".to_owned());
    } else if detected_kind == DetectedKind::GenotypeText {
        evidence.push("genotype-like sampled rows".to_owned());
    }
    let (has_index, index_path) = detect_index(path, detected_kind, options);

    FileInspection {
        path: path.to_path_buf(),
        container,
        detected_kind,
        confidence,
        source,
        assembly,
        phased,
        selected_entry: Some(selected_entry.to_owned()),
        has_index,
        index_path,
        reference_matches: None,
        evidence,
        warnings: Vec::new(),
        duration_ms: 0,
    }
}

fn read_plain_sample_lines(path: &Path) -> Result<Vec<String>, RuntimeError> {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open {}: {err}", path.display())))?;
    if lower.ends_with(".vcf.gz") {
        return read_sample_lines_from_reader(BufReader::new(bgzf::io::Reader::new(file)));
    }
    read_sample_lines_from_reader(BufReader::new(file))
}

fn read_zip_sample_lines(path: &Path, selected_entry: &str) -> Result<Vec<String>, RuntimeError> {
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open zip {}: {err}", path.display())))?;
    let mut archive = ZipArchive::new(file)
        .map_err(|err| RuntimeError::Io(format!("failed to read zip {}: {err}", path.display())))?;
    let mut entry = archive.by_name(selected_entry).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open zip entry {selected_entry} in {}: {err}",
            path.display()
        ))
    })?;

    if selected_entry.to_ascii_lowercase().ends_with(".vcf.gz") {
        let mut bytes = Vec::new();
        entry.read_to_end(&mut bytes).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read compressed zip entry {selected_entry} in {}: {err}",
                path.display()
            ))
        })?;
        let reader = bgzf::io::Reader::new(Cursor::new(bytes));
        return read_sample_lines_from_reader(BufReader::new(reader));
    }

    read_sample_lines_from_reader(BufReader::new(entry))
}

fn read_sample_lines_from_reader<R: BufRead>(mut reader: R) -> Result<Vec<String>, RuntimeError> {
    let mut out = Vec::new();
    let mut buf = String::new();
    for _ in 0..64 {
        buf.clear();
        let bytes = reader
            .read_line(&mut buf)
            .map_err(|err| RuntimeError::Io(format!("failed to read sample lines: {err}")))?;
        if bytes == 0 {
            break;
        }
        out.push(buf.trim_end_matches(['\n', '\r']).to_owned());
    }
    Ok(out)
}

fn select_zip_entry(path: &Path) -> Result<String, RuntimeError> {
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open zip {}: {err}", path.display())))?;
    let mut archive = ZipArchive::new(file)
        .map_err(|err| RuntimeError::Io(format!("failed to read zip {}: {err}", path.display())))?;
    let mut fallback = None;
    for idx in 0..archive.len() {
        let entry = archive.by_index(idx).map_err(|err| {
            RuntimeError::Io(format!("failed to inspect zip {}: {err}", path.display()))
        })?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        if name.starts_with("__MACOSX/") {
            continue;
        }
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
            || lower.ends_with(".txt")
            || lower.ends_with(".tsv")
            || lower.ends_with(".csv")
        {
            return Ok(name);
        }
        if fallback.is_none() {
            fallback = Some(name);
        }
    }
    fallback.ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "zip archive {} does not contain a supported file",
            path.display()
        ))
    })
}

fn looks_like_vcf_lines(lines: &[String]) -> bool {
    lines.iter().any(|line| {
        let trimmed = line.trim_start();
        trimmed.starts_with("##fileformat=VCF") || trimmed.starts_with("#CHROM\t")
    })
}

fn looks_like_genotype_text(lines: &[String]) -> bool {
    let mut checked = 0usize;
    let mut valid = 0usize;
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') || trimmed.starts_with("//") {
            continue;
        }
        let fields = split_fields(trimmed);
        checked += 1;
        if matches_genotype_shape(&fields) {
            valid += 1;
        }
    }
    checked > 0 && valid * 10 >= checked * 7
}

fn split_fields(line: &str) -> Vec<String> {
    if line.contains('\t') {
        return line
            .split('\t')
            .map(|field| field.trim().to_owned())
            .collect();
    }
    if line.contains(',') {
        return line
            .split(',')
            .map(|field| field.trim().trim_matches('"').to_owned())
            .collect();
    }
    line.split_whitespace().map(str::to_owned).collect()
}

fn matches_genotype_shape(fields: &[String]) -> bool {
    if fields.len() < 4 {
        return false;
    }
    let rsid_like = fields[0].starts_with("rs") || fields[0].starts_with('i');
    if !rsid_like {
        return false;
    }
    let chr_idx = fields.iter().position(|field| is_valid_chromosome(field));
    let Some(chr_idx) = chr_idx else {
        return false;
    };
    for pos_idx in (chr_idx + 1)..fields.len() {
        if fields[pos_idx].parse::<u64>().is_err() {
            continue;
        }
        for field in fields.iter().skip(pos_idx + 1) {
            if is_valid_genotype(field) {
                return true;
            }
        }
        if pos_idx + 2 < fields.len()
            && is_valid_allele(&fields[pos_idx + 1])
            && is_valid_allele(&fields[pos_idx + 2])
        {
            return true;
        }
    }
    false
}

fn is_valid_chromosome(value: &str) -> bool {
    let trimmed = value.trim().trim_start_matches("chr");
    if let Ok(n) = trimmed.parse::<u8>() {
        return (1..=26).contains(&n);
    }
    matches!(
        trimmed.to_ascii_uppercase().as_str(),
        "X" | "Y" | "M" | "MT" | "XY"
    )
}

fn is_valid_genotype(value: &str) -> bool {
    let trimmed = value.trim().to_ascii_uppercase();
    if trimmed.is_empty() || trimmed.len() > 4 {
        return false;
    }
    trimmed
        .chars()
        .all(|ch| matches!(ch, 'A' | 'C' | 'G' | 'T' | 'I' | 'D' | '-' | '0'))
}

fn is_valid_allele(value: &str) -> bool {
    let trimmed = value.trim().to_ascii_uppercase();
    matches!(
        trimmed.as_str(),
        "A" | "C" | "G" | "T" | "I" | "D" | "-" | "0"
    )
}

fn detect_source(
    lower_name: &str,
    sample_lines: &[String],
    kind: DetectedKind,
) -> Option<SourceMetadata> {
    let header = sample_lines
        .iter()
        .filter(|line| line.starts_with('#') || line.starts_with("//"))
        .map(|line| line.to_ascii_lowercase())
        .collect::<Vec<_>>()
        .join("\n");
    let combined = format!("{lower_name}\n{header}");
    let normalized = combined.replace(['_', '-', '.'], " ");
    let mut evidence = Vec::new();
    let mut vendor = None;
    let mut platform_version = None;
    let mut confidence = DetectionConfidence::Unknown;

    if normalized.contains("genes for good") || normalized.contains("geneforgood") {
        vendor = Some("Genes for Good".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("Genes for Good header".to_owned());
        if let Some(version) = extract_token_after_marker(&header, "genes for good ") {
            platform_version = Some(version);
            evidence.push("Genes for Good version header".to_owned());
        }
    } else if normalized.contains("23andme") || normalized.contains("23&me") {
        vendor = Some("23andMe".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("23andMe header/export name".to_owned());
        if normalized.contains(" v2 ") || lower_name.contains("/v2/") {
            platform_version = Some("v2".to_owned());
            evidence.push("v2 token".to_owned());
        } else if normalized.contains(" v3 ") || lower_name.contains("/v3/") {
            platform_version = Some("v3".to_owned());
            evidence.push("v3 token".to_owned());
        } else if normalized.contains(" v4 ") || lower_name.contains("/v4/") {
            platform_version = Some("v4".to_owned());
            evidence.push("v4 token".to_owned());
        } else if normalized.contains(" v5 ") || lower_name.contains("/v5/") {
            platform_version = Some("v5".to_owned());
            evidence.push("v5 token".to_owned());
        }
    } else if normalized.contains("ancestrydna") || normalized.contains("ancestry com dna") {
        vendor = Some("AncestryDNA".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("AncestryDNA header/export name".to_owned());
        if let Some(version) = extract_after_marker(&header, "array version:") {
            platform_version = Some(canonicalize_ancestry_version(&version));
            evidence.push("AncestryDNA array version header".to_owned());
        }
    } else if normalized.contains("family tree dna")
        || normalized.contains("familytreedna")
        || normalized.contains("ftdna")
    {
        vendor = Some("FamilyTreeDNA".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("FamilyTreeDNA header/export name".to_owned());
    } else if normalized.contains("dynamic dna")
        || normalized.contains("dynamicdnalabs")
        || normalized.contains("ddna laboratories")
        || normalized.contains("ddna")
    {
        vendor = Some("Dynamic DNA".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("Dynamic DNA header".to_owned());
        if normalized.contains("gsav3 dtc") {
            platform_version = Some("GSAv3-DTC".to_owned());
            evidence.push("GSAv3-DTC token".to_owned());
        } else if normalized.contains("gsav3") {
            platform_version = Some("GSAv3".to_owned());
            evidence.push("GSAv3 token".to_owned());
        }
    } else if normalized.contains("myheritage") {
        vendor = Some("MyHeritage".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("MyHeritage header/export name".to_owned());
    } else if normalized.contains("sequencing com") && kind == DetectedKind::Vcf {
        vendor = Some("Sequencing.com".to_owned());
        confidence = DetectionConfidence::WeakHeuristic;
        evidence.push("sequencing.com header text".to_owned());
    } else if normalized.contains("carigenetics") || normalized.contains("cari genetics") {
        vendor = Some("CariGenetics".to_owned());
        confidence = DetectionConfidence::StrongHeuristic;
        evidence.push("CariGenetics path/header text".to_owned());
    }

    vendor.map(|vendor| SourceMetadata {
        vendor: Some(vendor),
        platform_version,
        confidence,
        evidence,
    })
}

fn extract_after_marker(text: &str, marker: &str) -> Option<String> {
    text.lines().find_map(|line| {
        let trimmed = line.trim();
        let lower = trimmed.to_ascii_lowercase();
        lower.find(marker).map(|idx| {
            trimmed[idx + marker.len()..]
                .trim()
                .trim_end_matches('.')
                .to_owned()
        })
    })
}

fn extract_token_after_marker(text: &str, marker: &str) -> Option<String> {
    extract_after_marker(text, marker).map(|value| {
        value
            .split_whitespace()
            .next()
            .unwrap_or_default()
            .trim_end_matches(':')
            .to_owned()
    })
}

fn canonicalize_ancestry_version(value: &str) -> String {
    let trimmed = value.trim();
    if let Some(rest) = trimmed.strip_prefix('v') {
        return format!("V{rest}");
    }
    trimmed.to_owned()
}

fn detect_assembly(lower_name: &str, sample_lines: &[String]) -> Option<Assembly> {
    let header = sample_lines.join("\n").to_ascii_lowercase();
    let combined = format!("{lower_name}\n{header}");
    let looks_like_grch38 = combined.contains("build 38")
        || combined.contains("grch38")
        || combined.contains("hg38")
        || combined.contains("gca_000001405.15")
        || combined.contains("grch38_no_alt_analysis_set")
        || combined.contains("##contig=<id=chr1,length=248956422>");

    if looks_like_grch38 {
        Some(Assembly::Grch38)
    } else if combined.contains("build 37")
        || combined.contains("grch37")
        || combined.contains("hg19")
        || combined.contains("assembly=b37")
        || combined.contains("assembly=\"b37\"")
        || combined.contains("human_g1k_v37")
        || combined.contains("37.1")
    {
        Some(Assembly::Grch37)
    } else {
        None
    }
}

fn detect_vcf_phasing(lines: &[String]) -> Option<bool> {
    let mut saw_slash = false;
    for line in lines {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 10 {
            continue;
        }
        let gt = fields[9].split(':').next().unwrap_or_default().trim();
        if gt.contains('|') {
            return Some(true);
        }
        if gt.contains('/') {
            saw_slash = true;
        }
    }
    saw_slash.then_some(false)
}

fn detect_index(
    path: &Path,
    kind: DetectedKind,
    options: &InspectOptions,
) -> (Option<bool>, Option<PathBuf>) {
    if let Some(index) = options
        .input_index
        .as_ref()
        .or(options.reference_index.as_ref())
    {
        return (Some(index.exists()), Some(index.clone()));
    }

    match kind {
        DetectedKind::AlignmentCram => {
            let candidate = if path
                .to_string_lossy()
                .to_ascii_lowercase()
                .ends_with(".cram")
            {
                let first = path.with_extension("cram.crai");
                if first.exists() {
                    Some(first)
                } else {
                    Some(path.with_extension("crai"))
                }
            } else {
                None
            };
            match candidate {
                Some(candidate) => (Some(candidate.exists()), Some(candidate)),
                None => (Some(false), None),
            }
        }
        DetectedKind::AlignmentBam => {
            let first = path.with_extension("bam.bai");
            if first.exists() {
                return (Some(true), Some(first));
            }
            let second = path.with_extension("bai");
            (Some(second.exists()), Some(second))
        }
        DetectedKind::ReferenceFasta => {
            let candidate = if let Some(ext) = path.extension().and_then(|ext| ext.to_str()) {
                path.with_extension(format!("{ext}.fai"))
            } else {
                path.with_extension("fai")
            };
            (Some(candidate.exists()), Some(candidate))
        }
        _ => (None, None),
    }
}

fn is_reference_path(path: &Path) -> bool {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    lower.ends_with(".fa") || lower.ends_with(".fasta")
}

fn classify_confidence(
    kind: DetectedKind,
    sample_lines: &[String],
    source: Option<&SourceMetadata>,
) -> DetectionConfidence {
    match kind {
        DetectedKind::Vcf if looks_like_vcf_lines(sample_lines) => {
            DetectionConfidence::Authoritative
        }
        DetectedKind::AlignmentCram | DetectedKind::AlignmentBam | DetectedKind::ReferenceFasta => {
            DetectionConfidence::Authoritative
        }
        DetectedKind::GenotypeText if source.is_some() => DetectionConfidence::StrongHeuristic,
        DetectedKind::GenotypeText => DetectionConfidence::WeakHeuristic,
        DetectedKind::Unknown => DetectionConfidence::Unknown,
        DetectedKind::Vcf => DetectionConfidence::StrongHeuristic,
    }
}

fn render_container(value: FileContainer) -> &'static str {
    match value {
        FileContainer::Plain => "plain",
        FileContainer::Zip => "zip",
    }
}

fn render_kind(value: DetectedKind) -> &'static str {
    match value {
        DetectedKind::GenotypeText => "genotype_text",
        DetectedKind::Vcf => "vcf",
        DetectedKind::AlignmentCram => "alignment_cram",
        DetectedKind::AlignmentBam => "alignment_bam",
        DetectedKind::ReferenceFasta => "reference_fasta",
        DetectedKind::Unknown => "unknown",
    }
}

fn render_confidence(value: DetectionConfidence) -> &'static str {
    match value {
        DetectionConfidence::Authoritative => "authoritative",
        DetectionConfidence::StrongHeuristic => "strong_heuristic",
        DetectionConfidence::WeakHeuristic => "weak_heuristic",
        DetectionConfidence::Unknown => "unknown",
    }
}

fn render_assembly(value: Option<Assembly>) -> &'static str {
    match value {
        Some(Assembly::Grch37) => "grch37",
        Some(Assembly::Grch38) => "grch38",
        None => "",
    }
}

fn render_bool(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "true",
        Some(false) => "false",
        None => "",
    }
}
