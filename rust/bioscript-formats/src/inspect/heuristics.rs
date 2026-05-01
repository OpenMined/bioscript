use std::path::{Path, PathBuf};

use bioscript_core::Assembly;

use super::{DetectedKind, DetectionConfidence, InspectOptions, SourceMetadata};

pub(crate) fn looks_like_vcf_lines(lines: &[String]) -> bool {
    lines.iter().any(|line| {
        let trimmed = line.trim_start();
        trimmed.starts_with("##fileformat=VCF") || trimmed.starts_with("#CHROM\t")
    })
}

pub(crate) fn looks_like_genotype_text(lines: &[String]) -> bool {
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

pub(crate) fn split_fields(line: &str) -> Vec<String> {
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

pub(crate) fn matches_genotype_shape(fields: &[String]) -> bool {
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

pub(crate) fn is_valid_genotype(value: &str) -> bool {
    let trimmed = value.trim().to_ascii_uppercase();
    if trimmed.is_empty() || trimmed.len() > 4 {
        return false;
    }
    trimmed
        .chars()
        .all(|ch| matches!(ch, 'A' | 'C' | 'G' | 'T' | 'I' | 'D' | '-' | '0'))
}

pub(crate) fn is_valid_allele(value: &str) -> bool {
    let trimmed = value.trim().to_ascii_uppercase();
    matches!(
        trimmed.as_str(),
        "A" | "C" | "G" | "T" | "I" | "D" | "-" | "0"
    )
}

pub(crate) fn detect_source(
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

pub(crate) fn canonicalize_ancestry_version(value: &str) -> String {
    let trimmed = value.trim();
    if let Some(rest) = trimmed.strip_prefix('v') {
        return format!("V{rest}");
    }
    trimmed.to_owned()
}

pub(crate) fn detect_assembly(lower_name: &str, sample_lines: &[String]) -> Option<Assembly> {
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

pub(crate) fn detect_vcf_phasing(lines: &[String]) -> Option<bool> {
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

pub(crate) fn detect_index(
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

pub(crate) fn is_reference_path(path: &Path) -> bool {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    lower.ends_with(".fa") || lower.ends_with(".fasta")
}

pub(crate) fn classify_confidence(
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
