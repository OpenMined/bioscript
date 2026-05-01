use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use noodles::bgzf;
use noodles::csi;

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};

use crate::alignment;

use super::{
    describe_query, genotype_from_vcf_gt, is_bgzf_path, types::VcfBackend, variant_sort_key,
};

mod reader;

pub use reader::observe_vcf_snp_with_reader;

#[derive(Debug, Clone)]
pub(crate) struct ParsedVcfRow {
    pub(crate) rsid: Option<String>,
    pub(crate) chrom: String,
    pub(crate) position: i64,
    pub(crate) reference: String,
    pub(crate) alternates: Vec<String>,
    pub(crate) genotype: String,
}

pub(crate) fn scan_vcf_variants(
    backend: &VcfBackend,
    variants: &[VariantSpec],
) -> Result<Vec<VariantObservation>, RuntimeError> {
    let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
    indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

    let mut probe_lines = Vec::new();
    let detected_assembly = {
        let file = File::open(&backend.path).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open VCF file {}: {err}",
                backend.path.display()
            ))
        })?;
        let mut reader: Box<dyn BufRead> = if is_bgzf_path(&backend.path) {
            Box::new(BufReader::new(bgzf::io::Reader::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut buf = String::new();
        for _ in 0..256 {
            buf.clear();
            let bytes = reader.read_line(&mut buf).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read VCF file {}: {err}",
                    backend.path.display()
                ))
            })?;
            if bytes == 0 {
                break;
            }
            let line = buf.trim_end_matches(['\n', '\r']).to_owned();
            let stop = line.starts_with("#CHROM\t");
            probe_lines.push(line);
            if stop {
                break;
            }
        }

        detect_vcf_assembly(&backend.path, &probe_lines)
    };

    let mut rsid_targets: HashMap<String, Vec<usize>> = HashMap::new();
    let mut coord_targets: HashMap<(String, i64), Vec<usize>> = HashMap::new();
    let mut results = vec![VariantObservation::default(); variants.len()];
    let mut unresolved = variants.len();

    for (idx, variant) in &indexed {
        for rsid in &variant.rsids {
            rsid_targets.entry(rsid.clone()).or_default().push(*idx);
        }

        if let Some(locus) = choose_variant_locus_for_assembly(variant, detected_assembly) {
            let chrom = normalize_chromosome_name(&locus.chrom);
            coord_targets
                .entry((chrom.clone(), locus.start))
                .or_default()
                .push(*idx);
            if matches!(
                variant.kind,
                Some(VariantKind::Deletion | VariantKind::Insertion | VariantKind::Indel)
            ) {
                let anchor = locus.start.saturating_sub(1);
                coord_targets.entry((chrom, anchor)).or_default().push(*idx);
            }
        }
    }

    let targets = VcfResolutionTargets {
        variants,
        detected_assembly,
        rsid_targets: &rsid_targets,
        coord_targets: &coord_targets,
    };

    let file = File::open(&backend.path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open VCF file {}: {err}",
            backend.path.display()
        ))
    })?;
    let mut reader: Box<dyn BufRead> = if is_bgzf_path(&backend.path) {
        Box::new(BufReader::new(bgzf::io::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut buf = String::new();
    loop {
        buf.clear();
        let bytes = reader.read_line(&mut buf).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read VCF file {}: {err}",
                backend.path.display()
            ))
        })?;
        if bytes == 0 || unresolved == 0 {
            break;
        }
        if let Some(row) = parse_vcf_record(buf.trim_end_matches(['\n', '\r']))? {
            resolve_vcf_row(backend, &row, &targets, &mut results, &mut unresolved);
        }
    }

    for (idx, variant) in indexed {
        if results[idx].genotype.is_none() {
            results[idx] = VariantObservation {
                backend: backend.backend_name().to_owned(),
                assembly: detected_assembly,
                evidence: vec![format!(
                    "no matching rsid or locus found for {}",
                    describe_query(variant)
                )],
                ..VariantObservation::default()
            };
        }
    }

    Ok(results)
}

pub(crate) fn lookup_indexed_vcf_variants(
    backend: &VcfBackend,
    variants: &[VariantSpec],
) -> Result<Option<Vec<VariantObservation>>, RuntimeError> {
    let Some(input_index) = backend.options.input_index.as_ref() else {
        return Ok(None);
    };
    let detected_assembly = detect_vcf_assembly_from_path(&backend.path)?;
    let mut indexed_variants = Vec::with_capacity(variants.len());
    for (idx, variant) in variants.iter().enumerate() {
        let Some(locus) = choose_variant_locus_for_assembly(variant, detected_assembly) else {
            return Ok(None);
        };
        let Some(reference) = first_single_base_allele(variant.reference.as_deref()) else {
            return Ok(None);
        };
        let Some(alternate) = first_single_base_allele(variant.alternate.as_deref()) else {
            return Ok(None);
        };
        if !matches!(variant.kind, None | Some(VariantKind::Snp)) {
            return Ok(None);
        }
        indexed_variants.push((idx, variant, locus, reference, alternate));
    }

    let tabix_index = alignment::parse_tbi_bytes(&std::fs::read(input_index).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to read VCF index {}: {err}",
            input_index.display()
        ))
    })?)?;
    let mut indexed = csi::io::IndexedReader::new(
        File::open(&backend.path).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open VCF file {}: {err}",
                backend.path.display()
            ))
        })?,
        tabix_index,
    );

    let mut results = vec![VariantObservation::default(); variants.len()];
    for (idx, variant, locus, reference, alternate) in indexed_variants {
        results[idx] = observe_vcf_snp_with_reader(
            &mut indexed,
            &backend.path.display().to_string(),
            &locus,
            reference,
            alternate,
            variant.rsids.first().cloned(),
            detected_assembly,
        )?;
    }
    Ok(Some(results))
}

pub(crate) fn detect_vcf_assembly_from_path(path: &Path) -> Result<Option<Assembly>, RuntimeError> {
    let mut probe_lines = Vec::new();
    let file = File::open(path).map_err(|err| {
        RuntimeError::Io(format!("failed to open VCF file {}: {err}", path.display()))
    })?;
    let mut reader: Box<dyn BufRead> = if is_bgzf_path(path) {
        Box::new(BufReader::new(bgzf::io::Reader::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut buf = String::new();
    for _ in 0..256 {
        buf.clear();
        let bytes = reader.read_line(&mut buf).map_err(|err| {
            RuntimeError::Io(format!("failed to read VCF file {}: {err}", path.display()))
        })?;
        if bytes == 0 {
            break;
        }
        let line = buf.trim_end_matches(['\n', '\r']).to_owned();
        let stop = line.starts_with("#CHROM\t");
        probe_lines.push(line);
        if stop {
            break;
        }
    }
    Ok(detect_vcf_assembly(path, &probe_lines))
}

fn first_single_base_allele(value: Option<&str>) -> Option<char> {
    let value = value?;
    let mut chars = value.chars();
    let base = chars.next()?;
    chars.next().is_none().then_some(base)
}

struct VcfResolutionTargets<'a> {
    variants: &'a [VariantSpec],
    detected_assembly: Option<Assembly>,
    rsid_targets: &'a HashMap<String, Vec<usize>>,
    coord_targets: &'a HashMap<(String, i64), Vec<usize>>,
}

fn resolve_vcf_row(
    backend: &VcfBackend,
    row: &ParsedVcfRow,
    targets: &VcfResolutionTargets<'_>,
    results: &mut [VariantObservation],
    unresolved: &mut usize,
) {
    if let Some(rsid) = row.rsid.as_ref()
        && let Some(target_indexes) = targets.rsid_targets.get(rsid)
    {
        for &target_idx in target_indexes {
            if results[target_idx].genotype.is_none() {
                results[target_idx] = VariantObservation {
                    backend: backend.backend_name().to_owned(),
                    matched_rsid: Some(rsid.clone()),
                    assembly: targets.detected_assembly,
                    genotype: Some(row.genotype.clone()),
                    evidence: vec![format!("resolved by rsid {rsid}")],
                    ..VariantObservation::default()
                };
                *unresolved = (*unresolved).saturating_sub(1);
            }
        }
    }

    if *unresolved == 0 {
        return;
    }

    let key = (normalize_chromosome_name(&row.chrom), row.position);
    if let Some(target_indexes) = targets.coord_targets.get(&key) {
        for &target_idx in target_indexes {
            if results[target_idx].genotype.is_none()
                && vcf_row_matches_variant(
                    row,
                    &targets.variants[target_idx],
                    targets.detected_assembly,
                )
            {
                results[target_idx] = VariantObservation {
                    backend: backend.backend_name().to_owned(),
                    matched_rsid: row.rsid.clone(),
                    assembly: targets.detected_assembly,
                    genotype: Some(row.genotype.clone()),
                    evidence: vec![format!("resolved by locus {}:{}", row.chrom, row.position)],
                    ..VariantObservation::default()
                };
                *unresolved = (*unresolved).saturating_sub(1);
            }
        }
    }
}

pub(crate) fn parse_vcf_record(line: &str) -> Result<Option<ParsedVcfRow>, RuntimeError> {
    let trimmed = line.trim();
    if trimmed.is_empty() || trimmed.starts_with('#') {
        return Ok(None);
    }

    let fields: Vec<&str> = trimmed.split('\t').collect();
    if fields.len() < 10 {
        return Ok(None);
    }

    let chrom = fields[0].trim();
    let position = fields[1].trim().parse::<i64>().map_err(|err| {
        RuntimeError::Io(format!(
            "failed to parse VCF position '{}': {err}",
            fields[1].trim()
        ))
    })?;
    let rsid = {
        let value = fields[2].trim();
        (!value.is_empty() && value != ".").then(|| value.to_owned())
    };
    let reference = fields[3].trim();
    if reference.is_empty() || reference == "." {
        return Ok(None);
    }

    let alternates: Vec<String> = fields[4]
        .split(',')
        .map(str::trim)
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .map(ToOwned::to_owned)
        .collect();
    if alternates.is_empty() {
        return Ok(None);
    }

    let genotype = extract_vcf_sample_genotype(fields[8], fields[9], reference, &alternates)
        .unwrap_or_else(|| "--".to_owned());

    Ok(Some(ParsedVcfRow {
        rsid,
        chrom: chrom.to_owned(),
        position,
        reference: reference.to_owned(),
        alternates,
        genotype,
    }))
}

pub(crate) fn extract_vcf_sample_genotype(
    format_field: &str,
    sample_field: &str,
    reference: &str,
    alternates: &[String],
) -> Option<String> {
    let gt_index = format_field
        .split(':')
        .position(|field| field.eq_ignore_ascii_case("GT"))?;
    let sample_parts: Vec<&str> = sample_field.split(':').collect();
    let sample_gt = sample_parts.get(gt_index).copied().unwrap_or(".");
    let alternate_refs: Vec<&str> = alternates.iter().map(String::as_str).collect();
    genotype_from_vcf_gt(sample_gt, reference, &alternate_refs)
}

pub(crate) fn detect_vcf_assembly(path: &Path, probe_lines: &[String]) -> Option<Assembly> {
    let combined = probe_lines.join("\n").to_ascii_lowercase();
    if combined.contains("assembly=b37")
        || combined.contains("assembly=grch37")
        || combined.contains("assembly=hg19")
        || combined.contains("reference=grch37")
        || combined.contains("reference=hg19")
    {
        return Some(Assembly::Grch37);
    }
    if combined.contains("assembly=b38")
        || combined.contains("assembly=grch38")
        || combined.contains("assembly=hg38")
        || combined.contains("reference=grch38")
        || combined.contains("reference=hg38")
    {
        return Some(Assembly::Grch38);
    }

    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.contains("grch37") || lower.contains("hg19") || lower.contains("b37") {
        Some(Assembly::Grch37)
    } else if lower.contains("grch38") || lower.contains("hg38") || lower.contains("b38") {
        Some(Assembly::Grch38)
    } else {
        None
    }
}

pub(crate) fn choose_variant_locus_for_assembly(
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> Option<GenomicLocus> {
    match assembly {
        Some(Assembly::Grch37) => variant.grch37.clone().or_else(|| variant.grch38.clone()),
        Some(Assembly::Grch38) => variant.grch38.clone().or_else(|| variant.grch37.clone()),
        None => variant.grch37.clone().or_else(|| variant.grch38.clone()),
    }
}

pub(crate) fn normalize_chromosome_name(value: &str) -> String {
    value.trim().trim_start_matches("chr").to_ascii_lowercase()
}

pub(crate) fn vcf_row_matches_variant(
    row: &ParsedVcfRow,
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> bool {
    let Some(locus) = choose_variant_locus_for_assembly(variant, assembly) else {
        return false;
    };

    if normalize_chromosome_name(&row.chrom) != normalize_chromosome_name(&locus.chrom) {
        return false;
    }

    match variant.kind.unwrap_or(VariantKind::Other) {
        VariantKind::Snp => {
            row.position == locus.start
                && variant
                    .reference
                    .as_ref()
                    .is_none_or(|reference| reference.eq_ignore_ascii_case(&row.reference))
                && variant.alternate.as_ref().is_none_or(|alternate| {
                    row.alternates
                        .iter()
                        .any(|candidate| candidate.eq_ignore_ascii_case(alternate))
                })
        }
        VariantKind::Deletion => {
            let expected_len = variant.deletion_length.unwrap_or(0);
            row.position == locus.start.saturating_sub(1)
                && row.alternates.iter().any(|alternate| {
                    let actual_len = row.reference.len().saturating_sub(alternate.len());
                    (expected_len == 0 || actual_len == expected_len)
                        && alternate.len() < row.reference.len()
                })
        }
        VariantKind::Insertion | VariantKind::Indel => {
            row.position == locus.start.saturating_sub(1)
        }
        VariantKind::Other => row.position == locus.start,
    }
}
