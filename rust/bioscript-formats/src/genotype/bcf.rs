use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{Cursor, Read},
};

use bioscript_core::{Assembly, RuntimeError, VariantKind, VariantObservation, VariantSpec};
use noodles::bcf;
use zip::ZipArchive;

mod records;

use records::{bcf_record_to_vcf_row, read_bcf_header_lenient};

use super::{
    common::variant_sort_key,
    types::{BcfBackend, BcfSource},
    vcf::{
        ParsedVcfRow, choose_variant_locus_for_assembly, imputed_reference_observation,
        normalize_chromosome_name, vcf_row_genotype_for_variant, vcf_row_matches_variant,
    },
};

pub(crate) fn scan_bcf_variants(
    backend: &BcfBackend,
    variants: &[VariantSpec],
) -> Result<Vec<VariantObservation>, RuntimeError> {
    let assembly = backend.options.assembly.or(Some(Assembly::Grch38));
    let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
    indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

    let targets = BcfTargets::new(variants, assembly);
    let mut results = vec![VariantObservation::default(); variants.len()];
    let mut unresolved = variants.len();

    for shard in backend.shards()? {
        if unresolved == 0 {
            break;
        }
        if !targets.should_scan_shard(&shard.name) {
            continue;
        }
        scan_bcf_shard(
            backend,
            &shard,
            &targets,
            &mut results,
            &mut unresolved,
            assembly,
        )?;
    }

    for (idx, variant) in indexed {
        if results[idx].genotype.is_none() {
            let label = backend.label();
            results[idx] = if backend.options.impute_vcf_missing_as_reference {
                choose_variant_locus_for_assembly(variant, assembly)
                    .and_then(|locus| {
                        imputed_reference_observation(
                            backend.backend_name(),
                            &label,
                            variant,
                            &locus,
                            assembly,
                            backend.options.inferred_sex,
                            "no matching rsid or locus found in BCF",
                        )
                    })
                    .unwrap_or_else(|| missing_observation(backend, assembly))
            } else {
                missing_observation(backend, assembly)
            };
        }
    }

    Ok(results)
}

struct BcfTargets<'a> {
    variants: &'a [VariantSpec],
    assembly: Option<Assembly>,
    rsid_targets: HashMap<String, Vec<usize>>,
    coord_targets: HashMap<(String, i64), Vec<usize>>,
    chroms: HashSet<String>,
}

impl<'a> BcfTargets<'a> {
    fn new(variants: &'a [VariantSpec], assembly: Option<Assembly>) -> Self {
        let mut rsid_targets: HashMap<String, Vec<usize>> = HashMap::new();
        let mut coord_targets: HashMap<(String, i64), Vec<usize>> = HashMap::new();
        let mut chroms = HashSet::new();
        for (idx, variant) in variants.iter().enumerate() {
            for rsid in &variant.rsids {
                rsid_targets.entry(rsid.clone()).or_default().push(idx);
            }
            if let Some(locus) = choose_variant_locus_for_assembly(variant, assembly) {
                let chrom = normalize_chromosome_name(&locus.chrom);
                chroms.insert(chrom.clone());
                coord_targets
                    .entry((chrom.clone(), locus.start))
                    .or_default()
                    .push(idx);
                if matches!(
                    variant.kind,
                    Some(VariantKind::Deletion | VariantKind::Insertion | VariantKind::Indel)
                ) {
                    coord_targets
                        .entry((chrom, locus.start.saturating_sub(1)))
                        .or_default()
                        .push(idx);
                }
            }
        }
        Self {
            variants,
            assembly,
            rsid_targets,
            coord_targets,
            chroms,
        }
    }

    fn should_scan_shard(&self, name: &str) -> bool {
        if self.chroms.is_empty() {
            return true;
        }
        let Some(chrom) = chrom_from_bcf_shard_name(name) else {
            return true;
        };
        self.chroms.contains(&chrom)
    }
}

struct BcfShard {
    name: String,
    data: BcfShardData,
}

enum BcfShardData {
    File(std::path::PathBuf),
    Bytes(Vec<u8>),
    ZipEntry {
        zip_path: std::path::PathBuf,
        entry_name: String,
    },
}

impl BcfBackend {
    fn label(&self) -> String {
        match &self.source {
            BcfSource::File(path) | BcfSource::ZipFile { path, .. } => path.display().to_string(),
            BcfSource::Bytes { name, .. } | BcfSource::ZipBytes { name, .. } => name.clone(),
        }
    }

    fn shards(&self) -> Result<Vec<BcfShard>, RuntimeError> {
        match &self.source {
            BcfSource::File(path) => Ok(vec![BcfShard {
                name: path.display().to_string(),
                data: BcfShardData::File(path.clone()),
            }]),
            BcfSource::Bytes { name, data } => Ok(vec![BcfShard {
                name: name.clone(),
                data: BcfShardData::Bytes(data.clone()),
            }]),
            BcfSource::ZipFile { path, entries } => Ok(entries
                .iter()
                .map(|entry_name| BcfShard {
                    name: entry_name.clone(),
                    data: BcfShardData::ZipEntry {
                        zip_path: path.clone(),
                        entry_name: entry_name.clone(),
                    },
                })
                .collect()),
            BcfSource::ZipBytes { entries, .. } => Ok(entries
                .iter()
                .map(|(name, data)| BcfShard {
                    name: name.clone(),
                    data: BcfShardData::Bytes(data.clone()),
                })
                .collect()),
        }
    }
}

fn scan_bcf_shard(
    backend: &BcfBackend,
    shard: &BcfShard,
    targets: &BcfTargets<'_>,
    results: &mut [VariantObservation],
    unresolved: &mut usize,
    assembly: Option<Assembly>,
) -> Result<(), RuntimeError> {
    match &shard.data {
        BcfShardData::File(path) => {
            let file = File::open(path).map_err(|err| {
                RuntimeError::Io(format!("failed to open BCF file {}: {err}", path.display()))
            })?;
            scan_bcf_reader(
                backend,
                shard,
                bcf::io::Reader::new(file),
                targets,
                results,
                unresolved,
                assembly,
            )
        }
        BcfShardData::Bytes(data) => scan_bcf_reader(
            backend,
            shard,
            bcf::io::Reader::new(Cursor::new(data)),
            targets,
            results,
            unresolved,
            assembly,
        ),
        BcfShardData::ZipEntry {
            zip_path,
            entry_name,
        } => {
            let file = File::open(zip_path).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype zip {}: {err}",
                    zip_path.display()
                ))
            })?;
            let mut archive = ZipArchive::new(file).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype zip {}: {err}",
                    zip_path.display()
                ))
            })?;
            let mut entry = archive.by_name(entry_name).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype entry {entry_name} in {}: {err}",
                    zip_path.display()
                ))
            })?;
            let mut data = Vec::new();
            entry.read_to_end(&mut data).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype entry {entry_name} in {}: {err}",
                    zip_path.display()
                ))
            })?;
            scan_bcf_reader(
                backend,
                shard,
                bcf::io::Reader::new(Cursor::new(data)),
                targets,
                results,
                unresolved,
                assembly,
            )
        }
    }
}

fn scan_bcf_reader<R: Read>(
    backend: &BcfBackend,
    shard: &BcfShard,
    mut reader: bcf::io::Reader<R>,
    targets: &BcfTargets<'_>,
    results: &mut [VariantObservation],
    unresolved: &mut usize,
    assembly: Option<Assembly>,
) -> Result<(), RuntimeError> {
    let header = read_bcf_header_lenient(&mut reader, &shard.name)?;
    let string_maps = header.string_maps().clone();

    for record_result in reader.records() {
        if *unresolved == 0 {
            break;
        }
        let record = record_result.map_err(|err| {
            RuntimeError::Io(format!("failed to read BCF record {}: {err}", shard.name))
        })?;
        let row = bcf_record_to_vcf_row(&header, &string_maps, &record, &shard.name)?;
        resolve_bcf_row(backend, &row, targets, results, unresolved, assembly);
    }
    Ok(())
}

fn resolve_bcf_row(
    backend: &BcfBackend,
    row: &ParsedVcfRow,
    targets: &BcfTargets<'_>,
    results: &mut [VariantObservation],
    unresolved: &mut usize,
    assembly: Option<Assembly>,
) {
    if let Some(rsid) = row.rsid.as_ref()
        && let Some(target_indexes) = targets.rsid_targets.get(rsid)
    {
        for &target_idx in target_indexes {
            if results[target_idx].genotype.is_none() {
                results[target_idx] = observation_from_row(
                    backend,
                    row,
                    target_idx,
                    targets,
                    assembly,
                    Some(format!("resolved by rsid {rsid}")),
                );
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
                && vcf_row_matches_variant(row, &targets.variants[target_idx], targets.assembly)
            {
                results[target_idx] = observation_from_row(
                    backend,
                    row,
                    target_idx,
                    targets,
                    assembly,
                    Some(format!("resolved by locus {}:{}", row.chrom, row.position)),
                );
                *unresolved = (*unresolved).saturating_sub(1);
            }
        }
    }
}

fn observation_from_row(
    backend: &BcfBackend,
    row: &ParsedVcfRow,
    target_idx: usize,
    targets: &BcfTargets<'_>,
    assembly: Option<Assembly>,
    resolved: Option<String>,
) -> VariantObservation {
    let resolved = resolved.unwrap_or_else(|| "resolved by BCF record".to_owned());
    VariantObservation {
        backend: backend.backend_name().to_owned(),
        matched_rsid: row.rsid.clone(),
        assembly,
        genotype: Some(vcf_row_genotype_for_variant(
            row,
            &targets.variants[target_idx],
        )),
        evidence: vec![resolved, format!("source record: {}", row.raw_line)],
        ..VariantObservation::default()
    }
}

fn missing_observation(backend: &BcfBackend, assembly: Option<Assembly>) -> VariantObservation {
    VariantObservation {
        backend: backend.backend_name().to_owned(),
        assembly,
        evidence: vec!["no matching rsid or locus found in BCF".to_owned()],
        ..VariantObservation::default()
    }
}

fn chrom_from_bcf_shard_name(name: &str) -> Option<String> {
    let file = name.rsplit('/').next().unwrap_or(name);
    let lower = file.to_ascii_lowercase();
    let chrom = lower.strip_suffix(".bcf").unwrap_or(&lower);
    let chrom = chrom
        .strip_prefix("chr")
        .unwrap_or(chrom)
        .split('.')
        .next()
        .unwrap_or(chrom);
    (!chrom.is_empty()).then(|| chrom.to_owned())
}
