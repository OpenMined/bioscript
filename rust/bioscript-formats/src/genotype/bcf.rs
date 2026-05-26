use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{Cursor, Read},
};

use bioscript_core::{Assembly, RuntimeError, VariantKind, VariantObservation, VariantSpec};
use noodles::bcf;
use noodles::vcf::{
    self,
    variant::record::{
        AlternateBases as _, Ids as _, ReferenceBases as _,
        samples::Sample as _,
        samples::series::{Value as SampleValue, value::Array},
    },
};
use zip::ZipArchive;

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
            BcfSource::File(path) => path.display().to_string(),
            BcfSource::ZipFile { path, .. } => path.display().to_string(),
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

fn read_bcf_header_lenient<R: Read>(
    reader: &mut bcf::io::Reader<R>,
    label: &str,
) -> Result<vcf::Header, RuntimeError> {
    let mut header_reader = reader.header_reader();
    header_reader
        .read_magic_number()
        .map_err(|err| RuntimeError::Io(format!("failed to read BCF magic {label}: {err}")))?;
    header_reader
        .read_format_version()
        .map_err(|err| RuntimeError::Io(format!("failed to read BCF version {label}: {err}")))?;
    let mut raw_reader = header_reader
        .raw_vcf_header_reader()
        .map_err(|err| RuntimeError::Io(format!("failed to open BCF VCF header {label}: {err}")))?;
    let mut raw = String::new();
    raw_reader
        .read_to_string(&mut raw)
        .map_err(|err| RuntimeError::Io(format!("failed to read BCF VCF header {label}: {err}")))?;
    raw_reader.discard_to_end().map_err(|err| {
        RuntimeError::Io(format!(
            "failed to discard BCF VCF header padding {label}: {err}"
        ))
    })?;

    let sanitized = sanitize_bcf_vcf_header(&raw);
    sanitized
        .parse::<vcf::Header>()
        .map_err(|err| RuntimeError::Io(format!("failed to parse BCF VCF header {label}: {err}")))
}

fn sanitize_bcf_vcf_header(raw: &str) -> String {
    let is_23andme_imputed = raw.contains("##DISCLAIMER=") && raw.contains("23andMe");
    let mut inserted_23andme_contigs = false;
    raw.lines()
        .filter_map(|line| {
            if is_23andme_imputed && line.starts_with("##contig=<") {
                if inserted_23andme_contigs {
                    return None;
                }
                inserted_23andme_contigs = true;
                return Some(
                    (1..=22)
                        .map(|chrom| format!("##contig=<ID=chr{chrom}>"))
                        .chain(std::iter::once("##contig=<ID=chrX>".to_owned()))
                        .collect::<Vec<_>>()
                        .join("\n"),
                );
            }
            if line.starts_with("##FORMAT=<")
                && line.contains("ID=HDS")
                && !line.contains("Description=")
            {
                Some(
                    "##FORMAT=<ID=HDS,Number=.,Type=Float,Description=\"haploid alternate dosage\">"
                        .to_owned(),
                )
            } else {
                Some(line.to_owned())
            }
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn bcf_record_to_vcf_row(
    header: &vcf::Header,
    string_maps: &vcf::header::StringMaps,
    record: &bcf::Record,
    label: &str,
) -> Result<ParsedVcfRow, RuntimeError> {
    let chrom = record
        .reference_sequence_name(string_maps)
        .map(str::to_owned)
        .or_else(|_| {
            chrom_from_bcf_shard_name(label)
                .map(|chrom| {
                    if chrom == "mt" {
                        "chrM".to_owned()
                    } else {
                        format!("chr{chrom}")
                    }
                })
                .ok_or_else(|| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        "missing reference sequence name in contig string map",
                    )
                })
        })
        .map_err(|err| {
            RuntimeError::Io(format!("{label}: failed to read BCF chromosome: {err}"))
        })?;
    let position = record
        .variant_start()
        .transpose()
        .map_err(|err| RuntimeError::Io(format!("{label}: failed to read BCF position: {err}")))?
        .ok_or_else(|| RuntimeError::Io(format!("{label}: BCF record missing position")))?
        .get() as i64;
    let ids_buf = record.ids();
    let ids: Vec<&str> = ids_buf.iter().collect();
    let rsid = ids
        .iter()
        .find(|id| !id.is_empty() && **id != ".")
        .map(|id| (*id).to_owned());
    let reference = String::from_utf8(
        record
            .reference_bases()
            .iter()
            .collect::<Result<Vec<_>, _>>()
            .map_err(|err| {
                RuntimeError::Io(format!(
                    "{label}: failed to read BCF reference bases: {err}"
                ))
            })?,
    )
    .map_err(|err| RuntimeError::Io(format!("{label}: invalid BCF reference bases: {err}")))?;
    let alternates = record
        .alternate_bases()
        .iter()
        .map(|result| result.map(str::to_owned))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|err| {
            RuntimeError::Io(format!(
                "{label}: failed to read BCF alternate bases: {err}"
            ))
        })?;
    let hds = extract_hds(header, record, label)?;
    let genotype = genotype_from_hds(&reference, alternates.first().map(String::as_str), &hds)
        .unwrap_or_else(|| "--".to_owned());
    let raw_line = format!(
        "{}  {}  {}  {}  {}  FORMAT=HDS  HDS={}",
        chrom,
        position,
        rsid.as_deref().unwrap_or("."),
        reference,
        alternates.join(","),
        hds.iter()
            .map(|value| format_hds_value(*value))
            .collect::<Vec<_>>()
            .join(",")
    );

    Ok(ParsedVcfRow {
        rsid,
        chrom,
        position,
        reference,
        alternates,
        genotype,
        raw_line,
    })
}

fn extract_hds(
    header: &vcf::Header,
    record: &bcf::Record,
    label: &str,
) -> Result<Vec<f32>, RuntimeError> {
    let samples = record
        .samples()
        .map_err(|err| RuntimeError::Io(format!("{label}: failed to read BCF samples: {err}")))?;
    if let Some(values) = parse_first_bcf_float_series(samples.as_ref()) {
        return Ok(values);
    }
    let Some(sample) = samples.get_index(0) else {
        return Ok(Vec::new());
    };
    let Some(value) = sample
        .get(header, "HDS")
        .transpose()
        .map_err(|err| RuntimeError::Io(format!("{label}: failed to read BCF HDS: {err}")))?
        .flatten()
    else {
        return Ok(Vec::new());
    };
    match value {
        SampleValue::Array(Array::Float(values)) => values
            .iter()
            .filter_map(|result| {
                result
                    .map_err(|err| {
                        RuntimeError::Io(format!("{label}: failed to read HDS value: {err}"))
                    })
                    .transpose()
            })
            .collect(),
        SampleValue::Float(value) => Ok(vec![value]),
        _ => Ok(Vec::new()),
    }
}

fn parse_first_bcf_float_series(src: &[u8]) -> Option<Vec<f32>> {
    let mut offset = 0usize;
    skip_bcf_typed_value(src, &mut offset)?;
    let descriptor = *src.get(offset)?;
    offset += 1;
    let (len, ty) = bcf_descriptor_len_ty(src, descriptor, &mut offset)?;
    if ty != 5 {
        return None;
    }
    let mut values = Vec::with_capacity(len);
    for _ in 0..len {
        let bytes: [u8; 4] = src.get(offset..offset + 4)?.try_into().ok()?;
        offset += 4;
        let value = f32::from_le_bytes(bytes);
        if !value.is_nan() {
            values.push(value);
        }
    }
    Some(values)
}

fn skip_bcf_typed_value(src: &[u8], offset: &mut usize) -> Option<()> {
    let descriptor = *src.get(*offset)?;
    *offset += 1;
    let (len, ty) = bcf_descriptor_len_ty(src, descriptor, offset)?;
    let width = match ty {
        1 | 7 => 1,
        2 => 2,
        3 | 5 => 4,
        _ => return None,
    };
    *offset = offset.checked_add(len.checked_mul(width)?)?;
    (*offset <= src.len()).then_some(())
}

fn bcf_descriptor_len_ty(src: &[u8], descriptor: u8, offset: &mut usize) -> Option<(usize, u8)> {
    let ty = descriptor & 0x0f;
    let mut len = usize::from(descriptor >> 4);
    if len == 15 {
        let len_descriptor = *src.get(*offset)?;
        *offset += 1;
        let (len_len, len_ty) = bcf_descriptor_len_ty(src, len_descriptor, offset)?;
        if len_len != 1 {
            return None;
        }
        len = match len_ty {
            1 => usize::from(*src.get(*offset)?),
            2 => {
                let bytes: [u8; 2] = src.get(*offset..*offset + 2)?.try_into().ok()?;
                usize::try_from(i16::from_le_bytes(bytes)).ok()?
            }
            3 => {
                let bytes: [u8; 4] = src.get(*offset..*offset + 4)?.try_into().ok()?;
                usize::try_from(i32::from_le_bytes(bytes)).ok()?
            }
            _ => return None,
        };
        *offset += match len_ty {
            1 => 1,
            2 => 2,
            3 => 4,
            _ => return None,
        };
    }
    Some((len, ty))
}

fn genotype_from_hds(reference: &str, alternate: Option<&str>, hds: &[f32]) -> Option<String> {
    let alternate = alternate?;
    if reference.is_empty() || alternate.is_empty() || hds.is_empty() {
        return None;
    }
    let mut alleles = Vec::with_capacity(hds.len().max(2));
    for dosage in hds {
        if *dosage >= 0.5 {
            alleles.push(alternate.to_owned());
        } else {
            alleles.push(reference.to_owned());
        }
    }
    if alleles.len() == 1 {
        alleles.push(reference.to_owned());
    }
    if alleles.iter().any(|allele| allele.len() > 1) {
        return Some(alleles.join("/"));
    }
    Some(super::normalize_genotype(&alleles.join("/")))
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

fn format_hds_value(value: f32) -> String {
    if (value.fract()).abs() < f32::EPSILON {
        format!("{value:.0}")
    } else {
        format!("{value:.3}")
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_owned()
    }
}
