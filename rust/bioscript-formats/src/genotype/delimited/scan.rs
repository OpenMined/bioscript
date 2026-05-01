use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use zip::ZipArchive;

use bioscript_core::{RuntimeError, VariantObservation, VariantSpec};

use super::{
    super::{GenotypeSourceFormat, describe_query, types::DelimitedBackend, variant_sort_key},
    DelimitedColumnIndexes, detect_delimiter, parse_streaming_row,
};

pub(crate) fn scan_delimited_variants(
    backend: &DelimitedBackend,
    variants: &[VariantSpec],
) -> Result<Vec<VariantObservation>, RuntimeError> {
    let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
    indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

    let mut rsid_targets: HashMap<String, Vec<usize>> = HashMap::new();
    let mut coord_targets: HashMap<(String, i64), Vec<usize>> = HashMap::new();
    let mut results = vec![VariantObservation::default(); variants.len()];
    let mut unresolved = variants.len();

    for (idx, variant) in &indexed {
        for rsid in &variant.rsids {
            rsid_targets.entry(rsid.clone()).or_default().push(*idx);
        }
        if let Some(locus) = variant.grch38.as_ref().or(variant.grch37.as_ref()) {
            coord_targets
                .entry((
                    locus.chrom.trim_start_matches("chr").to_ascii_lowercase(),
                    locus.start,
                ))
                .or_default()
                .push(*idx);
        }
    }

    let mut scan_reader = |reader: &mut dyn BufRead| -> Result<(), RuntimeError> {
        let mut probe_lines = Vec::new();
        let mut buf = String::new();
        for _ in 0..8 {
            buf.clear();
            let bytes = reader.read_line(&mut buf).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype stream {}: {err}",
                    backend.path.display()
                ))
            })?;
            if bytes == 0 {
                break;
            }
            probe_lines.push(buf.trim_end_matches(['\n', '\r']).to_owned());
        }

        let delimiter = detect_delimiter(&probe_lines);
        let mut column_indexes: Option<DelimitedColumnIndexes> = None;
        let mut comment_header: Option<Vec<String>> = None;

        let mut process_line = |line: &str| -> Result<bool, RuntimeError> {
            let Some(row) =
                parse_streaming_row(line, delimiter, &mut column_indexes, &mut comment_header)?
            else {
                return Ok(unresolved == 0);
            };

            if let Some(rsid) = row.rsid.as_ref()
                && let Some(target_indexes) = rsid_targets.get(rsid)
            {
                for &target_idx in target_indexes {
                    if results[target_idx].genotype.is_none() {
                        results[target_idx] = VariantObservation {
                            backend: backend.backend_name().to_owned(),
                            matched_rsid: Some(rsid.clone()),
                            genotype: Some(row.genotype.clone()),
                            evidence: vec![format!("resolved by rsid {rsid}")],
                            ..VariantObservation::default()
                        };
                        unresolved = unresolved.saturating_sub(1);
                    }
                }
            }

            if unresolved == 0 {
                return Ok(true);
            }

            if let (Some(chrom), Some(position)) = (row.chrom.as_ref(), row.position) {
                let key = (
                    chrom.trim_start_matches("chr").to_ascii_lowercase(),
                    position,
                );
                if let Some(target_indexes) = coord_targets.get(&key) {
                    for &target_idx in target_indexes {
                        if results[target_idx].genotype.is_none() {
                            results[target_idx] = VariantObservation {
                                backend: backend.backend_name().to_owned(),
                                matched_rsid: row.rsid.clone(),
                                genotype: Some(row.genotype.clone()),
                                evidence: vec![format!("resolved by locus {}:{}", chrom, position)],
                                ..VariantObservation::default()
                            };
                            unresolved = unresolved.saturating_sub(1);
                        }
                    }
                }
            }
            Ok(unresolved == 0)
        };

        for line in &probe_lines {
            if process_line(line)? {
                return Ok(());
            }
        }

        loop {
            buf.clear();
            let bytes = reader.read_line(&mut buf).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype stream {}: {err}",
                    backend.path.display()
                ))
            })?;
            if bytes == 0 {
                break;
            }
            if process_line(buf.trim_end_matches(['\n', '\r']))? {
                break;
            }
        }
        Ok(())
    };

    match backend.format {
        GenotypeSourceFormat::Text => {
            let file = File::open(&backend.path).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype file {}: {err}",
                    backend.path.display()
                ))
            })?;
            let mut reader = BufReader::new(file);
            scan_reader(&mut reader)?;
        }
        GenotypeSourceFormat::Zip => {
            let entry_name = backend.zip_entry_name.as_ref().ok_or_else(|| {
                RuntimeError::Unsupported(format!(
                    "zip backend missing selected entry for {}",
                    backend.path.display()
                ))
            })?;
            let file = File::open(&backend.path).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype zip {}: {err}",
                    backend.path.display()
                ))
            })?;
            let mut archive = ZipArchive::new(file).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype zip {}: {err}",
                    backend.path.display()
                ))
            })?;
            let entry = archive.by_name(entry_name).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype entry {entry_name} in {}: {err}",
                    backend.path.display()
                ))
            })?;
            let mut reader = BufReader::new(entry);
            scan_reader(&mut reader)?;
        }
        _ => {
            return Err(RuntimeError::Unsupported(
                "streaming delimited backend only supports text and zip".to_owned(),
            ));
        }
    }

    for (idx, variant) in indexed {
        if results[idx].genotype.is_none() {
            results[idx] = VariantObservation {
                backend: backend.backend_name().to_owned(),
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
