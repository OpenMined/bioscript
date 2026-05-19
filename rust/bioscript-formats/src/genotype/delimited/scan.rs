use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
};

use zip::ZipArchive;

use crate::inspect::{AssemblyAnchorScorer, detect_assembly};
use bioscript_core::{RuntimeError, VariantObservation, VariantSpec};

use super::{
    super::{
        GenotypeSourceFormat, backends::delimited_locus_for_assembly, describe_query,
        types::DelimitedBackend, variant_sort_key,
    },
    DelimitedColumnIndexes, GsgtParser, detect_delimiter, is_no_call as gsgt_is_no_call,
    lines_look_like_gsgt, parse_streaming_row,
};
use bioscript_core::Assembly;

/// Apply a candidate genotype to a result slot.
///
/// Non-GSGT keeps the historical first-match-wins behaviour. GSGT merges
/// replicate probes in a single pass (spec §5): a real call replaces an
/// earlier no-call; two disagreeing real calls collapse to a no-call; equal
/// calls and no-call-after-call are ignored.
#[allow(clippy::too_many_arguments)]
fn apply_match(
    slot: &mut VariantObservation,
    unresolved: &mut usize,
    gsgt: bool,
    matched_rsid: Option<String>,
    genotype: &str,
    evidence_head: &str,
    raw_line: &str,
    backend_name: &str,
) {
    let make = |g: &str, rsid: Option<String>| VariantObservation {
        backend: backend_name.to_owned(),
        matched_rsid: rsid,
        genotype: Some(g.to_owned()),
        evidence: vec![evidence_head.to_owned(), format!("source line: {raw_line}")],
        ..VariantObservation::default()
    };

    if slot.genotype.is_none() {
        *slot = make(genotype, matched_rsid);
        *unresolved = unresolved.saturating_sub(1);
        return;
    }
    if !gsgt {
        return;
    }
    let existing = slot.genotype.clone().unwrap_or_default();
    let prev_rsid = slot.matched_rsid.clone();
    let existing_nc = gsgt_is_no_call(&existing);
    let new_nc = gsgt_is_no_call(genotype);
    if existing_nc && !new_nc {
        *slot = make(genotype, matched_rsid);
    } else if !existing_nc && !new_nc && existing != genotype {
        *slot = make("--", prev_rsid);
    }
}

/// One bounded pass over the file feeding the rsID/locus anchor vote, used
/// only when no build metadata was found. Mirrors the main scan's open
/// logic; reuses the GSGT / delimited row parsers.
fn prescan_assembly_anchors(
    backend: &DelimitedBackend,
    is_gsgt: bool,
) -> Result<Option<Assembly>, RuntimeError> {
    fn vote<R: BufRead>(mut reader: R, is_gsgt: bool) -> Result<Option<Assembly>, RuntimeError> {
        let mut scorer = AssemblyAnchorScorer::new();
        let mut probe = Vec::new();
        let mut buf = String::new();
        for _ in 0..8 {
            buf.clear();
            if reader.read_line(&mut buf).unwrap_or(0) == 0 {
                break;
            }
            probe.push(buf.trim_end_matches(['\n', '\r']).to_owned());
        }
        let delimiter = detect_delimiter(&probe);
        let mut ci: Option<DelimitedColumnIndexes> = None;
        let mut ch: Option<Vec<String>> = None;
        let mut gsgt = GsgtParser::new();
        let feed = |line: &str,
                    scorer: &mut AssemblyAnchorScorer,
                    gsgt: &mut GsgtParser,
                    ci: &mut Option<DelimitedColumnIndexes>,
                    ch: &mut Option<Vec<String>>|
         -> Result<(), RuntimeError> {
            let row = if is_gsgt {
                gsgt.consume(line)?
            } else {
                parse_streaming_row(line, delimiter, ci, ch)?
            };
            if let Some(row) = row
                && let (Some(c), Some(p)) = (row.chrom.as_ref(), row.position)
            {
                scorer.observe(row.rsid.as_deref().unwrap_or(""), c, p, &row.genotype);
            }
            Ok(())
        };
        for line in &probe {
            feed(line, &mut scorer, &mut gsgt, &mut ci, &mut ch)?;
        }
        loop {
            buf.clear();
            if reader
                .read_line(&mut buf)
                .map_err(|err| RuntimeError::Io(format!("anchor prescan read failed: {err}")))?
                == 0
            {
                break;
            }
            feed(
                buf.trim_end_matches(['\n', '\r']),
                &mut scorer,
                &mut gsgt,
                &mut ci,
                &mut ch,
            )?;
        }
        Ok(scorer.decide())
    }

    match backend.format {
        GenotypeSourceFormat::Text => {
            let file = File::open(&backend.path).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype file {}: {err}",
                    backend.path.display()
                ))
            })?;
            vote(BufReader::new(file), is_gsgt)
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
            vote(BufReader::new(entry), is_gsgt)
        }
        _ => Ok(None),
    }
}

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

    let has_coordinate_queries = variants
        .iter()
        .any(|variant| variant.has_coordinates() && !variant.has_rsids());
    let mut detected_assembly = backend.options.assembly;

    for (idx, variant) in &indexed {
        for rsid in &variant.rsids {
            rsid_targets.entry(rsid.clone()).or_default().push(*idx);
        }
        if let Some(locus) = delimited_locus_for_assembly(variant, detected_assembly) {
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
        let is_gsgt = lines_look_like_gsgt(&probe_lines);
        if detected_assembly.is_none() {
            let mut label = backend.path.to_string_lossy().to_ascii_lowercase();
            if let Some(entry_name) = backend.zip_entry_name.as_ref() {
                label.push('\n');
                label.push_str(&entry_name.to_ascii_lowercase());
            }
            detected_assembly = detect_assembly(&label, &probe_lines);
            // No declared build (e.g. a GSGT Final Report): resolve it from
            // the rsID/locus anchor vote over the whole file instead of
            // assuming. Only runs when metadata gave us nothing.
            if detected_assembly.is_none() {
                detected_assembly = prescan_assembly_anchors(backend, is_gsgt)?;
            }
            if let Some(assembly) = detected_assembly {
                for (idx, variant) in &indexed {
                    if let Some(locus) = delimited_locus_for_assembly(variant, Some(assembly)) {
                        coord_targets
                            .entry((
                                locus.chrom.trim_start_matches("chr").to_ascii_lowercase(),
                                locus.start,
                            ))
                            .or_default()
                            .push(*idx);
                    }
                }
            } else if has_coordinate_queries {
                return Err(RuntimeError::Unsupported(format!(
                    "delimited genotype input assembly is unknown for {}; refusing coordinate lookup",
                    backend.path.display()
                )));
            }
        }
        let mut column_indexes: Option<DelimitedColumnIndexes> = None;
        let mut comment_header: Option<Vec<String>> = None;
        let mut gsgt_parser = GsgtParser::new();
        let backend_name = backend.backend_name();

        let mut process_line = |line: &str| -> Result<bool, RuntimeError> {
            let row = if is_gsgt {
                gsgt_parser.consume(line)?
            } else {
                parse_streaming_row(line, delimiter, &mut column_indexes, &mut comment_header)?
            };
            let Some(row) = row else {
                // GSGT must scan the whole file so later replicate probes can
                // be merged; non-GSGT may stop once everything resolved.
                return Ok(!is_gsgt && unresolved == 0);
            };

            if let Some(rsid) = row.rsid.as_ref()
                && let Some(target_indexes) = rsid_targets.get(rsid)
            {
                for &target_idx in target_indexes {
                    apply_match(
                        &mut results[target_idx],
                        &mut unresolved,
                        is_gsgt,
                        Some(rsid.clone()),
                        &row.genotype,
                        &format!("resolved by rsid {rsid}"),
                        &row.raw_line,
                        backend_name,
                    );
                }
            }

            if !is_gsgt && unresolved == 0 {
                return Ok(true);
            }

            if let (Some(chrom), Some(position)) = (row.chrom.as_ref(), row.position) {
                let key = (
                    chrom.trim_start_matches("chr").to_ascii_lowercase(),
                    position,
                );
                if let Some(target_indexes) = coord_targets.get(&key) {
                    for &target_idx in target_indexes {
                        apply_match(
                            &mut results[target_idx],
                            &mut unresolved,
                            is_gsgt,
                            row.rsid.clone(),
                            &row.genotype,
                            &format!("resolved by locus {chrom}:{position}"),
                            &row.raw_line,
                            backend_name,
                        );
                    }
                }
            }
            Ok(!is_gsgt && unresolved == 0)
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
            let evidence = if variant.has_coordinates() {
                format!(
                    "no matching rsid or locus found for {}",
                    describe_query(variant)
                )
            } else {
                "no matching rsid found".to_owned()
            };
            results[idx] = VariantObservation {
                backend: backend.backend_name().to_owned(),
                evidence: vec![evidence],
                ..VariantObservation::default()
            };
        }
    }

    Ok(results)
}
