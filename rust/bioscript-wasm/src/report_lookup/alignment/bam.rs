use super::*;

pub(crate) struct BamReportLookup<R: std::io::Read + std::io::Seek> {
    pub(crate) reader: std::cell::RefCell<
        noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    >,
    pub(crate) label: String,
}

impl<R: std::io::Read + std::io::Seek> report_workspace::VariantLookup for BamReportLookup<R> {
    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(observe_bam_variant(&mut reader, &self.label, spec)?);
        }
        Ok(out)
    }
}

/// Build a minimal 23andMe-style text from the observations we already
/// computed. Format: `rsid\tchrom\tpos\tgenotype` per line. The runtime's
/// delimited-text loader reads this back as a `RsidMap`/`Delimited` backend
/// so analysis scripts can call `bioscript.load_genotypes(input_file)` and
/// have rsid lookups answered from the cached table.
#[allow(dead_code)]
fn synthesize_genotype_text_from_observations(observations: &[serde_json::Value]) -> String {
    let mut out = String::from("# rsid\tchromosome\tposition\tgenotype\n");
    for observation in observations {
        let rsid = observation
            .get("rsid")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("");
        if rsid.is_empty() {
            continue;
        }
        let chrom = observation
            .get("chrom")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("");
        let pos = observation
            .get("pos_start")
            .and_then(|v| {
                v.as_i64()
                    .or_else(|| v.as_str().and_then(|s| s.parse::<i64>().ok()))
            })
            .unwrap_or(0);
        let genotype = observation
            .get("genotype_display")
            .and_then(serde_json::Value::as_str)
            .filter(|s| !s.is_empty() && *s != "??")
            .unwrap_or("--");
        out.push_str(&format!("{rsid}\t{chrom}\t{pos}\t{genotype}\n"));
    }
    out
}

fn observe_bam_variant<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    label: &str,
    variant: &VariantSpec,
) -> Result<VariantObservation, RuntimeError> {
    let assembly = variant
        .grch38
        .as_ref()
        .map(|_| Assembly::Grch38)
        .or_else(|| variant.grch37.as_ref().map(|_| Assembly::Grch37));
    let locus = variant
        .grch38
        .as_ref()
        .or(variant.grch37.as_ref())
        .ok_or_else(|| {
            RuntimeError::Io(format!(
                "variant {} has no GRCh37/GRCh38 locus",
                variant
                    .rsids
                    .first()
                    .map(|s| s.as_str())
                    .unwrap_or("variant")
            ))
        })?;
    let locus = GenomicLocus {
        chrom: locus.chrom.clone(),
        start: locus.start,
        end: locus.end,
    };
    match variant.kind.unwrap_or(VariantKind::Snp) {
        VariantKind::Snp => {
            let ref_char = variant
                .reference
                .as_deref()
                .and_then(|s| s.chars().next())
                .ok_or_else(|| {
                    RuntimeError::Io(format!(
                        "variant {} missing reference allele",
                        variant
                            .rsids
                            .first()
                            .map(|s| s.as_str())
                            .unwrap_or("variant")
                    ))
                })?;
            let alt_char = variant
                .alternate
                .as_deref()
                .and_then(|s| s.chars().next())
                .ok_or_else(|| {
                    RuntimeError::Io(format!(
                        "variant {} missing alternate allele",
                        variant
                            .rsids
                            .first()
                            .map(|s| s.as_str())
                            .unwrap_or("variant")
                    ))
                })?;
            observe_bam_snp_with_reader(
                reader,
                label,
                &locus,
                ref_char,
                alt_char,
                variant.rsids.first().cloned(),
                assembly,
            )
        }
        VariantKind::Deletion => {
            observe_bam_deletion_with_reader(reader, label, &locus, variant, assembly)
        }
        VariantKind::Insertion | VariantKind::Indel => {
            let reference = variant.reference.as_deref().ok_or_else(|| {
                RuntimeError::Io(format!(
                    "variant {} missing reference allele",
                    variant
                        .rsids
                        .first()
                        .map(|s| s.as_str())
                        .unwrap_or("variant")
                ))
            })?;
            let alternate = variant.alternate.as_deref().ok_or_else(|| {
                RuntimeError::Io(format!(
                    "variant {} missing alternate allele",
                    variant
                        .rsids
                        .first()
                        .map(|s| s.as_str())
                        .unwrap_or("variant")
                ))
            })?;
            observe_bam_indel_with_reader(
                reader,
                label,
                &locus,
                reference,
                alternate,
                variant.rsids.first().cloned(),
                assembly,
            )
        }
        other => Err(RuntimeError::Io(format!(
            "variant {} kind {:?} not supported on BAM via wasm",
            variant
                .rsids
                .first()
                .map(|s| s.as_str())
                .unwrap_or("variant"),
            other
        ))),
    }
}

fn observe_bam_snp_with_reader<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    use noodles::core::Position;
    let mut counts = BamSnpPileupCounts::default();
    let header = read_bam_header(reader, label)?;
    let region = bam_region(&header, locus)?;
    let target_position = Position::try_from(usize::try_from(locus.start).map_err(|_| {
        RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned())
    })?)
    .map_err(|_| RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned()))?;

    let query = reader
        .query(&header, &region)
        .map_err(|err| RuntimeError::Io(format!("failed to query BAM {label}: {err}")))?;
    for result in query.records() {
        let record = result
            .map_err(|err| RuntimeError::Io(format!("failed to read BAM record {label}: {err}")))?;
        let flags = record.flags();
        if flags.is_unmapped() {
            counts.filtered_unmapped += 1;
            continue;
        }
        if flags.is_secondary() {
            counts.filtered_secondary += 1;
            continue;
        }
        if flags.is_qc_fail() {
            counts.filtered_qc_fail += 1;
            continue;
        }
        if flags.is_duplicate() {
            counts.filtered_duplicate += 1;
            continue;
        }
        if flags.is_segmented() && !flags.is_properly_segmented() {
            counts.filtered_improper_pair += 1;
            continue;
        }

        let Some((base, base_quality)) =
            bam_base_quality_at_reference_position(&record, target_position)?
        else {
            continue;
        };
        let normalized_base = normalize_pileup_base(base);
        let is_reverse = flags.is_reverse_complemented();
        if let Some(base) = normalized_base {
            counts.raw_depth += 1;
            *counts.raw_base_counts.entry(base.to_string()).or_insert(0) += 1;
            let strand_counts = if is_reverse {
                &mut counts.raw_reverse_counts
            } else {
                &mut counts.raw_forward_counts
            };
            *strand_counts.entry(base.to_string()).or_insert(0) += 1;
            if base == reference {
                counts.raw_ref_count += 1;
            } else if base == alternate {
                counts.raw_alt_count += 1;
            }
        }

        if base_quality < 13 {
            counts.filtered_low_base_quality += 1;
            continue;
        }

        let Some(base) = normalized_base else {
            counts.filtered_non_acgt += 1;
            continue;
        };

        counts.filtered_depth += 1;
        *counts
            .filtered_base_counts
            .entry(base.to_string())
            .or_insert(0) += 1;
        if base == reference {
            counts.filtered_ref_count += 1;
        } else if base == alternate {
            counts.filtered_alt_count += 1;
        }
    }

    let ref_count = counts.filtered_ref_count;
    let alt_count = counts.filtered_alt_count;
    let depth = counts.filtered_depth;
    let evidence = counts.evidence_lines(
        &format!("{}:{}-{}", locus.chrom, locus.start, locus.end),
        locus.start,
    );
    Ok(VariantObservation {
        backend: "bam".to_owned(),
        matched_rsid,
        assembly,
        genotype: infer_snp_genotype(reference, alternate, ref_count, alt_count, depth),
        ref_count: Some(ref_count),
        alt_count: Some(alt_count),
        depth: Some(depth),
        raw_counts: counts.raw_base_counts,
        decision: Some(describe_snp_decision_rule(
            reference, alternate, ref_count, alt_count, depth,
        )),
        evidence,
    })
}

fn observe_bam_deletion_with_reader<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    label: &str,
    locus: &GenomicLocus,
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    let deletion_length = variant.deletion_length.ok_or_else(|| {
        RuntimeError::InvalidArguments("deletion variant requires deletion_length".to_owned())
    })?;
    let reference = variant.reference.clone().unwrap_or_else(|| "I".to_owned());
    let alternate = variant.alternate.clone().unwrap_or_else(|| "D".to_owned());
    let anchor_pos = locus.start.saturating_sub(1);
    let anchor_locus = GenomicLocus {
        chrom: locus.chrom.clone(),
        start: anchor_pos,
        end: anchor_pos,
    };

    let mut alt_count = 0u32;
    let mut ref_count = 0u32;
    let mut depth = 0u32;

    let header = read_bam_header(reader, label)?;
    let region = bam_region(&header, &anchor_locus)?;
    let query = reader
        .query(&header, &region)
        .map_err(|err| RuntimeError::Io(format!("failed to query BAM {label}: {err}")))?;
    for result in query.records() {
        let record = result
            .map_err(|err| RuntimeError::Io(format!("failed to read BAM record {label}: {err}")))?;
        let alignment_record = bam_alignment_record(label, &record)?;
        if alignment_record.is_unmapped || !spans_position(&alignment_record, anchor_pos) {
            continue;
        }
        depth += 1;
        match indel_at_anchor(&alignment_record, anchor_pos) {
            Some((bioscript_formats::alignment::AlignmentOpKind::Deletion, len))
                if len == deletion_length =>
            {
                alt_count += 1;
            }
            _ => ref_count += 1,
        }
    }

    Ok(VariantObservation {
        backend: "bam".to_owned(),
        matched_rsid: variant.rsids.first().cloned(),
        assembly,
        genotype: infer_copy_number_genotype(&reference, &alternate, ref_count, alt_count, depth),
        ref_count: Some(ref_count),
        alt_count: Some(alt_count),
        depth: Some(depth),
        raw_counts: BTreeMap::new(),
        decision: Some(describe_copy_number_decision_rule(
            &reference, &alternate, ref_count, alt_count, depth,
        )),
        evidence: vec![format!(
            "observed BAM deletion anchor {}:{} len={} depth={} ref_count={} alt_count={}",
            locus.chrom, anchor_pos, deletion_length, depth, ref_count, alt_count
        )],
    })
}

fn observe_bam_indel_with_reader<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    label: &str,
    locus: &GenomicLocus,
    reference: &str,
    alternate: &str,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    let mut alt_count = 0u32;
    let mut ref_count = 0u32;
    let mut depth = 0u32;
    let mut matching_alt_lengths = std::collections::BTreeSet::new();

    let header = read_bam_header(reader, label)?;
    let region = bam_region(&header, locus)?;
    let query = reader
        .query(&header, &region)
        .map_err(|err| RuntimeError::Io(format!("failed to query BAM {label}: {err}")))?;
    for result in query.records() {
        let record = result
            .map_err(|err| RuntimeError::Io(format!("failed to read BAM record {label}: {err}")))?;
        let alignment_record = bam_alignment_record(label, &record)?;
        if alignment_record.is_unmapped || !record_overlaps_locus(&alignment_record, locus) {
            continue;
        }
        let classification =
            classify_expected_indel(&alignment_record, locus, reference.len(), alternate)?;
        if !classification.covering {
            continue;
        }
        depth += 1;
        if classification.matches_alt {
            alt_count += 1;
            matching_alt_lengths.insert(classification.observed_len);
        } else if classification.reference_like {
            ref_count += 1;
        }
    }

    let evidence_label = if matching_alt_lengths.is_empty() {
        "none".to_owned()
    } else {
        matching_alt_lengths
            .into_iter()
            .map(|len| len.to_string())
            .collect::<Vec<_>>()
            .join(",")
    };

    Ok(VariantObservation {
        backend: "bam".to_owned(),
        matched_rsid,
        assembly,
        genotype: infer_copy_number_genotype(reference, alternate, ref_count, alt_count, depth),
        ref_count: Some(ref_count),
        alt_count: Some(alt_count),
        depth: Some(depth),
        raw_counts: BTreeMap::new(),
        decision: Some(describe_copy_number_decision_rule(
            reference, alternate, ref_count, alt_count, depth,
        )),
        evidence: vec![format!(
            "observed BAM indel at {}:{}-{} depth={} ref_count={} alt_count={} matching_alt_lengths={}",
            locus.chrom, locus.start, locus.end, depth, ref_count, alt_count, evidence_label
        )],
    })
}

fn read_bam_header<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    label: &str,
) -> Result<noodles::sam::Header, RuntimeError> {
    reader
        .get_mut()
        .seek(noodles::bgzf::VirtualPosition::MIN)
        .map_err(|err| RuntimeError::Io(format!("failed to rewind BAM {label}: {err}")))?;
    reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header {label}: {err}")))
}

fn bam_region(
    header: &noodles::sam::Header,
    locus: &GenomicLocus,
) -> Result<noodles::core::Region, RuntimeError> {
    let chrom = resolve_bam_reference_name(header, &locus.chrom).ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "indexed BAM does not contain contig {} for {}:{}-{}",
            locus.chrom, locus.chrom, locus.start, locus.end
        ))
    })?;
    format!("{chrom}:{}-{}", locus.start, locus.end)
        .parse()
        .map_err(|err| RuntimeError::Io(format!("invalid BAM query region: {err}")))
}

fn resolve_bam_reference_name(header: &noodles::sam::Header, chrom: &str) -> Option<String> {
    let candidates = [
        chrom.to_owned(),
        format!("chr{chrom}"),
        chrom.trim_start_matches("chr").to_owned(),
    ];
    candidates.into_iter().find(|candidate| {
        header.reference_sequences().iter().any(|(name, _)| {
            let name_bytes: &[u8] = name.as_ref();
            name_bytes == candidate.as_bytes()
        })
    })
}

#[path = "bam/pileup.rs"]
mod pileup;

use pileup::{
    BamSnpPileupCounts, bam_alignment_record, bam_base_quality_at_reference_position,
    classify_expected_indel, describe_copy_number_decision_rule, describe_snp_decision_rule,
    indel_at_anchor, infer_copy_number_genotype, infer_snp_genotype, normalize_pileup_base,
    record_overlaps_locus, spans_position,
};
