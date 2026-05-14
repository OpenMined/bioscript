use super::*;

/// Per-variant CRAM lookup that satisfies the workspace's `VariantLookup`
/// trait. Holds the IndexedReader in a `RefCell` so &self lookup methods can
/// mutably read while still being object-safe.
pub(super) struct CramReportLookup<R: std::io::Read + std::io::Seek> {
    pub(super) reader: std::cell::RefCell<noodles::cram::io::indexed_reader::IndexedReader<R>>,
    pub(super) label: String,
}

pub(super) struct BamReportLookup<R: std::io::Read + std::io::Seek> {
    pub(super) reader: std::cell::RefCell<
        noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    >,
    pub(super) label: String,
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

impl<R: std::io::Read + std::io::Seek> report_workspace::VariantLookup for CramReportLookup<R> {
    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(observe_cram_variant(&mut reader, &self.label, spec)?);
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

fn bam_base_quality_at_reference_position(
    record: &noodles::bam::Record,
    target_position: noodles::core::Position,
) -> Result<Option<(u8, u8)>, RuntimeError> {
    use noodles::sam::alignment::record::cigar::op::Kind;

    let Some(alignment_start) = record
        .alignment_start()
        .transpose()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM alignment start: {err}")))?
    else {
        return Ok(None);
    };
    let mut reference_position = usize::from(alignment_start);
    let target = usize::from(target_position);
    let mut read_position = 0usize;
    let sequence = record.sequence();
    let qualities = record.quality_scores();

    for result in record.cigar().iter() {
        let op =
            result.map_err(|err| RuntimeError::Io(format!("failed to read BAM CIGAR: {err}")))?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                if target >= reference_position && target < reference_position + len {
                    let offset = target - reference_position;
                    let read_index = read_position + offset;
                    let Some(base) = sequence.get(read_index) else {
                        return Ok(None);
                    };
                    let quality = qualities
                        .as_ref()
                        .get(read_index)
                        .copied()
                        .unwrap_or(u8::MAX);
                    return Ok(Some((base, quality)));
                }
                reference_position += len;
                read_position += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                read_position += len;
            }
            Kind::Deletion | Kind::Skip => {
                if target >= reference_position && target < reference_position + len {
                    return Ok(None);
                }
                reference_position += len;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }

    Ok(None)
}

fn bam_alignment_record(
    label: &str,
    record: &noodles::bam::Record,
) -> Result<bioscript_formats::alignment::AlignmentRecord, RuntimeError> {
    use noodles::sam::alignment::Record as _;

    let flags = record.flags();
    let is_unmapped = flags.is_unmapped();
    let start = record
        .alignment_start()
        .transpose()
        .map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read BAM alignment start from {label}: {err}"
            ))
        })?
        .map(|pos| i64::try_from(usize::from(pos)))
        .transpose()
        .map_err(|_| {
            RuntimeError::Unsupported(format!(
                "record alignment start exceeds i64 range in {label}"
            ))
        })?
        .unwrap_or(0);
    let end = record
        .alignment_end()
        .transpose()
        .map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read BAM alignment end from {label}: {err}"
            ))
        })?
        .map(|pos| i64::try_from(usize::from(pos)))
        .transpose()
        .map_err(|_| {
            RuntimeError::Unsupported(format!("record alignment end exceeds i64 range in {label}"))
        })?
        .unwrap_or(start);
    let cigar = record
        .cigar()
        .iter()
        .map(|result| {
            result.map(map_bam_op).map_err(|err| {
                RuntimeError::Io(format!("failed to read BAM CIGAR from {label}: {err}"))
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(bioscript_formats::alignment::AlignmentRecord {
        start,
        end,
        is_unmapped,
        cigar,
    })
}

fn map_bam_op(
    op: noodles::sam::alignment::record::cigar::Op,
) -> bioscript_formats::alignment::AlignmentOp {
    use bioscript_formats::alignment::{AlignmentOp, AlignmentOpKind};
    use noodles::sam::alignment::record::cigar::op::Kind;

    let kind = match op.kind() {
        Kind::Match => AlignmentOpKind::Match,
        Kind::Insertion => AlignmentOpKind::Insertion,
        Kind::Deletion => AlignmentOpKind::Deletion,
        Kind::Skip => AlignmentOpKind::Skip,
        Kind::SoftClip => AlignmentOpKind::SoftClip,
        Kind::HardClip => AlignmentOpKind::HardClip,
        Kind::Pad => AlignmentOpKind::Pad,
        Kind::SequenceMatch => AlignmentOpKind::SequenceMatch,
        Kind::SequenceMismatch => AlignmentOpKind::SequenceMismatch,
    };

    AlignmentOp {
        kind,
        len: op.len(),
    }
}

fn normalize_pileup_base(base: u8) -> Option<char> {
    match (base as char).to_ascii_uppercase() {
        'A' | 'C' | 'G' | 'T' => Some((base as char).to_ascii_uppercase()),
        _ => None,
    }
}

struct IndelClassification {
    covering: bool,
    reference_like: bool,
    matches_alt: bool,
    observed_len: usize,
}

fn record_overlaps_locus(
    record: &bioscript_formats::alignment::AlignmentRecord,
    locus: &GenomicLocus,
) -> bool {
    record.end >= locus.start && record.start <= locus.end
}

fn spans_position(record: &bioscript_formats::alignment::AlignmentRecord, pos: i64) -> bool {
    pos >= record.start.saturating_sub(1) && pos <= record.end
}

fn indel_at_anchor(
    record: &bioscript_formats::alignment::AlignmentRecord,
    anchor_pos: i64,
) -> Option<(bioscript_formats::alignment::AlignmentOpKind, usize)> {
    let mut ref_pos = record.start;

    for op in &record.cigar {
        match op.kind {
            bioscript_formats::alignment::AlignmentOpKind::Match
            | bioscript_formats::alignment::AlignmentOpKind::SequenceMatch
            | bioscript_formats::alignment::AlignmentOpKind::SequenceMismatch
            | bioscript_formats::alignment::AlignmentOpKind::Skip => {
                ref_pos += i64::try_from(op.len).ok()?;
            }
            bioscript_formats::alignment::AlignmentOpKind::Insertion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((
                        bioscript_formats::alignment::AlignmentOpKind::Insertion,
                        op.len,
                    ));
                }
            }
            bioscript_formats::alignment::AlignmentOpKind::Deletion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((
                        bioscript_formats::alignment::AlignmentOpKind::Deletion,
                        op.len,
                    ));
                }
                ref_pos += i64::try_from(op.len).ok()?;
            }
            bioscript_formats::alignment::AlignmentOpKind::SoftClip
            | bioscript_formats::alignment::AlignmentOpKind::HardClip
            | bioscript_formats::alignment::AlignmentOpKind::Pad => {}
        }
    }

    None
}

fn classify_expected_indel(
    record: &bioscript_formats::alignment::AlignmentRecord,
    locus: &GenomicLocus,
    reference_len: usize,
    alternate: &str,
) -> Result<IndelClassification, RuntimeError> {
    let alt_len = alternate.len();
    let anchor_start = locus.start.saturating_sub(1);
    let anchor_end = locus.end;

    let covering = record.start <= locus.start && record.end >= locus.end;
    if !covering {
        return Ok(IndelClassification {
            covering: false,
            reference_like: false,
            matches_alt: false,
            observed_len: reference_len,
        });
    }

    let mut observed_len = reference_len;

    for anchor in anchor_start..=anchor_end {
        if let Some((kind, len)) = indel_at_anchor(record, anchor) {
            observed_len = match kind {
                bioscript_formats::alignment::AlignmentOpKind::Insertion => reference_len + len,
                bioscript_formats::alignment::AlignmentOpKind::Deletion => {
                    reference_len.saturating_sub(len)
                }
                _ => reference_len,
            };

            return Ok(IndelClassification {
                covering: true,
                reference_like: false,
                matches_alt: observed_len == alt_len,
                observed_len,
            });
        }
    }

    Ok(IndelClassification {
        covering: true,
        reference_like: true,
        matches_alt: false,
        observed_len,
    })
}

#[derive(Default)]
struct BamSnpPileupCounts {
    filtered_depth: u32,
    filtered_ref_count: u32,
    filtered_alt_count: u32,
    filtered_base_counts: BTreeMap<String, u32>,
    raw_depth: u32,
    raw_ref_count: u32,
    raw_alt_count: u32,
    raw_base_counts: BTreeMap<String, u32>,
    filtered_low_base_quality: u32,
    filtered_low_mapping_quality: u32,
    filtered_non_acgt: u32,
    filtered_unmapped: u32,
    filtered_secondary: u32,
    filtered_qc_fail: u32,
    filtered_duplicate: u32,
    filtered_improper_pair: u32,
    raw_forward_counts: BTreeMap<String, u32>,
    raw_reverse_counts: BTreeMap<String, u32>,
}

impl BamSnpPileupCounts {
    fn evidence_lines(&self, locus: &str, target_pos: i64) -> Vec<String> {
        vec![
            format!(
                "observed BAM SNP pileup at {locus} target_pos={target_pos} filtered_depth={} ref_count={} alt_count={}",
                self.filtered_depth, self.filtered_ref_count, self.filtered_alt_count
            ),
            format!(
                "raw pileup depth={} ref_count={} alt_count={} raw_counts={:?}",
                self.raw_depth, self.raw_ref_count, self.raw_alt_count, self.raw_base_counts
            ),
            format!(
                "raw strand counts: forward={:?} reverse={:?}",
                self.raw_forward_counts, self.raw_reverse_counts
            ),
            format!(
                "filters applied: min_base_quality=13 min_mapping_quality=0 filtered_low_base_quality={} filtered_low_mapping_quality={} filtered_non_acgt={} filtered_unmapped={} filtered_secondary={} filtered_qc_fail={} filtered_duplicate={} filtered_improper_pair={}",
                self.filtered_low_base_quality,
                self.filtered_low_mapping_quality,
                self.filtered_non_acgt,
                self.filtered_unmapped,
                self.filtered_secondary,
                self.filtered_qc_fail,
                self.filtered_duplicate,
                self.filtered_improper_pair
            ),
        ]
    }
}

fn infer_snp_genotype(
    reference: char,
    alternate: char,
    ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> Option<String> {
    if depth == 0 || ref_count + alt_count == 0 {
        return None;
    }
    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    if alt_fraction >= 0.8 {
        Some(format!("{alternate}{alternate}"))
    } else if alt_fraction <= 0.2 {
        Some(format!("{reference}{reference}"))
    } else {
        let mut alleles = [
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase(),
        ];
        alleles.sort_by_key(|allele| match allele {
            'A' => 0,
            'C' => 1,
            'G' => 2,
            'T' => 3,
            _ => 99,
        });
        Some(alleles.iter().collect())
    }
}

fn describe_snp_decision_rule(
    reference: char,
    alternate: char,
    ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> String {
    if depth == 0 {
        return format!(
            "no covering reads for SNP; genotype unresolved (ref={reference}, alt={alternate})"
        );
    }
    if ref_count + alt_count == 0 {
        return format!(
            "no reads matched the declared SNP alleles; genotype unresolved; counts ref={ref_count} alt={alt_count} depth={depth} for {reference}>{alternate}"
        );
    }

    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    format!(
        "SNP genotype rule: alt_fraction={alt_fraction:.3} with thresholds ref<=0.200, het=(0.200,0.800), alt>=0.800; counts ref={ref_count} alt={alt_count} depth={depth} for {reference}>{alternate}"
    )
}

fn infer_copy_number_genotype(
    reference: &str,
    alternate: &str,
    _ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> Option<String> {
    if depth == 0 {
        return None;
    }
    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    if alt_fraction >= 0.8 {
        Some(format!("{alternate}{alternate}"))
    } else if alt_fraction <= 0.2 {
        Some(format!("{reference}{reference}"))
    } else {
        let mut alleles = [
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase(),
        ];
        alleles.sort_by_key(|allele| {
            allele.chars().next().map_or(u8::MAX, |ch| match ch {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'T' => 3,
                'I' => 4,
                'D' => 5,
                _ => 99,
            })
        });
        Some(alleles.concat())
    }
}

fn describe_copy_number_decision_rule(
    reference: &str,
    alternate: &str,
    _ref_count: u32,
    alt_count: u32,
    depth: u32,
) -> String {
    if depth == 0 {
        return format!(
            "no covering reads for copy-number style variant; genotype unresolved (ref={reference}, alt={alternate})"
        );
    }

    let alt_fraction = f64::from(alt_count) / f64::from(depth);
    format!(
        "copy-number genotype rule: alt_fraction={alt_fraction:.3} with thresholds ref<=0.200, het=(0.200,0.800), alt>=0.800; counts alt={alt_count} depth={depth} for {reference}->{alternate}"
    )
}

fn observe_cram_variant<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::cram::io::indexed_reader::IndexedReader<R>,
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
    let kind = variant.kind.unwrap_or(VariantKind::Snp);
    match kind {
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
            bioscript_formats::observe_cram_snp_with_reader(
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
            let deletion_length = variant.deletion_length.ok_or_else(|| {
                RuntimeError::Io(format!(
                    "variant {} missing deletion_length",
                    variant
                        .rsids
                        .first()
                        .map(|s| s.as_str())
                        .unwrap_or("variant")
                ))
            })?;
            let reference = variant.reference.as_deref().unwrap_or("I");
            let alternate = variant.alternate.as_deref().unwrap_or("D");
            bioscript_formats::observe_cram_deletion_with_reader(
                reader,
                label,
                &locus,
                deletion_length,
                reference,
                alternate,
                variant.rsids.first().cloned(),
                assembly,
            )
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
            bioscript_formats::observe_cram_indel_with_reader(
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
            "variant {} kind {:?} not supported on CRAM via wasm",
            variant
                .rsids
                .first()
                .map(|s| s.as_str())
                .unwrap_or("variant"),
            other
        ))),
    }
}

pub(super) struct VcfReportLookup<R: std::io::Read + std::io::Seek> {
    pub(super) reader: std::cell::RefCell<
        noodles::csi::io::IndexedReader<noodles::bgzf::io::Reader<R>, noodles::tabix::Index>,
    >,
    pub(super) label: String,
    /// Assembly resolved from the VCF header (via `inspect_head_via_js_reader`).
    /// Matches the CLI's `lookup_indexed_vcf_variants` flow which calls
    /// `detect_vcf_assembly_from_path`. Without this the wasm picks GRCh38
    /// over GRCh37 for any panel variant that declares both loci, then
    /// misses the variant in a GRCh37-coded VCF (NA06985.clean.vcf.gz etc.)
    /// and falls through to "imputed reference".
    pub(super) detected_assembly: Option<Assembly>,
}

impl<R: std::io::Read + std::io::Seek> report_workspace::VariantLookup for VcfReportLookup<R> {
    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(observe_vcf_variant(
                &mut reader,
                &self.label,
                spec,
                self.detected_assembly,
            )?);
        }
        Ok(out)
    }
}

fn observe_vcf_variant<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::csi::io::IndexedReader<
        noodles::bgzf::io::Reader<R>,
        noodles::tabix::Index,
    >,
    label: &str,
    variant: &VariantSpec,
    detected_assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    // Use the existing CLI helper so wasm picks the same locus the path-based
    // path does: detected GRCh37 → grch37 first; detected GRCh38 → grch38
    // first; None → grch37 first (CLI default for variant-only VCFs).
    let raw_locus =
        bioscript_formats::choose_variant_locus_for_assembly(variant, detected_assembly)
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
    let assembly = detected_assembly.or_else(|| {
        if variant.grch37.as_ref().is_some_and(|l| l == &raw_locus) {
            Some(Assembly::Grch37)
        } else if variant.grch38.as_ref().is_some_and(|l| l == &raw_locus) {
            Some(Assembly::Grch38)
        } else {
            None
        }
    });
    let locus = GenomicLocus {
        chrom: raw_locus.chrom.clone(),
        start: raw_locus.start,
        end: raw_locus.end,
    };
    let observation = bioscript_formats::observe_vcf_variant_with_reader(
        reader,
        label,
        &locus,
        variant,
        variant.rsids.first().cloned(),
        assembly,
    )?;
    // Mirror the CLI report flow's `impute_vcf_missing_as_reference: true`
    // default. The full-file VCF scanner only marks a variant resolved when a
    // row actually matches the query. Unrelated rows in the indexed window
    // must not block reference imputation for absent variant-only calls.
    if observation.genotype.is_none()
        && !observation.evidence.iter().any(|line| {
            line.contains("tabix index has no contig")
                || line.contains("has no GRCh37/GRCh38 locus")
        })
    {
        if let Some(imputed) = bioscript_formats::imputed_reference_observation(
            "vcf",
            label,
            variant,
            &locus,
            assembly,
            None,
            &observation.evidence.join(" | "),
        ) {
            return Ok(imputed);
        }
    }
    Ok(observation)
}
