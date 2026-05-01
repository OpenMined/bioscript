use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt::Write as _,
    fs::File,
    io::{BufRead, BufReader, Cursor, Read, Seek},
    path::Path,
};

use noodles::bgzf;
use noodles::core::{Position, Region};
use noodles::cram;
use noodles::csi::{self, BinningIndex};
use noodles::sam::alignment::{
    Record as _,
    record::{Cigar as _, QualityScores as _, Sequence as _, cigar::op::Kind as CigarOpKind},
};
use noodles::tabix;
use zip::ZipArchive;

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};

use crate::alignment::{self, AlignmentOpKind, AlignmentRecord};

mod delimited;
mod io;
mod types;
mod vcf_tokens;

#[cfg(test)]
use delimited::{
    Delimiter, GENOTYPE_ALIASES, parse_streaming_row, split_csv_line, strip_bom,
    strip_inline_comment,
};
use delimited::{RowParser, detect_delimiter, scan_delimited_variants};
#[cfg(test)]
use delimited::{
    build_column_indexes, default_column_indexes, find_header_index, looks_like_header_fields,
    normalize_name,
};
#[cfg(test)]
use io::looks_like_vcf_lines;
use io::{
    detect_source_format, is_bgzf_path, read_lines_from_reader, read_zip_entry_limited,
    select_zip_entry,
};
pub use types::{
    BackendCapabilities, GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, QueryKind,
};
use types::{CramBackend, DelimitedBackend, QueryBackend, RsidMapBackend, VcfBackend};
use vcf_tokens::genotype_from_vcf_gt;
#[cfg(test)]
use vcf_tokens::{
    is_symbolic_vcf_alt, normalize_sequence_token, vcf_alt_token, vcf_reference_token,
};

const DEFAULT_MPILEUP_MIN_BASE_QUALITY: u8 = 13;
const DEFAULT_MPILEUP_MIN_MAPPING_QUALITY: u8 = 0;
const MAX_ZIP_ENTRY_BYTES: u64 = 128 * 1024 * 1024;

impl GenotypeStore {
    pub fn from_file(path: &Path) -> Result<Self, RuntimeError> {
        Self::from_file_with_options(path, &GenotypeLoadOptions::default())
    }

    pub fn from_file_with_options(
        path: &Path,
        options: &GenotypeLoadOptions,
    ) -> Result<Self, RuntimeError> {
        match detect_source_format(path, options.format)? {
            GenotypeSourceFormat::Text => Ok(Self::from_delimited_file(
                path,
                GenotypeSourceFormat::Text,
                None,
            )),
            GenotypeSourceFormat::Zip => Self::from_zip_file(path),
            GenotypeSourceFormat::Vcf => Ok(Self::from_vcf_file(path, options)),
            GenotypeSourceFormat::Cram => Self::from_cram_file(path, options),
        }
    }

    pub fn from_bytes(name: &str, bytes: &[u8]) -> Result<Self, RuntimeError> {
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".zip") {
            return Self::from_zip_bytes(name, bytes);
        }
        if lower.ends_with(".vcf") {
            let lines =
                read_lines_from_reader(BufReader::new(Cursor::new(bytes)), Path::new(name))?;
            return Self::from_vcf_lines(lines);
        }
        let lines = read_lines_from_reader(BufReader::new(Cursor::new(bytes)), Path::new(name))?;
        Self::from_delimited_lines(GenotypeSourceFormat::Text, lines)
    }

    fn from_zip_bytes(name: &str, bytes: &[u8]) -> Result<Self, RuntimeError> {
        let mut archive = ZipArchive::new(Cursor::new(bytes)).map_err(|err| {
            RuntimeError::Io(format!("failed to read genotype zip {name}: {err}"))
        })?;
        let mut selected = None;
        for idx in 0..archive.len() {
            let entry = archive.by_index(idx).map_err(|err| {
                RuntimeError::Io(format!("failed to inspect genotype zip {name}: {err}"))
            })?;
            if entry.is_dir() {
                continue;
            }
            let entry_name = entry.name().to_owned();
            let lower = entry_name.to_ascii_lowercase();
            if lower.ends_with(".vcf")
                || lower.ends_with(".txt")
                || lower.ends_with(".tsv")
                || lower.ends_with(".csv")
            {
                selected = Some(entry_name);
                break;
            }
        }
        let selected = selected.ok_or_else(|| {
            RuntimeError::Unsupported(format!(
                "zip archive {name} does not contain a supported genotype file"
            ))
        })?;
        let mut entry = archive.by_name(&selected).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open genotype entry {selected} in {name}: {err}"
            ))
        })?;
        let contents = read_zip_entry_limited(
            &mut entry,
            MAX_ZIP_ENTRY_BYTES,
            &format!("genotype entry {selected} in {name}"),
        )?;
        let lines =
            read_lines_from_reader(BufReader::new(Cursor::new(contents)), Path::new(&selected))?;
        if selected.to_ascii_lowercase().ends_with(".vcf") {
            return Self::from_vcf_lines(lines);
        }
        Self::from_delimited_lines(GenotypeSourceFormat::Zip, lines)
    }

    fn from_vcf_file(path: &Path, options: &GenotypeLoadOptions) -> Self {
        Self {
            backend: QueryBackend::Vcf(VcfBackend {
                path: path.to_path_buf(),
                options: options.clone(),
            }),
        }
    }

    fn from_zip_file(path: &Path) -> Result<Self, RuntimeError> {
        let selected = select_zip_entry(path)?;
        let lower = selected.to_ascii_lowercase();
        if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
            let file = File::open(path).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype zip {}: {err}",
                    path.display()
                ))
            })?;
            let mut archive = ZipArchive::new(file).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype zip {}: {err}",
                    path.display()
                ))
            })?;
            let entry = archive.by_name(&selected).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype entry {selected} in {}: {err}",
                    path.display()
                ))
            })?;
            let lines = read_lines_from_reader(BufReader::new(entry), path)?;
            return Self::from_vcf_lines(lines);
        }
        Ok(Self::from_delimited_file(
            path,
            GenotypeSourceFormat::Zip,
            Some(selected),
        ))
    }

    fn from_cram_file(path: &Path, options: &GenotypeLoadOptions) -> Result<Self, RuntimeError> {
        Ok(Self {
            backend: QueryBackend::Cram(CramBackend {
                path: path.to_path_buf(),
                options: options.clone(),
            }),
        })
    }

    fn from_vcf_lines(lines: Vec<String>) -> Result<Self, RuntimeError> {
        let mut values = HashMap::new();

        for line in lines {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with("##") || trimmed.starts_with("#CHROM") {
                continue;
            }

            let fields: Vec<&str> = trimmed.split('\t').collect();
            if fields.len() < 10 {
                continue;
            }

            let rsid = fields[2].trim();
            if rsid.is_empty() || rsid == "." {
                continue;
            }

            let reference = fields[3].trim();
            let alternates: Vec<&str> = fields[4]
                .split(',')
                .map(str::trim)
                .filter(|alt| !alt.is_empty() && *alt != ".")
                .collect();
            if reference.is_empty() || alternates.is_empty() {
                continue;
            }

            let sample_gt = fields[9].split(':').next().unwrap_or(".");
            if let Some(genotype) = genotype_from_vcf_gt(sample_gt, reference, &alternates) {
                values.insert(rsid.to_owned(), genotype);
            }
        }

        Ok(Self::from_rsid_map(GenotypeSourceFormat::Vcf, values))
    }

    fn from_delimited_lines(
        format: GenotypeSourceFormat,
        lines: Vec<String>,
    ) -> Result<Self, RuntimeError> {
        let delimiter = detect_delimiter(&lines);
        let mut parser = RowParser::new(delimiter);
        let mut values = HashMap::new();
        for line in lines {
            if let Some((rsid, genotype)) = parser.consume_line(&line)? {
                values.insert(rsid, genotype);
            }
        }
        Ok(Self::from_rsid_map(format, values))
    }

    fn from_rsid_map(format: GenotypeSourceFormat, values: HashMap<String, String>) -> Self {
        Self {
            backend: QueryBackend::RsidMap(RsidMapBackend { format, values }),
        }
    }

    fn from_delimited_file(
        path: &Path,
        format: GenotypeSourceFormat,
        zip_entry_name: Option<String>,
    ) -> Self {
        Self {
            backend: QueryBackend::Delimited(DelimitedBackend {
                format,
                path: path.to_path_buf(),
                zip_entry_name,
            }),
        }
    }

    pub fn capabilities(&self) -> BackendCapabilities {
        match &self.backend {
            QueryBackend::RsidMap(_) => BackendCapabilities {
                rsid_lookup: true,
                locus_lookup: false,
            },
            QueryBackend::Delimited(_) | QueryBackend::Vcf(_) => BackendCapabilities {
                rsid_lookup: true,
                locus_lookup: true,
            },
            QueryBackend::Cram(_) => BackendCapabilities {
                rsid_lookup: false,
                locus_lookup: true,
            },
        }
    }

    pub fn supports(&self, query: QueryKind) -> bool {
        let caps = self.capabilities();
        match query {
            QueryKind::GenotypeByRsid => caps.rsid_lookup,
            QueryKind::GenotypeByLocus => caps.locus_lookup,
        }
    }

    pub fn backend_name(&self) -> &'static str {
        match &self.backend {
            QueryBackend::RsidMap(map) => map.backend_name(),
            QueryBackend::Delimited(backend) => backend.backend_name(),
            QueryBackend::Vcf(backend) => backend.backend_name(),
            QueryBackend::Cram(backend) => backend.backend_name(),
        }
    }

    pub fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => Ok(map.values.get(rsid).cloned()),
            QueryBackend::Delimited(backend) => backend.get(rsid),
            QueryBackend::Vcf(backend) => backend.get(rsid),
            QueryBackend::Cram(backend) => backend
                .lookup_variant(&VariantSpec {
                    rsids: vec![rsid.to_owned()],
                    ..VariantSpec::default()
                })
                .map(|obs| obs.genotype),
        }
    }

    pub fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => map.lookup_variant(variant),
            QueryBackend::Delimited(backend) => backend.lookup_variant(variant),
            QueryBackend::Vcf(backend) => backend.lookup_variant(variant),
            QueryBackend::Cram(backend) => backend.lookup_variant(variant),
        }
    }

    pub fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        if let QueryBackend::Delimited(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::Vcf(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
        indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

        let mut results = vec![VariantObservation::default(); variants.len()];
        for (original_idx, variant) in indexed {
            results[original_idx] = self.lookup_variant(variant)?;
        }
        Ok(results)
    }
}

impl RsidMapBackend {
    fn backend_name(&self) -> &'static str {
        match self.format {
            GenotypeSourceFormat::Text => "text",
            GenotypeSourceFormat::Zip => "zip",
            GenotypeSourceFormat::Vcf => "vcf",
            GenotypeSourceFormat::Cram => "cram",
        }
    }

    fn lookup_variant(&self, variant: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        for rsid in &variant.rsids {
            if let Some(value) = self.values.get(rsid) {
                return Ok(VariantObservation {
                    backend: self.backend_name().to_owned(),
                    matched_rsid: Some(rsid.clone()),
                    genotype: Some(value.clone()),
                    evidence: vec![format!("resolved by rsid {rsid}")],
                    ..VariantObservation::default()
                });
            }
        }

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            evidence: vec!["no matching rsid found".to_owned()],
            ..VariantObservation::default()
        })
    }
}

impl DelimitedBackend {
    fn backend_name(&self) -> &'static str {
        match self.format {
            GenotypeSourceFormat::Text => "text",
            GenotypeSourceFormat::Zip => "zip",
            GenotypeSourceFormat::Vcf => "vcf",
            GenotypeSourceFormat::Cram => "cram",
        }
    }

    fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        let results = self.lookup_variants(&[VariantSpec {
            rsids: vec![rsid.to_owned()],
            ..VariantSpec::default()
        }])?;
        Ok(results.into_iter().next().and_then(|obs| obs.genotype))
    }

    fn lookup_variant(&self, variant: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        let mut results = self.lookup_variants(std::slice::from_ref(variant))?;
        Ok(results.pop().unwrap_or_default())
    }

    fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        scan_delimited_variants(self, variants)
    }
}

impl VcfBackend {
    fn backend_name(&self) -> &'static str {
        "vcf"
    }

    fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        let results = self.lookup_variants(&[VariantSpec {
            rsids: vec![rsid.to_owned()],
            ..VariantSpec::default()
        }])?;
        Ok(results.into_iter().next().and_then(|obs| obs.genotype))
    }

    fn lookup_variant(&self, variant: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        let mut results = self.lookup_variants(std::slice::from_ref(variant))?;
        Ok(results.pop().unwrap_or_default())
    }

    fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        if let Some(results) = lookup_indexed_vcf_variants(self, variants)? {
            return Ok(results);
        }
        scan_vcf_variants(self, variants)
    }
}

impl CramBackend {
    fn backend_name(&self) -> &'static str {
        "cram"
    }

    fn lookup_variant(&self, variant: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        let Some(reference_file) = self.options.reference_file.as_ref() else {
            return Err(RuntimeError::Unsupported(format!(
                "backend '{}' cannot satisfy query '{}' for {} without --reference-file",
                self.backend_name(),
                describe_query(variant),
                self.path.display()
            )));
        };

        let Some((assembly, locus)) = choose_variant_locus(variant, reference_file) else {
            let mut detail = format!(
                "backend '{}' cannot satisfy query '{}' for {} using reference {}",
                self.backend_name(),
                describe_query(variant),
                self.path.display(),
                reference_file.display()
            );
            detail.push_str(". This backend needs GRCh37/GRCh38 coordinates, not only rsIDs");
            if let Some(reference_index) = self.options.reference_index.as_ref() {
                let _ = write!(detail, " (reference index {})", reference_index.display());
            }
            if let Some(input_index) = self.options.input_index.as_ref() {
                let _ = write!(detail, " (input index {})", input_index.display());
            }
            return Err(RuntimeError::Unsupported(detail));
        };

        let observation = match variant.kind.unwrap_or(VariantKind::Other) {
            VariantKind::Snp => self.observe_snp(variant, assembly, &locus, reference_file)?,
            VariantKind::Deletion => {
                self.observe_deletion(variant, assembly, &locus, reference_file)?
            }
            VariantKind::Insertion | VariantKind::Indel => {
                self.observe_indel(variant, assembly, &locus, reference_file)?
            }
            VariantKind::Other => {
                return Err(RuntimeError::Unsupported(format!(
                    "backend '{}' does not yet support {:?} observation for {}",
                    self.backend_name(),
                    variant.kind.unwrap_or(VariantKind::Other),
                    self.path.display()
                )));
            }
        };

        Ok(observation)
    }

    fn observe_snp(
        &self,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
        reference_file: &Path,
    ) -> Result<VariantObservation, RuntimeError> {
        let reference = variant
            .reference
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| {
                RuntimeError::InvalidArguments("SNP variant requires ref/reference".to_owned())
            })?;
        let alternate = variant
            .alternate
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| {
                RuntimeError::InvalidArguments("SNP variant requires alt/alternate".to_owned())
            })?;

        let target_pos = locus.start;
        let pileup = observe_snp_pileup(
            &self.path,
            &self.options,
            reference_file,
            locus,
            reference,
            alternate,
        )?;
        let ref_count = pileup.filtered_ref_count;
        let alt_count = pileup.filtered_alt_count;
        let depth = pileup.filtered_depth;

        let evidence = pileup.evidence_lines(&describe_locus(locus), target_pos);

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_snp_genotype(reference, alternate, ref_count, alt_count, depth),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            raw_counts: pileup.raw_base_counts,
            decision: Some(describe_snp_decision_rule(
                reference, alternate, ref_count, alt_count, depth,
            )),
            evidence,
        })
    }

    fn observe_deletion(
        &self,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
        reference_file: &Path,
    ) -> Result<VariantObservation, RuntimeError> {
        let deletion_length = variant.deletion_length.ok_or_else(|| {
            RuntimeError::InvalidArguments("deletion variant requires deletion_length".to_owned())
        })?;
        let reference = variant.reference.clone().unwrap_or_else(|| "I".to_owned());
        let alternate = variant.alternate.clone().unwrap_or_else(|| "D".to_owned());
        let anchor_pos = locus.start.saturating_sub(1);

        let mut alt_count = 0u32;
        let mut ref_count = 0u32;
        let mut depth = 0u32;

        alignment::for_each_cram_record(
            &self.path,
            &self.options,
            reference_file,
            &anchor_window(locus),
            |record| {
                if record.is_unmapped || !spans_position(&record, anchor_pos) {
                    return Ok(true);
                }
                depth += 1;
                match indel_at_anchor(&record, anchor_pos) {
                    Some((AlignmentOpKind::Deletion, len)) if len == deletion_length => {
                        alt_count += 1;
                    }
                    _ => ref_count += 1,
                }
                Ok(true)
            },
        )?;

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_copy_number_genotype(
                &reference, &alternate, ref_count, alt_count, depth,
            ),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            raw_counts: BTreeMap::new(),
            decision: Some(describe_copy_number_decision_rule(
                &reference, &alternate, ref_count, alt_count, depth,
            )),
            evidence: vec![format!(
                "observed deletion anchor {}:{} len={} depth={} ref_count={} alt_count={}",
                locus.chrom, anchor_pos, deletion_length, depth, ref_count, alt_count
            )],
        })
    }

    fn observe_indel(
        &self,
        variant: &VariantSpec,
        assembly: Assembly,
        locus: &GenomicLocus,
        reference_file: &Path,
    ) -> Result<VariantObservation, RuntimeError> {
        let reference = variant.reference.clone().ok_or_else(|| {
            RuntimeError::InvalidArguments("indel variant requires ref/reference".to_owned())
        })?;
        let alternate = variant.alternate.clone().ok_or_else(|| {
            RuntimeError::InvalidArguments("indel variant requires alt/alternate".to_owned())
        })?;
        let records =
            alignment::query_cram_records(&self.path, &self.options, reference_file, locus)?;

        let mut alt_count = 0u32;
        let mut ref_count = 0u32;
        let mut depth = 0u32;
        let mut matching_alt_lengths = BTreeSet::new();

        for record in records {
            if record.is_unmapped {
                continue;
            }
            if !record_overlaps_locus(&record, locus) {
                continue;
            }
            let classification =
                classify_expected_indel(&record, locus, reference.len(), &alternate)?;
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
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_copy_number_genotype(
                &reference, &alternate, ref_count, alt_count, depth,
            ),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            raw_counts: BTreeMap::new(),
            decision: Some(describe_copy_number_decision_rule(
                &reference, &alternate, ref_count, alt_count, depth,
            )),
            evidence: vec![format!(
                "observed indel at {} depth={} ref_count={} alt_count={} matching_alt_lengths={}",
                describe_locus(locus),
                depth,
                ref_count,
                alt_count,
                evidence_label
            )],
        })
    }
}

fn choose_variant_locus(
    variant: &VariantSpec,
    reference_file: &Path,
) -> Option<(Assembly, GenomicLocus)> {
    match detect_reference_assembly(reference_file) {
        Some(Assembly::Grch38) => variant
            .grch38
            .clone()
            .map(|locus| (Assembly::Grch38, locus))
            .or_else(|| {
                variant
                    .grch37
                    .clone()
                    .map(|locus| (Assembly::Grch37, locus))
            }),
        Some(Assembly::Grch37) => variant
            .grch37
            .clone()
            .map(|locus| (Assembly::Grch37, locus))
            .or_else(|| {
                variant
                    .grch38
                    .clone()
                    .map(|locus| (Assembly::Grch38, locus))
            }),
        None => variant
            .grch38
            .clone()
            .map(|locus| (Assembly::Grch38, locus))
            .or_else(|| {
                variant
                    .grch37
                    .clone()
                    .map(|locus| (Assembly::Grch37, locus))
            }),
    }
}

fn detect_reference_assembly(reference_file: &Path) -> Option<Assembly> {
    let lower = reference_file.to_string_lossy().to_ascii_lowercase();
    if lower.contains("grch38") || lower.contains("hg38") || lower.contains("assembly38") {
        Some(Assembly::Grch38)
    } else if lower.contains("grch37") || lower.contains("hg19") || lower.contains("assembly37") {
        Some(Assembly::Grch37)
    } else {
        None
    }
}

fn describe_locus(locus: &GenomicLocus) -> String {
    format!("{}:{}-{}", locus.chrom, locus.start, locus.end)
}

fn anchor_window(locus: &GenomicLocus) -> GenomicLocus {
    let anchor = locus.start.saturating_sub(1);
    GenomicLocus {
        chrom: locus.chrom.clone(),
        start: anchor,
        end: anchor,
    }
}

fn first_base(value: &str) -> Option<char> {
    value
        .trim()
        .chars()
        .next()
        .map(|ch| ch.to_ascii_uppercase())
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
        Some(format!("{reference}{alternate}"))
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
        Some(format!("{reference}{alternate}"))
    }
}

fn describe_copy_number_decision_rule(
    reference: &str,
    alternate: &str,
    ref_count: u32,
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
        "copy-number genotype rule: alt_fraction={alt_fraction:.3} with thresholds ref<=0.200, het=(0.200,0.800), alt>=0.800; counts ref={ref_count} alt={alt_count} depth={depth} for {reference}->{alternate}"
    )
}

#[derive(Debug, Clone, Default)]
struct SnpPileupCounts {
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

impl SnpPileupCounts {
    fn evidence_lines(&self, locus: &str, target_pos: i64) -> Vec<String> {
        vec![
            format!(
                "observed SNP pileup at {locus} target_pos={target_pos} filtered_depth={} ref_count={} alt_count={}",
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
                "filters applied: min_base_quality={} min_mapping_quality={} filtered_low_base_quality={} filtered_low_mapping_quality={} filtered_non_acgt={} filtered_unmapped={} filtered_secondary={} filtered_qc_fail={} filtered_duplicate={} filtered_improper_pair={}",
                DEFAULT_MPILEUP_MIN_BASE_QUALITY,
                DEFAULT_MPILEUP_MIN_MAPPING_QUALITY,
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

fn observe_snp_pileup(
    cram_path: &Path,
    options: &GenotypeLoadOptions,
    reference_file: &Path,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
) -> Result<SnpPileupCounts, RuntimeError> {
    let repository = alignment::build_reference_repository(reference_file)?;
    let mut reader =
        alignment::build_cram_indexed_reader_from_path(cram_path, options, repository)?;
    let label = cram_path.display().to_string();
    snp_pileup_with_reader(
        &mut reader,
        &label,
        locus,
        reference,
        alternate,
        options.allow_reference_md5_mismatch,
    )
}

fn snp_pileup_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    allow_reference_md5_mismatch: bool,
) -> Result<SnpPileupCounts, RuntimeError> {
    let mut counts = SnpPileupCounts::default();
    let target_position = Position::try_from(usize::try_from(locus.start).map_err(|_| {
        RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned())
    })?)
    .map_err(|_| RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned()))?;
    let reference_base = reference as u8;

    alignment::for_each_raw_cram_record_with_reader_inner(
        reader,
        label,
        locus,
        allow_reference_md5_mismatch,
        |record| {
            let flags = record
                .flags()
                .map_err(|err| RuntimeError::Io(format!("failed to read CRAM flags: {err}")))?;
            if flags.is_unmapped() {
                counts.filtered_unmapped += 1;
                return Ok(true);
            }
            if flags.is_secondary() {
                counts.filtered_secondary += 1;
                return Ok(true);
            }
            if flags.is_qc_fail() {
                counts.filtered_qc_fail += 1;
                return Ok(true);
            }
            if flags.is_duplicate() {
                counts.filtered_duplicate += 1;
                return Ok(true);
            }
            if flags.is_segmented() && !flags.is_properly_segmented() {
                counts.filtered_improper_pair += 1;
                return Ok(true);
            }

            let Some((base, base_quality)) =
                cram_base_quality_at_reference_position(&record, target_position, reference_base)?
            else {
                return Ok(true);
            };

            let normalized_base = normalize_pileup_base(base);
            record.mapping_quality().transpose().map_err(|err| {
                RuntimeError::Io(format!("failed to read CRAM mapping quality: {err}"))
            })?;
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

            if base_quality < DEFAULT_MPILEUP_MIN_BASE_QUALITY {
                counts.filtered_low_base_quality += 1;
                return Ok(true);
            }

            let Some(base) = normalized_base else {
                counts.filtered_non_acgt += 1;
                return Ok(true);
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
            Ok(true)
        },
    )?;

    Ok(counts)
}

fn cram_base_quality_at_reference_position(
    record: &cram::Record<'_>,
    target_position: Position,
    reference_base: u8,
) -> Result<Option<(u8, u8)>, RuntimeError> {
    let Some(alignment_start) = record.alignment_start() else {
        return Ok(None);
    };
    let alignment_start = alignment_start
        .map_err(|err| RuntimeError::Io(format!("failed to read CRAM alignment start: {err}")))?;
    let mut reference_position = usize::from(alignment_start);
    let target = usize::from(target_position);
    let mut read_position = 0usize;
    let sequence = record.sequence();
    let qualities = record.quality_scores();

    for op in record.cigar().iter() {
        let op = op.map_err(|err| RuntimeError::Io(format!("failed to read CRAM CIGAR: {err}")))?;
        match op.kind() {
            CigarOpKind::Match | CigarOpKind::SequenceMatch | CigarOpKind::SequenceMismatch => {
                for offset in 0..op.len() {
                    if reference_position + offset == target {
                        let base = sequence
                            .get(read_position + offset)
                            .unwrap_or(reference_base);
                        let quality = qualities
                            .iter()
                            .nth(read_position + offset)
                            .transpose()
                            .map_err(|err| {
                                RuntimeError::Io(format!("failed to read CRAM base quality: {err}"))
                            })?
                            .unwrap_or(0);
                        return Ok(Some((base, quality)));
                    }
                }
                reference_position += op.len();
                read_position += op.len();
            }
            CigarOpKind::Insertion | CigarOpKind::SoftClip => {
                read_position += op.len();
            }
            CigarOpKind::Deletion | CigarOpKind::Skip => {
                if target >= reference_position && target < reference_position + op.len() {
                    return Ok(None);
                }
                reference_position += op.len();
            }
            CigarOpKind::HardClip | CigarOpKind::Pad => {}
        }
    }

    Ok(None)
}

/// Observe a SNP at `locus` over an already-built CRAM `IndexedReader` and
/// reference repository (held by the reader). Mirrors the internal
/// `CramBackend::observe_snp` but reader-based, so non-filesystem callers
/// (e.g. wasm with a JS-backed reader) don't need a `GenotypeStore` or paths.
///
/// `matched_rsid` and `assembly` are passed through to the returned
/// observation unchanged — callers that already know them (e.g. from
/// compiling a YAML variant) should supply them; otherwise `None`.
pub fn observe_cram_snp_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    let pileup = snp_pileup_with_reader(reader, label, locus, reference, alternate, false)?;
    let ref_count = pileup.filtered_ref_count;
    let alt_count = pileup.filtered_alt_count;
    let depth = pileup.filtered_depth;
    let evidence = pileup.evidence_lines(&describe_locus(locus), locus.start);

    Ok(VariantObservation {
        backend: "cram".to_owned(),
        matched_rsid,
        assembly,
        genotype: infer_snp_genotype(reference, alternate, ref_count, alt_count, depth),
        ref_count: Some(ref_count),
        alt_count: Some(alt_count),
        depth: Some(depth),
        raw_counts: pileup.raw_base_counts,
        decision: Some(describe_snp_decision_rule(
            reference, alternate, ref_count, alt_count, depth,
        )),
        evidence,
    })
}

/// Observe an insertion/indel-like variant at `locus` over an already-built
/// CRAM `IndexedReader`.
pub fn observe_cram_indel_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
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
    let mut matching_alt_lengths = BTreeSet::new();

    alignment::for_each_cram_record_with_reader(reader, label, locus, |record| {
        if record.is_unmapped || !record_overlaps_locus(&record, locus) {
            return Ok(true);
        }
        let classification = classify_expected_indel(&record, locus, reference.len(), alternate)?;
        if !classification.covering {
            return Ok(true);
        }
        depth += 1;
        if classification.matches_alt {
            alt_count += 1;
            matching_alt_lengths.insert(classification.observed_len);
        } else if classification.reference_like {
            ref_count += 1;
        }
        Ok(true)
    })?;

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
        backend: "cram".to_owned(),
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
            "observed indel at {} depth={} ref_count={} alt_count={} matching_alt_lengths={}",
            describe_locus(locus),
            depth,
            ref_count,
            alt_count,
            evidence_label
        )],
    })
}

/// Observe a SNP at `locus` over an already-built tabix-indexed bgzipped VCF
/// reader. Mirrors the CRAM variant for VCF: caller builds
/// `csi::io::IndexedReader::new(reader, tabix_index)` once and calls this per
/// variant. Non-filesystem callers (wasm with a JS-backed reader) go through
/// this path.
///
/// Only the locus-based match path is implemented — rsid-only variants would
/// need a linear scan, which is a follow-up. Pass `matched_rsid` through if
/// the caller already resolved it.
pub fn observe_vcf_snp_with_reader<R>(
    indexed: &mut csi::io::IndexedReader<bgzf::io::Reader<R>, tabix::Index>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
    matched_rsid: Option<String>,
    assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError>
where
    R: Read + Seek,
{
    let locus_label = format!("{}:{}", locus.chrom, locus.start);

    let Some(seq_name) = resolve_vcf_chrom_name(indexed.index(), &locus.chrom) else {
        return Ok(VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid,
            assembly,
            evidence: vec![format!(
                "{label}: tabix index has no contig matching {} (tried chr-prefixed and bare forms)",
                locus.chrom
            )],
            ..VariantObservation::default()
        });
    };

    let pos_usize = usize::try_from(locus.start).map_err(|err| {
        RuntimeError::Io(format!(
            "{label}: invalid VCF position {} for {locus_label}: {err}",
            locus.start
        ))
    })?;
    let position = Position::try_from(pos_usize).map_err(|err| {
        RuntimeError::Io(format!(
            "{label}: invalid VCF position {} for {locus_label}: {err}",
            locus.start
        ))
    })?;
    let region = Region::new(seq_name.as_str(), position..=position);

    let query = indexed.query(&region).map_err(|err| {
        RuntimeError::Io(format!("{label}: tabix query for {locus_label}: {err}"))
    })?;

    let reference_str = reference.to_ascii_uppercase().to_string();
    let alternate_str = alternate.to_ascii_uppercase().to_string();

    let mut saw_any = false;
    for record_result in query {
        let record = record_result
            .map_err(|err| RuntimeError::Io(format!("{label}: tabix record iter: {err}")))?;
        let line: &str = record.as_ref();
        let Some(row) = parse_vcf_record(line)? else {
            continue;
        };
        if row.position != locus.start {
            continue;
        }
        saw_any = true;
        if !row.reference.eq_ignore_ascii_case(&reference_str) {
            continue;
        }
        if !row
            .alternates
            .iter()
            .any(|alt| alt.eq_ignore_ascii_case(&alternate_str))
        {
            continue;
        }

        return Ok(VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid: matched_rsid.or_else(|| row.rsid.clone()),
            assembly,
            genotype: Some(row.genotype.clone()),
            evidence: vec![format!("{label}: resolved by locus {locus_label}")],
            ..VariantObservation::default()
        });
    }

    let evidence = if saw_any {
        vec![format!(
            "{label}: {locus_label} present but ref={reference}/alt={alternate} did not match any record"
        )]
    } else {
        vec![format!("{label}: no VCF record at {locus_label}")]
    };
    Ok(VariantObservation {
        backend: "vcf".to_owned(),
        matched_rsid,
        assembly,
        evidence,
        ..VariantObservation::default()
    })
}

/// Match the user-provided chromosome name against the tabix index's set of
/// reference sequence names. VCFs vary: some use `chr22`, others `22`. Try
/// the user's spelling verbatim, then toggle the `chr` prefix, then fall back
/// to a case-insensitive compare against the normalized suffix.
fn resolve_vcf_chrom_name(index: &tabix::Index, user_chrom: &str) -> Option<String> {
    let header = index.header()?;
    let names = header.reference_sequence_names();

    let trimmed = user_chrom.trim();
    let stripped = trimmed.strip_prefix("chr").unwrap_or(trimmed);

    let candidates = [
        trimmed.to_owned(),
        stripped.to_owned(),
        format!("chr{stripped}"),
    ];
    for cand in &candidates {
        if names.contains(cand.as_bytes()) {
            return Some(cand.clone());
        }
    }
    // Case-insensitive fallback against the full set.
    let target = stripped.to_ascii_lowercase();
    for name in names {
        let as_str = std::str::from_utf8(name.as_ref()).ok()?;
        let as_stripped = as_str.strip_prefix("chr").unwrap_or(as_str);
        if as_stripped.eq_ignore_ascii_case(&target) {
            return Some(as_str.to_owned());
        }
    }
    None
}

fn normalize_pileup_base(base: u8) -> Option<char> {
    match (base as char).to_ascii_uppercase() {
        'A' | 'C' | 'G' | 'T' => Some((base as char).to_ascii_uppercase()),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy)]
struct IndelClassification {
    covering: bool,
    reference_like: bool,
    matches_alt: bool,
    observed_len: usize,
}

fn len_as_i64(len: usize) -> Option<i64> {
    i64::try_from(len).ok()
}

fn spans_position(record: &AlignmentRecord, pos: i64) -> bool {
    pos >= record.start.saturating_sub(1) && pos <= record.end
}

fn record_overlaps_locus(record: &AlignmentRecord, locus: &GenomicLocus) -> bool {
    record.end >= locus.start && record.start <= locus.end
}

fn indel_at_anchor(record: &AlignmentRecord, anchor_pos: i64) -> Option<(AlignmentOpKind, usize)> {
    let mut ref_pos = record.start;

    for op in &record.cigar {
        match op.kind {
            AlignmentOpKind::Match
            | AlignmentOpKind::SequenceMatch
            | AlignmentOpKind::SequenceMismatch
            | AlignmentOpKind::Skip => {
                ref_pos += len_as_i64(op.len)?;
            }
            AlignmentOpKind::Insertion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((AlignmentOpKind::Insertion, op.len));
                }
            }
            AlignmentOpKind::Deletion => {
                let anchor = ref_pos.saturating_sub(1);
                if anchor == anchor_pos {
                    return Some((AlignmentOpKind::Deletion, op.len));
                }
                ref_pos += len_as_i64(op.len)?;
            }
            AlignmentOpKind::SoftClip | AlignmentOpKind::HardClip | AlignmentOpKind::Pad => {}
        }
    }

    None
}

fn classify_expected_indel(
    record: &AlignmentRecord,
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
                AlignmentOpKind::Insertion => reference_len + len,
                AlignmentOpKind::Deletion => reference_len.saturating_sub(len),
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

fn describe_query(variant: &VariantSpec) -> &'static str {
    if variant.has_coordinates() {
        "variant_by_locus"
    } else {
        "variant_by_rsid"
    }
}

fn variant_sort_key(variant: &VariantSpec) -> (u8, String, i64, i64, String) {
    if let Some(locus) = &variant.grch38 {
        return (
            0,
            chrom_sort_key(&locus.chrom),
            locus.start,
            locus.end,
            variant.rsids.first().cloned().unwrap_or_default(),
        );
    }
    if let Some(locus) = &variant.grch37 {
        return (
            1,
            chrom_sort_key(&locus.chrom),
            locus.start,
            locus.end,
            variant.rsids.first().cloned().unwrap_or_default(),
        );
    }
    (
        2,
        "~".to_owned(),
        i64::MAX,
        i64::MAX,
        variant.rsids.first().cloned().unwrap_or_default(),
    )
}

fn chrom_sort_key(raw: &str) -> String {
    let chrom = raw.trim().strip_prefix("chr").unwrap_or(raw.trim());
    if let Ok(value) = chrom.parse::<u32>() {
        return format!("{value:03}");
    }
    match chrom.to_ascii_uppercase().as_str() {
        "X" => "023".to_owned(),
        "Y" => "024".to_owned(),
        "M" | "MT" => "025".to_owned(),
        other => format!("999-{other}"),
    }
}

fn normalize_genotype(value: &str) -> String {
    let cleaned = value.trim().replace(' ', "").to_ascii_uppercase();
    if cleaned.is_empty() || matches!(cleaned.as_str(), "NA" | "N/A" | "#N/A" | "NONE") {
        return "--".to_owned();
    }
    if cleaned.contains('/') {
        let parts: Vec<&str> = cleaned.split('/').collect();
        if parts.iter().any(|part| part.is_empty() || *part == "-") {
            return "ID".to_owned();
        }
        return parts.concat();
    }
    cleaned
}

#[derive(Debug, Clone)]
struct ParsedVcfRow {
    rsid: Option<String>,
    chrom: String,
    position: i64,
    reference: String,
    alternates: Vec<String>,
    genotype: String,
}

fn scan_vcf_variants(
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

fn lookup_indexed_vcf_variants(
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

fn detect_vcf_assembly_from_path(path: &Path) -> Result<Option<Assembly>, RuntimeError> {
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

fn parse_vcf_record(line: &str) -> Result<Option<ParsedVcfRow>, RuntimeError> {
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

fn extract_vcf_sample_genotype(
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

fn detect_vcf_assembly(path: &Path, probe_lines: &[String]) -> Option<Assembly> {
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

fn choose_variant_locus_for_assembly(
    variant: &VariantSpec,
    assembly: Option<Assembly>,
) -> Option<GenomicLocus> {
    match assembly {
        Some(Assembly::Grch37) => variant.grch37.clone().or_else(|| variant.grch38.clone()),
        Some(Assembly::Grch38) => variant.grch38.clone().or_else(|| variant.grch37.clone()),
        None => variant.grch37.clone().or_else(|| variant.grch38.clone()),
    }
}

fn normalize_chromosome_name(value: &str) -> String {
    value.trim().trim_start_matches("chr").to_ascii_lowercase()
}

fn vcf_row_matches_variant(
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        io::Write,
        path::PathBuf,
        str::FromStr,
        time::{SystemTime, UNIX_EPOCH},
    };

    use zip::write::SimpleFileOptions;

    use crate::alignment::AlignmentOp;

    fn temp_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock drift")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "bioscript-genotype-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn locus(chrom: &str, start: i64, end: i64) -> GenomicLocus {
        GenomicLocus {
            chrom: chrom.to_owned(),
            start,
            end,
        }
    }

    fn mini_fixtures_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
    }

    fn variant_with_loci() -> VariantSpec {
        VariantSpec {
            rsids: vec!["rs1".to_owned()],
            grch37: Some(locus("1", 10, 10)),
            grch38: Some(locus("2", 20, 20)),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            deletion_length: None,
            motifs: Vec::new(),
        }
    }

    #[test]
    fn genotype_private_helpers_cover_assembly_sorting_and_decision_rules() {
        let variant = variant_with_loci();

        assert_eq!(
            choose_variant_locus(&variant, Path::new("ref/hg38.fa")),
            Some((Assembly::Grch38, locus("2", 20, 20)))
        );
        assert_eq!(
            choose_variant_locus(&variant, Path::new("ref/hg19.fa")),
            Some((Assembly::Grch37, locus("1", 10, 10)))
        );
        assert_eq!(
            choose_variant_locus(
                &VariantSpec {
                    grch38: Some(locus("3", 30, 30)),
                    ..VariantSpec::default()
                },
                Path::new("ref/hg19.fa")
            ),
            Some((Assembly::Grch38, locus("3", 30, 30)))
        );
        assert_eq!(
            choose_variant_locus(&variant, Path::new("ref/unknown.fa")),
            Some((Assembly::Grch38, locus("2", 20, 20)))
        );
        assert_eq!(
            choose_variant_locus_for_assembly(&variant, Some(Assembly::Grch37)),
            Some(locus("1", 10, 10))
        );
        assert_eq!(
            detect_reference_assembly(Path::new("assembly37.fa")),
            Some(Assembly::Grch37)
        );
        assert_eq!(
            detect_reference_assembly(Path::new("assembly38.fa")),
            Some(Assembly::Grch38)
        );
        assert_eq!(detect_reference_assembly(Path::new("other.fa")), None);

        assert_eq!(describe_locus(&locus("chr1", 7, 9)), "chr1:7-9");
        assert_eq!(anchor_window(&locus("1", 1, 4)), locus("1", 0, 0));
        assert_eq!(first_base(" tg"), Some('T'));
        assert_eq!(first_base(""), None);

        assert_eq!(infer_snp_genotype('A', 'G', 0, 0, 0), None);
        assert_eq!(
            infer_snp_genotype('A', 'G', 9, 1, 10).as_deref(),
            Some("AA")
        );
        assert_eq!(
            infer_snp_genotype('A', 'G', 1, 9, 10).as_deref(),
            Some("GG")
        );
        assert_eq!(
            infer_snp_genotype('A', 'G', 5, 5, 10).as_deref(),
            Some("AG")
        );
        assert!(describe_snp_decision_rule('A', 'G', 0, 0, 0).contains("no covering reads"));
        assert!(describe_snp_decision_rule('A', 'G', 0, 0, 3).contains("no reads matched"));
        assert!(describe_snp_decision_rule('A', 'G', 2, 8, 10).contains("alt_fraction=0.800"));

        assert_eq!(infer_copy_number_genotype("I", "D", 0, 0, 0), None);
        assert_eq!(
            infer_copy_number_genotype("I", "D", 9, 1, 10).as_deref(),
            Some("II")
        );
        assert_eq!(
            infer_copy_number_genotype("I", "D", 1, 9, 10).as_deref(),
            Some("DD")
        );
        assert_eq!(
            infer_copy_number_genotype("I", "D", 5, 5, 10).as_deref(),
            Some("ID")
        );
        assert!(
            describe_copy_number_decision_rule("I", "D", 0, 0, 0).contains("no covering reads")
        );

        assert_eq!(chrom_sort_key("chr2"), "002");
        assert_eq!(chrom_sort_key("X"), "023");
        assert_eq!(chrom_sort_key("MT"), "025");
        assert_eq!(chrom_sort_key("GL0001"), "999-GL0001");
        assert_eq!(variant_sort_key(&variant).0, 0);
        assert_eq!(describe_query(&variant), "variant_by_locus");
        assert_eq!(describe_query(&VariantSpec::default()), "variant_by_rsid");
    }

    #[test]
    fn genotype_private_helpers_cover_row_parsing_and_normalization() {
        assert!(matches!(
            detect_delimiter(&["# skip".to_owned(), "a,b".to_owned()]),
            Delimiter::Comma
        ));
        assert!(matches!(
            detect_delimiter(&["a b".to_owned()]),
            Delimiter::Space
        ));
        assert!(matches!(detect_delimiter(&Vec::new()), Delimiter::Tab));

        assert_eq!(strip_bom("\u{feff}rs1"), "rs1");
        assert_eq!(normalize_name("Base Pair-Position"), "basepairposition");
        assert_eq!(strip_inline_comment("AG # note"), "AG");
        assert_eq!(strip_inline_comment("AG // note"), "AG");
        assert_eq!(normalize_genotype("n/a"), "--");
        assert_eq!(normalize_genotype("a / g"), "AG");
        assert_eq!(normalize_genotype("A/-"), "ID");
        assert_eq!(split_csv_line(r#"rs1,"1,2",AG"#), vec!["rs1", "1,2", "AG"]);

        let mut parser = RowParser::new(Delimiter::Comma);
        assert!(
            parser
                .consume_record("# snpid,chr,pos,allele_a,allele_b")
                .unwrap()
                .is_none()
        );
        let row = parser.consume_record("rs1,1,10,A,G").unwrap().unwrap();
        assert_eq!(row.rsid.as_deref(), Some("rs1"));
        assert_eq!(row.chrom.as_deref(), Some("1"));
        assert_eq!(row.position, Some(10));
        assert_eq!(row.genotype, "AG");
        let short_row = parser.consume_record("bad,row").unwrap().unwrap();
        assert_eq!(short_row.rsid.as_deref(), Some("bad"));
        assert_eq!(short_row.genotype, "--");
        assert_eq!(
            parser.consume_line("rs4,4,40,tt").unwrap(),
            Some(("rs4".to_owned(), "TT".to_owned()))
        );
        assert!(parser.consume_record(",,").unwrap().is_none());
        assert_eq!(parser.default_header(2), vec!["rsid", "chromosome"]);
        assert_eq!(parser.default_header(6).len(), 6);

        let mut indexes = None;
        let mut comment_header = None;
        assert!(
            parse_streaming_row("", Delimiter::Space, &mut indexes, &mut comment_header)
                .unwrap()
                .is_none()
        );
        assert!(
            parse_streaming_row(
                "// marker chromosome position result",
                Delimiter::Space,
                &mut indexes,
                &mut comment_header
            )
            .unwrap()
            .is_none()
        );
        let row = parse_streaming_row(
            "rs2 chr2 20 ct",
            Delimiter::Space,
            &mut indexes,
            &mut comment_header,
        )
        .unwrap()
        .unwrap();
        assert_eq!(row.rsid.as_deref(), Some("rs2"));
        assert_eq!(row.genotype, "CT");

        let header = vec![
            "marker".to_owned(),
            "chrom".to_owned(),
            "base_pair_position".to_owned(),
            "allele1".to_owned(),
            "allele2".to_owned(),
        ];
        let cols = build_column_indexes(&header);
        assert_eq!(cols.rsid, Some(0));
        assert_eq!(cols.chrom, Some(1));
        assert_eq!(cols.position, Some(2));
        assert_eq!(cols.allele1, Some(3));
        assert_eq!(cols.allele2, Some(4));
        assert_eq!(default_column_indexes(2).position, None);
        assert_eq!(find_header_index(&header, GENOTYPE_ALIASES), None);
        assert!(looks_like_header_fields(&["rsid".to_owned()]));
        assert!(!looks_like_header_fields(&["sample".to_owned()]));
    }

    #[test]
    fn genotype_private_helpers_cover_vcf_parsing_and_matching() {
        assert!(parse_vcf_record("").unwrap().is_none());
        assert!(parse_vcf_record("#CHROM\tPOS").unwrap().is_none());
        assert!(parse_vcf_record("1\t10\trs1").unwrap().is_none());
        assert!(
            parse_vcf_record("1\t10\trs1\t.\tG\t.\tPASS\t.\tGT\t0/1")
                .unwrap()
                .is_none()
        );
        assert!(
            parse_vcf_record("1\t10\trs1\tA\t.\t.\tPASS\t.\tGT\t0/1")
                .unwrap()
                .is_none()
        );
        assert!(parse_vcf_record("1\tbad\trs1\tA\tG\t.\tPASS\t.\tGT\t0/1").is_err());

        let row = parse_vcf_record("chr1\t10\trs1\tA\tG,T\t.\tPASS\t.\tDP:GT\t8:1|2")
            .unwrap()
            .unwrap();
        assert_eq!(row.rsid.as_deref(), Some("rs1"));
        assert_eq!(row.genotype, "GT");
        assert_eq!(
            extract_vcf_sample_genotype("DP:AD", "8:1,2", "A", &["G".to_owned()]),
            None
        );
        assert_eq!(
            genotype_from_vcf_gt(".", "A", &["G"]).as_deref(),
            Some("--")
        );
        assert_eq!(
            genotype_from_vcf_gt("./1", "A", &["G"]).as_deref(),
            Some("--")
        );
        assert_eq!(
            genotype_from_vcf_gt("bad", "A", &["G"]).as_deref(),
            Some("--")
        );
        assert_eq!(genotype_from_vcf_gt("2/2", "A", &["G"]), None);
        assert_eq!(vcf_reference_token("AT", &["A"]), "I");
        assert_eq!(vcf_reference_token("A", &["AT"]), "D");
        assert_eq!(vcf_reference_token("A", &["<NON_REF>"]), "A");
        assert_eq!(vcf_alt_token("AT", "A"), "D");
        assert_eq!(vcf_alt_token("A", "AT"), "I");
        assert_eq!(vcf_alt_token("A", "<NON_REF>"), "--");
        assert!(is_symbolic_vcf_alt("<DEL>"));
        assert_eq!(normalize_sequence_token(" ag "), "AG");

        assert_eq!(
            detect_vcf_assembly(Path::new("sample.vcf"), &["##reference=hg19".to_owned()]),
            Some(Assembly::Grch37)
        );
        assert_eq!(
            detect_vcf_assembly(Path::new("sample.vcf"), &["##assembly=GRCh38".to_owned()]),
            Some(Assembly::Grch38)
        );
        assert_eq!(
            detect_vcf_assembly(Path::new("sample.b37.vcf"), &[]),
            Some(Assembly::Grch37)
        );
        assert_eq!(
            detect_vcf_assembly(Path::new("sample.b38.vcf"), &[]),
            Some(Assembly::Grch38)
        );
        assert_eq!(normalize_chromosome_name("chrX"), "x");

        let snp = VariantSpec {
            grch38: Some(locus("1", 10, 10)),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };
        assert!(vcf_row_matches_variant(&row, &snp, Some(Assembly::Grch38)));
        let deletion_row = parse_vcf_record("1\t9\trsdel\tATC\tA\t.\tPASS\t.\tGT\t0/1")
            .unwrap()
            .unwrap();
        let deletion = VariantSpec {
            grch38: Some(locus("1", 10, 12)),
            kind: Some(VariantKind::Deletion),
            deletion_length: Some(2),
            ..VariantSpec::default()
        };
        assert!(vcf_row_matches_variant(
            &deletion_row,
            &deletion,
            Some(Assembly::Grch38)
        ));
        let insertion = VariantSpec {
            grch38: Some(locus("1", 10, 10)),
            kind: Some(VariantKind::Insertion),
            ..VariantSpec::default()
        };
        assert!(vcf_row_matches_variant(
            &deletion_row,
            &insertion,
            Some(Assembly::Grch38)
        ));
        assert!(!vcf_row_matches_variant(
            &row,
            &VariantSpec::default(),
            None
        ));
        assert!(!vcf_row_matches_variant(
            &row,
            &VariantSpec {
                grch38: Some(locus("2", 10, 10)),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            },
            Some(Assembly::Grch38)
        ));
        assert!(vcf_row_matches_variant(
            &row,
            &VariantSpec {
                grch38: Some(locus("1", 10, 10)),
                kind: Some(VariantKind::Other),
                ..VariantSpec::default()
            },
            Some(Assembly::Grch38)
        ));
    }

    #[test]
    fn genotype_private_helpers_cover_indel_record_classification() {
        let record = AlignmentRecord {
            start: 10,
            end: 20,
            is_unmapped: false,
            cigar: vec![
                AlignmentOp {
                    kind: AlignmentOpKind::Match,
                    len: 2,
                },
                AlignmentOp {
                    kind: AlignmentOpKind::Insertion,
                    len: 3,
                },
                AlignmentOp {
                    kind: AlignmentOpKind::Deletion,
                    len: 2,
                },
                AlignmentOp {
                    kind: AlignmentOpKind::SoftClip,
                    len: 4,
                },
            ],
        };
        assert!(spans_position(&record, 9));
        assert!(record_overlaps_locus(&record, &locus("1", 15, 16)));
        assert_eq!(
            indel_at_anchor(&record, 11),
            Some((AlignmentOpKind::Insertion, 3))
        );
        let deletion_record = AlignmentRecord {
            start: 10,
            end: 20,
            is_unmapped: false,
            cigar: vec![
                AlignmentOp {
                    kind: AlignmentOpKind::Match,
                    len: 3,
                },
                AlignmentOp {
                    kind: AlignmentOpKind::Deletion,
                    len: 2,
                },
            ],
        };
        assert_eq!(
            indel_at_anchor(&deletion_record, 12),
            Some((AlignmentOpKind::Deletion, 2))
        );
        assert_eq!(indel_at_anchor(&record, 30), None);
        assert_eq!(len_as_i64(usize::MAX), None);

        let insertion = classify_expected_indel(&record, &locus("1", 12, 12), 1, "ATGC").unwrap();
        assert!(insertion.covering);
        assert!(insertion.matches_alt);
        assert_eq!(insertion.observed_len, 4);
        let deletion =
            classify_expected_indel(&deletion_record, &locus("1", 13, 13), 3, "A").unwrap();
        assert!(deletion.matches_alt);
        assert_eq!(deletion.observed_len, 1);
        let reference_like =
            classify_expected_indel(&record, &locus("1", 18, 18), 1, "AT").unwrap();
        assert!(reference_like.reference_like);
        let not_covering = classify_expected_indel(&record, &locus("1", 1, 2), 2, "A").unwrap();
        assert!(!not_covering.covering);

        assert_eq!(normalize_pileup_base(b'a'), Some('A'));
        assert_eq!(normalize_pileup_base(b'n'), None);
        let pileup = SnpPileupCounts {
            filtered_depth: 2,
            filtered_ref_count: 1,
            filtered_alt_count: 1,
            raw_depth: 3,
            raw_ref_count: 2,
            raw_alt_count: 1,
            filtered_low_base_quality: 1,
            filtered_non_acgt: 1,
            ..SnpPileupCounts::default()
        };
        let evidence = pileup.evidence_lines("1:10-10", 10);
        assert_eq!(evidence.len(), 4);
        assert!(evidence[0].contains("filtered_depth=2"));
    }

    #[test]
    fn genotype_private_helpers_cover_file_and_zip_scanning_paths() {
        let dir = temp_dir("file-zip-scanning");
        let text = dir.join("sample.txt");
        fs::write(
            &text,
            "# rsid chromosome position genotype\n\
             rs1 1 10 AG\n\
             rs2 2 20 CT\n",
        )
        .unwrap();
        assert!(matches!(
            detect_source_format(&text, None).unwrap(),
            GenotypeSourceFormat::Text
        ));
        assert!(matches!(
            detect_source_format(&text, Some(GenotypeSourceFormat::Cram)).unwrap(),
            GenotypeSourceFormat::Cram
        ));
        assert!(!looks_like_vcf_lines(&["rsid\tgenotype".to_owned()]));

        let backend = DelimitedBackend {
            format: GenotypeSourceFormat::Text,
            path: text.clone(),
            zip_entry_name: None,
        };
        let variants = vec![
            VariantSpec {
                rsids: vec!["rs2".to_owned()],
                ..VariantSpec::default()
            },
            VariantSpec {
                grch38: Some(locus("1", 10, 10)),
                ..VariantSpec::default()
            },
            VariantSpec {
                rsids: vec!["missing".to_owned()],
                ..VariantSpec::default()
            },
        ];
        let results = scan_delimited_variants(&backend, &variants).unwrap();
        assert_eq!(results[0].genotype.as_deref(), Some("CT"));
        assert_eq!(results[1].genotype.as_deref(), Some("AG"));
        assert!(results[2].evidence[0].contains("no matching rsid"));
        assert_eq!(backend.get("rs1").unwrap().as_deref(), Some("AG"));
        assert_eq!(
            backend
                .lookup_variant(&VariantSpec {
                    rsids: vec!["rs2".to_owned()],
                    ..VariantSpec::default()
                })
                .unwrap()
                .genotype
                .as_deref(),
            Some("CT")
        );

        let zip_path = dir.join("sample.zip");
        let cursor = std::io::Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .add_directory("nested/", SimpleFileOptions::default())
            .unwrap();
        writer
            .start_file("nested/sample.csv", SimpleFileOptions::default())
            .unwrap();
        writer
            .write_all(b"rsid,chromosome,position,genotype\nrs3,3,30,GG\n")
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        fs::write(&zip_path, bytes).unwrap();
        assert_eq!(select_zip_entry(&zip_path).unwrap(), "nested/sample.csv");
        let zip_backend = GenotypeStore::from_file(&zip_path).unwrap();
        assert_eq!(zip_backend.get("rs3").unwrap().as_deref(), Some("GG"));

        let unsupported_backend = DelimitedBackend {
            format: GenotypeSourceFormat::Vcf,
            path: text,
            zip_entry_name: None,
        };
        let err = scan_delimited_variants(&unsupported_backend, &variants).unwrap_err();
        assert!(
            err.to_string()
                .contains("streaming delimited backend only supports")
        );
    }

    #[test]
    fn genotype_public_entry_points_cover_in_memory_sources_and_fallbacks() {
        let text_store =
            GenotypeStore::from_bytes("sample.txt", b"rsid genotype\nrs1 AG\nrs2 CT\n").unwrap();
        assert_eq!(text_store.backend_name(), "text");
        assert!(text_store.supports(QueryKind::GenotypeByRsid));
        assert!(!text_store.supports(QueryKind::GenotypeByLocus));
        assert_eq!(text_store.get("rs1").unwrap().as_deref(), Some("AG"));
        let observations = text_store
            .lookup_variants(&[
                VariantSpec {
                    rsids: vec!["rs2".to_owned()],
                    ..VariantSpec::default()
                },
                VariantSpec {
                    rsids: vec!["missing".to_owned()],
                    ..VariantSpec::default()
                },
            ])
            .unwrap();
        assert_eq!(observations[0].genotype.as_deref(), Some("CT"));
        assert!(observations[1].evidence[0].contains("no matching rsid"));

        let vcf_store = GenotypeStore::from_bytes(
            "sample.vcf",
            b"##fileformat=VCFv4.3\n\
              #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
              1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
              1\t20\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
              1\t30\trs_bad\t.\tG\t.\tPASS\t.\tGT\t0/1\n",
        )
        .unwrap();
        assert_eq!(vcf_store.backend_name(), "vcf");
        assert_eq!(vcf_store.get("rs10").unwrap().as_deref(), Some("AG"));

        let cursor = std::io::Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .add_directory("nested/", SimpleFileOptions::default())
            .unwrap();
        writer
            .start_file("nested/sample.vcf", SimpleFileOptions::default())
            .unwrap();
        writer
            .write_all(
                b"##fileformat=VCFv4.3\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                  2\t20\trs20\tC\tT\t.\tPASS\t.\tGT\t1/1\n",
            )
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let zip_store = GenotypeStore::from_bytes("sample.zip", &bytes).unwrap();
        assert_eq!(zip_store.backend_name(), "vcf");
        assert_eq!(zip_store.get("rs20").unwrap().as_deref(), Some("TT"));

        let cursor = std::io::Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .add_directory("empty/", SimpleFileOptions::default())
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        let err = GenotypeStore::from_bytes("empty.zip", &bytes).unwrap_err();
        assert!(
            err.to_string()
                .contains("does not contain a supported genotype file"),
            "{err}"
        );

        assert!(
            GenotypeSourceFormat::from_str("unknown")
                .unwrap_err()
                .contains("unsupported input format")
        );
    }

    #[test]
    fn genotype_private_helpers_cover_vcf_file_zip_and_error_paths() {
        let dir = temp_dir("vcf-file-zip-errors");
        let vcf_path = dir.join("sample.grch38.vcf");
        fs::write(
            &vcf_path,
            "##fileformat=VCFv4.3\n\
             ##reference=GRCh38\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             chr1\t10\trs10\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:12\n\
             chr1\t19\trsDel\tAT\tA\t.\tPASS\t.\tGT\t1/1\n\
             chr2\t30\t.\tC\tT\t.\tPASS\t.\tDP:GT\t8:0/0\n",
        )
        .unwrap();

        let store = GenotypeStore::from_file(&vcf_path).unwrap();
        assert_eq!(store.get("rs10").unwrap().as_deref(), Some("AG"));
        let observations = store
            .lookup_variants(&[
                VariantSpec {
                    grch38: Some(locus("1", 10, 10)),
                    reference: Some("A".to_owned()),
                    alternate: Some("G".to_owned()),
                    kind: Some(VariantKind::Snp),
                    ..VariantSpec::default()
                },
                VariantSpec {
                    grch38: Some(locus("1", 20, 20)),
                    deletion_length: Some(1),
                    kind: Some(VariantKind::Deletion),
                    ..VariantSpec::default()
                },
                VariantSpec {
                    grch38: Some(locus("2", 31, 31)),
                    kind: Some(VariantKind::Other),
                    ..VariantSpec::default()
                },
                VariantSpec {
                    rsids: vec!["missing".to_owned()],
                    ..VariantSpec::default()
                },
            ])
            .unwrap();
        assert_eq!(observations[0].genotype.as_deref(), Some("AG"));
        assert_eq!(observations[1].genotype.as_deref(), Some("DD"));
        assert_eq!(observations[2].genotype.as_deref(), None);
        assert!(observations[3].evidence[0].contains("variant_by_rsid"));
        assert_eq!(observations[0].assembly, Some(Assembly::Grch38));

        let zip_path = dir.join("vcf.zip");
        let cursor = std::io::Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .start_file("nested/sample.vcf", SimpleFileOptions::default())
            .unwrap();
        writer
            .write_all(fs::read(&vcf_path).unwrap().as_slice())
            .unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        fs::write(&zip_path, bytes).unwrap();
        let zip_store = GenotypeStore::from_file(&zip_path).unwrap();
        assert_eq!(zip_store.get("rs10").unwrap().as_deref(), Some("AG"));

        let err = scan_vcf_variants(
            &VcfBackend {
                path: dir.join("missing.vcf"),
                options: GenotypeLoadOptions::default(),
            },
            &[VariantSpec::default()],
        )
        .unwrap_err();
        assert!(err.to_string().contains("failed to open VCF file"));

        let bad_zip_backend = DelimitedBackend {
            format: GenotypeSourceFormat::Zip,
            path: zip_path.clone(),
            zip_entry_name: None,
        };
        let err = scan_delimited_variants(&bad_zip_backend, &[VariantSpec::default()]).unwrap_err();
        assert!(
            err.to_string()
                .contains("zip backend missing selected entry"),
            "{err}"
        );

        let bad_entry_backend = DelimitedBackend {
            format: GenotypeSourceFormat::Zip,
            path: zip_path,
            zip_entry_name: Some("missing.csv".to_owned()),
        };
        let err =
            scan_delimited_variants(&bad_entry_backend, &[VariantSpec::default()]).unwrap_err();
        assert!(
            err.to_string().contains("failed to open genotype entry"),
            "{err}"
        );
    }

    #[test]
    fn genotype_private_helpers_cover_cram_backend_paths_with_mini_fixture() {
        let dir = mini_fixtures_dir();
        let cram = dir.join("mini.cram");
        let cram_index = dir.join("mini.cram.crai");
        let reference = dir.join("mini.fa");
        let options = GenotypeLoadOptions {
            input_index: Some(cram_index.clone()),
            reference_file: Some(reference.clone()),
            ..GenotypeLoadOptions::default()
        };
        let store = GenotypeStore::from_file_with_options(&cram, &options).unwrap();
        assert_eq!(store.backend_name(), "cram");
        assert!(store.supports(QueryKind::GenotypeByLocus));
        assert!(!store.supports(QueryKind::GenotypeByRsid));

        let snp = VariantSpec {
            rsids: vec!["mini_locus_1000".to_owned()],
            grch38: Some(locus("chr_test", 1000, 1000)),
            reference: Some("A".to_owned()),
            alternate: Some("C".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };
        let observation = store.lookup_variant(&snp).unwrap();
        assert_eq!(observation.backend, "cram");
        assert_eq!(observation.matched_rsid.as_deref(), Some("mini_locus_1000"));
        assert_eq!(observation.genotype.as_deref(), Some("AA"));
        assert_eq!(observation.depth, Some(50));

        let deletion = VariantSpec {
            rsids: vec!["mini_del".to_owned()],
            grch38: Some(locus("chr_test", 1000, 1000)),
            reference: Some("I".to_owned()),
            alternate: Some("D".to_owned()),
            kind: Some(VariantKind::Deletion),
            deletion_length: Some(1),
            ..VariantSpec::default()
        };
        let deletion_observation = store.lookup_variant(&deletion).unwrap();
        assert_eq!(deletion_observation.genotype.as_deref(), Some("II"));
        assert_eq!(deletion_observation.ref_count, Some(50));
        assert_eq!(deletion_observation.alt_count, Some(0));

        let indel = VariantSpec {
            rsids: vec!["mini_indel".to_owned()],
            grch38: Some(locus("chr_test", 1000, 1000)),
            reference: Some("A".to_owned()),
            alternate: Some("AT".to_owned()),
            kind: Some(VariantKind::Indel),
            ..VariantSpec::default()
        };
        let indel_observation = store.lookup_variant(&indel).unwrap();
        assert_eq!(indel_observation.genotype.as_deref(), Some("AA"));
        assert_eq!(indel_observation.ref_count, Some(50));
        assert_eq!(indel_observation.alt_count, Some(0));

        let missing_reference = GenotypeStore::from_file_with_options(
            &cram,
            &GenotypeLoadOptions {
                input_index: Some(cram_index.clone()),
                ..GenotypeLoadOptions::default()
            },
        )
        .unwrap();
        let err = missing_reference.lookup_variant(&snp).unwrap_err();
        assert!(err.to_string().contains("without --reference-file"));

        let err = store.get("rs-only").unwrap_err();
        assert!(err.to_string().contains("needs GRCh37/GRCh38 coordinates"));

        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(locus("chr_test", 1000, 1000)),
                kind: Some(VariantKind::Other),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("does not yet support"));

        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(locus("chr_test", 1000, 1000)),
                kind: Some(VariantKind::Snp),
                alternate: Some("C".to_owned()),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("SNP variant requires ref"));

        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(locus("chr_test", 1000, 1000)),
                kind: Some(VariantKind::Snp),
                reference: Some("A".to_owned()),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("SNP variant requires alt"));

        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(locus("chr_test", 1000, 1000)),
                kind: Some(VariantKind::Deletion),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("deletion_length"));

        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(locus("chr_test", 1000, 1000)),
                kind: Some(VariantKind::Indel),
                alternate: Some("AT".to_owned()),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("indel variant requires ref"));

        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(locus("chr_test", 1000, 1000)),
                kind: Some(VariantKind::Insertion),
                reference: Some("A".to_owned()),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("indel variant requires alt"));
    }

    #[test]
    fn genotype_public_cram_reader_snp_wrapper_uses_mini_fixture() {
        let dir = mini_fixtures_dir();
        let cram = dir.join("mini.cram");
        let cram_index = dir.join("mini.cram.crai");
        let reference = dir.join("mini.fa");
        let repository = crate::alignment::build_reference_repository(&reference).unwrap();
        let index = crate::alignment::parse_crai_bytes(&fs::read(cram_index).unwrap()).unwrap();
        let mut reader = crate::alignment::build_cram_indexed_reader_from_reader(
            fs::File::open(cram).unwrap(),
            index,
            repository,
        )
        .unwrap();

        let observation = observe_cram_snp_with_reader(
            &mut reader,
            "mini.cram",
            &locus("chr_test", 1000, 1000),
            'A',
            'C',
            Some("mini_locus_1000".to_owned()),
            Some(Assembly::Grch38),
        )
        .unwrap();

        assert_eq!(observation.genotype.as_deref(), Some("AA"));
        assert_eq!(observation.ref_count, Some(50));
        assert_eq!(observation.alt_count, Some(0));
        assert_eq!(observation.depth, Some(50));
        assert_eq!(observation.assembly, Some(Assembly::Grch38));
    }

    #[test]
    fn genotype_public_cram_reader_indel_wrapper_uses_mini_fixture() {
        let dir = mini_fixtures_dir();
        let cram = dir.join("mini.cram");
        let cram_index = dir.join("mini.cram.crai");
        let reference = dir.join("mini.fa");
        let repository = crate::alignment::build_reference_repository(&reference).unwrap();
        let index = crate::alignment::parse_crai_bytes(&fs::read(cram_index).unwrap()).unwrap();
        let mut reader = crate::alignment::build_cram_indexed_reader_from_reader(
            fs::File::open(cram).unwrap(),
            index,
            repository,
        )
        .unwrap();

        let observation = observe_cram_indel_with_reader(
            &mut reader,
            "mini.cram",
            &locus("chr_test", 1000, 1000),
            "A",
            "AT",
            Some("mini_indel".to_owned()),
            Some(Assembly::Grch38),
        )
        .unwrap();

        assert_eq!(observation.backend, "cram");
        assert_eq!(observation.matched_rsid.as_deref(), Some("mini_indel"));
        assert_eq!(observation.assembly, Some(Assembly::Grch38));
        assert_eq!(observation.genotype.as_deref(), Some("AA"));
        assert_eq!(observation.ref_count, Some(50));
        assert_eq!(observation.alt_count, Some(0));
        assert_eq!(observation.depth, Some(50));
        assert!(observation.evidence[0].contains("matching_alt_lengths=none"));
    }

    #[test]
    fn genotype_public_vcf_reader_wrapper_uses_tiny_tabix_fixture() {
        use noodles::vcf;

        let dir = temp_dir("vcf-reader-wrapper");
        let vcf_path = dir.join("sample.vcf.gz");
        let mut writer = bgzf::io::Writer::new(fs::File::create(&vcf_path).unwrap());
        writer
            .write_all(
                b"##fileformat=VCFv4.3\n\
                  #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
                  chr1\t10\trs10\tA\tG\t.\tPASS\t.\tGT\t0/1\n\
                  chr1\t20\trs20\tC\tT\t.\tPASS\t.\tGT\t1/1\n\
                  chr2\t30\trs30\tG\tA\t.\tPASS\t.\tGT\t0/0\n",
            )
            .unwrap();
        writer.finish().unwrap();

        let open_indexed = || {
            let index = vcf::fs::index(&vcf_path).unwrap();
            csi::io::IndexedReader::new(fs::File::open(&vcf_path).unwrap(), index)
        };

        let mut indexed = open_indexed();
        let observation = observe_vcf_snp_with_reader(
            &mut indexed,
            "tiny.vcf.gz",
            &locus("1", 10, 10),
            'A',
            'G',
            None,
            Some(Assembly::Grch38),
        )
        .unwrap();
        assert_eq!(observation.backend, "vcf");
        assert_eq!(observation.matched_rsid.as_deref(), Some("rs10"));
        assert_eq!(observation.genotype.as_deref(), Some("AG"));
        assert_eq!(observation.assembly, Some(Assembly::Grch38));

        let mut indexed = open_indexed();
        let observation = observe_vcf_snp_with_reader(
            &mut indexed,
            "tiny.vcf.gz",
            &locus("chr1", 10, 10),
            'A',
            'T',
            Some("requested".to_owned()),
            None,
        )
        .unwrap();
        assert_eq!(observation.matched_rsid.as_deref(), Some("requested"));
        assert!(observation.evidence[0].contains("did not match"));

        let mut indexed = open_indexed();
        let observation = observe_vcf_snp_with_reader(
            &mut indexed,
            "tiny.vcf.gz",
            &locus("1", 11, 11),
            'A',
            'G',
            None,
            None,
        )
        .unwrap();
        assert!(observation.evidence[0].contains("no VCF record"));

        let mut indexed = open_indexed();
        let observation = observe_vcf_snp_with_reader(
            &mut indexed,
            "tiny.vcf.gz",
            &locus("missing", 10, 10),
            'A',
            'G',
            Some("missing-rsid".to_owned()),
            Some(Assembly::Grch37),
        )
        .unwrap();
        assert_eq!(observation.matched_rsid.as_deref(), Some("missing-rsid"));
        assert_eq!(observation.assembly, Some(Assembly::Grch37));
        assert!(observation.evidence[0].contains("has no contig"));

        let mut indexed = open_indexed();
        let err = observe_vcf_snp_with_reader(
            &mut indexed,
            "tiny.vcf.gz",
            &locus("1", -1, -1),
            'A',
            'G',
            None,
            None,
        )
        .unwrap_err();
        assert!(err.to_string().contains("invalid VCF position"));
    }

    #[test]
    fn zip_entry_limited_reader_rejects_oversized_output() {
        let mut reader = std::io::Cursor::new(b"abcdef".to_vec());
        let err = read_zip_entry_limited(&mut reader, 5, "test zip entry").unwrap_err();
        assert!(
            err.to_string()
                .contains("test zip entry exceeds decompressed limit of 5 bytes"),
            "{err}"
        );
    }
}
