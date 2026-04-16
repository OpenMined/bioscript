use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    fmt::Write as _,
    fs::File,
    io::{BufRead, BufReader, Read, Seek},
    path::{Path, PathBuf},
    str::FromStr,
};

use noodles::bgzf;
use noodles::core::Position;
use noodles::cram;
use noodles::sam::alignment::Record as _;
use zip::ZipArchive;

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};

use crate::alignment::{self, AlignmentOpKind, AlignmentRecord};

const COMMENT_PREFIXES: [&str; 2] = ["#", "//"];
const DEFAULT_MPILEUP_MIN_BASE_QUALITY: u8 = 13;
const DEFAULT_MPILEUP_MIN_MAPPING_QUALITY: u8 = 0;

const RSID_ALIASES: &[&str] = &["rsid", "name", "snp", "marker", "id", "snpid"];
const CHROM_ALIASES: &[&str] = &["chromosome", "chr", "chrom"];
const POSITION_ALIASES: &[&str] = &[
    "position",
    "pos",
    "coordinate",
    "basepairposition",
    "basepair",
];
const GENOTYPE_ALIASES: &[&str] = &[
    "genotype",
    "gt",
    "result",
    "results",
    "result1",
    "call",
    "calls",
    "yourcode",
    "code",
    "genotypevalue",
    "variation",
];
const ALLELE1_ALIASES: &[&str] = &["allele1", "allelea", "allele_a", "allele1top"];
const ALLELE2_ALIASES: &[&str] = &["allele2", "alleleb", "allele_b", "allele2top"];

#[derive(Debug, Clone)]
pub struct GenotypeStore {
    backend: QueryBackend,
}

#[derive(Debug, Clone)]
enum QueryBackend {
    RsidMap(RsidMapBackend),
    Delimited(DelimitedBackend),
    Cram(CramBackend),
}

#[derive(Debug, Clone)]
struct RsidMapBackend {
    format: GenotypeSourceFormat,
    values: HashMap<String, String>,
}

#[derive(Debug, Clone)]
struct DelimitedBackend {
    format: GenotypeSourceFormat,
    path: PathBuf,
    zip_entry_name: Option<String>,
}

#[derive(Debug, Clone)]
struct CramBackend {
    path: PathBuf,
    options: GenotypeLoadOptions,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryKind {
    GenotypeByRsid,
    GenotypeByLocus,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BackendCapabilities {
    pub rsid_lookup: bool,
    pub locus_lookup: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GenotypeSourceFormat {
    Text,
    Zip,
    Vcf,
    Cram,
}

impl FromStr for GenotypeSourceFormat {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "txt" | "text" | "genotype" => Ok(Self::Text),
            "zip" => Ok(Self::Zip),
            "vcf" => Ok(Self::Vcf),
            "cram" => Ok(Self::Cram),
            other => Err(format!("unsupported input format: {other}")),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct GenotypeLoadOptions {
    pub format: Option<GenotypeSourceFormat>,
    pub input_index: Option<PathBuf>,
    pub reference_file: Option<PathBuf>,
    pub reference_index: Option<PathBuf>,
}

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
            GenotypeSourceFormat::Vcf => Self::from_vcf_file(path),
            GenotypeSourceFormat::Cram => Self::from_cram_file(path, options),
        }
    }

    fn from_vcf_file(path: &Path) -> Result<Self, RuntimeError> {
        let lines = if path
            .extension()
            .and_then(|ext| ext.to_str())
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            read_bgzf_lines(path)?
        } else {
            read_plain_lines(path)?
        };

        Self::from_vcf_lines(lines)
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
            QueryBackend::Delimited(_) => BackendCapabilities {
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
            QueryBackend::Cram(backend) => backend.backend_name(),
        }
    }

    pub fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => Ok(map.values.get(rsid).cloned()),
            QueryBackend::Delimited(backend) => backend.get(rsid),
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

#[derive(Debug, Clone)]
struct ParsedDelimitedRow {
    rsid: Option<String>,
    chrom: Option<String>,
    position: Option<i64>,
    genotype: String,
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
    snp_pileup_with_reader(&mut reader, &label, locus, reference, alternate)
}

fn snp_pileup_with_reader<R: Read + Seek>(
    reader: &mut cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    locus: &GenomicLocus,
    reference: char,
    alternate: char,
) -> Result<SnpPileupCounts, RuntimeError> {
    let mut counts = SnpPileupCounts::default();
    let target_position = Position::try_from(usize::try_from(locus.start).map_err(|_| {
        RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned())
    })?)
    .map_err(|_| RuntimeError::InvalidArguments("SNP locus start is out of range".to_owned()))?;
    let reference_base = reference as u8;

    alignment::for_each_raw_cram_record_with_reader(reader, label, locus, |record| {
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
            record.base_quality_at_reference_position(target_position, reference_base)
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
    })?;

    Ok(counts)
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
    let pileup = snp_pileup_with_reader(reader, label, locus, reference, alternate)?;
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

#[derive(Debug, Clone, Copy)]
enum Delimiter {
    Tab,
    Comma,
    Space,
}

fn detect_delimiter(lines: &[String]) -> Delimiter {
    for line in lines {
        let trimmed = line.trim();
        if trimmed.is_empty()
            || COMMENT_PREFIXES
                .iter()
                .any(|prefix| trimmed.starts_with(prefix))
        {
            continue;
        }
        if line.contains('\t') {
            return Delimiter::Tab;
        }
        if line.contains(',') {
            return Delimiter::Comma;
        }
        if trimmed.split_whitespace().count() > 1 {
            return Delimiter::Space;
        }
    }
    Delimiter::Tab
}

#[allow(dead_code)]
struct RowParser {
    delimiter: Delimiter,
    header: Option<Vec<String>>,
    comment_header: Option<Vec<String>>,
    alias_map: HashMap<&'static str, BTreeSet<&'static str>>,
}

#[allow(dead_code)]
impl RowParser {
    fn new(delimiter: Delimiter) -> Self {
        let mut alias_map = HashMap::new();
        alias_map.insert("rsid", RSID_ALIASES.iter().copied().collect());
        alias_map.insert("chromosome", CHROM_ALIASES.iter().copied().collect());
        alias_map.insert("position", POSITION_ALIASES.iter().copied().collect());
        alias_map.insert("genotype", GENOTYPE_ALIASES.iter().copied().collect());
        alias_map.insert("allele1", ALLELE1_ALIASES.iter().copied().collect());
        alias_map.insert("allele2", ALLELE2_ALIASES.iter().copied().collect());
        Self {
            delimiter,
            header: None,
            comment_header: None,
            alias_map,
        }
    }

    fn consume_line(&mut self, line: &str) -> Result<Option<(String, String)>, RuntimeError> {
        Ok(self
            .consume_record(line)?
            .and_then(|row| row.rsid.map(|rsid| (rsid, row.genotype))))
    }

    fn consume_record(&mut self, line: &str) -> Result<Option<ParsedDelimitedRow>, RuntimeError> {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }

        let trimmed = strip_bom(trimmed);
        if let Some(prefix) = COMMENT_PREFIXES
            .iter()
            .find(|prefix| trimmed.starts_with(**prefix))
        {
            let candidate = trimmed.trim_start_matches(prefix).trim();
            if !candidate.is_empty() {
                let fields = self.parse_fields(candidate);
                if self.looks_like_header(&fields) {
                    self.comment_header = Some(fields);
                }
            }
            return Ok(None);
        }

        let fields = self.parse_fields(strip_bom(line));
        if fields.is_empty() {
            return Ok(None);
        }

        if self.header.is_none() {
            if self.looks_like_header(&fields) {
                self.header = Some(fields);
                return Ok(None);
            }
            if let Some(header) = self.comment_header.take() {
                self.header = Some(header);
            } else {
                self.header = Some(self.default_header(fields.len()));
            }
        }

        let header = self.header.as_ref().expect("header initialized");
        let mut row_map = HashMap::new();
        for (idx, value) in fields.into_iter().enumerate() {
            if idx >= header.len() {
                continue;
            }
            row_map.insert(normalize_name(&header[idx]), strip_inline_comment(&value));
        }

        let rsid = self
            .lookup(&row_map, "rsid")
            .filter(|value| !value.is_empty());
        let chrom = self
            .lookup(&row_map, "chromosome")
            .filter(|value| !value.is_empty());
        let position = self
            .lookup(&row_map, "position")
            .and_then(|value| value.parse::<i64>().ok());
        if rsid.is_none() && (chrom.is_none() || position.is_none()) {
            return Ok(None);
        }

        let genotype = if let Some(gt) = self.lookup(&row_map, "genotype") {
            gt
        } else {
            let allele1 = self.lookup(&row_map, "allele1").unwrap_or_default();
            let allele2 = self.lookup(&row_map, "allele2").unwrap_or_default();
            format!("{allele1}{allele2}")
        };

        Ok(Some(ParsedDelimitedRow {
            rsid,
            chrom,
            position,
            genotype: normalize_genotype(&genotype),
        }))
    }

    fn parse_fields(&self, line: &str) -> Vec<String> {
        match self.delimiter {
            Delimiter::Tab => line
                .split('\t')
                .map(|field| field.trim().to_owned())
                .collect(),
            Delimiter::Space => line.split_whitespace().map(str::to_owned).collect(),
            Delimiter::Comma => split_csv_line(line),
        }
    }

    fn looks_like_header(&self, fields: &[String]) -> bool {
        fields.first().is_some_and(|first| {
            self.alias_map
                .get("rsid")
                .is_some_and(|aliases| aliases.contains(normalize_name(first).as_str()))
        })
    }

    fn lookup(&self, row_map: &HashMap<String, String>, key: &str) -> Option<String> {
        let aliases = self.alias_map.get(key)?;
        for alias in aliases {
            let key = normalize_name(alias);
            if let Some(value) = row_map.get(&key)
                && !value.is_empty()
            {
                return Some(value.clone());
            }
        }
        None
    }

    fn default_header(&self, field_count: usize) -> Vec<String> {
        let base = ["rsid", "chromosome", "position", "genotype"];
        if field_count <= base.len() {
            base[..field_count]
                .iter()
                .map(|s| (*s).to_owned())
                .collect()
        } else {
            let mut header: Vec<String> = base.iter().map(|s| (*s).to_owned()).collect();
            for idx in 0..(field_count - header.len()) {
                header.push(format!("extra_{idx}"));
            }
            header
        }
    }
}
fn strip_bom(value: &str) -> &str {
    value.strip_prefix('\u{feff}').unwrap_or(value)
}

fn normalize_name(name: &str) -> String {
    name.trim()
        .to_ascii_lowercase()
        .chars()
        .filter(|ch| !matches!(ch, ' ' | '_' | '-'))
        .collect()
}

fn strip_inline_comment(value: &str) -> String {
    for marker in ["#", "//"] {
        if let Some(idx) = value.find(marker) {
            return value[..idx].trim().to_owned();
        }
    }
    value.trim().to_owned()
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

fn split_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut current = String::new();
    let mut in_quotes = false;
    let chars = line.chars().peekable();

    for ch in chars {
        match ch {
            '"' => in_quotes = !in_quotes,
            ',' if !in_quotes => {
                fields.push(current.trim().to_owned());
                current.clear();
            }
            _ => current.push(ch),
        }
    }
    fields.push(current.trim().to_owned());
    fields
}

fn read_plain_lines(path: &Path) -> Result<Vec<String>, RuntimeError> {
    let file = File::open(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open genotype file {}: {err}",
            path.display()
        ))
    })?;
    read_lines_from_reader(BufReader::new(file), path)
}

fn select_zip_entry(path: &Path) -> Result<String, RuntimeError> {
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

    let mut selected_name: Option<String> = None;
    for idx in 0..archive.len() {
        let entry = archive.by_index(idx).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to inspect genotype zip {}: {err}",
                path.display()
            ))
        })?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".txt")
            || lower.ends_with(".csv")
            || lower.ends_with(".tsv")
            || lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
        {
            return Ok(name);
        }
        if selected_name.is_none() {
            selected_name = Some(name);
        }
    }

    selected_name.ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "zip archive {} does not contain a supported genotype file",
            path.display()
        ))
    })
}

fn scan_delimited_variants(
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

#[derive(Debug, Clone, Copy)]
struct DelimitedColumnIndexes {
    rsid: Option<usize>,
    chrom: Option<usize>,
    position: Option<usize>,
    genotype: Option<usize>,
    allele1: Option<usize>,
    allele2: Option<usize>,
}

fn parse_streaming_row(
    line: &str,
    delimiter: Delimiter,
    column_indexes: &mut Option<DelimitedColumnIndexes>,
    comment_header: &mut Option<Vec<String>>,
) -> Result<Option<ParsedDelimitedRow>, RuntimeError> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }

    let trimmed = strip_bom(trimmed);
    if let Some(prefix) = COMMENT_PREFIXES
        .iter()
        .find(|prefix| trimmed.starts_with(**prefix))
    {
        let candidate = trimmed.trim_start_matches(prefix).trim();
        if !candidate.is_empty() {
            let fields = parse_owned_fields(candidate, delimiter);
            if looks_like_header_fields(&fields) {
                *comment_header = Some(fields);
            }
        }
        return Ok(None);
    }

    let fields = parse_owned_fields(strip_bom(line), delimiter);
    if fields.is_empty() {
        return Ok(None);
    }

    if column_indexes.is_none() {
        if looks_like_header_fields(&fields) {
            *column_indexes = Some(build_column_indexes(&fields));
            return Ok(None);
        }
        if let Some(header) = comment_header.take() {
            *column_indexes = Some(build_column_indexes(&header));
        } else {
            *column_indexes = Some(default_column_indexes(fields.len()));
        }
    }

    let indexes = column_indexes.expect("streaming column indexes initialized");
    let rsid = indexes
        .rsid
        .and_then(|idx| fields.get(idx))
        .map(|value| strip_inline_comment(value).trim().to_owned())
        .filter(|value| !value.is_empty());
    let chrom = indexes
        .chrom
        .and_then(|idx| fields.get(idx))
        .map(|value| strip_inline_comment(value).trim().to_owned())
        .filter(|value| !value.is_empty());
    let position = indexes
        .position
        .and_then(|idx| fields.get(idx))
        .and_then(|value| strip_inline_comment(value).trim().parse::<i64>().ok());
    if rsid.is_none() && (chrom.is_none() || position.is_none()) {
        return Ok(None);
    }

    let genotype = if let Some(idx) = indexes.genotype {
        fields
            .get(idx)
            .map(|value| strip_inline_comment(value))
            .unwrap_or_default()
            .clone()
    } else {
        let allele1 = indexes
            .allele1
            .and_then(|idx| fields.get(idx))
            .map(|value| strip_inline_comment(value))
            .unwrap_or_default();
        let allele2 = indexes
            .allele2
            .and_then(|idx| fields.get(idx))
            .map(|value| strip_inline_comment(value))
            .unwrap_or_default();
        format!("{allele1}{allele2}")
    };

    Ok(Some(ParsedDelimitedRow {
        rsid,
        chrom,
        position,
        genotype: normalize_genotype(&genotype),
    }))
}

fn parse_owned_fields(line: &str, delimiter: Delimiter) -> Vec<String> {
    match delimiter {
        Delimiter::Tab => line
            .split('\t')
            .map(|field| field.trim().to_owned())
            .collect(),
        Delimiter::Space => line.split_whitespace().map(str::to_owned).collect(),
        Delimiter::Comma => split_csv_line(line),
    }
}

fn looks_like_header_fields(fields: &[String]) -> bool {
    fields
        .first()
        .is_some_and(|first| RSID_ALIASES.contains(&normalize_name(first).as_str()))
}

fn build_column_indexes(header: &[String]) -> DelimitedColumnIndexes {
    DelimitedColumnIndexes {
        rsid: find_header_index(header, RSID_ALIASES),
        chrom: find_header_index(header, CHROM_ALIASES),
        position: find_header_index(header, POSITION_ALIASES),
        genotype: find_header_index(header, GENOTYPE_ALIASES),
        allele1: find_header_index(header, ALLELE1_ALIASES),
        allele2: find_header_index(header, ALLELE2_ALIASES),
    }
}

fn default_column_indexes(field_count: usize) -> DelimitedColumnIndexes {
    DelimitedColumnIndexes {
        rsid: (field_count > 0).then_some(0),
        chrom: (field_count > 1).then_some(1),
        position: (field_count > 2).then_some(2),
        genotype: (field_count > 3).then_some(3),
        allele1: None,
        allele2: None,
    }
}

fn find_header_index(header: &[String], aliases: &[&str]) -> Option<usize> {
    header.iter().position(|field| {
        aliases
            .iter()
            .any(|alias| normalize_name(field) == normalize_name(alias))
    })
}

fn read_bgzf_lines(path: &Path) -> Result<Vec<String>, RuntimeError> {
    let file = File::open(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open genotype file {}: {err}",
            path.display()
        ))
    })?;
    let reader = bgzf::io::Reader::new(file);
    read_lines_from_reader(BufReader::new(reader), path)
}

fn read_lines_from_reader<R: BufRead>(
    mut reader: R,
    path: &Path,
) -> Result<Vec<String>, RuntimeError> {
    let mut lines = Vec::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        let bytes = reader.read_line(&mut buf).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read genotype file {}: {err}",
                path.display()
            ))
        })?;
        if bytes == 0 {
            break;
        }
        lines.push(buf.trim_end_matches(['\n', '\r']).to_owned());
    }
    Ok(lines)
}

fn detect_source_format(
    path: &Path,
    forced: Option<GenotypeSourceFormat>,
) -> Result<GenotypeSourceFormat, RuntimeError> {
    if let Some(format) = forced {
        return Ok(format);
    }

    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".zip") {
        return Ok(GenotypeSourceFormat::Zip);
    }
    if lower.ends_with(".cram") {
        return Ok(GenotypeSourceFormat::Cram);
    }
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
        return Ok(GenotypeSourceFormat::Vcf);
    }

    let lines = read_plain_lines(path)?;
    if looks_like_vcf_lines(&lines) {
        Ok(GenotypeSourceFormat::Vcf)
    } else {
        Ok(GenotypeSourceFormat::Text)
    }
}

fn looks_like_vcf_lines(lines: &[String]) -> bool {
    lines.iter().any(|line| {
        let trimmed = line.trim_start();
        trimmed.starts_with("##fileformat=VCF") || trimmed.starts_with("#CHROM\t")
    })
}

fn genotype_from_vcf_gt(gt: &str, reference: &str, alternates: &[&str]) -> Option<String> {
    if matches!(gt.trim(), "" | "." | "./." | ".|.") {
        return Some("--".to_owned());
    }

    let cleaned = gt.trim().replace('|', "/");
    let parts: Vec<&str> = cleaned.split('/').collect();
    if parts.len() != 2 || parts.contains(&".") {
        return Some("--".to_owned());
    }

    let ref_token = vcf_reference_token(reference, alternates);
    let mut out = String::new();
    for part in parts {
        let idx = part.parse::<usize>().ok()?;
        if idx == 0 {
            out.push_str(&ref_token);
        } else {
            let alt = alternates.get(idx - 1)?;
            out.push_str(&vcf_alt_token(reference, alt));
        }
    }

    Some(normalize_genotype(&out))
}

fn vcf_reference_token(reference: &str, alternates: &[&str]) -> String {
    let mut saw_shorter = false;
    let mut saw_longer = false;

    for alt in alternates {
        if is_symbolic_vcf_alt(alt) {
            continue;
        }
        match alt.len().cmp(&reference.len()) {
            std::cmp::Ordering::Less => saw_shorter = true,
            std::cmp::Ordering::Greater => saw_longer = true,
            std::cmp::Ordering::Equal => {}
        }
    }

    match (saw_shorter, saw_longer) {
        (true, false) => "I".to_owned(),
        (false, true) => "D".to_owned(),
        _ => normalize_sequence_token(reference),
    }
}

fn vcf_alt_token(reference: &str, alternate: &str) -> String {
    if is_symbolic_vcf_alt(alternate) {
        return "--".to_owned();
    }
    match alternate.len().cmp(&reference.len()) {
        std::cmp::Ordering::Less => "D".to_owned(),
        std::cmp::Ordering::Greater => "I".to_owned(),
        std::cmp::Ordering::Equal => normalize_sequence_token(alternate),
    }
}

fn is_symbolic_vcf_alt(alternate: &str) -> bool {
    let trimmed = alternate.trim();
    trimmed.starts_with('<') && trimmed.ends_with('>')
}

fn normalize_sequence_token(value: &str) -> String {
    value.trim().to_ascii_uppercase()
}
