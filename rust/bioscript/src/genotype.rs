use std::{
    collections::{BTreeSet, HashMap},
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
    str::FromStr,
};

use noodles::bgzf;
use rust_htslib::bam::{self, Read};
use zip::ZipArchive;

use crate::runtime::RuntimeError;
use crate::variant::{Assembly, VariantKind, VariantObservation, VariantSpec};

const COMMENT_PREFIXES: [&str; 2] = ["#", "//"];

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
    Cram(CramBackend),
}

#[derive(Debug, Clone)]
struct RsidMapBackend {
    format: GenotypeSourceFormat,
    values: HashMap<String, String>,
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicLocus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
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

    pub fn from_file_with_options(path: &Path, options: &GenotypeLoadOptions) -> Result<Self, RuntimeError> {
        match detect_source_format(path, options.format)? {
            GenotypeSourceFormat::Text => Self::from_text_lines(read_plain_lines(path)?),
            GenotypeSourceFormat::Zip => Self::from_zip_file(path),
            GenotypeSourceFormat::Vcf => Self::from_vcf_file(path),
            GenotypeSourceFormat::Cram => Self::from_cram_file(path, options),
        }
    }

    fn from_text_lines(lines: Vec<String>) -> Result<Self, RuntimeError> {
        if lines.is_empty() {
            return Ok(Self::from_rsid_map(GenotypeSourceFormat::Text, HashMap::new()));
        }

        let delimiter = detect_delimiter(&lines);
        let mut parser = RowParser::new(delimiter);
        let mut values = HashMap::new();

        for line in lines {
            if let Some((rsid, genotype)) = parser.consume_line(&line)? {
                values.insert(rsid, genotype);
            }
        }

        Ok(Self::from_rsid_map(GenotypeSourceFormat::Text, values))
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
        let file = File::open(path)
            .map_err(|err| RuntimeError::Io(format!("failed to open genotype zip {}: {err}", path.display())))?;
        let mut archive = ZipArchive::new(file)
            .map_err(|err| RuntimeError::Io(format!("failed to read genotype zip {}: {err}", path.display())))?;

        let mut selected_idx = None;
        for idx in 0..archive.len() {
            let entry = archive.by_index(idx).map_err(|err| {
                RuntimeError::Io(format!("failed to inspect genotype zip {}: {err}", path.display()))
            })?;
            if entry.is_dir() {
                continue;
            }
            let name = entry.name().to_ascii_lowercase();
            if name.ends_with(".txt")
                || name.ends_with(".csv")
                || name.ends_with(".tsv")
                || name.ends_with(".vcf")
                || name.ends_with(".vcf.gz")
            {
                selected_idx = Some(idx);
                break;
            }
            if selected_idx.is_none() {
                selected_idx = Some(idx);
            }
        }

        let Some(idx) = selected_idx else {
            return Err(RuntimeError::Unsupported(format!(
                "zip archive {} does not contain a supported genotype file",
                path.display()
            )));
        };

        let entry_name = {
            let entry = archive.by_index(idx).map_err(|err| {
                RuntimeError::Io(format!("failed to inspect genotype zip {}: {err}", path.display()))
            })?;
            entry.name().to_owned()
        };

        let lines = {
            let entry = archive.by_index(idx).map_err(|err| {
                RuntimeError::Io(format!("failed to open genotype entry in {}: {err}", path.display()))
            })?;
            read_lines_from_reader(BufReader::new(entry), path)?
        };

        if looks_like_vcf_lines(&lines) || entry_name.to_ascii_lowercase().ends_with(".vcf") {
            Self::from_vcf_lines(lines)
        } else {
            Self::from_text_lines(lines)
        }
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

    pub fn capabilities(&self) -> BackendCapabilities {
        match &self.backend {
            QueryBackend::RsidMap(_) => BackendCapabilities {
                rsid_lookup: true,
                locus_lookup: false,
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
            QueryBackend::Cram(backend) => backend.backend_name(),
        }
    }

    pub fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => Ok(map.values.get(rsid).cloned()),
            QueryBackend::Cram(backend) => backend.lookup_variant(&VariantSpec {
                rsids: vec![rsid.to_owned()],
                ..VariantSpec::default()
            }).map(|obs| obs.genotype),
        }
    }

    pub fn lookup_variant(&self, variant: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => map.lookup_variant(variant),
            QueryBackend::Cram(backend) => backend.lookup_variant(variant),
        }
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
                detail.push_str(&format!(" (reference index {})", reference_index.display()));
            }
            if let Some(input_index) = self.options.input_index.as_ref() {
                detail.push_str(&format!(" (input index {})", input_index.display()));
            }
            return Err(RuntimeError::Unsupported(detail));
        };

        let observation = match variant.kind.unwrap_or(VariantKind::Other) {
            VariantKind::Snp => self.observe_snp(variant, assembly, &locus, reference_file)?,
            VariantKind::Deletion => self.observe_deletion(variant, assembly, &locus, reference_file)?,
            VariantKind::Insertion | VariantKind::Indel | VariantKind::Other => {
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
            .ok_or_else(|| RuntimeError::InvalidArguments("SNP variant requires ref/reference".to_owned()))?;
        let alternate = variant
            .alternate
            .as_deref()
            .and_then(first_base)
            .ok_or_else(|| RuntimeError::InvalidArguments("SNP variant requires alt/alternate".to_owned()))?;

        let target_pos = locus.start;
        let mut alt_count = 0u32;
        let mut ref_count = 0u32;
        let mut depth = 0u32;

        self.with_pileups(reference_file, locus, |pileup| {
            let pos1 = i64::from(pileup.pos()) + 1;
            if pos1 != target_pos {
                return;
            }

            for alignment in pileup.alignments() {
                if alignment.is_del() || alignment.is_refskip() {
                    continue;
                }
                let Some(qpos) = alignment.qpos() else {
                    continue;
                };
                let record = alignment.record();
                let bases = record.seq().as_bytes();
                let Some(base) = bases.get(qpos).copied() else {
                    continue;
                };
                let base = (base as char).to_ascii_uppercase();
                depth += 1;
                if base == reference {
                    ref_count += 1;
                } else if base == alternate {
                    alt_count += 1;
                }
            }
        })?;

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_snp_genotype(reference, alternate, ref_count, alt_count, depth),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            evidence: vec![format!(
                "observed SNP at {}:{} depth={} ref_count={} alt_count={}",
                locus.chrom, target_pos, depth, ref_count, alt_count
            )],
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

        self.with_pileups(reference_file, &anchor_window(locus), |pileup| {
            let pos1 = i64::from(pileup.pos()) + 1;
            if pos1 != anchor_pos {
                return;
            }

            for alignment in pileup.alignments() {
                if alignment.is_refskip() {
                    continue;
                }
                depth += 1;
                match alignment.indel() {
                    bam::pileup::Indel::Del(len) if usize::try_from(len).ok() == Some(deletion_length) => {
                        alt_count += 1;
                    }
                    _ => {
                        ref_count += 1;
                    }
                }
            }
        })?;

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            assembly: Some(assembly),
            genotype: infer_copy_number_genotype(&reference, &alternate, ref_count, alt_count, depth),
            ref_count: Some(ref_count),
            alt_count: Some(alt_count),
            depth: Some(depth),
            evidence: vec![format!(
                "observed deletion anchor {}:{} len={} depth={} ref_count={} alt_count={}",
                locus.chrom, anchor_pos, deletion_length, depth, ref_count, alt_count
            )],
        })
    }

    fn with_pileups<F>(
        &self,
        reference_file: &Path,
        locus: &GenomicLocus,
        mut on_pileup: F,
    ) -> Result<(), RuntimeError>
    where
        F: FnMut(&bam::pileup::Pileup),
    {
        if let Some(index_path) = self.options.input_index.as_ref() {
            let mut reader = bam::IndexedReader::from_path_and_index(&self.path, index_path).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open indexed CRAM {} with index {}: {err}",
                    self.path.display(),
                    index_path.display()
                ))
            })?;
            reader.set_reference(reference_file).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to set CRAM reference {} for {}: {err}",
                    reference_file.display(),
                    self.path.display()
                ))
            })?;
            fetch_locus(&mut reader, locus)?;
            for pileup in reader.pileup() {
                let pileup = pileup.map_err(|err| {
                    RuntimeError::Io(format!("failed while piling up {}: {err}", self.path.display()))
                })?;
                on_pileup(&pileup);
            }
            return Ok(());
        }

        if self.path.with_extension("cram.crai").exists() || self.path.with_extension("crai").exists() {
            let mut reader = bam::IndexedReader::from_path(&self.path).map_err(|err| {
                RuntimeError::Io(format!("failed to open indexed CRAM {}: {err}", self.path.display()))
            })?;
            reader.set_reference(reference_file).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to set CRAM reference {} for {}: {err}",
                    reference_file.display(),
                    self.path.display()
                ))
            })?;
            fetch_locus(&mut reader, locus)?;
            for pileup in reader.pileup() {
                let pileup = pileup.map_err(|err| {
                    RuntimeError::Io(format!("failed while piling up {}: {err}", self.path.display()))
                })?;
                on_pileup(&pileup);
            }
            return Ok(());
        }

        let mut reader = bam::Reader::from_path(&self.path)
            .map_err(|err| RuntimeError::Io(format!("failed to open CRAM {}: {err}", self.path.display())))?;
        reader.set_reference(reference_file).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to set CRAM reference {} for {}: {err}",
                reference_file.display(),
                self.path.display()
            ))
        })?;

        let target_tid = header_tid(reader.header(), &locus.chrom).ok_or_else(|| {
            RuntimeError::Unsupported(format!(
                "reference {} does not contain contig {} for {}",
                self.path.display(),
                locus.chrom,
                describe_locus(locus)
            ))
        })?;

        for pileup in reader.pileup() {
            let pileup = pileup.map_err(|err| {
                RuntimeError::Io(format!("failed while piling up {}: {err}", self.path.display()))
            })?;
            if pileup.tid() != target_tid {
                continue;
            }
            let pos1 = i64::from(pileup.pos()) + 1;
            if pos1 < locus.start {
                continue;
            }
            if pos1 > locus.end {
                break;
            }
            on_pileup(&pileup);
        }

        Ok(())
    }
}

fn choose_variant_locus(variant: &VariantSpec, reference_file: &Path) -> Option<(Assembly, GenomicLocus)> {
    match detect_reference_assembly(reference_file) {
        Some(Assembly::Grch38) => variant
            .grch38
            .clone()
            .map(|locus| (Assembly::Grch38, locus))
            .or_else(|| variant.grch37.clone().map(|locus| (Assembly::Grch37, locus))),
        Some(Assembly::Grch37) => variant
            .grch37
            .clone()
            .map(|locus| (Assembly::Grch37, locus))
            .or_else(|| variant.grch38.clone().map(|locus| (Assembly::Grch38, locus))),
        None => variant
            .grch38
            .clone()
            .map(|locus| (Assembly::Grch38, locus))
            .or_else(|| variant.grch37.clone().map(|locus| (Assembly::Grch37, locus))),
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

fn fetch_locus(reader: &mut bam::IndexedReader, locus: &GenomicLocus) -> Result<(), RuntimeError> {
    let tid = header_tid(reader.header(), &locus.chrom).ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "indexed CRAM does not contain contig {} for {}",
            locus.chrom,
            describe_locus(locus)
        ))
    })?;

    let start = locus.start.saturating_sub(1);
    let end = locus.end;
    reader.fetch((tid as i32, start, end)).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to fetch {}:{}-{}: {err}",
            locus.chrom, locus.start, locus.end
        ))
    })
}

fn header_tid(header: &bam::HeaderView, chrom: &str) -> Option<u32> {
    let candidates = [chrom.to_owned(), format!("chr{chrom}"), chrom.trim_start_matches("chr").to_owned()];
    candidates.iter().find_map(|candidate| header.tid(candidate.as_bytes()))
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
    value.trim().chars().next().map(|ch| ch.to_ascii_uppercase())
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

fn describe_query(variant: &VariantSpec) -> &'static str {
    if variant.has_coordinates() {
        "variant_by_locus"
    } else {
        "variant_by_rsid"
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
        if trimmed.is_empty() || COMMENT_PREFIXES.iter().any(|prefix| trimmed.starts_with(prefix)) {
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

struct RowParser {
    delimiter: Delimiter,
    header: Option<Vec<String>>,
    comment_header: Option<Vec<String>>,
    alias_map: HashMap<&'static str, BTreeSet<&'static str>>,
}

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
        let trimmed = line.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }

        let trimmed = strip_bom(trimmed);
        if let Some(prefix) = COMMENT_PREFIXES.iter().find(|prefix| trimmed.starts_with(**prefix)) {
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

        let Some(rsid) = self.lookup(&row_map, "rsid").filter(|value| !value.is_empty()) else {
            return Ok(None);
        };
        let Some(_chrom) = self.lookup(&row_map, "chromosome").filter(|value| !value.is_empty()) else {
            return Ok(None);
        };
        let Some(_pos) = self
            .lookup(&row_map, "position")
            .and_then(|value| value.parse::<i64>().ok())
        else {
            return Ok(None);
        };

        let genotype = if let Some(gt) = self.lookup(&row_map, "genotype") {
            gt
        } else {
            let allele1 = self.lookup(&row_map, "allele1").unwrap_or_default();
            let allele2 = self.lookup(&row_map, "allele2").unwrap_or_default();
            format!("{allele1}{allele2}")
        };

        Ok(Some((rsid, normalize_genotype(&genotype))))
    }

    fn parse_fields(&self, line: &str) -> Vec<String> {
        match self.delimiter {
            Delimiter::Tab => line.split('\t').map(|field| field.trim().to_owned()).collect(),
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
            if let Some(value) = row_map.get(&key) && !value.is_empty() {
                return Some(value.clone());
            }
        }
        None
    }

    fn default_header(&self, field_count: usize) -> Vec<String> {
        let base = ["rsid", "chromosome", "position", "genotype"];
        if field_count <= base.len() {
            base[..field_count].iter().map(|s| (*s).to_owned()).collect()
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
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open genotype file {}: {err}", path.display())))?;
    read_lines_from_reader(BufReader::new(file), path)
}

fn read_bgzf_lines(path: &Path) -> Result<Vec<String>, RuntimeError> {
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open genotype file {}: {err}", path.display())))?;
    let reader = bgzf::io::Reader::new(file);
    read_lines_from_reader(BufReader::new(reader), path)
}

fn read_lines_from_reader<R: BufRead>(mut reader: R, path: &Path) -> Result<Vec<String>, RuntimeError> {
    let mut lines = Vec::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        let bytes = reader
            .read_line(&mut buf)
            .map_err(|err| RuntimeError::Io(format!("failed to read genotype file {}: {err}", path.display())))?;
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
    if parts.len() != 2 || parts.iter().any(|part| *part == ".") {
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
    match alternate.len().cmp(&reference.len()) {
        std::cmp::Ordering::Less => "D".to_owned(),
        std::cmp::Ordering::Greater => "I".to_owned(),
        std::cmp::Ordering::Equal => normalize_sequence_token(alternate),
    }
}

fn normalize_sequence_token(value: &str) -> String {
    value.trim().to_ascii_uppercase()
}
