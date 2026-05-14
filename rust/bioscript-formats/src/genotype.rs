use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Cursor},
    path::Path,
};

use zip::ZipArchive;

#[cfg(test)]
use bioscript_core::{Assembly, GenomicLocus, VariantKind};
use bioscript_core::{RuntimeError, VariantObservation, VariantSpec};

mod backends;
mod common;
mod cram_backend;
mod delimited;
mod io;
mod loaders;
mod types;
mod vcf;
mod vcf_tokens;

#[cfg(test)]
use common::chrom_sort_key;
pub(crate) use common::{describe_query, normalize_genotype, variant_sort_key};
#[cfg(test)]
use cram_backend::{
    SnpPileupCounts, anchor_window, choose_variant_locus, classify_expected_indel,
    describe_copy_number_decision_rule, describe_locus, describe_snp_decision_rule,
    detect_reference_assembly, first_base, indel_at_anchor, infer_copy_number_genotype,
    infer_snp_genotype, len_as_i64, normalize_pileup_base, record_overlaps_locus, spans_position,
};
pub use cram_backend::{
    observe_cram_deletion_with_reader, observe_cram_indel_with_reader, observe_cram_snp_with_reader,
};
pub(crate) use delimited::{
    COMMENT_PREFIXES, DelimitedColumnIndexes, Delimiter, detect_delimiter, parse_streaming_row,
};
#[cfg(test)]
use delimited::{GENOTYPE_ALIASES, split_csv_line, strip_bom, strip_inline_comment};
use delimited::{RowParser, scan_delimited_variants};
#[cfg(test)]
use delimited::{
    build_column_indexes, default_column_indexes, find_header_index, looks_like_header_fields,
    normalize_name,
};
#[cfg(test)]
use io::looks_like_vcf_lines;
use io::{detect_source_format, is_bgzf_path, read_lines_from_reader, select_zip_entry};
pub use types::{
    BackendCapabilities, GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, QueryKind,
};
use types::{CramBackend, DelimitedBackend, QueryBackend, RsidMapBackend, VcfBackend};
pub use vcf::{
    choose_variant_locus_for_assembly, imputed_reference_observation, observe_vcf_snp_with_reader,
    observe_vcf_variant_with_reader,
};
#[cfg(test)]
use vcf::{
    detect_vcf_assembly, extract_vcf_sample_genotype, normalize_chromosome_name, parse_vcf_record,
    vcf_row_matches_variant,
};
use vcf::{lookup_indexed_vcf_variants, scan_vcf_variants};
pub(crate) use vcf_tokens::genotype_from_vcf_gt;
#[cfg(test)]
use vcf_tokens::{
    is_symbolic_vcf_alt, normalize_sequence_token, vcf_alt_token, vcf_reference_token,
};

impl GenotypeStore {
    pub fn from_file(path: &Path) -> Result<Self, RuntimeError> {
        Self::from_file_with_options(path, &GenotypeLoadOptions::default())
    }

    /// Wrap an existing store with a cache of pre-resolved observations.
    /// `lookup_variant`/`lookup_variants` consult `observations` first; on
    /// miss they delegate to the underlying store. This is how the report
    /// pipeline avoids re-walking the genome inside analysis Python scripts —
    /// the variants the panel/assays declared have already been observed in
    /// `run_manifest_rows`, so the cache hits and the fallback only fires
    /// for novel rsids the script asks about on its own.
    pub fn with_cached_observations(
        observations: Vec<VariantObservation>,
        fallback: GenotypeStore,
    ) -> Self {
        Self {
            backend: QueryBackend::Cached {
                observations,
                fallback: Box::new(fallback.backend),
                require_hit: false,
            },
        }
    }

    /// Wrap an existing store with a required cache of pre-resolved
    /// observations. Unlike `with_cached_observations`, lookup misses are
    /// errors. Report analysis uses this to avoid silently producing analysis
    /// output from a different lookup path than the observations table.
    pub fn with_required_cached_observations(
        observations: Vec<VariantObservation>,
        fallback: GenotypeStore,
    ) -> Self {
        Self {
            backend: QueryBackend::Cached {
                observations,
                fallback: Box::new(fallback.backend),
                require_hit: true,
            },
        }
    }

    /// Empty store whose `lookup_*` always returns "no observation". Useful as
    /// the wasm-side fallback when paired with `with_cached_observations` —
    /// every rsid the panel cares about is already in the cache so the
    /// fallback never fires; novel rsids return None gracefully.
    pub fn empty() -> Self {
        Self {
            backend: QueryBackend::RsidMap(RsidMapBackend {
                format: GenotypeSourceFormat::Text,
                values: HashMap::new(),
                locus_values: HashMap::new(),
                source_lines: HashMap::new(),
            }),
        }
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
            GenotypeSourceFormat::Bam => Err(RuntimeError::Unsupported(format!(
                "BAM alignment lookup is not implemented yet for {}",
                path.display()
            ))),
        }
    }

    pub fn from_bytes(name: &str, bytes: &[u8]) -> Result<Self, RuntimeError> {
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".zip") {
            return Self::from_zip_bytes(name, bytes);
        }
        let reader = BufReader::new(Cursor::new(bytes));
        if lower.ends_with(".vcf") {
            return Self::from_vcf_reader(reader, name);
        }
        Self::from_delimited_reader(GenotypeSourceFormat::Text, reader, name)
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
        let entry = archive.by_name(&selected).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open genotype entry {selected} in {name}: {err}"
            ))
        })?;
        let label = format!("genotype entry {selected} in {name}");
        // Stream-decompress directly off the zip reader so we never have to
        // materialize the entire decompressed entry in memory. GenesForGood
        // exports decompress to >128MB which used to trip the old
        // `read_zip_entry_limited` cap; the cap is gone because the streaming
        // parser keeps memory bounded to the rsid map itself.
        let reader = BufReader::new(entry);
        if selected.to_ascii_lowercase().ends_with(".vcf") {
            return Self::from_vcf_reader(reader, &label);
        }
        Self::from_delimited_reader(GenotypeSourceFormat::Zip, reader, &label)
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

    fn from_vcf_reader<R: BufRead>(reader: R, label: &str) -> Result<Self, RuntimeError> {
        loaders::from_vcf_reader(reader, label)
    }

    fn from_delimited_reader<R: BufRead>(
        format: GenotypeSourceFormat,
        reader: R,
        label: &str,
    ) -> Result<Self, RuntimeError> {
        loaders::from_delimited_reader(format, reader, label)
    }

    fn from_vcf_lines(lines: Vec<String>) -> Result<Self, RuntimeError> {
        loaders::from_vcf_lines(lines)
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
            QueryBackend::Cached { .. } => {
                // The cache itself answers both rsid and locus queries (we
                // match by either), so unioning with the fallback gives the
                // strongest set.
                BackendCapabilities {
                    rsid_lookup: true,
                    locus_lookup: true,
                }
            }
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
            QueryBackend::Cached { .. } => "cached",
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
            QueryBackend::Cached {
                observations,
                fallback,
                require_hit,
            } => {
                if let Some(matched) = observations
                    .iter()
                    .find(|obs| obs.matched_rsid.as_deref().is_some_and(|r| r == rsid))
                {
                    return Ok(matched.genotype.clone());
                }
                if *require_hit {
                    return Err(required_cache_miss(&VariantSpec {
                        rsids: vec![rsid.to_owned()],
                        ..VariantSpec::default()
                    }));
                }
                let inner = GenotypeStore {
                    backend: (**fallback).clone(),
                };
                inner.get(rsid)
            }
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
            QueryBackend::Cached {
                observations,
                fallback,
                require_hit,
            } => {
                if let Some(hit) = match_cached_observation(observations, variant) {
                    return Ok(hit.clone());
                }
                if *require_hit {
                    return Err(required_cache_miss(variant));
                }
                let inner = GenotypeStore {
                    backend: (**fallback).clone(),
                };
                inner.lookup_variant(variant)
            }
        }
    }

    pub fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        if let QueryBackend::Cached {
            observations,
            fallback,
            require_hit,
        } = &self.backend
        {
            // Resolve cache hits up-front; only round-trip the fallback for
            // misses so we don't pay for re-opening a CRAM/VCF when the panel
            // already covered every variant the script needs.
            let mut results: Vec<Option<VariantObservation>> = vec![None; variants.len()];
            let mut miss_indices = Vec::new();
            let mut miss_specs = Vec::new();
            for (idx, spec) in variants.iter().enumerate() {
                if let Some(hit) = match_cached_observation(observations, spec) {
                    results[idx] = Some(hit.clone());
                } else {
                    if *require_hit {
                        return Err(required_cache_miss(spec));
                    }
                    miss_indices.push(idx);
                    miss_specs.push(spec.clone());
                }
            }
            if !miss_specs.is_empty() {
                let inner = GenotypeStore {
                    backend: (**fallback).clone(),
                };
                let resolved = inner.lookup_variants(&miss_specs)?;
                for (idx, observation) in miss_indices.into_iter().zip(resolved) {
                    results[idx] = Some(observation);
                }
            }
            return Ok(results.into_iter().map(Option::unwrap_or_default).collect());
        }
        if let QueryBackend::Delimited(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::Vcf(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::Cram(backend) = &self.backend {
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

/// Match a `VariantSpec` against a pre-resolved observation list. Tries rsid
/// equality first (most common case for `PGx` panels), then falls back to a
/// chrom+pos+ref+alt match against either `GRCh37` or `GRCh38` loci so cached
/// observations from a CRAM lookup (which may have been done on one assembly)
/// can satisfy a script that supplies the spec on the other.
fn match_cached_observation<'a>(
    observations: &'a [VariantObservation],
    spec: &VariantSpec,
) -> Option<&'a VariantObservation> {
    if let Some(matched) = observations.iter().find(|obs| {
        obs.matched_rsid
            .as_deref()
            .is_some_and(|rsid| spec.rsids.iter().any(|target| target == rsid))
    }) {
        return Some(matched);
    }
    let assembly_loci = [spec.grch37.as_ref(), spec.grch38.as_ref()]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let target_ref = spec.reference.as_deref();
    let target_alt = spec.alternate.as_deref();
    observations.iter().find(|obs| {
        let evidence_match = assembly_loci.iter().any(|loci| {
            obs.evidence
                .iter()
                .any(|line| line.contains(&loci.chrom) && line.contains(&loci.start.to_string()))
        });
        if !evidence_match {
            return false;
        }
        match (target_ref, target_alt) {
            (Some(r), Some(a)) => obs
                .evidence
                .iter()
                .any(|line| line.contains(r) && line.contains(a)),
            _ => true,
        }
    })
}

fn required_cache_miss(spec: &VariantSpec) -> RuntimeError {
    let rsids = if spec.rsids.is_empty() {
        "<none>".to_owned()
    } else {
        spec.rsids.join("|")
    };
    let loci = [
        spec.grch37
            .as_ref()
            .map(|locus| format!("grch37:{}:{}-{}", locus.chrom, locus.start, locus.end)),
        spec.grch38
            .as_ref()
            .map(|locus| format!("grch38:{}:{}-{}", locus.chrom, locus.start, locus.end)),
    ]
    .into_iter()
    .flatten()
    .collect::<Vec<_>>()
    .join(",");
    RuntimeError::InvalidArguments(format!(
        "required preloaded genotype observation missing for rsids={rsids} loci={loci}"
    ))
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

    use noodles::bgzf;
    use noodles::csi;
    use zip::write::SimpleFileOptions;

    use crate::alignment::{AlignmentOp, AlignmentOpKind, AlignmentRecord};
    use crate::genotype::{io::read_plain_lines, vcf::detect_vcf_assembly_from_path};

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
            grch37_assembly_ref: None,
            grch38_assembly_ref: None,
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
        assert_eq!(
            infer_snp_genotype('G', 'A', 5, 5, 10).as_deref(),
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
    fn cached_observation_match_checks_all_variant_assemblies() {
        let spec = VariantSpec {
            rsids: vec!["rs60910145".to_owned()],
            grch37: Some(locus("22", 36_662_034, 36_662_034)),
            grch38: Some(locus("22", 36_265_988, 36_265_988)),
            reference: Some("T".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };
        let observations = vec![VariantObservation {
            backend: "cram".to_owned(),
            matched_rsid: None,
            assembly: Some(Assembly::Grch38),
            genotype: Some("TT".to_owned()),
            evidence: vec!["observed SNP pileup at 22:36265988-36265988 ref=T alt=G".to_owned()],
            ..VariantObservation::default()
        }];

        let matched = match_cached_observation(&observations, &spec)
            .expect("GRCh38 observation should match even when GRCh37 is listed first");

        assert_eq!(matched.genotype.as_deref(), Some("TT"));
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
        assert_eq!(normalize_genotype("g / a"), "AG");
        assert_eq!(normalize_genotype("T|C"), "CT");
        assert_eq!(normalize_genotype("TC"), "CT");
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
        assert_eq!(genotype_from_vcf_gt("1", "C", &["T"]).as_deref(), Some("T"));
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
            detect_vcf_assembly(
                Path::new("sample.vcf"),
                &["##contig=<ID=chr1,length=248956422>".to_owned()]
            ),
            Some(Assembly::Grch38)
        );
        assert_eq!(
            detect_vcf_assembly(
                Path::new("sample.vcf"),
                &["##contig=<ID=chr1,length=249250621>".to_owned()]
            ),
            Some(Assembly::Grch37)
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
        let assembly_ref_snp = VariantSpec {
            grch37: Some(locus("1", 10, 10)),
            grch37_assembly_ref: Some("G".to_owned()),
            reference: Some("A".to_owned()),
            alternate: Some("G".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };
        let inverted_row = parse_vcf_record("1\t10\trs10\tG\tA\t.\tPASS\t.\tGT\t1/1")
            .unwrap()
            .unwrap();
        assert!(vcf_row_matches_variant(
            &inverted_row,
            &assembly_ref_snp,
            Some(Assembly::Grch37)
        ));
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
    fn vcf_scan_translates_assembly_ref_snv_rows_and_missing_reference() {
        let dir = temp_dir("vcf-assembly-ref");
        let vcf_path = dir.join("sample.grch37.vcf");
        fs::write(
            &vcf_path,
            "##fileformat=VCFv4.3\n\
             ##reference=GRCh37\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n\
             1\t10\trs10\tT\tC,A\t.\tPASS\t.\tGT\t1/1\n",
        )
        .unwrap();

        let backend = VcfBackend {
            path: vcf_path,
            options: GenotypeLoadOptions {
                assembly: Some(Assembly::Grch37),
                ..GenotypeLoadOptions::default()
            },
        };
        let present = VariantSpec {
            grch37: Some(locus("1", 10, 10)),
            grch37_assembly_ref: Some("T".to_owned()),
            reference: Some("C".to_owned()),
            alternate: Some("T".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };
        let missing = VariantSpec {
            grch37: Some(locus("1", 20, 20)),
            grch37_assembly_ref: Some("T".to_owned()),
            reference: Some("C".to_owned()),
            alternate: Some("T".to_owned()),
            kind: Some(VariantKind::Snp),
            ..VariantSpec::default()
        };

        let observations = scan_vcf_variants(&backend, &[present, missing]).unwrap();
        assert_eq!(observations[0].genotype.as_deref(), Some("CC"));
        assert!(
            observations[0].evidence[0].contains("resolved by locus"),
            "{:?}",
            observations[0].evidence
        );
        assert_eq!(observations[1].genotype.as_deref(), Some("TT"));
        assert!(
            observations[1].evidence[0].contains("imputed reference genotype"),
            "{:?}",
            observations[1].evidence
        );
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
        let vcf_like = dir.join("sample.dat");
        fs::write(&vcf_like, "##fileformat=VCFv4.3\n#CHROM\tPOS\n").unwrap();
        assert!(matches!(
            detect_source_format(&vcf_like, None).unwrap(),
            GenotypeSourceFormat::Vcf
        ));
        assert!(is_bgzf_path(Path::new("sample.VCF.GZ")));
        assert!(!is_bgzf_path(Path::new("sample.vcf")));
        let missing = read_plain_lines(&dir.join("missing.txt")).unwrap_err();
        assert!(missing.to_string().contains("failed to open genotype file"));
        assert!(matches!(
            detect_source_format(&text, Some(GenotypeSourceFormat::Cram)).unwrap(),
            GenotypeSourceFormat::Cram
        ));
        let bam = dir.join("sample.bam");
        fs::write(&bam, b"BAM").unwrap();
        assert!(matches!(
            detect_source_format(&bam, None).unwrap(),
            GenotypeSourceFormat::Bam
        ));
        assert!(
            GenotypeStore::from_file(&bam)
                .unwrap_err()
                .to_string()
                .contains("BAM alignment lookup is not implemented yet")
        );
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

        let fallback_zip = dir.join("fallback.zip");
        let cursor = std::io::Cursor::new(Vec::new());
        let mut writer = zip::ZipWriter::new(cursor);
        writer
            .add_directory("nested/", SimpleFileOptions::default())
            .unwrap();
        writer
            .start_file("nested/notes.bin", SimpleFileOptions::default())
            .unwrap();
        writer.write_all(b"notes\n").unwrap();
        let bytes = writer.finish().unwrap().into_inner();
        fs::write(&fallback_zip, bytes).unwrap();
        assert_eq!(select_zip_entry(&fallback_zip).unwrap(), "nested/notes.bin");

        let empty_zip = dir.join("empty.zip");
        let cursor = std::io::Cursor::new(Vec::new());
        let writer = zip::ZipWriter::new(cursor);
        fs::write(&empty_zip, writer.finish().unwrap().into_inner()).unwrap();
        let err = select_zip_entry(&empty_zip).unwrap_err();
        assert!(
            err.to_string()
                .contains("does not contain a supported genotype file")
        );

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

        assert_eq!(
            detect_vcf_assembly_from_path(&vcf_path).unwrap(),
            Some(Assembly::Grch38)
        );

        let no_index = VcfBackend {
            path: vcf_path.clone(),
            options: GenotypeLoadOptions::default(),
        };
        assert!(
            lookup_indexed_vcf_variants(&no_index, &[VariantSpec::default()])
                .unwrap()
                .is_none()
        );

        let indexed = VcfBackend {
            path: vcf_path.clone(),
            options: GenotypeLoadOptions {
                input_index: Some(dir.join("missing.tbi")),
                ..GenotypeLoadOptions::default()
            },
        };
        assert!(
            lookup_indexed_vcf_variants(&indexed, &[VariantSpec::default()])
                .unwrap()
                .is_none()
        );
        let err = lookup_indexed_vcf_variants(
            &indexed,
            &[VariantSpec {
                grch38: Some(locus("1", 10, 10)),
                reference: Some("A".to_owned()),
                alternate: Some("G".to_owned()),
                kind: Some(VariantKind::Snp),
                ..VariantSpec::default()
            }],
        )
        .unwrap_err();
        assert!(err.to_string().contains("failed to read VCF index"));

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

        let missing_text_backend = DelimitedBackend {
            format: GenotypeSourceFormat::Text,
            path: dir.join("missing.txt"),
            zip_entry_name: None,
        };
        let err =
            scan_delimited_variants(&missing_text_backend, &[VariantSpec::default()]).unwrap_err();
        assert!(err.to_string().contains("failed to open genotype file"));
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
}
