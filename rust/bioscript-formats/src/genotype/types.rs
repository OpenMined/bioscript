use std::{collections::HashMap, path::PathBuf, str::FromStr};

use bioscript_core::{Assembly, VariantObservation};

use crate::inspect::InferredSex;

#[derive(Debug, Clone)]
pub struct GenotypeStore {
    pub(crate) backend: QueryBackend,
}

#[derive(Debug, Clone)]
pub(crate) enum QueryBackend {
    RsidMap(RsidMapBackend),
    Delimited(DelimitedBackend),
    Vcf(VcfBackend),
    Cram(CramBackend),
    /// Pre-resolved observations layered on top of any other backend.
    /// Variant lookups consult `observations` first (matched by rsid OR by
    /// chrom+pos+ref+alt). On miss, falls back to `fallback`. This is the
    /// abstraction that lets the report pipeline collect every observation
    /// up-front in `run_manifest_rows` and have the analysis Python scripts'
    /// `genotypes.lookup_variants(plan)` calls resolve from the cache without
    /// re-walking the underlying genome — works identically on CLI (path
    /// fallback) and wasm (rsid-map empty fallback).
    Cached {
        observations: Vec<VariantObservation>,
        fallback: Box<QueryBackend>,
        require_hit: bool,
    },
}

#[derive(Debug, Clone)]
pub(crate) struct RsidMapBackend {
    pub(crate) format: GenotypeSourceFormat,
    pub(crate) values: HashMap<String, String>,
    pub(crate) locus_values: HashMap<(String, i64), (String, Option<String>, String)>,
    pub(crate) assembly: Option<Assembly>,
    /// Original input line per rsid, retained so wasm-side `from_bytes` loads
    /// can emit the same `| source line: …` evidence that the CLI's
    /// path-backed `DelimitedBackend` does on every lookup. Empty for
    /// in-memory maps that don't have a line representation.
    pub(crate) source_lines: HashMap<String, String>,
}

#[derive(Debug, Clone)]
pub(crate) struct DelimitedBackend {
    pub(crate) format: GenotypeSourceFormat,
    pub(crate) path: PathBuf,
    pub(crate) zip_entry_name: Option<String>,
    pub(crate) options: GenotypeLoadOptions,
}

#[derive(Debug, Clone)]
pub(crate) struct VcfBackend {
    pub(crate) path: PathBuf,
    pub(crate) options: GenotypeLoadOptions,
}

#[derive(Debug, Clone)]
pub(crate) struct CramBackend {
    pub(crate) path: PathBuf,
    pub(crate) options: GenotypeLoadOptions,
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
    Bam,
}

impl FromStr for GenotypeSourceFormat {
    type Err = String;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "txt" | "text" | "genotype" => Ok(Self::Text),
            "zip" => Ok(Self::Zip),
            "vcf" => Ok(Self::Vcf),
            "cram" => Ok(Self::Cram),
            "bam" => Ok(Self::Bam),
            other => Err(format!("unsupported input format: {other}")),
        }
    }
}

#[derive(Debug, Clone)]
pub struct GenotypeLoadOptions {
    pub format: Option<GenotypeSourceFormat>,
    pub input_index: Option<PathBuf>,
    pub reference_file: Option<PathBuf>,
    pub reference_index: Option<PathBuf>,
    pub assembly: Option<Assembly>,
    pub inferred_sex: Option<InferredSex>,
    pub impute_vcf_missing_as_reference: bool,
    pub allow_reference_md5_mismatch: bool,
}

impl Default for GenotypeLoadOptions {
    fn default() -> Self {
        Self {
            format: None,
            input_index: None,
            reference_file: None,
            reference_index: None,
            assembly: None,
            inferred_sex: None,
            impute_vcf_missing_as_reference: true,
            allow_reference_md5_mismatch: false,
        }
    }
}
