use std::{collections::HashMap, path::PathBuf, str::FromStr};

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
}

#[derive(Debug, Clone)]
pub(crate) struct RsidMapBackend {
    pub(crate) format: GenotypeSourceFormat,
    pub(crate) values: HashMap<String, String>,
}

#[derive(Debug, Clone)]
pub(crate) struct DelimitedBackend {
    pub(crate) format: GenotypeSourceFormat,
    pub(crate) path: PathBuf,
    pub(crate) zip_entry_name: Option<String>,
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
    pub allow_reference_md5_mismatch: bool,
}
