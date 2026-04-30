#![allow(
    clippy::case_sensitive_file_extension_comparisons,
    clippy::missing_errors_doc,
    clippy::must_use_candidate,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::unnecessary_wraps,
    clippy::unused_self
)]

pub mod alignment;
mod genotype;
mod inspect;
mod prepare;

pub use genotype::{
    BackendCapabilities, GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, QueryKind,
    observe_cram_indel_with_reader, observe_cram_snp_with_reader, observe_vcf_snp_with_reader,
};
pub use inspect::{
    DetectedKind, DetectionConfidence, FileContainer, FileInspection, InspectOptions,
    SourceMetadata, inspect_bytes, inspect_file,
};
pub use prepare::{PrepareRequest, PreparedPaths, prepare_indexes, shell_flags};
