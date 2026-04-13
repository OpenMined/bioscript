#![allow(
    clippy::case_sensitive_file_extension_comparisons,
    clippy::missing_errors_doc,
    clippy::must_use_candidate,
    clippy::too_many_lines,
    clippy::unnecessary_wraps,
    clippy::unused_self
)]

mod genotype;
mod prepare;

pub use genotype::{
    BackendCapabilities, GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, QueryKind,
};
pub use prepare::{PrepareRequest, PreparedPaths, prepare_indexes, shell_flags};
