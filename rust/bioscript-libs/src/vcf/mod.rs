use crate::{LibError, LibResult};

pub const MODULE: &str = "vcf";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VcfDirection {
    PysamVariantFile,
}

pub fn chosen_initial_surface() -> VcfDirection {
    VcfDirection::PysamVariantFile
}

pub fn open_variant_file() -> LibResult<()> {
    Err(LibError::unsupported_feature(
        MODULE,
        "VariantFile; planned as bioscript.pysam.VariantFile first",
    ))
}
