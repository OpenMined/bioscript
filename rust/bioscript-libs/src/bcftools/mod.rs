use std::path::Path;

use crate::{
    LibResult,
    tools::{CommandSpec, path_arg},
};

pub const MODULE: &str = "bcftools";

pub fn sort(input_vcf: &Path, output_vcf_gz: &Path) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "bcftools",
        vec![
            "sort".to_owned(),
            "-Oz".to_owned(),
            "-o".to_owned(),
            path_arg(output_vcf_gz)?,
            path_arg(input_vcf)?,
        ],
    )
}

pub fn index(vcf_gz: &Path) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "bcftools",
        vec!["index".to_owned(), "-t".to_owned(), path_arg(vcf_gz)?],
    )
}

pub fn view_filter(
    input_vcf: &Path,
    output_vcf_gz: &Path,
    include_expr: &str,
) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "bcftools",
        vec![
            "view".to_owned(),
            "-i".to_owned(),
            include_expr.to_owned(),
            "-Oz".to_owned(),
            "-o".to_owned(),
            path_arg(output_vcf_gz)?,
            path_arg(input_vcf)?,
        ],
    )
}

pub fn norm(
    input_vcf: &Path,
    reference_fasta: &Path,
    output_vcf_gz: &Path,
) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "bcftools",
        vec![
            "norm".to_owned(),
            "-f".to_owned(),
            path_arg(reference_fasta)?,
            "-Oz".to_owned(),
            "-o".to_owned(),
            path_arg(output_vcf_gz)?,
            path_arg(input_vcf)?,
        ],
    )
}
