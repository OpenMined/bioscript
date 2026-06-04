use std::path::Path;
use std::{ffi::OsString, process::ExitCode};

use crate::{
    LibError, LibResult,
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

pub fn view(input_vcf: &Path, output_vcf: &Path, output_type: &str) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "bcftools",
        vec![
            "view".to_owned(),
            "-O".to_owned(),
            output_type.to_owned(),
            "-o".to_owned(),
            path_arg(output_vcf)?,
            path_arg(input_vcf)?,
        ],
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

pub fn view_native(input_vcf: &Path, output_vcf: &Path, output_type: &str) -> LibResult<()> {
    let argv = [
        OsString::from("view"),
        OsString::from("--no-version"),
        OsString::from("-O"),
        OsString::from(output_type),
        OsString::from("-o"),
        output_vcf.as_os_str().to_owned(),
        input_vcf.as_os_str().to_owned(),
    ];
    run_bcftools("view", bcftools_rs::commands::view::main(&argv))
}

pub fn sort_native(
    input_vcf: &Path,
    output_vcf: &Path,
    output_type: &str,
    write_index: bool,
) -> LibResult<()> {
    let mut argv = vec![
        OsString::from("sort"),
        input_vcf.as_os_str().to_owned(),
        OsString::from("-o"),
        output_vcf.as_os_str().to_owned(),
        OsString::from("-O"),
        OsString::from(output_type),
    ];
    if write_index {
        argv.push(OsString::from("-W"));
    }
    run_bcftools("sort", bcftools_rs::commands::sort::main(&argv))
}

pub fn index_native(
    input_vcf: &Path,
    output_index: Option<&Path>,
    tbi: bool,
    force: bool,
) -> LibResult<()> {
    let mut argv = vec![OsString::from("index")];
    if tbi {
        argv.push(OsString::from("-t"));
    }
    if force {
        argv.push(OsString::from("-f"));
    }
    if let Some(path) = output_index {
        argv.push(OsString::from("-o"));
        argv.push(path.as_os_str().to_owned());
    }
    argv.push(input_vcf.as_os_str().to_owned());

    run_bcftools("index", bcftools_rs::commands::index::main(&argv))
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

pub fn view_header_native(input_vcf: &Path, output_vcf: &Path) -> LibResult<()> {
    let argv = [
        OsString::from("view"),
        OsString::from("--no-version"),
        OsString::from("-h"),
        OsString::from("-o"),
        output_vcf.as_os_str().to_owned(),
        input_vcf.as_os_str().to_owned(),
    ];
    run_bcftools(
        "view header extraction",
        bcftools_rs::commands::view::main(&argv),
    )
}

fn run_bcftools(operation: &str, status: ExitCode) -> LibResult<()> {
    match status {
        ExitCode::SUCCESS => Ok(()),
        status => Err(LibError::InvalidArguments(format!(
            "bcftools.{operation} failed with status {status:?}"
        ))),
    }
}
