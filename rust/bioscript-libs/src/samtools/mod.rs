use std::path::Path;

use crate::{
    LibResult,
    tools::{CommandSpec, path_arg},
};

pub const MODULE: &str = "samtools";

pub fn view_region(
    bam: &Path,
    region: &str,
    output_bam: &Path,
    include_unmapped: bool,
) -> LibResult<CommandSpec> {
    let mut args = vec![
        "view".to_owned(),
        "-b".to_owned(),
        path_arg(bam)?,
        region.to_owned(),
        "-o".to_owned(),
        path_arg(output_bam)?,
    ];
    if include_unmapped {
        args.push("-f".to_owned());
        args.push("4".to_owned());
    }
    CommandSpec::new("samtools", args)
}

pub fn fastq(bam: &Path, fastq_1: &Path, fastq_2: &Path) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "samtools",
        vec![
            "fastq".to_owned(),
            "-1".to_owned(),
            path_arg(fastq_1)?,
            "-2".to_owned(),
            path_arg(fastq_2)?,
            path_arg(bam)?,
        ],
    )
}

pub fn depth(bam: &Path, region: &str) -> LibResult<CommandSpec> {
    CommandSpec::new(
        "samtools",
        vec![
            "depth".to_owned(),
            "-r".to_owned(),
            region.to_owned(),
            path_arg(bam)?,
        ],
    )
}

pub fn index(bam: &Path) -> LibResult<CommandSpec> {
    CommandSpec::new("samtools", vec!["index".to_owned(), path_arg(bam)?])
}
