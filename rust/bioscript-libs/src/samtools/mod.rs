use std::path::Path;

use bioscript_core::GenomicLocus;
use bioscript_formats::{GenotypeLoadOptions, alignment};

use crate::{
    LibError, LibResult,
    tools::{CommandSpec, path_arg},
};

pub const MODULE: &str = "samtools";

pub use alignment::{DepthSummary, FastqPairSummary};

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

pub fn view_region_native(
    bam: &Path,
    index: Option<&Path>,
    region: &str,
    output_bam: &Path,
) -> LibResult<usize> {
    alignment::write_bam_region(bam, output_bam, &options(index), &parse_region(region)?)
        .map_err(|err| LibError::InvalidArguments(err.to_string()))
}

pub fn depth_native(bam: &Path, index: Option<&Path>, region: &str) -> LibResult<DepthSummary> {
    alignment::query_bam_depth_summary(bam, &options(index), &parse_region(region)?)
        .map_err(|err| LibError::InvalidArguments(err.to_string()))
}

pub fn fastq_native(
    bam: &Path,
    index: Option<&Path>,
    region: &str,
    fastq_1: &Path,
    fastq_2: &Path,
) -> LibResult<FastqPairSummary> {
    alignment::write_bam_region_fastq_pair(
        bam,
        fastq_1,
        fastq_2,
        &options(index),
        &parse_region(region)?,
    )
    .map_err(|err| LibError::InvalidArguments(err.to_string()))
}

fn options(index: Option<&Path>) -> GenotypeLoadOptions {
    GenotypeLoadOptions {
        input_index: index.map(Path::to_path_buf),
        ..GenotypeLoadOptions::default()
    }
}

fn parse_region(region: &str) -> LibResult<GenomicLocus> {
    let Some((chrom, coordinates)) = region.split_once(':') else {
        return Err(LibError::InvalidArguments(format!(
            "samtools region must be chrom:start-end, got {region:?}"
        )));
    };
    if chrom.is_empty() {
        return Err(LibError::InvalidArguments(
            "samtools region chromosome cannot be empty".to_owned(),
        ));
    }
    let Some((start, end)) = coordinates.split_once('-') else {
        return Err(LibError::InvalidArguments(format!(
            "samtools region must include start-end, got {region:?}"
        )));
    };
    let start = parse_position(start, "start")?;
    let end = parse_position(end, "end")?;
    if end < start {
        return Err(LibError::InvalidArguments(
            "samtools region end must be >= start".to_owned(),
        ));
    }
    Ok(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    })
}

fn parse_position(value: &str, label: &str) -> LibResult<i64> {
    let position = value.replace(',', "").parse::<i64>().map_err(|_| {
        LibError::InvalidArguments(format!(
            "samtools region {label} is not an integer: {value:?}"
        ))
    })?;
    if position < 1 {
        return Err(LibError::InvalidArguments(format!(
            "samtools region {label} must be >= 1"
        )));
    }
    Ok(position)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_region_parser_accepts_commas() {
        let locus = parse_region("chr1:155,160,500-155,162,000").unwrap();
        assert_eq!(locus.chrom, "chr1");
        assert_eq!(locus.start, 155_160_500);
        assert_eq!(locus.end, 155_162_000);
    }

    #[test]
    fn native_region_parser_rejects_bad_ranges() {
        assert!(parse_region("chr1").is_err());
        assert!(parse_region(":1-2").is_err());
        assert!(parse_region("chr1:0-2").is_err());
        assert!(parse_region("chr1:3-2").is_err());
    }
}
