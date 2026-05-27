use std::{io::Read, path::Path};

use bioscript_formats::alignment::{DepthSummary, FastqPairSummary};
use samtools_rs::native as samtools_native;

use crate::{
    LibError, LibResult,
    tools::{CommandSpec, path_arg},
};

pub const MODULE: &str = "samtools";

pub fn view(bam: &Path, region: &str, output_bam: &Path) -> LibResult<CommandSpec> {
    view_region(bam, region, output_bam, false)
}

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

pub fn sort(bam: &Path, output_bam: &Path, by_name: bool) -> LibResult<CommandSpec> {
    let mut args = vec!["sort".to_owned()];
    if by_name {
        args.push("-n".to_owned());
    }
    args.extend(["-o".to_owned(), path_arg(output_bam)?, path_arg(bam)?]);
    CommandSpec::new("samtools", args)
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

pub fn faidx(fasta: &Path) -> LibResult<CommandSpec> {
    CommandSpec::new("samtools", vec!["faidx".to_owned(), path_arg(fasta)?])
}

/// Make an explicitly-provided index discoverable by samtools-rs/HTSlib.
///
/// HTSlib's primary index discovery is co-location: for a data file `X` it
/// probes `X.csi` then `X.bai`. The runtime hands us the genome and its index
/// as two independent (often materialized-temp) paths, so when the caller
/// passes an index that is not already co-located we mirror it next to the
/// data file under HTSlib's expected name. This keeps the samtools port
/// faithful to upstream (which finds the index by HTSlib lookup) while still
/// honoring an explicit index argument.
fn colocate_index(bam: &Path, index: Option<&Path>) -> LibResult<()> {
    let Some(index) = index else {
        return Ok(());
    };
    let expected = std::path::PathBuf::from(format!("{}.bai", bam.display()));
    if expected == index || expected.exists() {
        return Ok(());
    }
    std::fs::copy(index, &expected).map_err(samtools_error)?;
    Ok(())
}

pub fn view_region_native(
    bam: &Path,
    index: Option<&Path>,
    region: &str,
    output_bam: &Path,
) -> LibResult<usize> {
    colocate_index(bam, index)?;
    samtools_native::view_region_native(bam, region, output_bam, None, None)
        .map_err(samtools_error)?;
    Ok(0)
}

pub fn depth_native(bam: &Path, index: Option<&Path>, region: &str) -> LibResult<DepthSummary> {
    colocate_index(bam, index)?;
    let depths = samtools_native::depth_native(bam, region, true, None).map_err(samtools_error)?;
    Ok(depth_summary(depths.iter().map(|entry| entry.depth)))
}

pub fn sort_native(bam: &Path, output_bam: &Path, by_name: bool) -> LibResult<()> {
    samtools_native::sort_native(bam, output_bam, by_name, None).map_err(samtools_error)
}

pub fn index_native(bam: &Path, output_bai: Option<&Path>) -> LibResult<std::path::PathBuf> {
    samtools_native::index_native(bam, output_bai, None).map_err(samtools_error)
}

pub fn fastq_native(
    bam: &Path,
    index: Option<&Path>,
    region: &str,
    fastq_1: &Path,
    fastq_2: &Path,
) -> LibResult<FastqPairSummary> {
    colocate_index(bam, index)?;
    let temp_dir = tempfile::tempdir().map_err(samtools_error)?;
    let sliced_bam = temp_dir.path().join("slice.bam");
    let other_fastq = temp_dir.path().join("other.fastq.gz");
    let singleton_fastq = temp_dir.path().join("singleton.fastq.gz");
    samtools_native::view_region_native(bam, region, &sliced_bam, None, None)
        .map_err(samtools_error)?;
    samtools_native::fastq_native(
        &sliced_bam,
        fastq_1,
        fastq_2,
        Some(&other_fastq),
        Some(&singleton_fastq),
        true,
        None,
    )
    .map_err(samtools_error)?;
    Ok(FastqPairSummary {
        read1_records: fastq_record_count(fastq_1)?,
        read2_records: fastq_record_count(fastq_2)?,
        skipped_records: 0,
    })
}

pub fn fastq_all_native(bam: &Path, fastq_1: &Path, fastq_2: &Path) -> LibResult<FastqPairSummary> {
    let temp_dir = tempfile::tempdir().map_err(samtools_error)?;
    let other_fastq = temp_dir.path().join("other.fastq.gz");
    let singleton_fastq = temp_dir.path().join("singleton.fastq.gz");
    samtools_native::fastq_native(
        bam,
        fastq_1,
        fastq_2,
        Some(&other_fastq),
        Some(&singleton_fastq),
        true,
        None,
    )
    .map_err(samtools_error)?;
    Ok(FastqPairSummary {
        read1_records: fastq_record_count(fastq_1)?,
        read2_records: fastq_record_count(fastq_2)?,
        skipped_records: fastq_record_count(&singleton_fastq)?,
    })
}

fn depth_summary(depths: impl IntoIterator<Item = u32>) -> DepthSummary {
    let mut depths = depths.into_iter().collect::<Vec<_>>();
    if depths.is_empty() {
        return DepthSummary {
            mean: 0.0,
            median: 0.0,
            stdev: 0.0,
            min: 0,
            max: 0,
            region_length: 0,
            uncovered_bases: 0,
            percent_uncovered: 0.0,
        };
    }
    let region_length = depths.len();
    let uncovered_bases = depths.iter().filter(|depth| **depth == 0).count();
    let sum = depths.iter().map(|depth| f64::from(*depth)).sum::<f64>();
    let mean = sum / region_length as f64;
    let stdev = (depths
        .iter()
        .map(|depth| {
            let delta = f64::from(*depth) - mean;
            delta * delta
        })
        .sum::<f64>()
        / region_length as f64)
        .sqrt();
    let min = depths.iter().copied().min().unwrap_or(0);
    let max = depths.iter().copied().max().unwrap_or(0);
    depths.sort_unstable();
    let median = if region_length % 2 == 0 {
        let upper = region_length / 2;
        (f64::from(depths[upper - 1]) + f64::from(depths[upper])) / 2.0
    } else {
        f64::from(depths[region_length / 2])
    };
    DepthSummary {
        mean,
        median,
        stdev,
        min,
        max,
        region_length,
        uncovered_bases,
        percent_uncovered: uncovered_bases as f64 / region_length as f64 * 100.0,
    }
}

fn fastq_record_count(path: &Path) -> LibResult<usize> {
    let mut bytes = Vec::new();
    if path.extension().is_some_and(|extension| extension == "gz") {
        let file = std::fs::File::open(path).map_err(samtools_error)?;
        flate2::read::GzDecoder::new(file)
            .read_to_end(&mut bytes)
            .map_err(samtools_error)?;
    } else {
        bytes = std::fs::read(path).map_err(samtools_error)?;
    }
    let content = String::from_utf8(bytes)
        .map_err(|err| LibError::InvalidArguments(format!("FASTQ output is not UTF-8: {err}")))?;
    Ok(content.lines().step_by(4).count())
}

fn samtools_error(err: std::io::Error) -> LibError {
    LibError::InvalidArguments(err.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_depth_summary_matches_bioscript_shape() {
        let summary = depth_summary([10, 0, 20]);
        assert_eq!(summary.mean, 10.0);
        assert_eq!(summary.median, 10.0);
        assert_eq!(summary.min, 0);
        assert_eq!(summary.max, 20);
        assert_eq!(summary.region_length, 3);
        assert_eq!(summary.uncovered_bases, 1);
    }
}
