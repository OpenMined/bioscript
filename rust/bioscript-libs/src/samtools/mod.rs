use std::{io::Read, path::Path};

#[cfg(any(target_arch = "wasm32", test))]
use bioscript_core::GenomicLocus;
use bioscript_formats::alignment::{DepthSummary, FastqPairSummary};
#[cfg(target_arch = "wasm32")]
use bioscript_formats::{
    GenotypeLoadOptions, GenotypeSourceFormat,
    alignment::{
        query_bam_depth_summary, write_bam_region, write_bam_region_bytes,
        write_bam_region_fastq_pair, write_bam_region_fastq_pair_bytes,
        write_cram_region_as_bam_bytes, write_cram_region_fastq_pair_bytes,
    },
};
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
/// `HTSlib`'s primary index discovery is co-location: for a data file `X` it
/// probes `X.csi` then `X.bai`. The runtime hands us the genome and its index
/// as two independent (often materialized-temp) paths, so when the caller
/// passes an index that is not already co-located we mirror it next to the
/// data file under `HTSlib`'s expected name. This keeps the samtools port
/// faithful to upstream (which finds the index by `HTSlib` lookup) while still
/// honoring an explicit index argument.
fn colocate_index(bam: &Path, index: Option<&Path>) -> LibResult<()> {
    let Some(index) = index else {
        return Ok(());
    };
    let index_extension = if bam
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("cram"))
    {
        "crai"
    } else {
        "bai"
    };
    let expected = std::path::PathBuf::from(format!("{}.{}", bam.display(), index_extension));
    if expected == index || expected.exists() {
        return Ok(());
    }
    std::fs::copy(index, &expected).map_err(samtools_error)?;
    Ok(())
}

fn colocate_reference_index(
    reference_fasta: Option<&Path>,
    reference_index: Option<&Path>,
) -> LibResult<()> {
    let (Some(reference_fasta), Some(reference_index)) = (reference_fasta, reference_index) else {
        return Ok(());
    };
    let expected = std::path::PathBuf::from(format!("{}.fai", reference_fasta.display()));
    if expected == reference_index || expected.exists() {
        return Ok(());
    }
    std::fs::copy(reference_index, &expected).map_err(samtools_error)?;
    Ok(())
}

pub fn view_region_native(
    bam: &Path,
    index: Option<&Path>,
    reference_fasta: Option<&Path>,
    reference_index: Option<&Path>,
    region: &str,
    output_bam: &Path,
) -> LibResult<usize> {
    #[cfg(target_arch = "wasm32")]
    {
        if is_bam_path(bam) {
            let locus = parse_region(region)?;
            return write_bam_region(bam, output_bam, &bam_load_options(index), &locus)
                .map_err(samtools_error);
        }
    }

    colocate_index(bam, index)?;
    colocate_reference_index(reference_fasta, reference_index)?;
    samtools_native::view_region_native(bam, region, output_bam, None, reference_fasta)
        .map_err(samtools_error)?;
    Ok(0)
}

#[cfg(target_arch = "wasm32")]
pub fn view_region_native_bytes(
    alignment_bytes: &[u8],
    index_bytes: &[u8],
    reference_fasta_bytes: Option<&[u8]>,
    reference_index_bytes: Option<&[u8]>,
    region: &str,
    input_format: GenotypeSourceFormat,
) -> LibResult<(Vec<u8>, usize)> {
    let locus = parse_region(region)?;
    match input_format {
        GenotypeSourceFormat::Cram => {
            let reference_fasta_bytes = reference_fasta_bytes.ok_or_else(|| {
                LibError::InvalidArguments("CRAM region extraction requires reference FASTA bytes".to_owned())
            })?;
            let reference_index_bytes = reference_index_bytes.ok_or_else(|| {
                LibError::InvalidArguments("CRAM region extraction requires reference FASTA index bytes".to_owned())
            })?;
            write_cram_region_as_bam_bytes(
                alignment_bytes,
                index_bytes,
                reference_fasta_bytes,
                reference_index_bytes,
                &locus,
            )
            .map_err(samtools_error)
        }
        _ => write_bam_region_bytes(alignment_bytes, index_bytes, &locus).map_err(samtools_error),
    }
}

pub fn depth_native(bam: &Path, index: Option<&Path>, region: &str) -> LibResult<DepthSummary> {
    #[cfg(target_arch = "wasm32")]
    {
        if is_bam_path(bam) {
            let locus = parse_region(region)?;
            return query_bam_depth_summary(bam, &bam_load_options(index), &locus)
                .map_err(samtools_error);
        }
    }

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
    reference_fasta: Option<&Path>,
    reference_index: Option<&Path>,
    region: &str,
    fastq_1: &Path,
    fastq_2: &Path,
) -> LibResult<FastqPairSummary> {
    #[cfg(target_arch = "wasm32")]
    {
        if is_bam_path(bam) {
            let locus = parse_region(region)?;
            return write_bam_region_fastq_pair(
                bam,
                fastq_1,
                fastq_2,
                &bam_load_options(index),
                &locus,
            )
            .map_err(samtools_error);
        }
    }

    colocate_index(bam, index)?;
    colocate_reference_index(reference_fasta, reference_index)?;
    let temp_dir = tempfile::tempdir().map_err(samtools_error)?;
    let sliced_bam = temp_dir.path().join("slice.bam");
    let other_fastq = temp_dir.path().join("other.fastq.gz");
    let singleton_fastq = temp_dir.path().join("singleton.fastq.gz");
    samtools_native::view_region_native(bam, region, &sliced_bam, None, reference_fasta)
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

#[cfg(target_arch = "wasm32")]
pub fn fastq_native_bytes(
    alignment_bytes: &[u8],
    index_bytes: &[u8],
    reference_fasta_bytes: Option<&[u8]>,
    reference_index_bytes: Option<&[u8]>,
    region: &str,
    gzip: bool,
    input_format: GenotypeSourceFormat,
) -> LibResult<(Vec<u8>, Vec<u8>, FastqPairSummary)> {
    let locus = parse_region(region)?;
    match input_format {
        GenotypeSourceFormat::Cram => {
            let reference_fasta_bytes = reference_fasta_bytes.ok_or_else(|| {
                LibError::InvalidArguments("CRAM FASTQ extraction requires reference FASTA bytes".to_owned())
            })?;
            let reference_index_bytes = reference_index_bytes.ok_or_else(|| {
                LibError::InvalidArguments("CRAM FASTQ extraction requires reference FASTA index bytes".to_owned())
            })?;
            write_cram_region_fastq_pair_bytes(
                alignment_bytes,
                index_bytes,
                reference_fasta_bytes,
                reference_index_bytes,
                &locus,
                gzip,
            )
            .map_err(samtools_error)
        }
        _ => write_bam_region_fastq_pair_bytes(alignment_bytes, index_bytes, &locus, gzip)
            .map_err(samtools_error),
    }
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
    let region_length_f64 = usize_to_f64(region_length);
    let uncovered_bases = depths.iter().filter(|depth| **depth == 0).count();
    let sum = depths.iter().map(|depth| f64::from(*depth)).sum::<f64>();
    let mean = sum / region_length_f64;
    let stdev = (depths
        .iter()
        .map(|depth| {
            let delta = f64::from(*depth) - mean;
            delta * delta
        })
        .sum::<f64>()
        / region_length_f64)
        .sqrt();
    let min = depths.iter().copied().min().unwrap_or(0);
    let max = depths.iter().copied().max().unwrap_or(0);
    depths.sort_unstable();
    let median = if region_length % 2 == 0 {
        let upper = region_length / 2;
        f64::midpoint(f64::from(depths[upper - 1]), f64::from(depths[upper]))
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
        percent_uncovered: usize_to_f64(uncovered_bases) / region_length_f64 * 100.0,
    }
}

fn usize_to_f64(value: usize) -> f64 {
    f64::from(u32::try_from(value).expect("samtools depth region length must fit in u32"))
}

#[cfg(target_arch = "wasm32")]
fn is_bam_path(path: &Path) -> bool {
    path.extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("bam"))
}

#[cfg(target_arch = "wasm32")]
fn bam_load_options(index: Option<&Path>) -> GenotypeLoadOptions {
    GenotypeLoadOptions {
        format: Some(GenotypeSourceFormat::Bam),
        input_index: index.map(Path::to_path_buf),
        ..GenotypeLoadOptions::default()
    }
}

#[cfg(any(target_arch = "wasm32", test))]
fn parse_region(region: &str) -> LibResult<GenomicLocus> {
    let (chrom, range) = region
        .split_once(':')
        .ok_or_else(|| LibError::InvalidArguments(format!("invalid samtools region: {region}")))?;
    let (start, end) = range
        .split_once('-')
        .ok_or_else(|| LibError::InvalidArguments(format!("invalid samtools region: {region}")))?;
    let start = start.replace(',', "").parse::<i64>().map_err(|err| {
        LibError::InvalidArguments(format!("invalid samtools region start {start}: {err}"))
    })?;
    let end = end.replace(',', "").parse::<i64>().map_err(|err| {
        LibError::InvalidArguments(format!("invalid samtools region end {end}: {err}"))
    })?;
    if chrom.is_empty() || start < 1 || end < start {
        return Err(LibError::InvalidArguments(format!(
            "invalid samtools region: {region}"
        )));
    }
    Ok(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    })
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

fn samtools_error(err: impl std::fmt::Display) -> LibError {
    LibError::InvalidArguments(err.to_string())
}

#[cfg(test)]
#[allow(clippy::float_cmp)]
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

    #[test]
    fn parse_region_accepts_samtools_coordinate_shape() {
        let locus = parse_region("chr1:155,158,000-155,163,000").unwrap();
        assert_eq!(locus.chrom, "chr1");
        assert_eq!(locus.start, 155_158_000);
        assert_eq!(locus.end, 155_163_000);

        let err = parse_region("chr1:10-1").unwrap_err();
        assert!(err.to_string().contains("invalid samtools region"));
    }

    #[test]
    fn colocate_index_uses_alignment_specific_extension() {
        let dir = tempfile::tempdir().unwrap();
        let bam = dir.path().join("input.bam");
        let cram = dir.path().join("input.cram");
        let index = dir.path().join("provided.index");
        let reference = dir.path().join("reference.fa");
        let reference_index = dir.path().join("provided.fai");
        std::fs::write(&bam, b"bam").unwrap();
        std::fs::write(&cram, b"cram").unwrap();
        std::fs::write(&index, b"index").unwrap();
        std::fs::write(&reference, b">chr1\nACGT\n").unwrap();
        std::fs::write(&reference_index, b"fai").unwrap();

        colocate_index(&bam, Some(&index)).unwrap();
        colocate_index(&cram, Some(&index)).unwrap();
        colocate_reference_index(Some(&reference), Some(&reference_index)).unwrap();

        assert_eq!(
            std::fs::read(dir.path().join("input.bam.bai")).unwrap(),
            b"index"
        );
        assert_eq!(
            std::fs::read(dir.path().join("input.cram.crai")).unwrap(),
            b"index"
        );
        assert_eq!(
            std::fs::read(dir.path().join("reference.fa.fai")).unwrap(),
            b"fai"
        );
        assert!(!dir.path().join("input.cram.bai").exists());
    }
}
