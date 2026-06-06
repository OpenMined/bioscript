#![allow(clippy::missing_errors_doc)]

use std::{collections::HashMap, path::PathBuf};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
fn supported_modules() -> Vec<&'static str> {
    bioscript_libs::supported_modules()
        .iter()
        .map(|module| module.name.as_str())
        .collect()
}

#[pyfunction]
fn samtools_view_region_native(
    bam: &str,
    index: Option<&str>,
    reference_fasta: Option<&str>,
    reference_index: Option<&str>,
    region: &str,
    output_bam: &str,
) -> PyResult<usize> {
    bioscript_libs::samtools::view_region_native(
        PathBuf::from(bam).as_path(),
        index.map(PathBuf::from).as_deref(),
        reference_fasta.map(PathBuf::from).as_deref(),
        reference_index.map(PathBuf::from).as_deref(),
        region,
        PathBuf::from(output_bam).as_path(),
    )
    .map_err(to_py_value_error)
}

#[pyfunction]
fn samtools_depth_native(
    bam: &str,
    index: Option<&str>,
    region: &str,
) -> PyResult<HashMap<&'static str, f64>> {
    let summary = bioscript_libs::samtools::depth_native(
        PathBuf::from(bam).as_path(),
        index.map(PathBuf::from).as_deref(),
        region,
    )
    .map_err(to_py_value_error)?;
    Ok(HashMap::from([
        ("mean", summary.mean),
        ("median", summary.median),
        ("stdev", summary.stdev),
        ("min", f64::from(summary.min)),
        ("max", f64::from(summary.max)),
        ("region_length", summary.region_length as f64),
        ("uncovered_bases", summary.uncovered_bases as f64),
        ("percent_uncovered", summary.percent_uncovered),
    ]))
}

#[pyfunction]
fn samtools_fastq_native(
    bam: &str,
    index: Option<&str>,
    reference_fasta: Option<&str>,
    reference_index: Option<&str>,
    region: &str,
    fastq_1: &str,
    fastq_2: &str,
) -> PyResult<HashMap<&'static str, usize>> {
    let summary = bioscript_libs::samtools::fastq_native(
        PathBuf::from(bam).as_path(),
        index.map(PathBuf::from).as_deref(),
        reference_fasta.map(PathBuf::from).as_deref(),
        reference_index.map(PathBuf::from).as_deref(),
        region,
        PathBuf::from(fastq_1).as_path(),
        PathBuf::from(fastq_2).as_path(),
    )
    .map_err(to_py_value_error)?;
    Ok(HashMap::from([
        ("read1_records", summary.read1_records),
        ("read2_records", summary.read2_records),
        ("skipped_records", summary.skipped_records),
    ]))
}

#[pyfunction]
fn bcftools_view_header_native(input_vcf: &str, output_vcf: &str) -> PyResult<()> {
    bioscript_libs::bcftools::view_header_native(
        PathBuf::from(input_vcf).as_path(),
        PathBuf::from(output_vcf).as_path(),
    )
    .map_err(to_py_value_error)
}

#[pyfunction]
fn bcftools_view_native(input_vcf: &str, output_vcf: &str, output_type: &str) -> PyResult<()> {
    bioscript_libs::bcftools::view_native(
        PathBuf::from(input_vcf).as_path(),
        PathBuf::from(output_vcf).as_path(),
        output_type,
    )
    .map_err(to_py_value_error)
}

#[pyfunction]
fn bcftools_sort_native(
    input_vcf: &str,
    output_vcf: &str,
    output_type: &str,
    write_index: bool,
) -> PyResult<()> {
    bioscript_libs::bcftools::sort_native(
        PathBuf::from(input_vcf).as_path(),
        PathBuf::from(output_vcf).as_path(),
        output_type,
        write_index,
    )
    .map_err(to_py_value_error)
}

#[pyfunction]
fn bcftools_index_native(
    input_vcf: &str,
    output_index: Option<&str>,
    tbi: bool,
    force: bool,
) -> PyResult<()> {
    let input = PathBuf::from(input_vcf);
    let output = output_index.map(PathBuf::from);
    bioscript_libs::bcftools::index_native(&input, output.as_deref(), tbi, force)
        .map_err(to_py_value_error)
}

#[pyfunction]
fn pyfaidx_fetch_native(path: &str, contig: &str, start: usize, stop: usize) -> PyResult<String> {
    let fasta = bioscript_libs::pyfaidx::Fasta::from_path(PathBuf::from(path))
        .map_err(to_py_value_error)?;
    let record = fasta.get(contig).map_err(to_py_value_error)?;
    record.slice(start, stop).map_err(to_py_value_error)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn kestrel_call_sequences_native(
    reference_name: &str,
    reference_sequence: &str,
    read_sequences: Vec<String>,
    kmer_size: usize,
    sample_name: &str,
    source_version: Option<&str>,
    reference_md5: Option<&str>,
    minimum_difference: Option<u32>,
    difference_quantile: Option<f32>,
    anchor_both_ends: Option<bool>,
    decay_min: Option<f32>,
    decay_alpha: Option<f32>,
    peak_scan_length: Option<usize>,
    scan_limit_factor: Option<f32>,
    max_gap_size: Option<usize>,
    recover_right_anchor: Option<bool>,
    call_ambiguous_regions: Option<bool>,
    min_kmer_count: Option<u32>,
    max_haplotypes: Option<usize>,
    max_bases: Option<usize>,
    max_repeat_count: Option<usize>,
    max_saved_states: Option<usize>,
    locus_depth: Option<u32>,
) -> PyResult<String> {
    let _ = (
        source_version,
        reference_md5,
        max_gap_size,
        recover_right_anchor,
        max_bases,
        locus_depth,
    );
    let options = kestrel_options(
        sample_name,
        minimum_difference,
        difference_quantile,
        anchor_both_ends,
        decay_min,
        decay_alpha,
        peak_scan_length,
        scan_limit_factor,
        call_ambiguous_regions,
        min_kmer_count,
        max_haplotypes,
        max_repeat_count,
        max_saved_states,
    );
    bioscript_libs::kestrel::native::call_sequences_to_vcf(
        reference_name,
        reference_sequence,
        read_sequences.iter().map(String::as_str),
        kmer_size,
        &options,
    )
    .map_err(to_py_value_error)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn kestrel_call_fastq_native(
    reference_name: &str,
    reference_sequence: &str,
    fastq_paths: Vec<String>,
    kmer_size: usize,
    sample_name: &str,
    source_version: Option<&str>,
    reference_md5: Option<&str>,
    minimum_difference: Option<u32>,
    difference_quantile: Option<f32>,
    anchor_both_ends: Option<bool>,
    decay_min: Option<f32>,
    decay_alpha: Option<f32>,
    peak_scan_length: Option<usize>,
    scan_limit_factor: Option<f32>,
    max_gap_size: Option<usize>,
    recover_right_anchor: Option<bool>,
    call_ambiguous_regions: Option<bool>,
    min_kmer_count: Option<u32>,
    max_haplotypes: Option<usize>,
    max_bases: Option<usize>,
    max_repeat_count: Option<usize>,
    max_saved_states: Option<usize>,
    locus_depth: Option<u32>,
) -> PyResult<String> {
    let _ = (
        source_version,
        reference_md5,
        max_gap_size,
        recover_right_anchor,
        max_bases,
        locus_depth,
    );
    let options = kestrel_options(
        sample_name,
        minimum_difference,
        difference_quantile,
        anchor_both_ends,
        decay_min,
        decay_alpha,
        peak_scan_length,
        scan_limit_factor,
        call_ambiguous_regions,
        min_kmer_count,
        max_haplotypes,
        max_repeat_count,
        max_saved_states,
    );
    let paths: Vec<PathBuf> = fastq_paths.into_iter().map(PathBuf::from).collect();
    bioscript_libs::kestrel::native::call_fastq_paths_to_vcf(
        reference_name,
        reference_sequence,
        paths.iter().map(PathBuf::as_path),
        kmer_size,
        &options,
    )
    .map_err(to_py_value_error)
}

#[allow(clippy::too_many_arguments)]
#[pyfunction]
fn kestrel_call_fastq_references_native(
    references: Vec<(String, String, String)>,
    fastq_paths: Vec<String>,
    kmer_size: usize,
    sample_name: &str,
    source_version: Option<&str>,
    minimum_difference: Option<u32>,
    difference_quantile: Option<f32>,
    anchor_both_ends: Option<bool>,
    decay_min: Option<f32>,
    decay_alpha: Option<f32>,
    peak_scan_length: Option<usize>,
    scan_limit_factor: Option<f32>,
    max_gap_size: Option<usize>,
    recover_right_anchor: Option<bool>,
    call_ambiguous_regions: Option<bool>,
    min_kmer_count: Option<u32>,
    max_haplotypes: Option<usize>,
    max_bases: Option<usize>,
    max_repeat_count: Option<usize>,
    max_saved_states: Option<usize>,
    locus_depth: Option<u32>,
) -> PyResult<String> {
    let references: Vec<bioscript_libs::kestrel::native::NativeReferenceRegion> = references
        .into_iter()
        .map(|(name, sequence, md5)| {
            bioscript_libs::kestrel::native::NativeReferenceRegion::new(name, sequence, md5)
        })
        .collect();
    let _ = (
        source_version,
        max_gap_size,
        recover_right_anchor,
        max_bases,
        locus_depth,
    );
    let options = kestrel_options(
        sample_name,
        minimum_difference,
        difference_quantile,
        anchor_both_ends,
        decay_min,
        decay_alpha,
        peak_scan_length,
        scan_limit_factor,
        call_ambiguous_regions,
        min_kmer_count,
        max_haplotypes,
        max_repeat_count,
        max_saved_states,
    );
    let paths: Vec<PathBuf> = fastq_paths.into_iter().map(PathBuf::from).collect();
    bioscript_libs::kestrel::native::call_fastq_paths_to_vcf_references(
        &references,
        paths.iter().map(PathBuf::as_path),
        kmer_size,
        &options,
    )
    .map_err(to_py_value_error)
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(supported_modules, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_view_region_native, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_depth_native, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_fastq_native, module)?)?;
    module.add_function(wrap_pyfunction!(bcftools_view_header_native, module)?)?;
    module.add_function(wrap_pyfunction!(bcftools_view_native, module)?)?;
    module.add_function(wrap_pyfunction!(bcftools_sort_native, module)?)?;
    module.add_function(wrap_pyfunction!(bcftools_index_native, module)?)?;
    module.add_function(wrap_pyfunction!(pyfaidx_fetch_native, module)?)?;
    module.add_function(wrap_pyfunction!(kestrel_call_sequences_native, module)?)?;
    module.add_function(wrap_pyfunction!(kestrel_call_fastq_native, module)?)?;
    module.add_function(wrap_pyfunction!(
        kestrel_call_fastq_references_native,
        module
    )?)?;
    Ok(())
}

fn to_py_value_error(err: bioscript_libs::LibError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

#[allow(clippy::too_many_arguments)]
fn kestrel_options(
    sample_name: &str,
    minimum_difference: Option<u32>,
    difference_quantile: Option<f32>,
    anchor_both_ends: Option<bool>,
    decay_min: Option<f32>,
    decay_alpha: Option<f32>,
    peak_scan_length: Option<usize>,
    scan_limit_factor: Option<f32>,
    call_ambiguous_regions: Option<bool>,
    min_kmer_count: Option<u32>,
    max_haplotypes: Option<usize>,
    max_repeat_count: Option<usize>,
    max_saved_states: Option<usize>,
) -> bioscript_libs::kestrel::native::NativeKestrelRunOptions {
    let mut options = bioscript_libs::kestrel::native::NativeKestrelRunOptions::new(sample_name);
    options.minimum_difference = minimum_difference.unwrap_or(options.minimum_difference);
    options.difference_quantile = difference_quantile.unwrap_or(options.difference_quantile);
    options.anchor_both_ends = anchor_both_ends.unwrap_or(options.anchor_both_ends);
    options.decay_min = decay_min.unwrap_or(options.decay_min);
    options.decay_alpha = decay_alpha.unwrap_or(options.decay_alpha);
    options.peak_scan_length = peak_scan_length.unwrap_or(options.peak_scan_length);
    options.scan_limit_factor = scan_limit_factor.unwrap_or(options.scan_limit_factor);
    options.call_ambiguous_regions =
        call_ambiguous_regions.unwrap_or(options.call_ambiguous_regions);
    options.min_kmer_count = min_kmer_count.unwrap_or(options.min_kmer_count);
    options.max_haplotypes = max_haplotypes.unwrap_or(options.max_haplotypes);
    options.max_repeat_count = max_repeat_count.unwrap_or(options.max_repeat_count);
    options.max_saved_states = max_saved_states.unwrap_or(options.max_saved_states);
    options
}
