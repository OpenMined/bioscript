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
    region: &str,
    output_bam: &str,
) -> PyResult<usize> {
    bioscript_libs::samtools::view_region_native(
        PathBuf::from(bam).as_path(),
        index.map(PathBuf::from).as_deref(),
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
    region: &str,
    fastq_1: &str,
    fastq_2: &str,
) -> PyResult<HashMap<&'static str, usize>> {
    let summary = bioscript_libs::samtools::fastq_native(
        PathBuf::from(bam).as_path(),
        index.map(PathBuf::from).as_deref(),
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
    let region = bioscript_libs::kestrel::native::ReferenceRegion {
        reference_name: reference_name.to_owned(),
        sequence: reference_sequence.to_owned(),
    };
    let detector_config = bioscript_libs::kestrel::native::ActiveRegionDetectorConfig {
        minimum_difference: minimum_difference.unwrap_or(5),
        difference_quantile: difference_quantile.unwrap_or(0.90),
        count_reverse_kmers: true,
        anchor_both_ends: anchor_both_ends.unwrap_or(true),
        decay_min: decay_min.unwrap_or(0.55),
        decay_alpha: decay_alpha.unwrap_or(0.80),
        peak_scan_length: peak_scan_length.unwrap_or(7),
        scan_limit_factor: scan_limit_factor.unwrap_or(7.0),
        max_gap_size: max_gap_size.unwrap_or(0),
        recover_right_anchor: recover_right_anchor.unwrap_or(true),
        call_ambiguous_regions: call_ambiguous_regions.unwrap_or(true),
    };
    let assembly_config = bioscript_libs::kestrel::native::HaplotypeAssemblyConfig {
        min_kmer_count: min_kmer_count.unwrap_or(1),
        max_haplotypes: max_haplotypes.unwrap_or(40),
        max_bases: max_bases.unwrap_or(500),
        max_repeat_count: max_repeat_count.unwrap_or(0),
        max_saved_states: max_saved_states.unwrap_or(40),
        locus_depth: locus_depth.unwrap_or(1),
    };
    let call_config = bioscript_libs::kestrel::native::NativeKestrelCallConfig::new(
        source_version.unwrap_or("native"),
        sample_name,
        reference_md5.unwrap_or("."),
    );
    bioscript_libs::kestrel::native::call_sequences_to_vcf(
        &region,
        read_sequences.iter().map(String::as_str),
        kmer_size,
        &detector_config,
        &assembly_config,
        &call_config,
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
    let region = bioscript_libs::kestrel::native::ReferenceRegion {
        reference_name: reference_name.to_owned(),
        sequence: reference_sequence.to_owned(),
    };
    let detector_config = bioscript_libs::kestrel::native::ActiveRegionDetectorConfig {
        minimum_difference: minimum_difference.unwrap_or(5),
        difference_quantile: difference_quantile.unwrap_or(0.90),
        count_reverse_kmers: true,
        anchor_both_ends: anchor_both_ends.unwrap_or(true),
        decay_min: decay_min.unwrap_or(0.55),
        decay_alpha: decay_alpha.unwrap_or(0.80),
        peak_scan_length: peak_scan_length.unwrap_or(7),
        scan_limit_factor: scan_limit_factor.unwrap_or(7.0),
        max_gap_size: max_gap_size.unwrap_or(0),
        recover_right_anchor: recover_right_anchor.unwrap_or(true),
        call_ambiguous_regions: call_ambiguous_regions.unwrap_or(true),
    };
    let assembly_config = bioscript_libs::kestrel::native::HaplotypeAssemblyConfig {
        min_kmer_count: min_kmer_count.unwrap_or(1),
        max_haplotypes: max_haplotypes.unwrap_or(40),
        max_bases: max_bases.unwrap_or(500),
        max_repeat_count: max_repeat_count.unwrap_or(0),
        max_saved_states: max_saved_states.unwrap_or(40),
        locus_depth: locus_depth.unwrap_or(1),
    };
    let call_config = bioscript_libs::kestrel::native::NativeKestrelCallConfig::new(
        source_version.unwrap_or("native"),
        sample_name,
        reference_md5.unwrap_or("."),
    );
    let paths: Vec<PathBuf> = fastq_paths.into_iter().map(PathBuf::from).collect();
    bioscript_libs::kestrel::native::call_fastq_paths_to_vcf(
        &region,
        paths.iter().map(PathBuf::as_path),
        kmer_size,
        &detector_config,
        &assembly_config,
        &call_config,
    )
    .map_err(to_py_value_error)
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(supported_modules, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_view_region_native, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_depth_native, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_fastq_native, module)?)?;
    module.add_function(wrap_pyfunction!(kestrel_call_sequences_native, module)?)?;
    module.add_function(wrap_pyfunction!(kestrel_call_fastq_native, module)?)?;
    Ok(())
}

fn to_py_value_error(err: bioscript_libs::LibError) -> PyErr {
    PyValueError::new_err(err.to_string())
}
