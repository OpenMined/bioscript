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

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(supported_modules, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_view_region_native, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_depth_native, module)?)?;
    module.add_function(wrap_pyfunction!(samtools_fastq_native, module)?)?;
    Ok(())
}

fn to_py_value_error(err: bioscript_libs::LibError) -> PyErr {
    PyValueError::new_err(err.to_string())
}
