#![allow(clippy::missing_errors_doc)]

use pyo3::prelude::*;

#[pyfunction]
fn supported_modules() -> Vec<&'static str> {
    bioscript_libs::supported_modules()
        .iter()
        .map(|module| module.name.as_str())
        .collect()
}

#[pymodule]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(supported_modules, module)?)?;
    Ok(())
}
