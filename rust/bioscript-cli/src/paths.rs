use std::{
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
};

use bioscript_formats::GenotypeLoadOptions;
use bioscript_runtime::StageTiming;

pub(crate) fn write_timing_report(path: &PathBuf, timings: &[StageTiming]) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|err| {
            format!(
                "failed to create timing report dir {}: {err}",
                parent.display()
            )
        })?;
    }
    let mut output = String::from("stage\tduration_ms\tdetail\n");
    for timing in timings {
        let _ = writeln!(
            output,
            "{}\t{}\t{}",
            timing.stage,
            timing.duration_ms,
            timing.detail.replace('\t', " ")
        );
    }
    fs::write(path, output)
        .map_err(|err| format!("failed to write timing report {}: {err}", path.display()))
}

pub(crate) fn resolve_cli_path(root: &Path, value: &str) -> String {
    resolve_cli_path_buf(root, Path::new(value))
        .display()
        .to_string()
}

pub(crate) fn resolve_cli_path_buf(root: &Path, value: &Path) -> PathBuf {
    if value.is_absolute() {
        value.to_path_buf()
    } else {
        root.join(value)
    }
}

pub(crate) fn normalize_loader_paths(root: &Path, loader: &mut GenotypeLoadOptions) {
    if let Some(path) = loader.input_index.take() {
        loader.input_index = Some(resolve_cli_path_buf(root, &path));
    }
    if let Some(path) = loader.reference_file.take() {
        loader.reference_file = Some(resolve_cli_path_buf(root, &path));
    }
    if let Some(path) = loader.reference_index.take() {
        loader.reference_index = Some(resolve_cli_path_buf(root, &path));
    }
}
