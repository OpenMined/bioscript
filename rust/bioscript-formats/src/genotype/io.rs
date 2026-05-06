use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use zip::ZipArchive;

use bioscript_core::RuntimeError;

use super::GenotypeSourceFormat;

pub(crate) fn is_bgzf_path(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
}

pub(crate) fn read_plain_lines(path: &Path) -> Result<Vec<String>, RuntimeError> {
    let file = File::open(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open genotype file {}: {err}",
            path.display()
        ))
    })?;
    read_lines_from_reader(BufReader::new(file), path)
}

pub(crate) fn select_zip_entry(path: &Path) -> Result<String, RuntimeError> {
    let file = File::open(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open genotype zip {}: {err}",
            path.display()
        ))
    })?;
    let mut archive = ZipArchive::new(file).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to read genotype zip {}: {err}",
            path.display()
        ))
    })?;

    let mut selected_name: Option<String> = None;
    for idx in 0..archive.len() {
        let entry = archive.by_index(idx).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to inspect genotype zip {}: {err}",
                path.display()
            ))
        })?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".txt")
            || lower.ends_with(".csv")
            || lower.ends_with(".tsv")
            || lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
        {
            return Ok(name);
        }
        if selected_name.is_none() {
            selected_name = Some(name);
        }
    }

    selected_name.ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "zip archive {} does not contain a supported genotype file",
            path.display()
        ))
    })
}

pub(crate) fn read_lines_from_reader<R: BufRead>(
    mut reader: R,
    path: &Path,
) -> Result<Vec<String>, RuntimeError> {
    let mut lines = Vec::new();
    let mut buf = String::new();
    loop {
        buf.clear();
        let bytes = reader.read_line(&mut buf).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read genotype file {}: {err}",
                path.display()
            ))
        })?;
        if bytes == 0 {
            break;
        }
        lines.push(buf.trim_end_matches(['\n', '\r']).to_owned());
    }
    Ok(lines)
}

pub(crate) fn read_zip_entry_limited<R: Read>(
    reader: &mut R,
    max_bytes: u64,
    label: &str,
) -> Result<Vec<u8>, RuntimeError> {
    let mut contents = Vec::new();
    reader
        .take(max_bytes.saturating_add(1))
        .read_to_end(&mut contents)
        .map_err(|err| RuntimeError::Io(format!("failed to read {label}: {err}")))?;
    if u64::try_from(contents.len()).unwrap_or(u64::MAX) > max_bytes {
        return Err(RuntimeError::InvalidArguments(format!(
            "{label} exceeds decompressed limit of {max_bytes} bytes"
        )));
    }
    Ok(contents)
}

pub(crate) fn detect_source_format(
    path: &Path,
    forced: Option<GenotypeSourceFormat>,
) -> Result<GenotypeSourceFormat, RuntimeError> {
    if let Some(format) = forced {
        return Ok(format);
    }

    let lower = path.to_string_lossy().to_ascii_lowercase();
    if lower.ends_with(".zip") {
        return Ok(GenotypeSourceFormat::Zip);
    }
    if lower.ends_with(".cram") {
        return Ok(GenotypeSourceFormat::Cram);
    }
    if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
        return Ok(GenotypeSourceFormat::Vcf);
    }

    let lines = read_plain_lines(path)?;
    if looks_like_vcf_lines(&lines) {
        Ok(GenotypeSourceFormat::Vcf)
    } else {
        Ok(GenotypeSourceFormat::Text)
    }
}

pub(crate) fn looks_like_vcf_lines(lines: &[String]) -> bool {
    lines.iter().any(|line| {
        let trimmed = line.trim_start();
        trimmed.starts_with("##fileformat=VCF") || trimmed.starts_with("#CHROM\t")
    })
}
