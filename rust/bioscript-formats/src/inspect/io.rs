use std::{
    fs::File,
    io::{BufRead, BufReader, Cursor, Read},
    path::Path,
};

use bioscript_core::RuntimeError;
use noodles::bgzf;
use zip::ZipArchive;

const MAX_ZIP_SAMPLE_ENTRY_BYTES: u64 = 128 * 1024 * 1024;

pub(crate) fn read_plain_sample_lines_from_bytes(
    lower_name: &str,
    bytes: &[u8],
) -> Result<Vec<String>, RuntimeError> {
    if lower_name.ends_with(".vcf.gz") {
        return read_sample_lines_from_reader(BufReader::new(bgzf::io::Reader::new(Cursor::new(
            bytes,
        ))));
    }
    read_sample_lines_from_reader(BufReader::new(Cursor::new(bytes)))
}

pub(crate) fn read_zip_sample_lines_from_bytes(
    bytes: &[u8],
    selected_entry: &str,
) -> Result<Vec<String>, RuntimeError> {
    let mut archive = ZipArchive::new(Cursor::new(bytes))
        .map_err(|err| RuntimeError::Io(format!("failed to read zip bytes: {err}")))?;
    let mut entry = archive.by_name(selected_entry).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open zip entry {selected_entry} from bytes: {err}"
        ))
    })?;
    if selected_entry.to_ascii_lowercase().ends_with(".vcf.gz") {
        let inner = read_entry_limited(
            &mut entry,
            MAX_ZIP_SAMPLE_ENTRY_BYTES,
            &format!("compressed zip entry {selected_entry}"),
        )?;
        let reader = bgzf::io::Reader::new(Cursor::new(inner));
        return read_sample_lines_from_reader(BufReader::new(reader));
    }
    read_sample_lines_from_reader(BufReader::new(entry))
}

pub(crate) fn select_zip_entry_from_bytes(bytes: &[u8]) -> Result<String, RuntimeError> {
    let mut archive = ZipArchive::new(Cursor::new(bytes))
        .map_err(|err| RuntimeError::Io(format!("failed to read zip bytes: {err}")))?;
    let mut fallback = None;
    for idx in 0..archive.len() {
        let entry = archive
            .by_index(idx)
            .map_err(|err| RuntimeError::Io(format!("failed to inspect zip bytes: {err}")))?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        if name.starts_with("__MACOSX/") {
            continue;
        }
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
            || lower.ends_with(".txt")
            || lower.ends_with(".tsv")
            || lower.ends_with(".csv")
        {
            return Ok(name);
        }
        if fallback.is_none() {
            fallback = Some(name);
        }
    }
    fallback.ok_or_else(|| {
        RuntimeError::Unsupported("zip archive does not contain a supported file".to_owned())
    })
}

pub(crate) fn read_plain_sample_lines(path: &Path) -> Result<Vec<String>, RuntimeError> {
    let lower = path.to_string_lossy().to_ascii_lowercase();
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open {}: {err}", path.display())))?;
    if lower.ends_with(".vcf.gz") {
        return read_sample_lines_from_reader(BufReader::new(bgzf::io::Reader::new(file)));
    }
    read_sample_lines_from_reader(BufReader::new(file))
}

pub(crate) fn read_zip_sample_lines(
    path: &Path,
    selected_entry: &str,
) -> Result<Vec<String>, RuntimeError> {
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open zip {}: {err}", path.display())))?;
    let mut archive = ZipArchive::new(file)
        .map_err(|err| RuntimeError::Io(format!("failed to read zip {}: {err}", path.display())))?;
    let mut entry = archive.by_name(selected_entry).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open zip entry {selected_entry} in {}: {err}",
            path.display()
        ))
    })?;

    if selected_entry.to_ascii_lowercase().ends_with(".vcf.gz") {
        let bytes = read_entry_limited(
            &mut entry,
            MAX_ZIP_SAMPLE_ENTRY_BYTES,
            &format!(
                "compressed zip entry {selected_entry} in {}",
                path.display()
            ),
        )?;
        let reader = bgzf::io::Reader::new(Cursor::new(bytes));
        return read_sample_lines_from_reader(BufReader::new(reader));
    }

    read_sample_lines_from_reader(BufReader::new(entry))
}

pub(crate) fn read_entry_limited<R: Read>(
    reader: &mut R,
    max_bytes: u64,
    label: &str,
) -> Result<Vec<u8>, RuntimeError> {
    let mut bytes = Vec::new();
    reader
        .take(max_bytes.saturating_add(1))
        .read_to_end(&mut bytes)
        .map_err(|err| RuntimeError::Io(format!("failed to read {label}: {err}")))?;
    if u64::try_from(bytes.len()).unwrap_or(u64::MAX) > max_bytes {
        return Err(RuntimeError::InvalidArguments(format!(
            "{label} exceeds decompressed limit of {max_bytes} bytes"
        )));
    }
    Ok(bytes)
}

fn read_sample_lines_from_reader<R: BufRead>(mut reader: R) -> Result<Vec<String>, RuntimeError> {
    let mut out = Vec::new();
    let mut buf = String::new();
    for _ in 0..64 {
        buf.clear();
        let bytes = reader
            .read_line(&mut buf)
            .map_err(|err| RuntimeError::Io(format!("failed to read sample lines: {err}")))?;
        if bytes == 0 {
            break;
        }
        out.push(buf.trim_end_matches(['\n', '\r']).to_owned());
    }
    Ok(out)
}

pub(crate) fn select_zip_entry(path: &Path) -> Result<String, RuntimeError> {
    let file = File::open(path)
        .map_err(|err| RuntimeError::Io(format!("failed to open zip {}: {err}", path.display())))?;
    let mut archive = ZipArchive::new(file)
        .map_err(|err| RuntimeError::Io(format!("failed to read zip {}: {err}", path.display())))?;
    let mut fallback = None;
    for idx in 0..archive.len() {
        let entry = archive.by_index(idx).map_err(|err| {
            RuntimeError::Io(format!("failed to inspect zip {}: {err}", path.display()))
        })?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        if name.starts_with("__MACOSX/") {
            continue;
        }
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".vcf")
            || lower.ends_with(".vcf.gz")
            || lower.ends_with(".txt")
            || lower.ends_with(".tsv")
            || lower.ends_with(".csv")
        {
            return Ok(name);
        }
        if fallback.is_none() {
            fallback = Some(name);
        }
    }
    fallback.ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "zip archive {} does not contain a supported file",
            path.display()
        ))
    })
}
