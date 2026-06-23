use std::{
    fs::File,
    io::{BufRead, BufReader, Cursor, Read},
    path::Path,
};

use flate2::read::MultiGzDecoder;
use zip::ZipArchive;

use bioscript_core::RuntimeError;

use super::{
    io::{detect_source_format, looks_like_vcf_lines, read_lines_from_reader, select_zip_entry},
    loaders,
    types::{
        AlignmentBytesBackend, BamBackend, BcfBackend, BcfSource, CramBackend, DelimitedBackend,
        GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore, QueryBackend, VcfBackend,
    },
};

impl GenotypeStore {
    pub fn from_file(path: &Path) -> Result<Self, RuntimeError> {
        Self::from_file_with_options(path, &GenotypeLoadOptions::default())
    }

    pub fn from_file_with_options(
        path: &Path,
        options: &GenotypeLoadOptions,
    ) -> Result<Self, RuntimeError> {
        match detect_source_format(path, options.format)? {
            GenotypeSourceFormat::Text => Ok(Self::from_delimited_file(
                path,
                GenotypeSourceFormat::Text,
                None,
                options,
            )),
            GenotypeSourceFormat::Zip => Self::from_zip_file(path, options),
            GenotypeSourceFormat::Vcf => Ok(Self::from_vcf_file(path, options)),
            GenotypeSourceFormat::Bcf => {
                if path
                    .to_string_lossy()
                    .to_ascii_lowercase()
                    .ends_with(".zip")
                {
                    Self::from_zip_file(path, options)
                } else {
                    Ok(Self::from_bcf_file(path, options))
                }
            }
            GenotypeSourceFormat::Cram => Self::from_cram_file(path, options),
            GenotypeSourceFormat::Bam => Ok(Self::from_bam_file(path, options)),
        }
    }

    pub fn from_bytes(name: &str, bytes: &[u8]) -> Result<Self, RuntimeError> {
        Self::from_bytes_with_options(name, bytes, &GenotypeLoadOptions::default())
    }

    pub fn from_bytes_with_options(
        name: &str,
        bytes: &[u8],
        options: &GenotypeLoadOptions,
    ) -> Result<Self, RuntimeError> {
        if let Some(format) = options.format {
            return match format {
                GenotypeSourceFormat::Text => {
                    let reader = BufReader::new(Cursor::new(bytes));
                    Self::from_delimited_reader(GenotypeSourceFormat::Text, reader, name)
                }
                GenotypeSourceFormat::Zip => Self::from_zip_bytes(name, bytes),
                GenotypeSourceFormat::Vcf => Self::from_vcf_bytes(name, bytes),
                GenotypeSourceFormat::Bcf => Ok(Self::from_bcf_bytes(name, bytes, options)),
                GenotypeSourceFormat::Cram | GenotypeSourceFormat::Bam => {
                    Err(RuntimeError::Unsupported(format!(
                        "{format:?} input requires alignment byte loading"
                    )))
                }
            };
        }

        // The report pipeline hands us a fixed virtual path (`/input/genotypes`)
        // with no extension, so we cannot rely on `name` alone for format
        // detection the way `from_file_with_options` can. Sniff the leading
        // bytes so a zip/VCF payload is recognised regardless of the name.
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".zip") || bytes_look_like_zip(bytes) {
            return Self::from_zip_bytes(name, bytes);
        }
        if lower.ends_with(".bcf") {
            return Ok(Self::from_bcf_bytes(
                name,
                bytes,
                &GenotypeLoadOptions::default(),
            ));
        }
        if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") || bytes_look_like_vcf(bytes) {
            return Self::from_vcf_bytes(name, bytes);
        }
        if bytes_look_like_gzip(bytes) {
            if gzip_bytes_look_like_vcf(bytes)? {
                return Self::from_vcf_bytes(name, bytes);
            }
            return Self::from_delimited_reader(
                GenotypeSourceFormat::Text,
                BufReader::new(MultiGzDecoder::new(Cursor::new(bytes))),
                name,
            );
        }
        let reader = BufReader::new(Cursor::new(bytes));
        Self::from_delimited_reader(GenotypeSourceFormat::Text, reader, name)
    }

    fn from_zip_bytes(name: &str, bytes: &[u8]) -> Result<Self, RuntimeError> {
        let mut archive = ZipArchive::new(Cursor::new(bytes)).map_err(|err| {
            RuntimeError::Io(format!("failed to read genotype zip {name}: {err}"))
        })?;
        let bcf_entries = collect_bcf_zip_entries_from_archive(&mut archive, name)?;
        if !bcf_entries.is_empty() {
            let mut entries = Vec::with_capacity(bcf_entries.len());
            for entry_name in bcf_entries {
                let mut entry = archive.by_name(&entry_name).map_err(|err| {
                    RuntimeError::Io(format!(
                        "failed to open genotype entry {entry_name} in {name}: {err}"
                    ))
                })?;
                let mut data = Vec::new();
                std::io::Read::read_to_end(&mut entry, &mut data).map_err(|err| {
                    RuntimeError::Io(format!(
                        "failed to read genotype entry {entry_name} in {name}: {err}"
                    ))
                })?;
                entries.push((entry_name, data));
            }
            return Ok(Self {
                backend: QueryBackend::Bcf(BcfBackend {
                    source: BcfSource::ZipBytes {
                        name: name.to_owned(),
                        entries,
                    },
                    options: GenotypeLoadOptions::default(),
                }),
            });
        }
        let mut selected = None;
        for idx in 0..archive.len() {
            let entry = archive.by_index(idx).map_err(|err| {
                RuntimeError::Io(format!("failed to inspect genotype zip {name}: {err}"))
            })?;
            if entry.is_dir() {
                continue;
            }
            let entry_name = entry.name().to_owned();
            let lower = entry_name.to_ascii_lowercase();
            if lower.ends_with(".vcf")
                || lower.ends_with(".vcf.gz")
                || lower.ends_with(".bcf")
                || lower.ends_with(".txt")
                || lower.ends_with(".txt.gz")
                || lower.ends_with(".txt.bgz")
                || lower.ends_with(".tsv")
                || lower.ends_with(".tsv.gz")
                || lower.ends_with(".tsv.bgz")
                || lower.ends_with(".csv")
                || lower.ends_with(".csv.gz")
                || lower.ends_with(".csv.bgz")
            {
                selected = Some(entry_name);
                break;
            }
        }
        let selected = selected.ok_or_else(|| {
            RuntimeError::Unsupported(format!(
                "zip archive {name} does not contain a supported genotype file"
            ))
        })?;
        let entry = archive.by_name(&selected).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open genotype entry {selected} in {name}: {err}"
            ))
        })?;
        let label = format!("genotype entry {selected} in {name}");
        // Stream-decompress directly off the zip reader so we never have to
        // materialize the entire decompressed entry in memory. GenesForGood
        // exports decompress to >128MB which used to trip the old
        // `read_zip_entry_limited` cap; the cap is gone because the streaming
        // parser keeps memory bounded to the rsid map itself.
        if selected.to_ascii_lowercase().ends_with(".vcf.gz") {
            return Self::from_vcf_zip_entry(entry, &label);
        }
        if is_gzip_text_name(&selected.to_ascii_lowercase()) {
            return Self::from_delimited_reader(
                GenotypeSourceFormat::Zip,
                BufReader::new(MultiGzDecoder::new(entry)),
                &label,
            );
        }
        let reader = BufReader::new(entry);
        if selected.to_ascii_lowercase().ends_with(".vcf") {
            return Self::from_vcf_reader(reader, &label);
        }
        if selected.to_ascii_lowercase().ends_with(".bcf") {
            let mut entry = reader.into_inner();
            let mut data = Vec::new();
            std::io::Read::read_to_end(&mut entry, &mut data).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to read genotype entry {selected} in {name}: {err}"
                ))
            })?;
            return Ok(Self::from_bcf_bytes(
                &label,
                &data,
                &GenotypeLoadOptions::default(),
            ));
        }
        Self::from_delimited_reader(GenotypeSourceFormat::Zip, reader, &label)
    }

    fn from_vcf_file(path: &Path, options: &GenotypeLoadOptions) -> Self {
        Self {
            backend: QueryBackend::Vcf(VcfBackend {
                path: path.to_path_buf(),
                options: options.clone(),
            }),
        }
    }

    fn from_bcf_file(path: &Path, options: &GenotypeLoadOptions) -> Self {
        Self {
            backend: QueryBackend::Bcf(BcfBackend {
                source: BcfSource::File(path.to_path_buf()),
                options: options.clone(),
            }),
        }
    }

    fn from_bcf_bytes(name: &str, bytes: &[u8], options: &GenotypeLoadOptions) -> Self {
        Self {
            backend: QueryBackend::Bcf(BcfBackend {
                source: BcfSource::Bytes {
                    name: name.to_owned(),
                    data: bytes.to_vec(),
                },
                options: options.clone(),
            }),
        }
    }

    fn from_vcf_bytes(name: &str, bytes: &[u8]) -> Result<Self, RuntimeError> {
        let lower = name.to_ascii_lowercase();
        if lower.ends_with(".vcf.gz") || bytes_look_like_gzip(bytes) {
            return Self::from_vcf_reader(
                BufReader::new(MultiGzDecoder::new(Cursor::new(bytes))),
                name,
            );
        }
        Self::from_vcf_reader(BufReader::new(Cursor::new(bytes)), name)
    }

    fn from_zip_file(path: &Path, options: &GenotypeLoadOptions) -> Result<Self, RuntimeError> {
        let bcf_entries = collect_bcf_zip_entries(path)?;
        if !bcf_entries.is_empty() {
            return Ok(Self {
                backend: QueryBackend::Bcf(BcfBackend {
                    source: BcfSource::ZipFile {
                        path: path.to_path_buf(),
                        entries: bcf_entries,
                    },
                    options: options.clone(),
                }),
            });
        }
        let selected = select_zip_entry(path)?;
        let lower = selected.to_ascii_lowercase();
        if lower.ends_with(".vcf.gz") {
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
            let entry = archive.by_name(&selected).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype entry {selected} in {}: {err}",
                    path.display()
                ))
            })?;
            let label = format!("genotype entry {selected} in {}", path.display());
            return Self::from_vcf_zip_entry(entry, &label);
        }
        if lower.ends_with(".vcf") {
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
            let entry = archive.by_name(&selected).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to open genotype entry {selected} in {}: {err}",
                    path.display()
                ))
            })?;
            let lines = read_lines_from_reader(BufReader::new(entry), path)?;
            return Self::from_vcf_lines(lines);
        }
        Ok(Self::from_delimited_file(
            path,
            GenotypeSourceFormat::Zip,
            Some(selected),
            options,
        ))
    }

    fn from_cram_file(path: &Path, options: &GenotypeLoadOptions) -> Result<Self, RuntimeError> {
        Ok(Self {
            backend: QueryBackend::Cram(CramBackend {
                path: path.to_path_buf(),
                options: options.clone(),
            }),
        })
    }

    fn from_bam_file(path: &Path, options: &GenotypeLoadOptions) -> Self {
        Self {
            backend: QueryBackend::Bam(BamBackend {
                path: path.to_path_buf(),
                options: options.clone(),
            }),
        }
    }

    /// Build an in-memory CRAM/BAM store. `kind` is `Cram` or `Bam`; `index`
    /// is the `.crai`/`.bai` bytes; `reference`/`reference_index` are the
    /// FASTA + `.fai` bytes (CRAM only -- pass empty for BAM). Used by the
    /// report pipeline, which virtualizes the genotype input.
    #[must_use]
    pub fn from_alignment_bytes(
        kind: GenotypeSourceFormat,
        data: Vec<u8>,
        index: Vec<u8>,
        reference: Vec<u8>,
        reference_index: Vec<u8>,
        options: &GenotypeLoadOptions,
    ) -> Self {
        Self {
            backend: QueryBackend::AlignmentBytes(AlignmentBytesBackend {
                kind,
                data,
                index,
                reference,
                reference_index,
                options: options.clone(),
            }),
        }
    }

    fn from_vcf_reader<R: BufRead>(reader: R, label: &str) -> Result<Self, RuntimeError> {
        loaders::from_vcf_reader(reader, label)
    }

    fn from_vcf_zip_entry<R: Read>(entry: R, label: &str) -> Result<Self, RuntimeError> {
        let mut reader = BufReader::new(entry);
        let is_gzip = bytes_look_like_gzip(
            reader
                .fill_buf()
                .map_err(|err| RuntimeError::Io(format!("failed to inspect {label}: {err}")))?,
        );
        if is_gzip {
            return Self::from_vcf_reader(BufReader::new(MultiGzDecoder::new(reader)), label);
        }
        Self::from_vcf_reader(reader, label)
    }

    fn from_delimited_reader<R: BufRead>(
        format: GenotypeSourceFormat,
        reader: R,
        label: &str,
    ) -> Result<Self, RuntimeError> {
        loaders::from_delimited_reader(format, reader, label)
    }

    fn from_vcf_lines(lines: Vec<String>) -> Result<Self, RuntimeError> {
        loaders::from_vcf_lines(lines)
    }

    fn from_delimited_file(
        path: &Path,
        format: GenotypeSourceFormat,
        zip_entry_name: Option<String>,
        options: &GenotypeLoadOptions,
    ) -> Self {
        Self {
            backend: QueryBackend::Delimited(DelimitedBackend {
                format,
                path: path.to_path_buf(),
                zip_entry_name,
                options: options.clone(),
            }),
        }
    }
}

fn collect_bcf_zip_entries(path: &Path) -> Result<Vec<String>, RuntimeError> {
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
    collect_bcf_zip_entries_from_archive(&mut archive, &path.display().to_string())
}

fn collect_bcf_zip_entries_from_archive<R: std::io::Read + std::io::Seek>(
    archive: &mut ZipArchive<R>,
    label: &str,
) -> Result<Vec<String>, RuntimeError> {
    let mut entries = Vec::new();
    for idx in 0..archive.len() {
        let entry = archive.by_index(idx).map_err(|err| {
            RuntimeError::Io(format!("failed to inspect genotype zip {label}: {err}"))
        })?;
        if entry.is_dir() {
            continue;
        }
        let name = entry.name().to_owned();
        if name.to_ascii_lowercase().ends_with(".bcf") {
            entries.push(name);
        }
    }
    entries.sort();
    Ok(entries)
}

fn bytes_look_like_zip(bytes: &[u8]) -> bool {
    bytes.starts_with(b"PK\x03\x04")
        || bytes.starts_with(b"PK\x05\x06")
        || bytes.starts_with(b"PK\x07\x08")
}

fn bytes_look_like_gzip(bytes: &[u8]) -> bool {
    bytes.starts_with(&[0x1f, 0x8b])
}

fn gzip_bytes_look_like_vcf(bytes: &[u8]) -> Result<bool, RuntimeError> {
    let mut decoder = MultiGzDecoder::new(Cursor::new(bytes));
    let mut prefix = Vec::new();
    decoder
        .by_ref()
        .take(8192)
        .read_to_end(&mut prefix)
        .map_err(|err| RuntimeError::Io(format!("failed to read gzip genotype bytes: {err}")))?;
    Ok(bytes_look_like_vcf(&prefix))
}

fn bytes_look_like_vcf(bytes: &[u8]) -> bool {
    let prefix = &bytes[..bytes.len().min(8192)];
    let text = String::from_utf8_lossy(prefix);
    let lines: Vec<String> = text.lines().map(str::to_owned).collect();
    looks_like_vcf_lines(&lines)
}

fn is_gzip_text_name(lower_name: &str) -> bool {
    lower_name.ends_with(".txt.gz")
        || lower_name.ends_with(".txt.bgz")
        || lower_name.ends_with(".tsv.gz")
        || lower_name.ends_with(".tsv.bgz")
        || lower_name.ends_with(".csv.gz")
        || lower_name.ends_with(".csv.bgz")
}
