use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use flate2::read::MultiGzDecoder;

use crate::{LibError, LibResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceRecord {
    pub name: String,
    pub sequence: String,
    pub line: usize,
}

pub fn read_reference_records(path: &Path) -> LibResult<Vec<ReferenceRecord>> {
    match reference_format(path)? {
        ReferenceFormat::Fasta => read_fasta_records(path),
        ReferenceFormat::Fastq => read_fastq_records(path),
        ReferenceFormat::Raw => read_raw_records(path),
    }
}

pub fn reference_kmers(sequence: &str, kmer_size: usize) -> LibResult<Vec<String>> {
    if kmer_size == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel reference k-mer size must be greater than zero".to_owned(),
        ));
    }
    let normalized = normalize_reference_sequence(sequence);
    if normalized.len() < kmer_size {
        return Ok(Vec::new());
    }
    Ok((0..=normalized.len() - kmer_size)
        .map(|index| normalized[index..index + kmer_size].to_owned())
        .collect())
}

fn read_fasta_records(path: &Path) -> LibResult<Vec<ReferenceRecord>> {
    let mut reader = open_reader(path)?;
    let mut records = Vec::new();
    let mut line = String::new();
    let mut line_number = 0usize;
    let mut pending_name: Option<(String, usize)> = None;
    let mut sequence = String::new();

    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(read_error(path))? == 0 {
            break;
        }
        line_number += 1;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(name) = trimmed.strip_prefix('>') {
            push_pending_record(&mut records, &mut pending_name, &mut sequence)?;
            pending_name = Some((required_name(name, path, line_number)?, line_number));
        } else {
            if pending_name.is_none() {
                return Err(LibError::InvalidArguments(format!(
                    "FASTA {} line {line_number} is missing a > header",
                    path.display()
                )));
            }
            sequence.push_str(trimmed);
        }
    }
    push_pending_record(&mut records, &mut pending_name, &mut sequence)?;
    Ok(records)
}

fn read_fastq_records(path: &Path) -> LibResult<Vec<ReferenceRecord>> {
    let mut reader = open_reader(path)?;
    let mut records = Vec::new();
    let mut line = String::new();
    let mut line_number = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(read_error(path))? == 0 {
            break;
        }
        line_number += 1;
        let header = line.trim();
        if header.is_empty() {
            continue;
        }
        let Some(name) = header.strip_prefix('@') else {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ {} line {line_number} is missing an @ header",
                path.display()
            )));
        };
        let name = required_name(name, path, line_number)?;
        let sequence_line = read_required_line(&mut *reader, path, &mut line_number, "sequence")?;
        let separator = read_required_line(&mut *reader, path, &mut line_number, "separator")?;
        if !separator.trim().starts_with('+') {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ {} line {line_number} is missing a + separator",
                path.display()
            )));
        }
        let quality = read_required_line(&mut *reader, path, &mut line_number, "quality")?;
        let sequence = sequence_line.trim().to_owned();
        if quality.trim().len() != sequence.len() {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ {} record {name} has mismatched sequence and quality lengths",
                path.display()
            )));
        }
        records.push(ReferenceRecord {
            name,
            sequence,
            line: line_number - 3,
        });
    }
    Ok(records)
}

fn read_raw_records(path: &Path) -> LibResult<Vec<ReferenceRecord>> {
    let mut reader = open_reader(path)?;
    let mut records = Vec::new();
    let mut line = String::new();
    let mut sequence = String::new();
    let mut record_number = 0usize;
    let mut start_line = 0usize;
    let mut line_number = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(read_error(path))? == 0 {
            break;
        }
        line_number += 1;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            push_raw_record(&mut records, &mut sequence, &mut record_number, start_line)?;
            continue;
        }
        if sequence.is_empty() {
            start_line = line_number;
        }
        sequence.push_str(trimmed);
    }
    push_raw_record(&mut records, &mut sequence, &mut record_number, start_line)?;
    Ok(records)
}

fn push_pending_record(
    records: &mut Vec<ReferenceRecord>,
    pending_name: &mut Option<(String, usize)>,
    sequence: &mut String,
) -> LibResult<()> {
    if let Some((name, line)) = pending_name.take() {
        if sequence.is_empty() {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel reference record {name} has no sequence"
            )));
        }
        records.push(ReferenceRecord {
            name,
            sequence: std::mem::take(sequence),
            line,
        });
    }
    Ok(())
}

fn push_raw_record(
    records: &mut Vec<ReferenceRecord>,
    sequence: &mut String,
    record_number: &mut usize,
    line: usize,
) -> LibResult<()> {
    if sequence.is_empty() {
        return Ok(());
    }
    *record_number += 1;
    records.push(ReferenceRecord {
        name: format!("Sequence{record_number}"),
        sequence: std::mem::take(sequence),
        line,
    });
    Ok(())
}

fn required_name(name: &str, path: &Path, line_number: usize) -> LibResult<String> {
    let name = name.trim();
    if name.is_empty() {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel reference {} line {line_number} has an empty record name",
            path.display()
        )));
    }
    Ok(name.to_owned())
}

fn read_required_line(
    reader: &mut dyn BufRead,
    path: &Path,
    line_number: &mut usize,
    field: &str,
) -> LibResult<String> {
    let mut line = String::new();
    if reader.read_line(&mut line).map_err(read_error(path))? == 0 {
        return Err(LibError::InvalidArguments(format!(
            "FASTQ {} is missing {field}",
            path.display()
        )));
    }
    *line_number += 1;
    Ok(line)
}

fn normalize_reference_sequence(sequence: &str) -> String {
    let mut ambiguous_index = 0usize;
    sequence
        .chars()
        .map(|base| match base.to_ascii_uppercase() {
            'A' | 'C' | 'G' | 'T' => base.to_ascii_uppercase(),
            'U' => 'T',
            _ => {
                let base = ['A', 'C', 'G', 'T'][ambiguous_index % 4];
                ambiguous_index += 1;
                base
            }
        })
        .collect()
}

fn open_reader(path: &Path) -> LibResult<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|err| {
        LibError::InvalidArguments(format!(
            "failed to open Kestrel reference {}: {err}",
            path.display()
        ))
    })?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        return Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))));
    }
    Ok(Box::new(BufReader::new(file)))
}

fn read_error(path: &Path) -> impl Fn(std::io::Error) -> LibError + '_ {
    move |err| {
        LibError::InvalidArguments(format!(
            "failed to read Kestrel reference {}: {err}",
            path.display()
        ))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReferenceFormat {
    Fasta,
    Fastq,
    Raw,
}

fn reference_format(path: &Path) -> LibResult<ReferenceFormat> {
    let file_name = path
        .file_name()
        .and_then(|file_name| file_name.to_str())
        .ok_or_else(|| {
            LibError::InvalidArguments(format!(
                "Kestrel reference path has no valid file name: {}",
                path.display()
            ))
        })?
        .to_ascii_lowercase();
    let uncompressed = file_name.strip_suffix(".gz").unwrap_or(&file_name);
    if uncompressed.ends_with(".fasta") || uncompressed.ends_with(".fa") {
        Ok(ReferenceFormat::Fasta)
    } else if uncompressed.ends_with(".fastq") || uncompressed.ends_with(".fq") {
        Ok(ReferenceFormat::Fastq)
    } else if uncompressed.ends_with(".raw") {
        Ok(ReferenceFormat::Raw)
    } else {
        Err(LibError::InvalidArguments(format!(
            "unsupported Kestrel reference format: {}",
            path.display()
        )))
    }
}
