use std::collections::{BTreeMap, HashMap};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use crate::{LibError, LibResult};
use flate2::read::MultiGzDecoder;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KmerCountMap {
    kmer_size: usize,
    counts: HashMap<String, u32>,
    transitions: HashMap<(String, String), u32>,
}

impl KmerCountMap {
    pub fn from_sequences<'a>(
        sequences: impl IntoIterator<Item = &'a str>,
        kmer_size: usize,
    ) -> LibResult<Self> {
        validate_kmer_size(kmer_size)?;
        let mut counts = HashMap::new();
        let mut transitions = HashMap::new();
        for sequence in sequences {
            count_into(&mut counts, &mut transitions, sequence, kmer_size)?;
        }
        Ok(Self {
            kmer_size,
            counts,
            transitions,
        })
    }

    pub fn from_fastq_paths<'a>(
        paths: impl IntoIterator<Item = &'a Path>,
        kmer_size: usize,
    ) -> LibResult<Self> {
        validate_kmer_size(kmer_size)?;
        let mut counts = HashMap::new();
        let mut transitions = HashMap::new();
        for path in paths {
            count_fastq_path_into(&mut counts, &mut transitions, path, kmer_size)?;
        }
        Ok(Self {
            kmer_size,
            counts,
            transitions,
        })
    }

    pub fn kmer_size(&self) -> usize {
        self.kmer_size
    }

    pub fn get(&self, kmer: &str) -> LibResult<u32> {
        validate_kmer_size(self.kmer_size)?;
        let normalized = normalize_kmer(kmer, self.kmer_size)?;
        Ok(*self.counts.get(&normalized).unwrap_or(&0))
    }

    pub fn counts(&self) -> &HashMap<String, u32> {
        &self.counts
    }

    pub fn has_transition_counts(&self) -> bool {
        !self.transitions.is_empty()
    }

    pub fn transition_count(&self, from: &str, to: &str) -> LibResult<u32> {
        let from = normalize_kmer(from, self.kmer_size)?;
        let to = normalize_kmer(to, self.kmer_size)?;
        Ok(*self.transitions.get(&(from, to)).unwrap_or(&0))
    }

    pub fn reference_counts(
        &self,
        sequence: &str,
        count_reverse_kmers: bool,
    ) -> LibResult<Vec<u32>> {
        validate_kmer_size(self.kmer_size)?;
        let bases = normalize_sequence(sequence)?;
        if bases.len() < self.kmer_size {
            return Ok(Vec::new());
        }
        let mut counts = Vec::with_capacity(bases.len() - self.kmer_size + 1);
        for window in bases.windows(self.kmer_size) {
            if window.iter().any(|base| *base == b'N') {
                counts.push(0);
                continue;
            }
            let kmer = String::from_utf8(window.to_vec()).map_err(|err| {
                LibError::InvalidArguments(format!("Kestrel k-mer is not valid UTF-8: {err}"))
            })?;
            let mut count = *self.counts.get(&kmer).unwrap_or(&0);
            if count_reverse_kmers {
                let revcomp = reverse_complement(window);
                count += *self.counts.get(&revcomp).unwrap_or(&0);
            }
            counts.push(count);
        }
        Ok(counts)
    }
}

pub fn count_sequence_kmers(sequence: &str, kmer_size: usize) -> LibResult<BTreeMap<String, u32>> {
    Ok(KmerCountMap::from_sequences([sequence], kmer_size)?
        .counts
        .into_iter()
        .collect())
}

pub fn count_fastq_kmers(path: &Path, kmer_size: usize) -> LibResult<BTreeMap<String, u32>> {
    Ok(KmerCountMap::from_fastq_paths([path], kmer_size)?
        .counts
        .into_iter()
        .collect())
}

fn count_fastq_path_into(
    counts: &mut HashMap<String, u32>,
    transitions: &mut HashMap<(String, String), u32>,
    path: &Path,
    kmer_size: usize,
) -> LibResult<()> {
    let mut reader = open_fastq_reader(path)?;
    let mut header = String::new();
    let mut sequence = String::new();
    let mut separator = String::new();
    let mut quality = String::new();
    let mut record_number = 0usize;

    loop {
        header.clear();
        if reader.read_line(&mut header).map_err(|err| {
            LibError::InvalidArguments(format!("failed to read FASTQ header: {err}"))
        })? == 0
        {
            break;
        }
        record_number += 1;
        sequence.clear();
        separator.clear();
        quality.clear();
        read_required_fastq_line(&mut reader, &mut sequence, path, record_number, "sequence")?;
        read_required_fastq_line(
            &mut reader,
            &mut separator,
            path,
            record_number,
            "separator",
        )?;
        read_required_fastq_line(&mut reader, &mut quality, path, record_number, "quality")?;

        if !header.starts_with('@') {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ record {record_number} in {} does not start with @",
                path.display()
            )));
        }
        if !separator.starts_with('+') {
            return Err(LibError::InvalidArguments(format!(
                "FASTQ record {record_number} in {} has no + separator",
                path.display()
            )));
        }
        count_into(counts, transitions, sequence.trim_end(), kmer_size)?;
    }
    Ok(())
}

fn count_into(
    counts: &mut HashMap<String, u32>,
    transitions: &mut HashMap<(String, String), u32>,
    sequence: &str,
    kmer_size: usize,
) -> LibResult<()> {
    let bases = normalize_sequence(sequence)?;
    if bases.len() < kmer_size {
        return Ok(());
    }

    let mut previous_kmer: Option<String> = None;
    for window in bases.windows(kmer_size) {
        if window.iter().any(|base| *base == b'N') {
            previous_kmer = None;
            continue;
        }
        let current_kmer = String::from_utf8(window.to_vec()).map_err(|err| {
            LibError::InvalidArguments(format!("Kestrel k-mer is not valid UTF-8: {err}"))
        })?;
        *counts.entry(current_kmer.clone()).or_insert(0) += 1;
        if let Some(previous) = previous_kmer.replace(current_kmer.clone()) {
            *transitions.entry((previous, current_kmer)).or_insert(0) += 1;
        }
    }
    Ok(())
}

fn open_fastq_reader(path: &Path) -> LibResult<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|err| {
        LibError::InvalidArguments(format!("failed to open FASTQ {}: {err}", path.display()))
    })?;
    if path.extension().is_some_and(|extension| extension == "gz") {
        return Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))));
    }
    Ok(Box::new(BufReader::new(file)))
}

fn read_required_fastq_line(
    reader: &mut dyn BufRead,
    buffer: &mut String,
    path: &Path,
    record_number: usize,
    field: &str,
) -> LibResult<()> {
    if reader
        .read_line(buffer)
        .map_err(|err| LibError::InvalidArguments(format!("failed to read FASTQ {field}: {err}")))?
        == 0
    {
        return Err(LibError::InvalidArguments(format!(
            "FASTQ record {record_number} in {} is missing {field}",
            path.display()
        )));
    }
    Ok(())
}

fn normalize_kmer(kmer: &str, kmer_size: usize) -> LibResult<String> {
    let bases = normalize_sequence(kmer)?;
    if bases.len() != kmer_size {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel k-mer length must be {kmer_size}: {kmer:?}"
        )));
    }
    if bases.iter().any(|base| *base == b'N') {
        return Err(LibError::InvalidArguments(
            "Kestrel k-mer cannot contain ambiguous bases".to_owned(),
        ));
    }
    String::from_utf8(bases).map_err(|err| {
        LibError::InvalidArguments(format!("Kestrel k-mer is not valid UTF-8: {err}"))
    })
}

fn normalize_sequence(sequence: &str) -> LibResult<Vec<u8>> {
    let mut bases = Vec::with_capacity(sequence.len());
    for base in sequence.bytes() {
        let normalized = match base {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' | b'U' | b'u' => b'T',
            b'N' | b'n' | b'R' | b'r' | b'Y' | b'y' | b'S' | b's' | b'W' | b'w' | b'K' | b'k'
            | b'M' | b'm' | b'B' | b'b' | b'D' | b'd' | b'H' | b'h' | b'V' | b'v' | b'.' | b'-' => {
                b'N'
            }
            b'\n' | b'\r' | b'\t' | b' ' => continue,
            _ => {
                return Err(LibError::InvalidArguments(format!(
                    "Kestrel sequence contains unsupported base: {}",
                    char::from(base)
                )));
            }
        };
        bases.push(normalized);
    }
    Ok(bases)
}

fn reverse_complement(kmer: &[u8]) -> String {
    kmer.iter()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => 'N',
        })
        .collect()
}

fn validate_kmer_size(kmer_size: usize) -> LibResult<()> {
    if kmer_size == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel k-mer size must be greater than zero".to_owned(),
        ));
    }
    Ok(())
}
