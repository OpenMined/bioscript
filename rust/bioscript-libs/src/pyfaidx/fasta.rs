use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use crate::{LibError, LibResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fasta {
    path: PathBuf,
    records: BTreeMap<String, String>,
}

impl Fasta {
    pub fn open(path: impl Into<PathBuf>) -> Self {
        Self {
            path: path.into(),
            records: BTreeMap::new(),
        }
    }

    pub fn from_path(path: impl Into<PathBuf>) -> LibResult<Self> {
        let path = path.into();
        let contents = fs::read_to_string(&path).map_err(|err| {
            LibError::InvalidArguments(format!("failed to read FASTA {}: {err}", path.display()))
        })?;
        let records = parse_fasta_records(&contents)?;
        Ok(Self { path, records })
    }

    pub fn get(&self, contig: &str) -> LibResult<FastaRecord> {
        if contig.trim().is_empty() {
            return Err(LibError::InvalidArguments(
                "pyfaidx.Fasta contig name cannot be empty".to_owned(),
            ));
        }
        let sequence = self.records.get(contig).ok_or_else(|| {
            LibError::InvalidArguments(format!(
                "pyfaidx.Fasta record {contig:?} was not found in {}",
                self.path.display()
            ))
        })?;
        Ok(FastaRecord {
            name: contig.to_owned(),
            sequence: sequence.clone(),
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

fn parse_fasta_records(contents: &str) -> LibResult<BTreeMap<String, String>> {
    let mut records = BTreeMap::new();
    let mut current_name: Option<String> = None;
    let mut current_sequence = String::new();

    for line in contents.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            flush_record(&mut records, &mut current_name, &mut current_sequence)?;
            let name = rest
                .split_whitespace()
                .next()
                .filter(|value| !value.is_empty())
                .ok_or_else(|| LibError::InvalidArguments("FASTA header is empty".to_owned()))?;
            current_name = Some(name.to_owned());
        } else if current_name.is_none() {
            return Err(LibError::InvalidArguments(
                "FASTA sequence appeared before first header".to_owned(),
            ));
        } else {
            current_sequence.push_str(trimmed);
        }
    }

    flush_record(&mut records, &mut current_name, &mut current_sequence)?;
    if records.is_empty() {
        return Err(LibError::InvalidArguments(
            "FASTA did not contain any records".to_owned(),
        ));
    }
    Ok(records)
}

fn flush_record(
    records: &mut BTreeMap<String, String>,
    current_name: &mut Option<String>,
    current_sequence: &mut String,
) -> LibResult<()> {
    let Some(name) = current_name.take() else {
        return Ok(());
    };
    if records.contains_key(&name) {
        return Err(LibError::InvalidArguments(format!(
            "duplicate FASTA record {name:?}"
        )));
    }
    records.insert(name, std::mem::take(current_sequence));
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub name: String,
    pub sequence: String,
}

impl FastaRecord {
    pub fn slice(&self, start: usize, stop: usize) -> LibResult<String> {
        if stop < start {
            return Err(LibError::InvalidArguments(
                "pyfaidx slice stop must be >= start".to_owned(),
            ));
        }
        self.sequence
            .get(start..stop)
            .map(str::to_owned)
            .ok_or_else(|| LibError::InvalidArguments("pyfaidx slice is out of bounds".to_owned()))
    }
}
