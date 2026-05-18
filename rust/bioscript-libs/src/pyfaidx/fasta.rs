use std::{
    fs::File,
    io::{BufReader, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use crate::{LibError, LibResult};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fasta {
    path: PathBuf,
    index: Option<htslib_rs::faidx_compat::Index>,
}

impl Fasta {
    pub fn open(path: impl Into<PathBuf>) -> Self {
        Self {
            path: path.into(),
            index: None,
        }
    }

    pub fn from_path(path: impl Into<PathBuf>) -> LibResult<Self> {
        let path = path.into();
        let file = File::open(&path).map_err(|err| {
            LibError::InvalidArguments(format!("failed to open FASTA {}: {err}", path.display()))
        })?;
        let index = htslib_rs::faidx_compat::build_index(BufReader::new(file)).map_err(|err| {
            LibError::InvalidArguments(format!("failed to index FASTA {}: {err}", path.display()))
        })?;
        Ok(Self {
            path,
            index: Some(index),
        })
    }

    pub fn get(&self, contig: &str) -> LibResult<FastaRecord> {
        if contig.trim().is_empty() {
            return Err(LibError::InvalidArguments(
                "pyfaidx.Fasta contig name cannot be empty".to_owned(),
            ));
        }
        let Some(index) = self.index.as_ref() else {
            return Err(LibError::InvalidArguments(format!(
                "pyfaidx.Fasta record {contig:?} was not loaded from {}",
                self.path.display()
            )));
        };
        let mut file = File::open(&self.path).map_err(|err| {
            LibError::InvalidArguments(format!(
                "failed to open FASTA {}: {err}",
                self.path.display()
            ))
        })?;
        file.seek(SeekFrom::Start(0)).map_err(|err| {
            LibError::InvalidArguments(format!(
                "failed to seek FASTA {}: {err}",
                self.path.display()
            ))
        })?;
        let sequence = htslib_rs::faidx_compat::fetch_region_sequence(&mut file, index, contig)
            .map_err(|err| {
                LibError::InvalidArguments(format!(
                    "pyfaidx.Fasta record {contig:?} was not found in {}: {err}",
                    self.path.display()
                ))
            })?;
        let sequence = String::from_utf8(sequence).map_err(|err| {
            LibError::InvalidArguments(format!(
                "pyfaidx.Fasta record {contig:?} in {} is not UTF-8: {err}",
                self.path.display()
            ))
        })?;
        Ok(FastaRecord {
            name: contig.to_owned(),
            sequence,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
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
