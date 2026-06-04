use std::path::{Path, PathBuf};

use htslib_rs::{alignment_compat, core::Region};

use super::AlignedSegment;
use crate::{LibError, LibResult};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    Read,
    ReadCram,
    ReadBam,
}

impl AlignmentMode {
    pub fn parse(mode: &str) -> LibResult<Self> {
        match mode {
            "r" | "rb" => Ok(Self::ReadBam),
            "rc" => Ok(Self::ReadCram),
            "" => Ok(Self::Read),
            other if other.contains('w') || other.contains('a') => Err(LibError::UnsupportedMode {
                module: super::MODULE,
                object: "AlignmentFile",
                mode: other.to_owned(),
            }),
            other => Err(LibError::InvalidArguments(format!(
                "pysam.AlignmentFile mode {other:?} is not recognized"
            ))),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentFile {
    path: PathBuf,
    mode: AlignmentMode,
    reference_filename: Option<PathBuf>,
    index_filename: Option<PathBuf>,
}

impl AlignmentFile {
    pub fn open(
        path: impl Into<PathBuf>,
        mode: &str,
        reference_filename: Option<PathBuf>,
        index_filename: Option<PathBuf>,
    ) -> LibResult<Self> {
        let path = path.into();
        if is_remote_path(&path) {
            return Err(LibError::unsupported_feature(
                super::MODULE,
                "remote alignment files",
            ));
        }
        Ok(Self {
            path,
            mode: AlignmentMode::parse(mode)?,
            reference_filename,
            index_filename,
        })
    }

    pub fn fetch(
        &self,
        contig: &str,
        start: Option<u64>,
        stop: Option<u64>,
    ) -> LibResult<AlignmentFetch> {
        if contig.trim().is_empty() {
            return Err(LibError::InvalidArguments(
                "pysam.AlignmentFile.fetch requires a contig".to_owned(),
            ));
        }
        if matches!((start, stop), (Some(start), Some(stop)) if stop < start) {
            return Err(LibError::InvalidArguments(
                "pysam.AlignmentFile.fetch stop must be >= start".to_owned(),
            ));
        }
        let (Some(start), Some(stop)) = (start, stop) else {
            return Err(LibError::unsupported_feature(
                super::MODULE,
                "AlignmentFile.fetch without explicit start and stop",
            ));
        };
        let region = fetch_region(contig, start, stop)?;
        let records = match self.mode {
            AlignmentMode::ReadBam | AlignmentMode::Read => {
                alignment_compat::query_bam_records_from_path(&self.path, &region)
                    .map_err(|err| LibError::InvalidArguments(err.to_string()))?
                    .into_iter()
                    .map(|record| AlignedSegment::from_hts_record(contig, &record))
                    .collect::<LibResult<Vec<_>>>()
            }
            AlignmentMode::ReadCram => {
                let Some(reference_file) = self.reference_filename.as_ref() else {
                    return Err(LibError::InvalidArguments(
                        "pysam.AlignmentFile.fetch for CRAM requires reference_filename".to_owned(),
                    ));
                };
                alignment_compat::query_cram_records_from_path_with_reference(
                    &self.path,
                    &region,
                    reference_file,
                )
                .map_err(|err| LibError::InvalidArguments(err.to_string()))?
                .into_iter()
                .map(|record| AlignedSegment::from_hts_record(contig, &record))
                .collect::<LibResult<Vec<_>>>()
            }
        }?;
        Ok(AlignmentFetch {
            contig: contig.to_owned(),
            start: Some(start),
            stop: Some(stop),
            records,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn mode(&self) -> AlignmentMode {
        self.mode
    }

    pub fn reference_filename(&self) -> Option<&Path> {
        self.reference_filename.as_deref()
    }

    pub fn index_filename(&self) -> Option<&Path> {
        self.index_filename.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentFetch {
    pub contig: String,
    pub start: Option<u64>,
    pub stop: Option<u64>,
    pub records: Vec<AlignedSegment>,
}

fn is_remote_path(path: &Path) -> bool {
    let text = path.to_string_lossy();
    text.starts_with("http://") || text.starts_with("https://") || text.starts_with("s3://")
}

fn fetch_region(contig: &str, start: u64, stop: u64) -> LibResult<Region> {
    let one_based_start = start.saturating_add(1);
    format!("{contig}:{one_based_start}-{stop}")
        .parse()
        .map_err(|err| {
            LibError::InvalidArguments(format!("pysam.AlignmentFile.fetch region: {err}"))
        })
}
