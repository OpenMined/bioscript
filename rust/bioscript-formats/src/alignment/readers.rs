use std::{
    io::{BufRead, Read, Seek},
    path::Path,
};

use noodles::{
    cram::{self, crai},
    fasta::{self, repository::adapters::IndexedReader as FastaIndexedReader},
    tabix,
};

use bioscript_core::RuntimeError;

use crate::genotype::GenotypeLoadOptions;

/// Build a CRAM `IndexedReader` over any `Read + Seek` source given a parsed
/// CRAI index and a reference repository. Mirrors `build_from_path` but with
/// an externally-provided reader — the wasm path uses this with a JS-backed
/// reader; native paths still go through the path-based helper below.
pub fn build_cram_indexed_reader_from_reader<R>(
    reader: R,
    crai_index: crai::Index,
    repository: fasta::Repository,
) -> Result<cram::io::indexed_reader::IndexedReader<R>, RuntimeError>
where
    R: Read,
{
    cram::io::indexed_reader::Builder::default()
        .set_reference_sequence_repository(repository)
        .set_index(crai_index)
        .build_from_reader(reader)
        .map_err(|err| RuntimeError::Io(format!("failed to build indexed CRAM reader: {err}")))
}

/// Build a FASTA `Repository` over any `BufRead + Seek + Send + Sync` source
/// given a parsed FAI index. The `Send + Sync + 'static` bounds come from
/// `fasta::Repository`'s internal `Arc<RwLock<dyn Adapter + Send + Sync>>`
/// cache — on single-threaded wasm32 these can be met via `unsafe impl`.
pub fn build_reference_repository_from_readers<R>(
    reader: R,
    fai_index: fasta::fai::Index,
) -> fasta::Repository
where
    R: BufRead + Seek + Send + Sync + 'static,
{
    let indexed = fasta::io::IndexedReader::new(reader, fai_index);
    fasta::Repository::new(FastaIndexedReader::new(indexed))
}

/// Parse a CRAM index (`.crai`) from an in-memory byte buffer. Used by wasm
/// callers that receive the small index inline while the big CRAM stays on a
/// JS-backed reader.
pub fn parse_crai_bytes(bytes: &[u8]) -> Result<crai::Index, RuntimeError> {
    crai::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse CRAM index bytes: {err}")))
}

/// Parse a FASTA index (`.fai`) from an in-memory byte buffer.
pub fn parse_fai_bytes(bytes: &[u8]) -> Result<fasta::fai::Index, RuntimeError> {
    fasta::fai::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse FASTA index bytes: {err}")))
}

/// Parse a tabix index (`.tbi`) from an in-memory byte buffer. Used by wasm
/// callers that pass the small index inline while the bgzipped VCF stays on
/// a JS-backed `Read + Seek` reader.
pub fn parse_tbi_bytes(bytes: &[u8]) -> Result<tabix::Index, RuntimeError> {
    tabix::io::Reader::new(std::io::Cursor::new(bytes))
        .read_index()
        .map_err(|err| RuntimeError::Io(format!("failed to parse tabix index bytes: {err}")))
}

pub(crate) fn build_cram_indexed_reader_from_path(
    path: &Path,
    options: &GenotypeLoadOptions,
    repository: fasta::Repository,
) -> Result<cram::io::indexed_reader::IndexedReader<std::fs::File>, RuntimeError> {
    let mut builder =
        cram::io::indexed_reader::Builder::default().set_reference_sequence_repository(repository);

    if let Some(index_path) = options.input_index.as_ref() {
        let index = crai::fs::read(index_path).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read CRAM index {} for {}: {err}",
                index_path.display(),
                path.display()
            ))
        })?;
        builder = builder.set_index(index);
    }

    builder.build_from_path(path).map_err(|err| {
        RuntimeError::Io(format!(
            "failed to open indexed CRAM {}: {err}",
            path.display()
        ))
    })
}

pub(crate) fn build_reference_repository(
    reference_file: &Path,
) -> Result<fasta::Repository, RuntimeError> {
    let reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_file)
        .map_err(|err| {
            RuntimeError::Io(format!(
                "failed to open indexed FASTA {}: {err}",
                reference_file.display()
            ))
        })?;

    Ok(fasta::Repository::new(FastaIndexedReader::new(reader)))
}
