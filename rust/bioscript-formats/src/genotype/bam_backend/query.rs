use std::io::{Read, Seek};

use bioscript_core::{GenomicLocus, RuntimeError};

pub(super) fn read_bam_header<R: Read + Seek>(
    reader: &mut noodles::bam::io::indexed_reader::IndexedReader<noodles::bgzf::io::Reader<R>>,
    label: &str,
) -> Result<noodles::sam::Header, RuntimeError> {
    reader
        .get_mut()
        .seek(noodles::bgzf::VirtualPosition::MIN)
        .map_err(|err| RuntimeError::Io(format!("failed to rewind BAM {label}: {err}")))?;
    reader
        .read_header()
        .map_err(|err| RuntimeError::Io(format!("failed to read BAM header {label}: {err}")))
}

pub(super) fn bam_region(
    header: &noodles::sam::Header,
    locus: &GenomicLocus,
) -> Result<noodles::core::Region, RuntimeError> {
    let chrom = resolve_bam_reference_name(header, &locus.chrom).ok_or_else(|| {
        RuntimeError::Unsupported(format!(
            "indexed BAM does not contain contig {} for {}:{}-{}",
            locus.chrom, locus.chrom, locus.start, locus.end
        ))
    })?;
    format!("{chrom}:{}-{}", locus.start, locus.end)
        .parse()
        .map_err(|err| RuntimeError::Io(format!("invalid BAM query region: {err}")))
}

fn resolve_bam_reference_name(header: &noodles::sam::Header, chrom: &str) -> Option<String> {
    let candidates = [
        chrom.to_owned(),
        format!("chr{chrom}"),
        chrom.trim_start_matches("chr").to_owned(),
    ];
    candidates.into_iter().find(|candidate| {
        header.reference_sequences().iter().any(|(name, _)| {
            let name_bytes: &[u8] = name.as_ref();
            name_bytes == candidate.as_bytes()
        })
    })
}
