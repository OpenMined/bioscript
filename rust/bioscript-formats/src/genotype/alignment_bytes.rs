use std::io::Cursor;

use bioscript_core::{
    Assembly, GenomicLocus, RuntimeError, VariantKind, VariantObservation, VariantSpec,
};
use noodles::cram;

use crate::alignment;

use super::cram_backend::choose_variant_locus;
use super::types::{AlignmentBytesBackend, GenotypeSourceFormat};
use super::{
    bam_backend::observe_bam_variant, observe_cram_deletion_with_reader,
    observe_cram_indel_with_reader, observe_cram_snp_with_reader, variant_sort_key,
};

const LABEL: &str = "/input/genotypes";

impl AlignmentBytesBackend {
    pub(crate) fn backend_name(&self) -> &'static str {
        match self.kind {
            GenotypeSourceFormat::Bam => "bam",
            _ => "cram",
        }
    }

    pub(crate) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let mut results = self.lookup_variants(std::slice::from_ref(variant))?;
        Ok(results.pop().unwrap_or_default())
    }

    pub(crate) fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        match self.kind {
            GenotypeSourceFormat::Cram => self.lookup_cram(variants),
            GenotypeSourceFormat::Bam => self.lookup_bam(variants),
            other => Err(RuntimeError::Unsupported(format!(
                "alignment-bytes backend does not support {other:?}"
            ))),
        }
    }

    fn lookup_bam(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let bai = alignment::parse_bai_bytes(&self.index)?;
        let mut reader =
            alignment::build_bam_indexed_reader_from_reader(Cursor::new(self.data.clone()), bai)?;

        let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
        indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

        let mut results = vec![VariantObservation::default(); variants.len()];
        for (idx, variant) in indexed {
            results[idx] = observe_bam_variant(&mut reader, LABEL, variant)?;
        }
        Ok(results)
    }

    fn lookup_cram(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let fai = alignment::parse_fai_bytes(&self.reference_index)?;
        let repository = alignment::build_reference_repository_from_readers(
            Cursor::new(self.reference.clone()),
            fai,
        );
        let crai = alignment::parse_crai_bytes(&self.index)?;
        let mut reader = alignment::build_cram_indexed_reader_from_reader(
            Cursor::new(self.data.clone()),
            crai,
            repository,
        )?;

        let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
        indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

        let mut results = vec![VariantObservation::default(); variants.len()];
        for (idx, variant) in indexed {
            results[idx] = self.observe_cram(&mut reader, variant)?;
        }
        Ok(results)
    }

    fn observe_cram(
        &self,
        reader: &mut cram::io::indexed_reader::IndexedReader<Cursor<Vec<u8>>>,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let Some((assembly, locus)) = self.choose_locus(variant) else {
            return Ok(self.unsupported_locus_observation(variant));
        };
        let matched_rsid = variant.rsids.first().cloned();
        match variant.kind.unwrap_or(VariantKind::Other) {
            VariantKind::Snp => {
                let reference = first_base(variant.reference.as_deref()).ok_or_else(|| {
                    RuntimeError::InvalidArguments("SNP variant requires ref/reference".to_owned())
                })?;
                let alternate = first_base(variant.alternate.as_deref()).ok_or_else(|| {
                    RuntimeError::InvalidArguments("SNP variant requires alt/alternate".to_owned())
                })?;
                observe_cram_snp_with_reader(
                    reader,
                    LABEL,
                    &locus,
                    reference,
                    alternate,
                    matched_rsid,
                    Some(assembly),
                )
            }
            VariantKind::Deletion => {
                let deletion_length = variant.deletion_length.ok_or_else(|| {
                    RuntimeError::InvalidArguments(
                        "deletion variant requires deletion_length".to_owned(),
                    )
                })?;
                let reference = variant.reference.clone().unwrap_or_else(|| "I".to_owned());
                let alternate = variant.alternate.clone().unwrap_or_else(|| "D".to_owned());
                observe_cram_deletion_with_reader(
                    reader,
                    LABEL,
                    &locus,
                    deletion_length,
                    &reference,
                    &alternate,
                    matched_rsid,
                    Some(assembly),
                )
            }
            VariantKind::Insertion | VariantKind::Indel => {
                let reference = variant.reference.clone().unwrap_or_else(|| "I".to_owned());
                let alternate = variant.alternate.clone().unwrap_or_else(|| "D".to_owned());
                let alternate_lengths = indel_alternate_lengths(variant, &alternate);
                observe_cram_indel_with_reader(
                    reader,
                    LABEL,
                    &locus,
                    &reference,
                    &alternate,
                    &alternate_lengths,
                    matched_rsid,
                    Some(assembly),
                )
            }
            VariantKind::Other => Err(RuntimeError::Unsupported(format!(
                "backend '{}' does not yet support {:?} observation",
                self.backend_name(),
                variant.kind.unwrap_or(VariantKind::Other)
            ))),
        }
    }

    /// Mirror `CramBackend`'s reference-filename assembly heuristic: prefer an
    /// explicit `options.assembly`, otherwise fall back to the same
    /// GRCh38-first default `choose_variant_locus` uses when the reference
    /// filename carries no assembly hint.
    fn choose_locus(&self, variant: &VariantSpec) -> Option<(Assembly, GenomicLocus)> {
        if let Some(assembly) = self.options.assembly {
            let locus = match assembly {
                Assembly::Grch38 => variant.grch38.clone().or_else(|| variant.grch37.clone()),
                Assembly::Grch37 => variant.grch37.clone().or_else(|| variant.grch38.clone()),
            }?;
            return Some((assembly, locus));
        }
        choose_variant_locus(variant, std::path::Path::new(""))
    }

    fn unsupported_locus_observation(&self, variant: &VariantSpec) -> VariantObservation {
        VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            evidence: vec![format!(
                "backend '{}' cannot satisfy query for {} without GRCh37/GRCh38 coordinates",
                self.backend_name(),
                super::describe_query(variant)
            )],
            ..VariantObservation::default()
        }
    }
}

fn indel_alternate_lengths(variant: &VariantSpec, fallback_alternate: &str) -> Vec<usize> {
    let mut lengths = variant
        .observed_alternates
        .iter()
        .map(String::len)
        .filter(|len| *len > 0)
        .collect::<Vec<_>>();
    if lengths.is_empty() {
        lengths.push(fallback_alternate.len());
    }
    lengths.sort_unstable();
    lengths.dedup();
    lengths
}

fn first_base(value: Option<&str>) -> Option<char> {
    value.and_then(|s| s.chars().next())
}
