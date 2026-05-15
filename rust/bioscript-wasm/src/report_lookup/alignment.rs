use super::*;

/// Per-variant CRAM lookup that satisfies the reporting lookup
/// trait. Holds the IndexedReader in a `RefCell` so &self lookup methods can
/// mutably read while still being object-safe.
pub(crate) struct CramReportLookup<R: std::io::Read + std::io::Seek> {
    pub(crate) reader: std::cell::RefCell<noodles::cram::io::indexed_reader::IndexedReader<R>>,
    pub(crate) label: String,
}

#[path = "alignment/bam.rs"]
mod bam;

pub(crate) use bam::BamReportLookup;

impl<R: std::io::Read + std::io::Seek> bioscript_reporting::ReportVariantLookup
    for CramReportLookup<R>
{
    fn lookup_variants(&self, specs: &[VariantSpec]) -> Result<Vec<VariantObservation>, String> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(
                observe_cram_variant(&mut reader, &self.label, spec)
                    .map_err(|err| err.to_string())?,
            );
        }
        Ok(out)
    }
}

fn observe_cram_variant<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::cram::io::indexed_reader::IndexedReader<R>,
    label: &str,
    variant: &VariantSpec,
) -> Result<VariantObservation, RuntimeError> {
    let assembly = variant
        .grch38
        .as_ref()
        .map(|_| Assembly::Grch38)
        .or_else(|| variant.grch37.as_ref().map(|_| Assembly::Grch37));
    let locus = variant
        .grch38
        .as_ref()
        .or(variant.grch37.as_ref())
        .ok_or_else(|| {
            RuntimeError::Io(format!(
                "variant {} has no GRCh37/GRCh38 locus",
                variant
                    .rsids
                    .first()
                    .map(|s| s.as_str())
                    .unwrap_or("variant")
            ))
        })?;
    let locus = GenomicLocus {
        chrom: locus.chrom.clone(),
        start: locus.start,
        end: locus.end,
    };
    let kind = variant.kind.unwrap_or(VariantKind::Snp);
    match kind {
        VariantKind::Snp => {
            let ref_char = variant
                .reference
                .as_deref()
                .and_then(|s| s.chars().next())
                .ok_or_else(|| {
                    RuntimeError::Io(format!(
                        "variant {} missing reference allele",
                        variant
                            .rsids
                            .first()
                            .map(|s| s.as_str())
                            .unwrap_or("variant")
                    ))
                })?;
            let alt_char = variant
                .alternate
                .as_deref()
                .and_then(|s| s.chars().next())
                .ok_or_else(|| {
                    RuntimeError::Io(format!(
                        "variant {} missing alternate allele",
                        variant
                            .rsids
                            .first()
                            .map(|s| s.as_str())
                            .unwrap_or("variant")
                    ))
                })?;
            bioscript_formats::observe_cram_snp_with_reader(
                reader,
                label,
                &locus,
                ref_char,
                alt_char,
                variant.rsids.first().cloned(),
                assembly,
            )
        }
        VariantKind::Deletion => {
            let deletion_length = variant.deletion_length.ok_or_else(|| {
                RuntimeError::Io(format!(
                    "variant {} missing deletion_length",
                    variant
                        .rsids
                        .first()
                        .map(|s| s.as_str())
                        .unwrap_or("variant")
                ))
            })?;
            let reference = variant.reference.as_deref().unwrap_or("I");
            let alternate = variant.alternate.as_deref().unwrap_or("D");
            bioscript_formats::observe_cram_deletion_with_reader(
                reader,
                label,
                &locus,
                deletion_length,
                reference,
                alternate,
                variant.rsids.first().cloned(),
                assembly,
            )
        }
        VariantKind::Insertion | VariantKind::Indel => {
            let reference = variant.reference.as_deref().ok_or_else(|| {
                RuntimeError::Io(format!(
                    "variant {} missing reference allele",
                    variant
                        .rsids
                        .first()
                        .map(|s| s.as_str())
                        .unwrap_or("variant")
                ))
            })?;
            let alternate = variant.alternate.as_deref().ok_or_else(|| {
                RuntimeError::Io(format!(
                    "variant {} missing alternate allele",
                    variant
                        .rsids
                        .first()
                        .map(|s| s.as_str())
                        .unwrap_or("variant")
                ))
            })?;
            bioscript_formats::observe_cram_indel_with_reader(
                reader,
                label,
                &locus,
                reference,
                alternate,
                variant.rsids.first().cloned(),
                assembly,
            )
        }
        other => Err(RuntimeError::Io(format!(
            "variant {} kind {:?} not supported on CRAM via wasm",
            variant
                .rsids
                .first()
                .map(|s| s.as_str())
                .unwrap_or("variant"),
            other
        ))),
    }
}
