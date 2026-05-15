use super::*;

#[path = "report_lookup/alignment.rs"]
mod alignment;

pub(crate) use alignment::{BamReportLookup, CramReportLookup};

pub(super) struct VcfReportLookup<R: std::io::Read + std::io::Seek> {
    pub(super) reader: std::cell::RefCell<
        noodles::csi::io::IndexedReader<noodles::bgzf::io::Reader<R>, noodles::tabix::Index>,
    >,
    pub(super) label: String,
    /// Assembly resolved from the VCF header (via `inspect_head_via_js_reader`).
    /// Matches the CLI's `lookup_indexed_vcf_variants` flow which calls
    /// `detect_vcf_assembly_from_path`. Without this the wasm picks GRCh38
    /// over GRCh37 for any panel variant that declares both loci, then
    /// misses the variant in a GRCh37-coded VCF (NA06985.clean.vcf.gz etc.)
    /// and falls through to "imputed reference".
    pub(super) detected_assembly: Option<Assembly>,
}

impl<R: std::io::Read + std::io::Seek> bioscript_reporting::ReportVariantLookup
    for VcfReportLookup<R>
{
    fn lookup_variants(&self, specs: &[VariantSpec]) -> Result<Vec<VariantObservation>, String> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(
                observe_vcf_variant(&mut reader, &self.label, spec, self.detected_assembly)
                    .map_err(|err| err.to_string())?,
            );
        }
        Ok(out)
    }
}

fn observe_vcf_variant<R: std::io::Read + std::io::Seek>(
    reader: &mut noodles::csi::io::IndexedReader<
        noodles::bgzf::io::Reader<R>,
        noodles::tabix::Index,
    >,
    label: &str,
    variant: &VariantSpec,
    detected_assembly: Option<Assembly>,
) -> Result<VariantObservation, RuntimeError> {
    // Use the existing CLI helper so wasm picks the same locus the path-based
    // path does: detected GRCh37 → grch37 first; detected GRCh38 → grch38
    // first; None → grch37 first (CLI default for variant-only VCFs).
    let raw_locus =
        bioscript_formats::choose_variant_locus_for_assembly(variant, detected_assembly)
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
    let assembly = detected_assembly.or_else(|| {
        if variant.grch37.as_ref().is_some_and(|l| l == &raw_locus) {
            Some(Assembly::Grch37)
        } else if variant.grch38.as_ref().is_some_and(|l| l == &raw_locus) {
            Some(Assembly::Grch38)
        } else {
            None
        }
    });
    let locus = GenomicLocus {
        chrom: raw_locus.chrom.clone(),
        start: raw_locus.start,
        end: raw_locus.end,
    };
    let observation = bioscript_formats::observe_vcf_variant_with_reader(
        reader,
        label,
        &locus,
        variant,
        variant.rsids.first().cloned(),
        assembly,
    )?;
    // Mirror the CLI report flow's `impute_vcf_missing_as_reference: true`
    // default. The full-file VCF scanner only marks a variant resolved when a
    // row actually matches the query. Unrelated rows in the indexed window
    // must not block reference imputation for absent variant-only calls.
    if observation.genotype.is_none()
        && !observation.evidence.iter().any(|line| {
            line.contains("tabix index has no contig")
                || line.contains("has no GRCh37/GRCh38 locus")
        })
        && let Some(imputed) = bioscript_formats::imputed_reference_observation(
            "vcf",
            label,
            variant,
            &locus,
            assembly,
            None,
            &observation.evidence.join(" | "),
        )
    {
        return Ok(imputed);
    }
    Ok(observation)
}
