use super::*;

/// Per-variant CRAM lookup that satisfies the workspace's `VariantLookup`
/// trait. Holds the IndexedReader in a `RefCell` so &self lookup methods can
/// mutably read while still being object-safe.
pub(super) struct CramReportLookup<R: std::io::Read + std::io::Seek> {
    pub(super) reader: std::cell::RefCell<noodles::cram::io::indexed_reader::IndexedReader<R>>,
    pub(super) label: String,
}

impl<R: std::io::Read + std::io::Seek> report_workspace::VariantLookup for CramReportLookup<R> {
    fn lookup_variant(&self, spec: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        observe_cram_variant(&mut reader, &self.label, spec)
    }

    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(observe_cram_variant(&mut reader, &self.label, spec)?);
        }
        Ok(out)
    }
}

/// Build a minimal 23andMe-style text from the observations we already
/// computed. Format: `rsid\tchrom\tpos\tgenotype` per line. The runtime's
/// delimited-text loader reads this back as a `RsidMap`/`Delimited` backend
/// so analysis scripts can call `bioscript.load_genotypes(input_file)` and
/// have rsid lookups answered from the cached table.
#[allow(dead_code)]
fn synthesize_genotype_text_from_observations(observations: &[serde_json::Value]) -> String {
    let mut out = String::from("# rsid\tchromosome\tposition\tgenotype\n");
    for observation in observations {
        let rsid = observation
            .get("rsid")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("");
        if rsid.is_empty() {
            continue;
        }
        let chrom = observation
            .get("chrom")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("");
        let pos = observation
            .get("pos_start")
            .and_then(|v| {
                v.as_i64()
                    .or_else(|| v.as_str().and_then(|s| s.parse::<i64>().ok()))
            })
            .unwrap_or(0);
        let genotype = observation
            .get("genotype_display")
            .and_then(serde_json::Value::as_str)
            .filter(|s| !s.is_empty() && *s != "??")
            .unwrap_or("--");
        out.push_str(&format!("{rsid}\t{chrom}\t{pos}\t{genotype}\n"));
    }
    out
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

pub(super) struct VcfReportLookup<R: std::io::Read + std::io::Seek> {
    pub(super) reader: std::cell::RefCell<
        noodles::csi::io::IndexedReader<noodles::bgzf::io::Reader<R>, noodles::tabix::Index>,
    >,
    pub(super) label: String,
}

impl<R: std::io::Read + std::io::Seek> report_workspace::VariantLookup for VcfReportLookup<R> {
    fn lookup_variant(&self, spec: &VariantSpec) -> Result<VariantObservation, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        observe_vcf_variant(&mut reader, &self.label, spec)
    }

    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let mut reader = self.reader.borrow_mut();
        let mut out = Vec::with_capacity(specs.len());
        for spec in specs {
            out.push(observe_vcf_variant(&mut reader, &self.label, spec)?);
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
) -> Result<VariantObservation, RuntimeError> {
    let assembly = variant
        .grch38
        .as_ref()
        .map(|_| Assembly::Grch38)
        .or_else(|| variant.grch37.as_ref().map(|_| Assembly::Grch37));
    let raw_locus = variant
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
    // default: when the VCF has no record at this locus, treat the genotype
    // as homozygous reference.
    if observation.genotype.is_none()
        && observation
            .evidence
            .iter()
            .any(|line| line.contains("no VCF record at"))
    {
        if let Some(imputed) = bioscript_formats::imputed_reference_observation(
            "vcf",
            label,
            variant,
            &locus,
            assembly,
            None,
            &observation.evidence.join(" | "),
        ) {
            return Ok(imputed);
        }
    }
    Ok(observation)
}
