use bioscript_core::{RuntimeError, VariantObservation, VariantSpec};

use super::{
    lookup_indexed_vcf_variants, scan_delimited_variants, scan_vcf_variants,
    types::{DelimitedBackend, GenotypeSourceFormat, RsidMapBackend, VcfBackend},
};

impl RsidMapBackend {
    pub(super) fn backend_name(&self) -> &'static str {
        match self.format {
            GenotypeSourceFormat::Text => "text",
            GenotypeSourceFormat::Zip => "zip",
            GenotypeSourceFormat::Vcf => "vcf",
            GenotypeSourceFormat::Cram => "cram",
            GenotypeSourceFormat::Bam => "bam",
        }
    }

    pub(super) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        for rsid in &variant.rsids {
            if let Some(value) = self.values.get(rsid) {
                let mut evidence = vec![format!("resolved by rsid {rsid}")];
                // Mirror DelimitedBackend's `| source line: …` evidence so
                // wasm-side from_bytes loads produce byte-identical reports
                // to the CLI's path-backed DelimitedBackend.
                if let Some(source) = self.source_lines.get(rsid) {
                    evidence.push(format!("source line: {source}"));
                }
                return Ok(VariantObservation {
                    backend: self.backend_name().to_owned(),
                    matched_rsid: Some(rsid.clone()),
                    genotype: Some(value.clone()),
                    evidence,
                    ..VariantObservation::default()
                });
            }
        }

        Ok(VariantObservation {
            backend: self.backend_name().to_owned(),
            evidence: vec!["no matching rsid found".to_owned()],
            ..VariantObservation::default()
        })
    }
}

impl DelimitedBackend {
    pub(super) fn backend_name(&self) -> &'static str {
        match self.format {
            GenotypeSourceFormat::Text => "text",
            GenotypeSourceFormat::Zip => "zip",
            GenotypeSourceFormat::Vcf => "vcf",
            GenotypeSourceFormat::Cram => "cram",
            GenotypeSourceFormat::Bam => "bam",
        }
    }

    pub(super) fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        let results = self.lookup_variants(&[VariantSpec {
            rsids: vec![rsid.to_owned()],
            ..VariantSpec::default()
        }])?;
        Ok(results.into_iter().next().and_then(|obs| obs.genotype))
    }

    pub(super) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let mut results = self.lookup_variants(std::slice::from_ref(variant))?;
        Ok(results.pop().unwrap_or_default())
    }

    pub(super) fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        scan_delimited_variants(self, variants)
    }
}

impl VcfBackend {
    pub(super) fn backend_name(&self) -> &'static str {
        "vcf"
    }

    pub(super) fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        let results = self.lookup_variants(&[VariantSpec {
            rsids: vec![rsid.to_owned()],
            ..VariantSpec::default()
        }])?;
        Ok(results.into_iter().next().and_then(|obs| obs.genotype))
    }

    pub(super) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let mut results = self.lookup_variants(std::slice::from_ref(variant))?;
        Ok(results.pop().unwrap_or_default())
    }

    pub(super) fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        if let Some(results) = lookup_indexed_vcf_variants(self, variants)? {
            return Ok(results);
        }
        scan_vcf_variants(self, variants)
    }
}
