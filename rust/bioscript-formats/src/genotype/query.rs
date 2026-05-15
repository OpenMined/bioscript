use bioscript_core::{RuntimeError, VariantObservation, VariantSpec};

use super::types::QueryBackend;
use super::{
    BackendCapabilities, GenotypeStore, QueryKind, match_cached_observation, required_cache_miss,
    variant_sort_key,
};

impl GenotypeStore {
    pub fn capabilities(&self) -> BackendCapabilities {
        match &self.backend {
            QueryBackend::RsidMap(_) => BackendCapabilities {
                rsid_lookup: true,
                locus_lookup: false,
            },
            QueryBackend::Delimited(_) | QueryBackend::Vcf(_) => BackendCapabilities {
                rsid_lookup: true,
                locus_lookup: true,
            },
            QueryBackend::Cram(_) | QueryBackend::Bam(_) | QueryBackend::AlignmentBytes(_) => {
                BackendCapabilities {
                    rsid_lookup: false,
                    locus_lookup: true,
                }
            }
            QueryBackend::Cached { .. } => {
                // The cache itself answers both rsid and locus queries (we
                // match by either), so unioning with the fallback gives the
                // strongest set.
                BackendCapabilities {
                    rsid_lookup: true,
                    locus_lookup: true,
                }
            }
        }
    }

    pub fn supports(&self, query: QueryKind) -> bool {
        let caps = self.capabilities();
        match query {
            QueryKind::GenotypeByRsid => caps.rsid_lookup,
            QueryKind::GenotypeByLocus => caps.locus_lookup,
        }
    }

    pub fn backend_name(&self) -> &'static str {
        match &self.backend {
            QueryBackend::RsidMap(map) => map.backend_name(),
            QueryBackend::Delimited(backend) => backend.backend_name(),
            QueryBackend::Vcf(backend) => backend.backend_name(),
            QueryBackend::Cram(backend) => backend.backend_name(),
            QueryBackend::Bam(backend) => backend.backend_name(),
            QueryBackend::AlignmentBytes(backend) => backend.backend_name(),
            QueryBackend::Cached { .. } => "cached",
        }
    }

    pub fn get(&self, rsid: &str) -> Result<Option<String>, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => Ok(map.values.get(rsid).cloned()),
            QueryBackend::Delimited(backend) => backend.get(rsid),
            QueryBackend::Vcf(backend) => backend.get(rsid),
            QueryBackend::Cram(backend) => backend
                .lookup_variant(&VariantSpec {
                    rsids: vec![rsid.to_owned()],
                    ..VariantSpec::default()
                })
                .map(|obs| obs.genotype),
            QueryBackend::Bam(backend) => backend
                .lookup_variant(&VariantSpec {
                    rsids: vec![rsid.to_owned()],
                    ..VariantSpec::default()
                })
                .map(|obs| obs.genotype),
            QueryBackend::AlignmentBytes(backend) => backend
                .lookup_variant(&VariantSpec {
                    rsids: vec![rsid.to_owned()],
                    ..VariantSpec::default()
                })
                .map(|obs| obs.genotype),
            QueryBackend::Cached {
                observations,
                fallback,
                require_hit,
            } => {
                if let Some(matched) = observations
                    .iter()
                    .find(|obs| obs.matched_rsid.as_deref().is_some_and(|r| r == rsid))
                {
                    return Ok(matched.genotype.clone());
                }
                if *require_hit {
                    return Err(required_cache_miss(&VariantSpec {
                        rsids: vec![rsid.to_owned()],
                        ..VariantSpec::default()
                    }));
                }
                let inner = GenotypeStore {
                    backend: (**fallback).clone(),
                };
                inner.get(rsid)
            }
        }
    }

    pub fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        match &self.backend {
            QueryBackend::RsidMap(map) => map.lookup_variant(variant),
            QueryBackend::Delimited(backend) => backend.lookup_variant(variant),
            QueryBackend::Vcf(backend) => backend.lookup_variant(variant),
            QueryBackend::Cram(backend) => backend.lookup_variant(variant),
            QueryBackend::Bam(backend) => backend.lookup_variant(variant),
            QueryBackend::AlignmentBytes(backend) => backend.lookup_variant(variant),
            QueryBackend::Cached {
                observations,
                fallback,
                require_hit,
            } => {
                if let Some(hit) = match_cached_observation(observations, variant) {
                    return Ok(hit.clone());
                }
                if *require_hit {
                    return Err(required_cache_miss(variant));
                }
                let inner = GenotypeStore {
                    backend: (**fallback).clone(),
                };
                inner.lookup_variant(variant)
            }
        }
    }

    pub fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        if let QueryBackend::Cached {
            observations,
            fallback,
            require_hit,
        } = &self.backend
        {
            // Resolve cache hits up-front; only round-trip the fallback for
            // misses so we don't pay for re-opening a CRAM/VCF when the panel
            // already covered every variant the script needs.
            let mut results: Vec<Option<VariantObservation>> = vec![None; variants.len()];
            let mut miss_indices = Vec::new();
            let mut miss_specs = Vec::new();
            for (idx, spec) in variants.iter().enumerate() {
                if let Some(hit) = match_cached_observation(observations, spec) {
                    results[idx] = Some(hit.clone());
                } else {
                    if *require_hit {
                        return Err(required_cache_miss(spec));
                    }
                    miss_indices.push(idx);
                    miss_specs.push(spec.clone());
                }
            }
            if !miss_specs.is_empty() {
                let inner = GenotypeStore {
                    backend: (**fallback).clone(),
                };
                let resolved = inner.lookup_variants(&miss_specs)?;
                for (idx, observation) in miss_indices.into_iter().zip(resolved) {
                    results[idx] = Some(observation);
                }
            }
            return Ok(results.into_iter().map(Option::unwrap_or_default).collect());
        }
        if let QueryBackend::Delimited(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::Vcf(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::Cram(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::Bam(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        if let QueryBackend::AlignmentBytes(backend) = &self.backend {
            return backend.lookup_variants(variants);
        }
        let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
        indexed.sort_by_cached_key(|(_, variant)| variant_sort_key(variant));

        let mut results = vec![VariantObservation::default(); variants.len()];
        for (original_idx, variant) in indexed {
            results[original_idx] = self.lookup_variant(variant)?;
        }
        Ok(results)
    }
}
