use bioscript_core::{RuntimeError, VariantObservation, VariantSpec};

/// Match a `VariantSpec` against a pre-resolved observation list. Tries rsid
/// equality first (most common case for `PGx` panels), then falls back to a
/// chrom+pos+ref+alt match against either `GRCh37` or `GRCh38` loci so cached
/// observations from a CRAM lookup (which may have been done on one assembly)
/// can satisfy a script that supplies the spec on the other.
pub(crate) fn match_cached_observation<'a>(
    observations: &'a [VariantObservation],
    spec: &VariantSpec,
) -> Option<&'a VariantObservation> {
    if let Some(matched) = observations.iter().find(|obs| {
        obs.matched_rsid
            .as_deref()
            .is_some_and(|rsid| spec.rsids.iter().any(|target| target == rsid))
    }) {
        return Some(matched);
    }
    let assembly_loci = [spec.grch37.as_ref(), spec.grch38.as_ref()]
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let target_ref = spec.reference.as_deref();
    let target_alt = spec.alternate.as_deref();
    observations.iter().find(|obs| {
        let evidence_match = assembly_loci.iter().any(|loci| {
            obs.evidence
                .iter()
                .any(|line| line.contains(&loci.chrom) && line.contains(&loci.start.to_string()))
        });
        if !evidence_match {
            return false;
        }
        match (target_ref, target_alt) {
            (Some(r), Some(a)) => obs
                .evidence
                .iter()
                .any(|line| line.contains(r) && line.contains(a)),
            _ => true,
        }
    })
}

pub(crate) fn required_cache_miss(spec: &VariantSpec) -> RuntimeError {
    let rsids = if spec.rsids.is_empty() {
        "<none>".to_owned()
    } else {
        spec.rsids.join("|")
    };
    let loci = [
        spec.grch37
            .as_ref()
            .map(|locus| format!("grch37:{}:{}-{}", locus.chrom, locus.start, locus.end)),
        spec.grch38
            .as_ref()
            .map(|locus| format!("grch38:{}:{}-{}", locus.chrom, locus.start, locus.end)),
    ]
    .into_iter()
    .flatten()
    .collect::<Vec<_>>()
    .join(",");
    RuntimeError::InvalidArguments(format!(
        "required preloaded genotype observation missing for rsids={rsids} loci={loci}"
    ))
}
