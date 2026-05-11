use crate::LibResult;

use super::{
    alignment::{align_haplotype, call_alignment_variants},
    variant::ReferenceRegion,
    vcf::{KestrelVcfWriter, ReferenceSequence},
};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HaplotypeEvidence {
    pub sequence: String,
    pub variant_depth: u32,
    pub locus_depth: u32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NativeKestrelCallConfig {
    pub source_version: String,
    pub sample_name: String,
    pub reference_md5: String,
}

impl NativeKestrelCallConfig {
    pub fn new(
        source_version: impl Into<String>,
        sample_name: impl Into<String>,
        reference_md5: impl Into<String>,
    ) -> Self {
        Self {
            source_version: source_version.into(),
            sample_name: sample_name.into(),
            reference_md5: reference_md5.into(),
        }
    }
}

pub fn call_explicit_haplotypes_to_vcf(
    region: &ReferenceRegion,
    haplotypes: &[HaplotypeEvidence],
    config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let mut writer = KestrelVcfWriter::new(
        &config.source_version,
        vec![ReferenceSequence {
            name: region.reference_name.clone(),
            length: region.sequence.len(),
            md5: config.reference_md5.clone(),
        }],
    );
    writer.add_sample(&config.sample_name)?;

    for haplotype in haplotypes {
        let alignment = align_haplotype(&region.sequence, &haplotype.sequence)?;
        for variant in call_alignment_variants(
            &config.sample_name,
            &alignment,
            1,
            haplotype.variant_depth,
            haplotype.locus_depth,
        )? {
            writer.add_native_variant(&variant, region)?;
        }
    }
    Ok(writer.to_vcf_string())
}
