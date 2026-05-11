use crate::LibResult;

use super::{
    active_region::ActiveRegion,
    alignment::{align_haplotype, call_alignment_variants},
    detector::{ActiveRegionDetectorConfig, detect_active_regions},
    haplotype::{HaplotypeAssemblyConfig, assemble_haplotypes},
    kmer::KmerCountMap,
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
    let mut writer = new_writer(region, config)?;
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

pub fn call_assembled_haplotypes_to_vcf(
    region: &ReferenceRegion,
    active_region: &ActiveRegion,
    counts: &KmerCountMap,
    assembly_config: &HaplotypeAssemblyConfig,
    call_config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let mut writer = new_writer(region, call_config)?;
    let haplotypes = assemble_haplotypes(active_region, counts, assembly_config)?;
    add_active_region_haplotypes(
        &mut writer,
        region,
        active_region,
        &haplotypes,
        &call_config.sample_name,
    )?;
    Ok(writer.to_vcf_string())
}

pub fn call_sequences_to_vcf<'a>(
    region: &ReferenceRegion,
    read_sequences: impl IntoIterator<Item = &'a str>,
    kmer_size: usize,
    detector_config: &ActiveRegionDetectorConfig,
    assembly_config: &HaplotypeAssemblyConfig,
    call_config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let counts = KmerCountMap::from_sequences(read_sequences, kmer_size)?;
    let detection = detect_active_regions(region, &counts, detector_config)?;
    let mut writer = new_writer(region, call_config)?;
    for active_region in &detection.regions {
        let haplotypes = assemble_haplotypes(active_region, &counts, assembly_config)?;
        add_active_region_haplotypes(
            &mut writer,
            region,
            active_region,
            &haplotypes,
            &call_config.sample_name,
        )?;
    }
    Ok(writer.to_vcf_string())
}

fn new_writer(
    region: &ReferenceRegion,
    config: &NativeKestrelCallConfig,
) -> LibResult<KestrelVcfWriter> {
    let mut writer = KestrelVcfWriter::new(
        &config.source_version,
        vec![ReferenceSequence {
            name: region.reference_name.clone(),
            length: region.sequence.len(),
            md5: config.reference_md5.clone(),
        }],
    );
    writer.add_sample(&config.sample_name)?;
    Ok(writer)
}

fn add_active_region_haplotypes(
    writer: &mut KestrelVcfWriter,
    region: &ReferenceRegion,
    active_region: &ActiveRegion,
    haplotypes: &[HaplotypeEvidence],
    sample_name: &str,
) -> LibResult<()> {
    let active_reference = active_reference_sequence(region, active_region);
    let reference_start = u32::try_from(active_region.start_index + 1).unwrap_or(u32::MAX);
    for haplotype in haplotypes {
        let alignment = align_haplotype(&active_reference, &haplotype.sequence)?;
        for variant in call_alignment_variants(
            sample_name,
            &alignment,
            reference_start,
            haplotype.variant_depth,
            haplotype.locus_depth,
        )? {
            writer.add_native_variant(&variant, region)?;
        }
    }
    Ok(())
}

fn active_reference_sequence(region: &ReferenceRegion, active_region: &ActiveRegion) -> String {
    region.sequence[active_region.start_index..=active_region.end_index].to_owned()
}
