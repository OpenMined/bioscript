use std::path::Path;

use crate::LibResult;

use super::{
    active_region::ActiveRegion,
    alignment::{NativeAlignment, align_haplotype, call_alignment_variants, score_alignment},
    alignment_weight::AlignmentWeight,
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
pub struct NativeReferenceRegion {
    pub reference_name: String,
    pub sequence: String,
    pub md5: String,
}

impl NativeReferenceRegion {
    pub fn new(
        reference_name: impl Into<String>,
        sequence: impl Into<String>,
        md5: impl Into<String>,
    ) -> Self {
        Self {
            reference_name: reference_name.into(),
            sequence: sequence.into(),
            md5: md5.into(),
        }
    }

    fn region(&self) -> ReferenceRegion {
        ReferenceRegion {
            reference_name: self.reference_name.clone(),
            sequence: self.sequence.clone(),
        }
    }
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
    call_counted_kmers_to_vcf(
        region,
        &counts,
        detector_config,
        assembly_config,
        call_config,
    )
}

pub fn call_fastq_paths_to_vcf<'a>(
    region: &ReferenceRegion,
    fastq_paths: impl IntoIterator<Item = &'a Path>,
    kmer_size: usize,
    detector_config: &ActiveRegionDetectorConfig,
    assembly_config: &HaplotypeAssemblyConfig,
    call_config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let counts = KmerCountMap::from_fastq_paths(fastq_paths, kmer_size)?;
    call_counted_kmers_to_vcf(
        region,
        &counts,
        detector_config,
        assembly_config,
        call_config,
    )
}

pub fn call_fastq_paths_to_vcf_references<'a>(
    references: &[NativeReferenceRegion],
    fastq_paths: impl IntoIterator<Item = &'a Path>,
    kmer_size: usize,
    detector_config: &ActiveRegionDetectorConfig,
    assembly_config: &HaplotypeAssemblyConfig,
    call_config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let counts = KmerCountMap::from_fastq_paths(fastq_paths, kmer_size)?;
    call_counted_kmers_to_vcf_references(
        references,
        &counts,
        detector_config,
        assembly_config,
        call_config,
    )
}

pub fn call_counted_kmers_to_vcf(
    region: &ReferenceRegion,
    counts: &KmerCountMap,
    detector_config: &ActiveRegionDetectorConfig,
    assembly_config: &HaplotypeAssemblyConfig,
    call_config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let detection = detect_active_regions(region, counts, detector_config)?;
    let mut writer = new_writer(region, call_config)?;
    for active_region in &detection.regions {
        let haplotypes = assemble_haplotypes(active_region, counts, assembly_config)?;
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

pub fn call_counted_kmers_to_vcf_references(
    references: &[NativeReferenceRegion],
    counts: &KmerCountMap,
    detector_config: &ActiveRegionDetectorConfig,
    assembly_config: &HaplotypeAssemblyConfig,
    call_config: &NativeKestrelCallConfig,
) -> LibResult<String> {
    let mut writer = new_writer_for_references(references, call_config)?;
    for reference in references {
        let region = reference.region();
        let detection = detect_active_regions(&region, counts, detector_config)?;
        for active_region in &detection.regions {
            let haplotypes = assemble_haplotypes(active_region, counts, assembly_config)?;
            add_active_region_haplotypes(
                &mut writer,
                &region,
                active_region,
                &haplotypes,
                &call_config.sample_name,
            )?;
        }
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

fn new_writer_for_references(
    references: &[NativeReferenceRegion],
    config: &NativeKestrelCallConfig,
) -> LibResult<KestrelVcfWriter> {
    let reference_sequences = references
        .iter()
        .map(|reference| ReferenceSequence {
            name: reference.reference_name.clone(),
            length: reference.sequence.len(),
            md5: reference.md5.clone(),
        })
        .collect();
    let mut writer = KestrelVcfWriter::new(&config.source_version, reference_sequences);
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
    for (haplotype, alignment) in max_scoring_haplotypes(&active_reference, haplotypes)? {
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

fn max_scoring_haplotypes<'a>(
    active_reference: &str,
    haplotypes: &'a [HaplotypeEvidence],
) -> LibResult<Vec<(&'a HaplotypeEvidence, NativeAlignment)>> {
    let weight = AlignmentWeight::default();
    let mut scored = Vec::new();
    let mut max_score = f32::NEG_INFINITY;
    for haplotype in haplotypes {
        let alignment = align_haplotype(active_reference, &haplotype.sequence)?;
        let score = score_alignment(&alignment, &weight);
        if score > max_score {
            max_score = score;
        }
        scored.push((score, haplotype, alignment));
    }
    if scored
        .iter()
        .any(|(_, haplotype, _)| haplotype.sequence != active_reference)
    {
        scored.retain(|(_, haplotype, _)| haplotype.sequence != active_reference);
        max_score = scored
            .iter()
            .map(|(score, _, _)| *score)
            .fold(f32::NEG_INFINITY, f32::max);
    }
    Ok(scored
        .into_iter()
        .filter(|(score, _, _)| (*score - max_score).abs() <= f32::EPSILON)
        .map(|(_, haplotype, alignment)| (haplotype, alignment))
        .collect())
}

fn active_reference_sequence(region: &ReferenceRegion, active_region: &ActiveRegion) -> String {
    region.sequence[active_region.start_index..=active_region.end_index].to_owned()
}
