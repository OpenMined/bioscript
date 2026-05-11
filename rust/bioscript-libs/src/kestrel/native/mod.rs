mod active_region;
mod alignment;
mod alignment_weight;
mod detector;
mod engine;
mod haplotype;
mod kmer;
mod refreader;
mod variant;
mod vcf;

pub use active_region::{ActiveRegion, RegionStats};
pub use alignment::{
    AlignmentOp, NativeAlignment, align_haplotype, call_alignment_variants, score_alignment,
    score_haplotype_alignment,
};
pub use alignment_weight::AlignmentWeight;
pub use detector::{
    ActiveRegionDetection, ActiveRegionDetectorConfig, detect_active_regions, difference_threshold,
    recovery_threshold, scan_limit_length,
};
pub use engine::{
    HaplotypeEvidence, NativeKestrelCallConfig, NativeReferenceRegion,
    call_assembled_haplotypes_to_vcf, call_counted_kmers_to_vcf,
    call_counted_kmers_to_vcf_references, call_explicit_haplotypes_to_vcf, call_fastq_paths_to_vcf,
    call_fastq_paths_to_vcf_references, call_sequences_to_vcf,
};
pub use haplotype::{HaplotypeAssemblyConfig, assemble_haplotypes};
pub use kmer::{KmerCountMap, count_fastq_kmers, count_sequence_kmers};
pub use refreader::{ReferenceRecord, read_reference_records, reference_kmers};
pub use variant::{NativeVariantCall, ReferenceRegion, VariantKind};
pub use vcf::{KestrelVcfWriter, ReferenceSequence, VariantCall};
