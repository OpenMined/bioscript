mod active_region;
mod alignment;
mod detector;
mod engine;
mod kmer;
mod variant;
mod vcf;

pub use active_region::{ActiveRegion, RegionStats};
pub use alignment::{AlignmentOp, NativeAlignment, align_haplotype, call_alignment_variants};
pub use detector::{
    ActiveRegionDetection, ActiveRegionDetectorConfig, detect_active_regions, difference_threshold,
};
pub use engine::{HaplotypeEvidence, NativeKestrelCallConfig, call_explicit_haplotypes_to_vcf};
pub use kmer::{KmerCountMap, count_fastq_kmers, count_sequence_kmers};
pub use variant::{NativeVariantCall, ReferenceRegion, VariantKind};
pub use vcf::{KestrelVcfWriter, ReferenceSequence, VariantCall};
