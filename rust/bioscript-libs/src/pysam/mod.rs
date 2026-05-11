mod aligned_segment;
mod alignment_file;

pub use aligned_segment::AlignedSegment;
pub use alignment_file::{AlignmentFetch, AlignmentFile, AlignmentMode};

pub const MODULE: &str = "pysam";
