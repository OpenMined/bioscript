use crate::{LibError, LibResult};

use super::variant::ReferenceRegion;

#[derive(Debug, Clone, PartialEq)]
pub struct RegionStats {
    pub min: u32,
    pub pct25: f32,
    pub pct50: f32,
    pub pct75: f32,
    pub max: u32,
    pub n: usize,
}

impl RegionStats {
    pub fn from_counts(counts: &[u32], start: usize, end: usize) -> LibResult<Self> {
        if start > end {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel region stats start {start} is after end {end}"
            )));
        }
        if end > counts.len() || end == start {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel region stats range [{start}, {end}) is empty or outside {} counts",
                counts.len()
            )));
        }

        let mut slice = counts[start..end].to_vec();
        slice.sort_unstable();
        let n = slice.len();
        if n == 1 {
            let count = slice[0];
            return Ok(Self {
                min: count,
                pct25: count as f32,
                pct50: count as f32,
                pct75: count as f32,
                max: count,
                n,
            });
        }

        Ok(Self {
            min: slice[0],
            pct25: percentile(&slice, 0.25),
            pct50: percentile(&slice, 0.50),
            pct75: percentile(&slice, 0.75),
            max: slice[n - 1],
            n,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ActiveRegion {
    pub reference_name: String,
    pub start_index: usize,
    pub end_index: usize,
    pub start_kmer_index: usize,
    pub end_kmer_index: usize,
    pub left_end: bool,
    pub right_end: bool,
    pub left_end_kmer: Option<String>,
    pub right_end_kmer: Option<String>,
    pub stats: RegionStats,
}

impl ActiveRegion {
    pub fn new(
        region: &ReferenceRegion,
        start_kmer_index: Option<usize>,
        end_kmer_index: Option<usize>,
        counts: &[u32],
        kmer_size: usize,
    ) -> LibResult<Self> {
        validate_region_args(region, start_kmer_index, end_kmer_index, counts, kmer_size)?;
        let sequence_len = region.sequence.len();
        let left_end = start_kmer_index.is_none();
        let right_end = end_kmer_index.is_none();
        let start_kmer_index = start_kmer_index.unwrap_or(0);
        let end_kmer_index = end_kmer_index.unwrap_or(counts.len() - 1);
        let start_index = if left_end { 0 } else { start_kmer_index };
        let end_index = if right_end {
            sequence_len - 1
        } else {
            end_kmer_index + kmer_size - 1
        };
        let left_end_kmer = if left_end {
            None
        } else {
            Some(reference_kmer(region, start_kmer_index, kmer_size)?)
        };
        let right_end_kmer = if right_end {
            None
        } else {
            Some(reference_kmer(region, end_kmer_index, kmer_size)?)
        };

        Ok(Self {
            reference_name: region.reference_name.clone(),
            start_index,
            end_index,
            start_kmer_index,
            end_kmer_index,
            left_end,
            right_end,
            left_end_kmer,
            right_end_kmer,
            stats: RegionStats::from_counts(counts, start_kmer_index, end_kmer_index)?,
        })
    }

    pub fn matches_left_end(&self, kmer: &str) -> bool {
        self.left_end_kmer.as_deref() == Some(kmer)
    }

    pub fn matches_right_end(&self, kmer: &str) -> bool {
        self.right_end_kmer.as_deref() == Some(kmer)
    }
}

fn percentile(sorted_counts: &[u32], quantile: f32) -> f32 {
    let n_less_one = (sorted_counts.len() - 1) as f32;
    let position = n_less_one * quantile;
    let loc = position as usize;
    let offset = position - loc as f32;
    sorted_counts[loc] as f32 * (1.0 - offset) + sorted_counts[loc + 1] as f32 * offset
}

fn validate_region_args(
    region: &ReferenceRegion,
    start_kmer_index: Option<usize>,
    end_kmer_index: Option<usize>,
    counts: &[u32],
    kmer_size: usize,
) -> LibResult<()> {
    if kmer_size == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel active-region k-mer size must be greater than zero".to_owned(),
        ));
    }
    if counts.is_empty() {
        return Err(LibError::InvalidArguments(
            "Kestrel active-region counts cannot be empty".to_owned(),
        ));
    }
    if start_kmer_index.is_none() && end_kmer_index.is_none() {
        return Err(LibError::InvalidArguments(
            "Kestrel active region may not span the entire reference".to_owned(),
        ));
    }
    if let Some(end) = end_kmer_index {
        if end >= counts.len() || end + kmer_size > region.sequence.len() {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel active-region end k-mer index {end} is outside {} counts",
                counts.len()
            )));
        }
    }
    if let (Some(start), Some(end)) = (start_kmer_index, end_kmer_index) {
        if start >= end {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel active-region start {start} must come before end {end}"
            )));
        }
    }
    Ok(())
}

fn reference_kmer(region: &ReferenceRegion, start: usize, kmer_size: usize) -> LibResult<String> {
    let end = start + kmer_size;
    let kmer = region
        .sequence
        .get(start..end)
        .ok_or_else(|| {
            LibError::InvalidArguments(format!(
                "Kestrel reference k-mer [{start}, {end}) is outside {}",
                region.reference_name
            ))
        })?
        .to_ascii_uppercase();
    if kmer
        .bytes()
        .any(|base| !matches!(base, b'A' | b'C' | b'G' | b'T'))
    {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel reference k-mer contains ambiguous bases: {kmer}"
        )));
    }
    Ok(kmer)
}
