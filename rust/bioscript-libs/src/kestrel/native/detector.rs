use crate::{LibError, LibResult};

use super::{active_region::ActiveRegion, kmer::KmerCountMap, variant::ReferenceRegion};

#[derive(Debug, Clone, PartialEq)]
pub struct ActiveRegionDetectorConfig {
    pub minimum_difference: u32,
    pub difference_quantile: f32,
    pub count_reverse_kmers: bool,
}

impl Default for ActiveRegionDetectorConfig {
    fn default() -> Self {
        Self {
            minimum_difference: 5,
            difference_quantile: 0.90,
            count_reverse_kmers: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ActiveRegionDetection {
    pub reference_counts: Vec<u32>,
    pub difference_threshold: u32,
    pub regions: Vec<ActiveRegion>,
}

pub fn detect_active_regions(
    region: &ReferenceRegion,
    counts: &KmerCountMap,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<ActiveRegionDetection> {
    validate_config(config)?;
    let reference_counts = counts.reference_counts(&region.sequence, config.count_reverse_kmers)?;
    let difference_threshold = difference_threshold(
        &reference_counts,
        config.minimum_difference,
        config.difference_quantile,
    )?;
    let regions = candidate_regions(
        region,
        &reference_counts,
        counts.kmer_size(),
        difference_threshold,
    )?;
    Ok(ActiveRegionDetection {
        reference_counts,
        difference_threshold,
        regions,
    })
}

pub fn difference_threshold(
    counts: &[u32],
    minimum_difference: u32,
    difference_quantile: f32,
) -> LibResult<u32> {
    validate_difference_quantile(difference_quantile)?;
    if counts.len() < 3 {
        return Ok(minimum_difference);
    }

    let mut diffs = Vec::with_capacity(counts.len() - 1);
    let mut last_count = counts[0];
    for count in counts.iter().take(counts.len() - 1) {
        diffs.push(last_count.abs_diff(*count));
        last_count = *count;
    }
    diffs.sort_unstable();

    let threshold = if difference_quantile > 0.0 {
        let n_less_one = (diffs.len() - 1) as f32;
        let position = n_less_one * difference_quantile;
        let loc = position as usize;
        let offset = position - loc as f32;
        (diffs[loc] as f32 * (1.0 - offset) + diffs[loc + 1] as f32 * offset) as u32
    } else {
        minimum_difference
    };
    Ok(threshold.max(minimum_difference))
}

fn candidate_regions(
    region: &ReferenceRegion,
    counts: &[u32],
    kmer_size: usize,
    difference_threshold: u32,
) -> LibResult<Vec<ActiveRegion>> {
    if counts.len() < 2 {
        return Ok(Vec::new());
    }

    let mut regions = Vec::new();
    let mut index = 1usize;
    while index < counts.len() {
        let left = counts[index - 1];
        let right = counts[index];
        if left > right && left - right >= difference_threshold {
            let recovery_value = left.saturating_sub(difference_threshold).max(1);
            let mut end = index + 1;
            while end < counts.len() && counts[end] < recovery_value {
                end += 1;
            }
            if end < counts.len() && end.saturating_sub(index) >= kmer_size.saturating_sub(1) {
                regions.push(ActiveRegion::new(
                    region,
                    Some(index - 1),
                    Some(end),
                    counts,
                    kmer_size,
                )?);
                index = end + 1;
                continue;
            }
            if end == counts.len() && end.saturating_sub(index) >= kmer_size.saturating_sub(1) {
                regions.push(ActiveRegion::new(
                    region,
                    Some(index - 1),
                    None,
                    counts,
                    kmer_size,
                )?);
                break;
            }
        }
        index += 1;
    }
    Ok(regions)
}

fn validate_config(config: &ActiveRegionDetectorConfig) -> LibResult<()> {
    if config.minimum_difference == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel active-region minimum difference must be at least 1".to_owned(),
        ));
    }
    validate_difference_quantile(config.difference_quantile)
}

fn validate_difference_quantile(difference_quantile: f32) -> LibResult<()> {
    if !(0.0..1.0).contains(&difference_quantile) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel active-region difference quantile must be in [0.0, 1.0): {difference_quantile}"
        )));
    }
    Ok(())
}
