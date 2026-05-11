use crate::{LibError, LibResult};

use super::{active_region::ActiveRegion, kmer::KmerCountMap, variant::ReferenceRegion};

mod left_scan;
mod right_scan;

use left_scan::{scan_left_start, skip_left_peak};
use right_scan::scan_right_end;

#[derive(Debug, Clone, PartialEq)]
pub struct ActiveRegionDetectorConfig {
    pub minimum_difference: u32,
    pub difference_quantile: f32,
    pub count_reverse_kmers: bool,
    pub anchor_both_ends: bool,
    pub decay_min: f32,
    pub decay_alpha: f32,
    pub peak_scan_length: usize,
    pub scan_limit_factor: f32,
    pub max_gap_size: usize,
    pub recover_right_anchor: bool,
    pub call_ambiguous_regions: bool,
}

impl Default for ActiveRegionDetectorConfig {
    fn default() -> Self {
        Self {
            minimum_difference: 5,
            difference_quantile: 0.90,
            count_reverse_kmers: true,
            anchor_both_ends: true,
            decay_min: 0.55,
            decay_alpha: 0.80,
            peak_scan_length: 7,
            scan_limit_factor: 7.0,
            max_gap_size: 0,
            recover_right_anchor: true,
            call_ambiguous_regions: true,
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
        config,
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

pub fn recovery_threshold(
    anchor_count: u32,
    difference_threshold: u32,
    distance: usize,
    kmer_size: usize,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<f32> {
    validate_decay(config)?;
    if config.decay_min == 1.0 {
        return Ok(anchor_count.saturating_sub(difference_threshold).max(1) as f32);
    }

    let min_value = (anchor_count as f32 * config.decay_min).max(1.0);
    let range = anchor_count as f32 - min_value;
    let lambda = -config.decay_alpha.ln() / kmer_size as f32;
    Ok(range * (-(distance as f32) * lambda).exp() + min_value)
}

fn candidate_regions(
    region: &ReferenceRegion,
    counts: &[u32],
    kmer_size: usize,
    difference_threshold: u32,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<Vec<ActiveRegion>> {
    if counts.len() < 2 {
        return Ok(Vec::new());
    }

    let mut regions = Vec::new();
    let mut index = 1usize;
    let mut last_region_end = 0usize;
    while index < counts.len() {
        let left = counts[index - 1];
        let right = counts[index];
        if left > right && left - right >= difference_threshold {
            let Some(end) =
                scan_right_end(counts, index, left, kmer_size, difference_threshold, config)?
            else {
                index += 1;
                continue;
            };
            if end < counts.len() && end.saturating_sub(index) >= kmer_size.saturating_sub(1) {
                if !config.call_ambiguous_regions
                    && contains_ambiguous_region_base(region, index, end + kmer_size)
                {
                    index += 1;
                    continue;
                }
                regions.push(ActiveRegion::new(
                    region,
                    Some(index - 1),
                    Some(end),
                    counts,
                    kmer_size,
                )?);
                last_region_end = end;
                index = end + 1;
                continue;
            }
            if !config.anchor_both_ends
                && end == counts.len()
                && end.saturating_sub(index) >= kmer_size.saturating_sub(1)
            {
                if !config.call_ambiguous_regions
                    && contains_ambiguous_region_base(region, index, region.sequence.len())
                {
                    break;
                }
                regions.push(ActiveRegion::new(
                    region,
                    Some(index - 1),
                    None,
                    counts,
                    kmer_size,
                )?);
                break;
            }
        } else if right > left && right - left >= difference_threshold {
            if let Some(next_index) =
                skip_left_peak(counts, index, left, right, difference_threshold, config)
            {
                index = next_index;
                continue;
            }
            let Some(start) = scan_left_start(
                counts,
                index,
                right,
                kmer_size,
                difference_threshold,
                config,
            )?
            else {
                index += 1;
                continue;
            };
            if start.is_none() && (config.anchor_both_ends || index < kmer_size.saturating_sub(1)) {
                index += 1;
                continue;
            }
            let start_base = start.unwrap_or(0);
            if last_region_end > 0 && start_base < last_region_end {
                index += 1;
                continue;
            }
            if !config.call_ambiguous_regions
                && contains_ambiguous_region_base(region, start_base, index + kmer_size)
            {
                index += 1;
                continue;
            }
            regions.push(ActiveRegion::new(
                region,
                start,
                Some(index),
                counts,
                kmer_size,
            )?);
            last_region_end = index;
            index += 1;
            continue;
        }
        index += 1;
    }
    Ok(regions)
}

fn contains_ambiguous_region_base(region: &ReferenceRegion, start: usize, end: usize) -> bool {
    region.sequence[start.min(region.sequence.len())..end.min(region.sequence.len())]
        .bytes()
        .any(|base| !matches!(base, b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't'))
}

pub fn scan_limit_length(
    kmer_size: usize,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<usize> {
    validate_scan_limit(config)?;
    let scaled = (config.scan_limit_factor * kmer_size as f32) as usize;
    Ok(kmer_size.max(config.max_gap_size.saturating_add(scaled)))
}

fn validate_config(config: &ActiveRegionDetectorConfig) -> LibResult<()> {
    if config.minimum_difference == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel active-region minimum difference must be at least 1".to_owned(),
        ));
    }
    validate_difference_quantile(config.difference_quantile)?;
    validate_decay(config)?;
    validate_scan_limit(config)
}

fn validate_difference_quantile(difference_quantile: f32) -> LibResult<()> {
    if !(0.0..1.0).contains(&difference_quantile) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel active-region difference quantile must be in [0.0, 1.0): {difference_quantile}"
        )));
    }
    Ok(())
}

fn validate_scan_limit(config: &ActiveRegionDetectorConfig) -> LibResult<()> {
    if config.scan_limit_factor < 0.0 || !config.scan_limit_factor.is_finite() {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel active-region scan limit factor must be finite and nonnegative: {}",
            config.scan_limit_factor
        )));
    }
    Ok(())
}

fn validate_decay(config: &ActiveRegionDetectorConfig) -> LibResult<()> {
    if !(0.0..=1.0).contains(&config.decay_min) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel active-region decay minimum must be in [0.0, 1.0]: {}",
            config.decay_min
        )));
    }
    if !(0.0..1.0).contains(&config.decay_alpha) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel active-region decay alpha must be in (0.0, 1.0): {}",
            config.decay_alpha
        )));
    }
    Ok(())
}
