use crate::{LibError, LibResult};

use super::{active_region::ActiveRegion, kmer::KmerCountMap, variant::ReferenceRegion};

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
    pub recover_right_anchor: bool,
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
            recover_right_anchor: true,
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
            if !config.anchor_both_ends
                && end == counts.len()
                && end.saturating_sub(index) >= kmer_size.saturating_sub(1)
            {
                regions.push(ActiveRegion::new(
                    region,
                    Some(index - 1),
                    None,
                    counts,
                    kmer_size,
                )?);
                break;
            }
        } else if right > left
            && right - left >= difference_threshold
            && !config.anchor_both_ends
            && index >= kmer_size.saturating_sub(1)
        {
            regions.push(ActiveRegion::new(
                region,
                None,
                Some(index),
                counts,
                kmer_size,
            )?);
            index += 1;
            continue;
        }
        index += 1;
    }
    Ok(regions)
}

fn scan_right_end(
    counts: &[u32],
    start_index: usize,
    anchor_count: u32,
    kmer_size: usize,
    difference_threshold: u32,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<Option<usize>> {
    let mut end = start_index + 1;
    let mut peak_count = 0usize;
    let mut peak_scan_index = 0usize;
    let mut last_valley_index = 0usize;
    let scan_limit = scan_limit_length(kmer_size, config)?;

    'scan_loop: loop {
        while end < counts.len()
            && end.saturating_sub(start_index) <= scan_limit
            && (counts[end] as f32)
                < recovery_threshold(
                    anchor_count,
                    difference_threshold,
                    end - start_index,
                    kmer_size,
                    config,
                )?
        {
            end += 1;
        }
        if end.saturating_sub(start_index) > scan_limit {
            return Ok(None);
        }

        if config.peak_scan_length == 0 {
            if end == counts.len() && config.recover_right_anchor {
                if let Some(anchor) =
                    recover_right_anchor_index(counts, start_index, kmer_size, difference_threshold)
                {
                    return Ok(Some(anchor));
                }
            }
            return Ok(Some(end));
        }

        if peak_scan_index > 0 && end.saturating_sub(peak_scan_index) >= kmer_size {
            last_valley_index = end;
        } else if peak_scan_index == 0 && end.saturating_sub(start_index) >= kmer_size {
            last_valley_index = end;
        }

        let recovery_value = recovery_threshold(
            anchor_count,
            difference_threshold,
            end.saturating_sub(start_index),
            kmer_size,
            config,
        )?;
        peak_scan_index = end;
        let peak_scan_limit = end
            .saturating_add(config.peak_scan_length)
            .min(counts.len());

        while peak_scan_index < peak_scan_limit {
            if (counts[peak_scan_index] as f32) < recovery_value {
                peak_count += 1;
                end = peak_scan_index;
                if peak_count > 3 && end.saturating_sub(start_index) / peak_count < kmer_size {
                    return Ok(Some(last_valley_index.max(start_index + 1)));
                }
                continue 'scan_loop;
            }
            peak_scan_index += 1;
        }

        if peak_scan_index == counts.len() && last_valley_index > 0 {
            return Ok(Some(last_valley_index));
        }

        if end == counts.len() && config.recover_right_anchor {
            if let Some(anchor) =
                recover_right_anchor_index(counts, start_index, kmer_size, difference_threshold)
            {
                return Ok(Some(anchor));
            }
        }

        return Ok(Some(end));
    }
}

fn recover_right_anchor_index(
    counts: &[u32],
    start_index: usize,
    kmer_size: usize,
    difference_threshold: u32,
) -> Option<usize> {
    let mut index = start_index + kmer_size;
    while index < counts.len() {
        if counts[index] > counts[index - 1]
            && counts[index] - counts[index - 1] >= difference_threshold
        {
            return Some(index);
        }
        index += 1;
    }
    None
}

pub fn scan_limit_length(
    kmer_size: usize,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<usize> {
    validate_scan_limit(config)?;
    let scaled = (config.scan_limit_factor * kmer_size as f32) as usize;
    Ok(kmer_size.max(scaled))
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
