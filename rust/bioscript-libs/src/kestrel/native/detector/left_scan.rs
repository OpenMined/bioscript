use crate::LibResult;

use super::{ActiveRegionDetectorConfig, recovery_threshold, scan_limit_length};

pub(super) fn scan_left_start(
    counts: &[u32],
    index: usize,
    anchor_count: u32,
    kmer_size: usize,
    difference_threshold: u32,
    config: &ActiveRegionDetectorConfig,
) -> LibResult<Option<Option<usize>>> {
    let scan_limit = scan_limit_length(kmer_size, config)?;
    if index > scan_limit {
        return Ok(None);
    }

    let mut scan_end = index as isize - 1;
    while scan_end >= 0
        && (counts[scan_end as usize] as f32)
            < recovery_threshold(
                anchor_count,
                difference_threshold,
                index - scan_end as usize,
                kmer_size,
                config,
            )?
    {
        scan_end -= 1;
    }
    if scan_end > 0 {
        return Ok(None);
    }

    if config.recover_right_anchor && index < scan_limit {
        if let Some(anchor) =
            recover_left_anchor_index(counts, index, kmer_size, difference_threshold)
        {
            return Ok(Some(Some(anchor)));
        }
    }
    Ok(Some(None))
}

pub(super) fn skip_left_peak(
    counts: &[u32],
    index: usize,
    left: u32,
    right: u32,
    difference_threshold: u32,
    config: &ActiveRegionDetectorConfig,
) -> Option<usize> {
    if config.peak_scan_length == 0 {
        return None;
    }

    let java_difference_threshold = difference_threshold.saturating_sub(1);
    let recovery_value = left + java_difference_threshold;
    let scan_limit = index
        .saturating_add(config.peak_scan_length)
        .min(counts.len());
    let mut scan_index = index + 1;
    while scan_index < scan_limit {
        if counts[scan_index] <= recovery_value
            && right.saturating_sub(counts[scan_index]) < java_difference_threshold
        {
            return Some(scan_index + 1);
        }
        scan_index += 1;
    }
    None
}

fn recover_left_anchor_index(
    counts: &[u32],
    index: usize,
    kmer_size: usize,
    difference_threshold: u32,
) -> Option<usize> {
    let mut scan_index = index.saturating_sub(kmer_size);
    while scan_index > 0 {
        if counts[scan_index - 1] > counts[scan_index]
            && counts[scan_index - 1] - counts[scan_index] >= difference_threshold
        {
            return Some(scan_index);
        }
        scan_index -= 1;
    }
    None
}
