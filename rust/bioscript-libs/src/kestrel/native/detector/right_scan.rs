use crate::LibResult;

use super::{ActiveRegionDetectorConfig, recovery_threshold, scan_limit_length};

pub(super) fn scan_right_end(
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
