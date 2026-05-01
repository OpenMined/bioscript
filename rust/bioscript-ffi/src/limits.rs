use std::time::Duration;

use monty::ResourceLimits;

use crate::RunFileRequest;

pub(crate) const DEFAULT_MAX_DURATION_MS: u64 = 100;
pub(crate) const DEFAULT_MAX_MEMORY_BYTES: usize = 8 * 1024 * 1024;
pub(crate) const DEFAULT_MAX_ALLOCATIONS: usize = 200_000;
pub(crate) const DEFAULT_MAX_RECURSION_DEPTH: usize = 200;
pub(crate) const HARD_MAX_DURATION_MS: u64 = 60_000;
pub(crate) const HARD_MAX_MEMORY_BYTES: usize = 256 * 1024 * 1024;
pub(crate) const HARD_MAX_ALLOCATIONS: usize = 10_000_000;
pub(crate) const HARD_MAX_RECURSION_DEPTH: usize = 10_000;

pub(crate) fn build_limits(request: &RunFileRequest) -> ResourceLimits {
    let mut limits = ResourceLimits::new()
        .max_duration(Duration::from_millis(DEFAULT_MAX_DURATION_MS))
        .max_memory(DEFAULT_MAX_MEMORY_BYTES)
        .max_allocations(DEFAULT_MAX_ALLOCATIONS)
        .gc_interval(1000)
        .max_recursion_depth(Some(DEFAULT_MAX_RECURSION_DEPTH));

    if let Some(value) = request.max_duration_ms {
        limits = limits.max_duration(Duration::from_millis(value.min(HARD_MAX_DURATION_MS)));
    }
    if let Some(value) = request.max_memory_bytes {
        limits = limits.max_memory(value.min(HARD_MAX_MEMORY_BYTES));
    }
    if let Some(value) = request.max_allocations {
        limits = limits.max_allocations(value.min(HARD_MAX_ALLOCATIONS));
    }
    if let Some(value) = request.max_recursion_depth {
        limits = limits.max_recursion_depth(Some(value.min(HARD_MAX_RECURSION_DEPTH)));
    }

    limits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request_with_limits() -> RunFileRequest {
        RunFileRequest {
            script_path: "script.py".to_owned(),
            root: None,
            input_file: None,
            output_file: None,
            participant_id: None,
            trace_report_path: None,
            timing_report_path: None,
            input_format: None,
            input_index: None,
            reference_file: None,
            reference_index: None,
            allow_md5_mismatch: None,
            auto_index: None,
            cache_dir: None,
            max_duration_ms: Some(u64::MAX),
            max_memory_bytes: Some(usize::MAX),
            max_allocations: Some(usize::MAX),
            max_recursion_depth: Some(usize::MAX),
        }
    }

    #[test]
    fn ffi_resource_limits_are_clamped_to_hard_ceilings() {
        let limits = build_limits(&request_with_limits());

        assert_eq!(
            limits.max_duration,
            Some(Duration::from_millis(HARD_MAX_DURATION_MS))
        );
        assert_eq!(limits.max_memory, Some(HARD_MAX_MEMORY_BYTES));
        assert_eq!(limits.max_allocations, Some(HARD_MAX_ALLOCATIONS));
        assert_eq!(limits.max_recursion_depth, Some(HARD_MAX_RECURSION_DEPTH));
    }
}
