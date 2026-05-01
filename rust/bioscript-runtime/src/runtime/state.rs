use std::{
    collections::HashMap,
    sync::{
        Mutex,
        atomic::{AtomicU64, Ordering},
    },
    time::Duration,
};

use bioscript_formats::{GenotypeLoadOptions, GenotypeStore};
use monty::{MontyException, ResourceLimits};

use bioscript_core::RuntimeError;

#[derive(Debug, Clone)]
pub struct RuntimeConfig {
    pub limits: ResourceLimits,
    pub loader: GenotypeLoadOptions,
}

impl Default for RuntimeConfig {
    fn default() -> Self {
        let limits = ResourceLimits::new()
            .max_duration(Duration::from_millis(100))
            .max_memory(8 * 1024 * 1024)
            .max_allocations(200_000)
            .gc_interval(1000)
            .max_recursion_depth(Some(200));
        Self {
            limits,
            loader: GenotypeLoadOptions::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageTiming {
    pub stage: String,
    pub duration_ms: u128,
    pub detail: String,
}

pub(crate) fn monty_error(value: MontyException) -> RuntimeError {
    RuntimeError::Monty(value.to_string())
}

pub(crate) struct RuntimeState {
    pub(crate) next_handle: AtomicU64,
    pub(crate) genotype_files: Mutex<HashMap<u64, GenotypeStore>>,
    pub(crate) trace_lines: Mutex<Vec<usize>>,
    pub(crate) timings: Mutex<Vec<StageTiming>>,
}

impl RuntimeState {
    pub(crate) fn new() -> Self {
        Self {
            next_handle: AtomicU64::new(1),
            genotype_files: Mutex::new(HashMap::new()),
            trace_lines: Mutex::new(Vec::new()),
            timings: Mutex::new(Vec::new()),
        }
    }

    pub(crate) fn next_handle(&self) -> u64 {
        self.next_handle.fetch_add(1, Ordering::Relaxed)
    }
}
