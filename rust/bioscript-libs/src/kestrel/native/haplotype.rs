use std::{cmp::Reverse, collections::BTreeSet};

use crate::{LibError, LibResult};

use super::{active_region::ActiveRegion, engine::HaplotypeEvidence, kmer::KmerCountMap};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HaplotypeAssemblyConfig {
    pub min_kmer_count: u32,
    pub max_haplotypes: usize,
    pub max_bases: usize,
    pub max_repeat_count: usize,
    pub max_saved_states: usize,
    pub locus_depth: u32,
}

impl Default for HaplotypeAssemblyConfig {
    fn default() -> Self {
        Self {
            min_kmer_count: 1,
            max_haplotypes: 40,
            max_bases: 500,
            max_repeat_count: 0,
            max_saved_states: 40,
            locus_depth: 1,
        }
    }
}

pub fn assemble_haplotypes(
    active_region: &ActiveRegion,
    counts: &KmerCountMap,
    config: &HaplotypeAssemblyConfig,
) -> LibResult<Vec<HaplotypeEvidence>> {
    validate_config(config)?;
    let Some(left_anchor) = active_region.left_end_kmer.as_deref() else {
        return Ok(Vec::new());
    };
    let Some(right_anchor) = active_region.right_end_kmer.as_deref() else {
        return Ok(Vec::new());
    };
    if left_anchor.len() != counts.kmer_size() || right_anchor.len() != counts.kmer_size() {
        return Err(LibError::InvalidArguments(
            "Kestrel haplotype anchors must match k-mer size".to_owned(),
        ));
    }

    let mut stack = vec![AssemblyState {
        sequence: left_anchor.to_owned(),
        min_depth: counts.get(left_anchor)?,
        seen_kmers: BTreeSet::from([left_anchor.to_owned()]),
        repeat_count: 0,
    }];
    let mut haplotypes = Vec::new();

    while let Some(state) = stack.pop() {
        let current_kmer = &state.sequence[state.sequence.len() - counts.kmer_size()..];
        if state.sequence.len() > counts.kmer_size() && current_kmer == right_anchor {
            haplotypes.push(HaplotypeEvidence {
                sequence: state.sequence,
                variant_depth: state.min_depth,
                locus_depth: config.locus_depth.max(state.min_depth),
            });
            if haplotypes.len() == config.max_haplotypes {
                break;
            }
            continue;
        }
        if state.sequence.len() >= config.max_bases {
            continue;
        }

        let mut next = next_states(&state, current_kmer, counts, config.min_kmer_count)?;
        next.retain(|candidate| candidate.repeat_count <= config.max_repeat_count);
        next.sort_by_key(|candidate| Reverse(candidate.min_depth));
        stack.extend(next.into_iter().rev());
        trim_saved_states(&mut stack, config.max_saved_states);
    }

    Ok(haplotypes)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AssemblyState {
    sequence: String,
    min_depth: u32,
    seen_kmers: BTreeSet<String>,
    repeat_count: usize,
}

fn next_states(
    state: &AssemblyState,
    current_kmer: &str,
    counts: &KmerCountMap,
    min_kmer_count: u32,
) -> LibResult<Vec<AssemblyState>> {
    let suffix_start = current_kmer.len() - counts.kmer_size() + 1;
    let suffix = &current_kmer[suffix_start..];
    let mut states = Vec::new();
    for base in ['A', 'C', 'G', 'T'] {
        let next_kmer = format!("{suffix}{base}");
        let depth = counts.get(&next_kmer)?;
        if depth < min_kmer_count {
            continue;
        }
        let mut sequence = state.sequence.clone();
        sequence.push(base);
        let mut seen_kmers = state.seen_kmers.clone();
        let is_repeat = !seen_kmers.insert(next_kmer);
        states.push(AssemblyState {
            sequence,
            min_depth: state.min_depth.min(depth),
            seen_kmers,
            repeat_count: state.repeat_count + usize::from(is_repeat),
        });
    }
    Ok(states)
}

fn trim_saved_states(stack: &mut Vec<AssemblyState>, max_saved_states: usize) {
    if stack.len() <= max_saved_states {
        return;
    }
    stack.sort_by_key(|state| Reverse(state.min_depth));
    stack.truncate(max_saved_states);
}

fn validate_config(config: &HaplotypeAssemblyConfig) -> LibResult<()> {
    if config.min_kmer_count == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel haplotype minimum k-mer count must be at least 1".to_owned(),
        ));
    }
    if config.max_haplotypes == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel haplotype max_haplotypes must be at least 1".to_owned(),
        ));
    }
    if config.max_bases == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel haplotype max_bases must be at least 1".to_owned(),
        ));
    }
    if config.max_saved_states == 0 {
        return Err(LibError::InvalidArguments(
            "Kestrel haplotype max_saved_states must be at least 1".to_owned(),
        ));
    }
    Ok(())
}
