use crate::{LibError, LibResult};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AlignmentWeight {
    pub match_weight: f32,
    pub mismatch: f32,
    pub gap_open: f32,
    pub gap_extend: f32,
    pub init_score: f32,
}

impl AlignmentWeight {
    pub const DEFAULT_MATCH: f32 = 10.0;
    pub const DEFAULT_MISMATCH: f32 = -10.0;
    pub const DEFAULT_GAP_OPEN: f32 = -40.0;
    pub const DEFAULT_GAP_EXTEND: f32 = -4.0;
    pub const DEFAULT_INIT_SCORE: f32 = 0.0;

    pub fn new(
        match_weight: f32,
        mismatch: f32,
        gap_open: f32,
        gap_extend: f32,
        init_score: f32,
    ) -> LibResult<Self> {
        Ok(Self {
            match_weight: normalize_nonzero_positive("matching bases", match_weight)?,
            mismatch: normalize_nonzero_negative("mismatched bases", mismatch)?,
            gap_open: -gap_open.abs(),
            gap_extend: normalize_nonzero_negative("gap extension", gap_extend)?,
            init_score: init_score.abs(),
        })
    }

    pub fn initial_score(&self, kmer_size: usize) -> LibResult<f32> {
        if kmer_size == 0 {
            return Err(LibError::InvalidArguments(
                "Kestrel alignment weight requires k-mer size at least 1".to_owned(),
            ));
        }
        if is_zero(self.init_score) {
            return Ok(self.match_weight * kmer_size as f32);
        }
        Ok(self.init_score)
    }

    pub fn max_exclusive_gap_size(&self, kmer_size: usize) -> LibResult<usize> {
        let init_score = self.initial_score(kmer_size)? as i32 as f32;
        if init_score > self.gap_open {
            return Ok(((init_score + self.gap_open) / -self.gap_extend) as usize);
        }
        Ok(0)
    }
}

impl Default for AlignmentWeight {
    fn default() -> Self {
        Self {
            match_weight: Self::DEFAULT_MATCH,
            mismatch: Self::DEFAULT_MISMATCH,
            gap_open: Self::DEFAULT_GAP_OPEN,
            gap_extend: Self::DEFAULT_GAP_EXTEND,
            init_score: Self::DEFAULT_INIT_SCORE,
        }
    }
}

fn normalize_nonzero_positive(label: &str, value: f32) -> LibResult<f32> {
    if !value.is_finite() || is_zero(value) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel alignment weight for {label} must be finite and nonzero: {value}"
        )));
    }
    Ok(value.abs())
}

fn normalize_nonzero_negative(label: &str, value: f32) -> LibResult<f32> {
    if !value.is_finite() || is_zero(value) {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel alignment weight for {label} must be finite and nonzero: {value}"
        )));
    }
    Ok(-value.abs())
}

fn is_zero(value: f32) -> bool {
    value.abs() <= f32::EPSILON
}
