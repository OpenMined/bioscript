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

    pub fn parse(weight_string: Option<&str>) -> LibResult<Self> {
        let Some(mut value) = weight_string.map(str::trim) else {
            return Ok(Self::default());
        };
        if value.is_empty() {
            return Ok(Self::default());
        }
        value = strip_matching_bounds(value)?;

        let tokens: Vec<&str> = value.split(',').map(str::trim).collect();
        if tokens.len() > 5 {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel alignment weight vector has more than 5 comma-separated values: {}",
                tokens.len()
            )));
        }

        let mut weights = Self::default();
        if let Some(token) = tokens.first().filter(|token| !token.is_empty()) {
            weights.match_weight =
                normalize_nonzero_positive("matching bases", parse_number(token)?)?;
        }
        if let Some(token) = tokens.get(1).filter(|token| !token.is_empty()) {
            weights.mismatch =
                normalize_nonzero_negative("mismatched bases", parse_number(token)?)?;
        }
        if let Some(token) = tokens.get(2).filter(|token| !token.is_empty()) {
            weights.gap_open = -parse_number(token)?.abs();
        }
        if let Some(token) = tokens.get(3).filter(|token| !token.is_empty()) {
            weights.gap_extend = normalize_nonzero_negative("gap extension", parse_number(token)?)?;
        }
        if let Some(token) = tokens.get(4).filter(|token| !token.is_empty()) {
            weights.init_score = parse_number(token)?.abs();
        }
        Ok(weights)
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

fn strip_matching_bounds(value: &str) -> LibResult<&str> {
    let mut chars = value.chars();
    let Some(first) = chars.next() else {
        return Ok(value);
    };
    let Some(last) = value.chars().next_back() else {
        return Ok(value);
    };

    let expected = match first {
        '(' => Some(')'),
        '<' => Some('>'),
        '[' => Some(']'),
        '{' => Some('}'),
        _ => None,
    };
    if let Some(expected) = expected {
        if last != expected {
            return Err(LibError::InvalidArguments(format!(
                "Kestrel alignment weight vector has mismatched bounds: {value}"
            )));
        }
        return Ok(&value[first.len_utf8()..value.len() - last.len_utf8()]);
    }
    if matches!(last, ')' | '>' | ']' | '}') {
        return Err(LibError::InvalidArguments(format!(
            "Kestrel alignment weight vector has a closing bound without an opening bound: {value}"
        )));
    }
    Ok(value)
}

fn parse_number(value: &str) -> LibResult<f32> {
    value
        .parse::<f32>()
        .or_else(|_| parse_java_integer(value).map(|number| number as f32))
        .map_err(|_| {
            LibError::InvalidArguments(format!(
                "Kestrel alignment weight is not a valid number: {value}"
            ))
        })
}

fn parse_java_integer(value: &str) -> Result<i32, std::num::ParseIntError> {
    let (negative, unsigned) = value
        .strip_prefix('-')
        .map(|value| (true, value))
        .or_else(|| value.strip_prefix('+').map(|value| (false, value)))
        .unwrap_or((false, value));
    let (radix, digits) = if let Some(digits) = unsigned
        .strip_prefix("0x")
        .or_else(|| unsigned.strip_prefix("0X"))
    {
        (16, digits)
    } else if let Some(digits) = unsigned.strip_prefix('#') {
        (16, digits)
    } else if unsigned.len() > 1 && unsigned.starts_with('0') {
        (8, &unsigned[1..])
    } else {
        (10, unsigned)
    };
    let parsed = i32::from_str_radix(digits, radix)?;
    Ok(if negative { -parsed } else { parsed })
}
