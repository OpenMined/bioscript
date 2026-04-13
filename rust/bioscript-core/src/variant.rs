#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomicLocus {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Assembly {
    Grch37,
    Grch38,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantKind {
    Snp,
    Insertion,
    Deletion,
    Indel,
    Other,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VariantSpec {
    pub rsids: Vec<String>,
    pub grch37: Option<GenomicLocus>,
    pub grch38: Option<GenomicLocus>,
    pub reference: Option<String>,
    pub alternate: Option<String>,
    pub kind: Option<VariantKind>,
    pub deletion_length: Option<usize>,
    pub motifs: Vec<String>,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VariantObservation {
    pub backend: String,
    pub matched_rsid: Option<String>,
    pub assembly: Option<Assembly>,
    pub genotype: Option<String>,
    pub ref_count: Option<u32>,
    pub alt_count: Option<u32>,
    pub depth: Option<u32>,
    pub evidence: Vec<String>,
}

impl VariantSpec {
    #[must_use]
    pub fn has_rsids(&self) -> bool {
        !self.rsids.is_empty()
    }

    #[must_use]
    pub fn has_coordinates(&self) -> bool {
        self.grch37.is_some() || self.grch38.is_some()
    }
}
