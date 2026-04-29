#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationKind {
    Snp,
    Insertion,
    Deletion,
    Mnv,
    StructuralVariant,
    CopyNumberVariant,
    Fusion,
    Hybrid,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchStatus {
    Found,
    NotFound,
    Ambiguous,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoverageStatus {
    Covered,
    NotCovered,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CallStatus {
    Called,
    NoCall,
    LowQuality,
    Filtered,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    HomRef,
    Het,
    HomAlt,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObservationOutcome {
    Variant,
    Reference,
    NoCall,
    NotCovered,
    Unknown,
    Partial,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceType {
    VcfRecord,
    GvcfBlock,
    Mpileup,
    CramSlice,
    SplitReads,
    ReadPairs,
    Imputed,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ObservationRecord {
    pub participant_id: String,
    pub assay_id: String,
    pub assay_version: String,
    pub variant_key: String,
    pub rsid: Option<String>,
    pub assembly: String,
    pub chrom: String,
    pub pos_start: i64,
    pub pos_end: i64,
    pub reference: String,
    pub alternate: String,
    pub kind: ObservationKind,
    pub match_status: MatchStatus,
    pub coverage_status: CoverageStatus,
    pub call_status: CallStatus,
    pub genotype: String,
    pub genotype_display: String,
    pub zygosity: Zygosity,
    pub ref_count: Option<u32>,
    pub alt_count: Option<u32>,
    pub depth: Option<u32>,
    pub genotype_quality: Option<u32>,
    pub allele_balance: Option<f64>,
    pub outcome: ObservationOutcome,
    pub evidence_type: Option<EvidenceType>,
    pub evidence_raw: Option<String>,
    pub facets: Option<String>,
}

pub const OBSERVATION_TSV_HEADERS: [&str; 27] = [
    "participant_id",
    "assay_id",
    "assay_version",
    "variant_key",
    "rsid",
    "assembly",
    "chrom",
    "pos_start",
    "pos_end",
    "ref",
    "alt",
    "kind",
    "match_status",
    "coverage_status",
    "call_status",
    "genotype",
    "genotype_display",
    "zygosity",
    "ref_count",
    "alt_count",
    "depth",
    "genotype_quality",
    "allele_balance",
    "outcome",
    "evidence_type",
    "evidence_raw",
    "facets",
];
