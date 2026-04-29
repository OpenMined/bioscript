#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportStatus {
    Complete,
    Partial,
    Failed,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReportTiming {
    pub started_at: String,
    pub finished_at: String,
    pub duration_ms: u64,
    pub breakdown: ReportTimingBreakdown,
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ReportTimingBreakdown {
    pub load_ms: Option<u64>,
    pub query_ms: Option<u64>,
    pub compute_ms: Option<u64>,
    pub write_ms: Option<u64>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ReportValue {
    String(String),
    Number(f64),
    Boolean(bool),
    Json(String),
    Null,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReportField {
    pub key: String,
    pub label: String,
    pub value: ReportValue,
    pub value_type: String,
    pub format: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReportRecord {
    pub participant_id: String,
    pub assay_id: String,
    pub assay_version: String,
    pub report_status: ReportStatus,
    pub derived_from: Vec<String>,
    pub facets: Option<String>,
    pub timing: ReportTiming,
    pub fields: Vec<ReportField>,
    pub metrics_json: Option<String>,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}
