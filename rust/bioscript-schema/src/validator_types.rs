use std::{
    collections::BTreeSet,
    fmt::{self, Write as _},
    fs,
    path::{Path, PathBuf},
};

use bioscript_core::{GenomicLocus, VariantKind, VariantSpec};
use serde_yaml::{Mapping, Value};
use url::Url;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    Error,
    Warning,
}

impl fmt::Display for Severity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Error => f.write_str("error"),
            Self::Warning => f.write_str("warning"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Issue {
    pub severity: Severity,
    pub path: String,
    pub message: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FileReport {
    pub file: PathBuf,
    pub issues: Vec<Issue>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationReport {
    pub files_scanned: usize,
    pub reports: Vec<FileReport>,
}

impl ValidationReport {
    #[must_use]
    pub fn total_issues(&self) -> usize {
        self.reports.iter().map(|report| report.issues.len()).sum()
    }

    #[must_use]
    pub fn total_errors(&self) -> usize {
        self.reports
            .iter()
            .flat_map(|report| &report.issues)
            .filter(|issue| issue.severity == Severity::Error)
            .count()
    }

    #[must_use]
    pub fn total_warnings(&self) -> usize {
        self.reports
            .iter()
            .flat_map(|report| &report.issues)
            .filter(|issue| issue.severity == Severity::Warning)
            .count()
    }

    #[must_use]
    pub fn has_errors(&self) -> bool {
        self.total_errors() > 0
    }

    #[must_use]
    pub fn render_text(&self) -> String {
        let mut out = String::new();
        let _ = write!(
            out,
            "files_scanned: {}\nerrors: {}\nwarnings: {}\n",
            self.files_scanned,
            self.total_errors(),
            self.total_warnings()
        );
        for report in &self.reports {
            out.push('\n');
            let _ = writeln!(out, "file: {}", report.file.display());
            for issue in &report.issues {
                let _ = writeln!(
                    out,
                    "  - [{}] {}: {}",
                    issue.severity, issue.path, issue.message
                );
            }
        }
        out
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VariantManifest {
    pub path: PathBuf,
    pub name: String,
    pub tags: Vec<String>,
    pub spec: VariantSpec,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelManifest {
    pub path: PathBuf,
    pub name: String,
    pub tags: Vec<String>,
    pub permissions: Permissions,
    pub downloads: Vec<Download>,
    pub members: Vec<PanelMember>,
    pub interpretations: Vec<PanelInterpretation>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssayManifest {
    pub path: PathBuf,
    pub name: String,
    pub tags: Vec<String>,
    pub members: Vec<PanelMember>,
    pub interpretations: Vec<PanelInterpretation>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Permissions {
    pub domains: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Download {
    pub id: String,
    pub url: String,
    pub origin: String,
    pub sha256: String,
    pub version: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelMember {
    pub kind: String,
    pub path: Option<String>,
    pub download: Option<String>,
    pub sha256: Option<String>,
    pub version: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelInterpretation {
    pub id: String,
    pub kind: String,
    pub path: String,
    pub output_format: Option<String>,
    pub derived_from: Vec<String>,
    pub emits: Vec<PanelInterpretationEmit>,
    pub logic: Option<PanelInterpretationLogic>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelInterpretationLogic {
    pub source: Option<PanelInterpretationLogicSource>,
    pub description: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelInterpretationLogicSource {
    pub name: Option<String>,
    pub url: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PanelInterpretationEmit {
    pub key: String,
    pub label: Option<String>,
    pub value_type: Option<String>,
    pub format: Option<String>,
}
