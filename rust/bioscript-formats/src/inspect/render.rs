use bioscript_core::Assembly;

use super::{DetectedKind, DetectionConfidence, FileContainer, FileInspection};

impl FileInspection {
    #[must_use]
    pub fn render_text(&self) -> String {
        let mut lines = Vec::new();
        lines.push(format!("path\t{}", self.path.display()));
        lines.push(format!("container\t{}", render_container(self.container)));
        lines.push(format!("kind\t{}", render_kind(self.detected_kind)));
        lines.push(format!(
            "confidence\t{}",
            render_confidence(self.confidence)
        ));
        lines.push(format!("assembly\t{}", render_assembly(self.assembly)));
        lines.push(format!("phased\t{}", render_bool(self.phased)));
        lines.push(format!(
            "selected_entry\t{}",
            self.selected_entry.as_deref().unwrap_or("")
        ));
        lines.push(format!("has_index\t{}", render_bool(self.has_index)));
        lines.push(format!(
            "index_path\t{}",
            self.index_path
                .as_ref()
                .map(|path| path.display().to_string())
                .unwrap_or_default()
        ));
        lines.push(format!(
            "reference_matches\t{}",
            render_bool(self.reference_matches)
        ));
        if let Some(source) = &self.source {
            lines.push(format!(
                "vendor\t{}",
                source.vendor.as_deref().unwrap_or_default()
            ));
            lines.push(format!(
                "platform_version\t{}",
                source.platform_version.as_deref().unwrap_or_default()
            ));
            lines.push(format!(
                "source_confidence\t{}",
                render_confidence(source.confidence)
            ));
            lines.push(format!("source_evidence\t{}", source.evidence.join(" | ")));
        } else {
            lines.push("vendor\t".to_owned());
            lines.push("platform_version\t".to_owned());
            lines.push("source_confidence\t".to_owned());
            lines.push("source_evidence\t".to_owned());
        }
        lines.push(format!("evidence\t{}", self.evidence.join(" | ")));
        lines.push(format!("warnings\t{}", self.warnings.join(" | ")));
        lines.push(format!("duration_ms\t{}", self.duration_ms));
        lines.join("\n")
    }
}

pub(crate) fn render_container(value: FileContainer) -> &'static str {
    match value {
        FileContainer::Plain => "plain",
        FileContainer::Zip => "zip",
    }
}

pub(crate) fn render_kind(value: DetectedKind) -> &'static str {
    match value {
        DetectedKind::GenotypeText => "genotype_text",
        DetectedKind::Vcf => "vcf",
        DetectedKind::AlignmentCram => "alignment_cram",
        DetectedKind::AlignmentBam => "alignment_bam",
        DetectedKind::ReferenceFasta => "reference_fasta",
        DetectedKind::Unknown => "unknown",
    }
}

pub(crate) fn render_confidence(value: DetectionConfidence) -> &'static str {
    match value {
        DetectionConfidence::Authoritative => "authoritative",
        DetectionConfidence::StrongHeuristic => "strong_heuristic",
        DetectionConfidence::WeakHeuristic => "weak_heuristic",
        DetectionConfidence::Unknown => "unknown",
    }
}

pub(crate) fn render_assembly(value: Option<Assembly>) -> &'static str {
    match value {
        Some(Assembly::Grch37) => "grch37",
        Some(Assembly::Grch38) => "grch38",
        None => "",
    }
}

pub(crate) fn render_bool(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "true",
        Some(false) => "false",
        None => "",
    }
}
