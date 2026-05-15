//! Faithful port of upstream VNtyper `motif_correction_and_annotation`.
//!
//! This is a whole-set operation (left/right split by position,
//! frameshift/depth-priority dedupe per genomic locus, the legacy GG
//! `.any()` guard, then the exclude lists). Keeping it in its own module
//! keeps `vntyper.rs` focused and under the production line limit.

use std::collections::{HashMap, HashSet};

use super::VcfRecord;

pub(super) const MOTIF_POSITION_THRESHOLD: i64 = 60;
pub(super) const EXCLUDE_MOTIFS_RIGHT: &[&str] = &["8", "9", "7", "6p", "6"];
pub(super) const ALT_FOR_MOTIF_RIGHT_GG: &str = "GG";
pub(super) const MOTIFS_FOR_ALT_GG: &[&str] = &[];
pub(super) const EXCLUDE_ALTS_COMBINED: &[&str] = &["CCGCC", "CGGCG", "CGGCC"];
pub(super) const EXCLUDE_MOTIFS_COMBINED: &[&str] = &["6", "6p", "7"];

pub(super) struct MotifCorrection {
    pub(super) surviving: HashSet<usize>,
    pub(super) motif_by_index: HashMap<usize, String>,
}

fn row_pos(row: &VcfRecord) -> i64 {
    row.get("POS")
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0) as i64
}

fn row_depth_score(row: &VcfRecord) -> f64 {
    row.get("Depth_Score")
        .and_then(|value| value.parse::<f64>().ok())
        .unwrap_or(0.0)
}

fn row_is_valid_frameshift(row: &VcfRecord) -> bool {
    row.get("is_valid_frameshift").map(String::as_str) == Some("True")
}

/// Mirror upstream's pandas `str.contains(r"\bGG\b")`: with `[ACGT]+` ALTs the
/// only match is an ALT that is exactly `GG`.
fn gg_word_match(alt: &str, gg: &str) -> bool {
    alt == gg
}

/// Sort by (is_valid_frameshift, Depth_Score, POS) descending, then keep the
/// first row per (POS, REF, ALT) locus. Port of upstream
/// `_prioritize_frameshift_and_dedupe`.
fn prioritize_and_dedupe(rows: &[VcfRecord], mut items: Vec<usize>) -> Vec<usize> {
    items.sort_by(|&a, &b| {
        let fa = i32::from(row_is_valid_frameshift(&rows[a]));
        let fb = i32::from(row_is_valid_frameshift(&rows[b]));
        fb.cmp(&fa)
            .then(
                row_depth_score(&rows[b])
                    .partial_cmp(&row_depth_score(&rows[a]))
                    .unwrap_or(std::cmp::Ordering::Equal),
            )
            .then(row_pos(&rows[b]).cmp(&row_pos(&rows[a])))
    });
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for idx in items {
        let key = (
            row_pos(&rows[idx]),
            rows[idx].get("REF").cloned().unwrap_or_default(),
            rows[idx].get("ALT").cloned().unwrap_or_default(),
        );
        if seen.insert(key) {
            out.push(idx);
        }
    }
    out
}

pub(super) fn motif_correction(rows: &[VcfRecord]) -> MotifCorrection {
    let motifs_of = |row: &VcfRecord| -> String {
        row.get("Motifs")
            .or_else(|| row.get("CHROM"))
            .cloned()
            .unwrap_or_default()
    };

    // Upstream guard: every Motifs must contain exactly one dash, otherwise
    // the split fails and nothing passes (empty combined set).
    let max_dash = rows
        .iter()
        .map(|row| motifs_of(row).matches('-').count())
        .max()
        .unwrap_or(0);

    let mut motif_by_index = HashMap::new();
    let mut surviving = HashSet::new();
    if rows.is_empty() || max_dash != 1 {
        return MotifCorrection {
            surviving,
            motif_by_index,
        };
    }

    let mut motif_left = Vec::new();
    let mut motif_right = Vec::new();
    for (idx, row) in rows.iter().enumerate() {
        let motifs = motifs_of(row);
        let parts: Vec<&str> = motifs.split('-').collect();
        if parts.len() != 2 {
            continue;
        }
        let pos = row_pos(row);
        if pos < MOTIF_POSITION_THRESHOLD {
            // left motif: Motif = right token
            motif_by_index.insert(idx, parts[1].to_owned());
            motif_left.push(idx);
        } else {
            // right motif: Motif = left token
            motif_by_index.insert(idx, parts[0].to_owned());
            motif_right.push(idx);
        }
    }

    let motif_left = prioritize_and_dedupe(rows, motif_left);

    // Legacy GG branch (use_uniform_filtering = false).
    let has_gg = motif_right.iter().any(|&idx| {
        gg_word_match(
            rows[idx].get("ALT").map(String::as_str).unwrap_or_default(),
            ALT_FOR_MOTIF_RIGHT_GG,
        )
    });
    let motif_right = if has_gg {
        let kept: Vec<usize> = motif_right
            .into_iter()
            .filter(|&idx| {
                !EXCLUDE_MOTIFS_RIGHT
                    .contains(&motif_by_index.get(&idx).map(String::as_str).unwrap_or(""))
            })
            .collect();
        let kept = prioritize_and_dedupe(rows, kept);
        let any_allowed = kept.iter().any(|&idx| {
            MOTIFS_FOR_ALT_GG.contains(&motif_by_index.get(&idx).map(String::as_str).unwrap_or(""))
        });
        if any_allowed {
            kept.into_iter()
                .filter(|&idx| {
                    MOTIFS_FOR_ALT_GG
                        .contains(&motif_by_index.get(&idx).map(String::as_str).unwrap_or(""))
                })
                .collect()
        } else {
            kept
        }
    } else {
        motif_right
    };

    for idx in motif_right.into_iter().chain(motif_left.into_iter()) {
        let alt = rows[idx].get("ALT").map(String::as_str).unwrap_or_default();
        let motif = motif_by_index.get(&idx).map(String::as_str).unwrap_or("");
        if EXCLUDE_ALTS_COMBINED.contains(&alt) {
            continue;
        }
        if EXCLUDE_MOTIFS_COMBINED.contains(&motif) {
            continue;
        }
        surviving.insert(idx);
    }

    MotifCorrection {
        surviving,
        motif_by_index,
    }
}
