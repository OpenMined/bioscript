use std::collections::BTreeSet;

use serde_yaml::{Mapping, Value};

use super::{
    Issue, Severity,
    common::{
        i64_at_mapping, mapping_at, require_path, scalar_at, seq_of_strings,
        validate_optional_strings, validate_schema_and_identity, validate_tags,
        validate_url_string, value_at,
    },
};

pub(crate) fn validate_variant_root(root: &Value, issues: &mut Vec<Issue>) {
    validate_schema_and_identity(
        root,
        "bioscript:variant:1.0",
        Some("bioscript:variant"),
        issues,
    );
    validate_optional_strings(root, &["name", "label", "gene", "summary"], issues);
    validate_tags(root, issues);
    validate_identifiers(root, issues);
    validate_coordinates(root, issues);
    validate_alleles(root, issues);
    validate_findings(root, issues);
    validate_provenance(root, issues);

    let has_identifiers = value_at(root, &["identifiers"])
        .and_then(Value::as_mapping)
        .is_some_and(|mapping| !mapping.is_empty());
    let has_coordinates = ["grch37", "grch38"]
        .iter()
        .any(|assembly| value_at(root, &["coordinates", assembly]).is_some());
    if !has_identifiers && !has_coordinates {
        issues.push(Issue {
            severity: Severity::Error,
            path: "identifiers/coordinates".to_owned(),
            message: "expected at least one identifier block or one coordinate block".to_owned(),
        });
    }
}

pub(crate) fn validate_identifiers(root: &Value, issues: &mut Vec<Issue>) {
    for field in ["rsids", "aliases"] {
        let Some(values) = value_at(root, &["identifiers", field]) else {
            continue;
        };
        let Some(items) = values.as_sequence() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("identifiers.{field}"),
                message: "expected a sequence of strings".to_owned(),
            });
            continue;
        };
        let mut seen = BTreeSet::new();
        for (idx, item) in items.iter().enumerate() {
            let Some(value) = item.as_str() else {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: "expected string".to_owned(),
                });
                continue;
            };
            if !is_rsid(value) {
                issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: format!("expected rsid like rs123, found '{value}'"),
                });
            }
            if !seen.insert(value.to_owned()) {
                issues.push(Issue {
                    severity: Severity::Warning,
                    path: format!("identifiers.{field}[{idx}]"),
                    message: format!("duplicate identifier '{value}'"),
                });
            }
        }
    }
}

pub(crate) fn validate_coordinates(root: &Value, issues: &mut Vec<Issue>) {
    for assembly in ["grch37", "grch38"] {
        let Some(coord) = mapping_at(root, &["coordinates", assembly]) else {
            continue;
        };

        let Some(chrom) = coord
            .get(Value::String("chrom".to_owned()))
            .and_then(Value::as_str)
        else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.chrom"),
                message: "missing chrom".to_owned(),
            });
            continue;
        };
        if !is_allowed_chromosome(chrom) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.chrom"),
                message: format!("invalid chromosome '{chrom}'; expected 1-22, X, Y, or MT"),
            });
        }

        let has_pos = coord.contains_key(Value::String("pos".to_owned()));
        let has_start = coord.contains_key(Value::String("start".to_owned()));
        let has_end = coord.contains_key(Value::String("end".to_owned()));
        if has_pos && (has_start || has_end) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}"),
                message: "use either pos or start/end, not both".to_owned(),
            });
            continue;
        }
        if !(has_pos || has_start && has_end) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}"),
                message: "expected either pos or start/end".to_owned(),
            });
            continue;
        }

        if has_pos {
            validate_coordinate_pos(coord, assembly, issues);
        } else {
            validate_coordinate_range(coord, assembly, issues);
        }
    }
}

fn validate_coordinate_pos(coord: &Mapping, assembly: &str, issues: &mut Vec<Issue>) {
    if let Some(pos) = i64_at_mapping(coord, "pos") {
        if pos < 1 {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("coordinates.{assembly}.pos"),
                message: "expected integer >= 1".to_owned(),
            });
        }
    } else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.pos"),
            message: "expected integer".to_owned(),
        });
    }
}

fn validate_coordinate_range(coord: &Mapping, assembly: &str, issues: &mut Vec<Issue>) {
    let start = i64_at_mapping(coord, "start");
    let end = i64_at_mapping(coord, "end");
    match (start, end) {
        (Some(start), Some(end)) => validate_coordinate_range_values(start, end, assembly, issues),
        _ => issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}"),
            message: "expected integer start/end".to_owned(),
        }),
    }
}

fn validate_coordinate_range_values(start: i64, end: i64, assembly: &str, issues: &mut Vec<Issue>) {
    if start < 1 {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.start"),
            message: "expected integer >= 1".to_owned(),
        });
    }
    if end < 1 {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.end"),
            message: "expected integer >= 1".to_owned(),
        });
    }
    if end < start {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("coordinates.{assembly}.end"),
            message: "expected end >= start".to_owned(),
        });
    }
    if start == end {
        issues.push(Issue {
            severity: Severity::Warning,
            path: format!("coordinates.{assembly}"),
            message: "single-position coordinate uses start/end; prefer pos".to_owned(),
        });
    }
}

pub(crate) fn validate_alleles(root: &Value, issues: &mut Vec<Issue>) {
    require_path(root, &["alleles"], issues);
    require_path(root, &["alleles", "kind"], issues);
    require_path(root, &["alleles", "ref"], issues);
    require_path(root, &["alleles", "alts"], issues);

    let Some(kind) = scalar_at(root, &["alleles", "kind"]) else {
        return;
    };
    if !matches!(kind.as_str(), "snv" | "deletion" | "insertion" | "indel") {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.kind".to_owned(),
            message: "expected one of snv, deletion, insertion, indel".to_owned(),
        });
    }

    if value_at(root, &["alleles", "canonical_alt"]).is_some() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.canonical_alt".to_owned(),
            message: "canonical_alt is not part of the current schema".to_owned(),
        });
    }

    let Some(reference) = scalar_at(root, &["alleles", "ref"]) else {
        return;
    };
    if reference.trim().is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.ref".to_owned(),
            message: "empty string".to_owned(),
        });
    }

    let Some(alts_value) = value_at(root, &["alleles", "alts"]) else {
        return;
    };
    let Some(alts_seq) = alts_value.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.alts".to_owned(),
            message: "expected a non-empty sequence of strings".to_owned(),
        });
        return;
    };
    if alts_seq.is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.alts".to_owned(),
            message: "expected at least one alternate allele".to_owned(),
        });
        return;
    }

    let mut alts = Vec::new();
    for (idx, item) in alts_seq.iter().enumerate() {
        let Some(alt) = item.as_str() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message: "expected string".to_owned(),
            });
            continue;
        };
        if alt.trim().is_empty() {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message: "empty string".to_owned(),
            });
            continue;
        }
        alts.push(alt.to_owned());
    }
    validate_symbolic_alleles(&reference, &alts, issues);
    validate_snv_alleles(&kind, &reference, &alts, issues);
}

fn validate_symbolic_alleles(reference: &str, alts: &[String], issues: &mut Vec<Issue>) {
    if reference == "I" || reference == "D" {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.ref".to_owned(),
            message: "symbolic I/D alleles are not allowed in stored YAML; use biological alleles"
                .to_owned(),
        });
    }
    for (idx, alt) in alts.iter().enumerate() {
        if alt == "I" || alt == "D" {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message:
                    "symbolic I/D alleles are not allowed in stored YAML; use biological alleles"
                        .to_owned(),
            });
        }
    }
}

fn validate_snv_alleles(kind: &str, reference: &str, alts: &[String], issues: &mut Vec<Issue>) {
    if kind != "snv" {
        return;
    }
    if !is_base_allele(reference) {
        issues.push(Issue {
            severity: Severity::Error,
            path: "alleles.ref".to_owned(),
            message: "snv ref must be one of A/C/G/T".to_owned(),
        });
    }
    for (idx, alt) in alts.iter().enumerate() {
        if !is_base_allele(alt) {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("alleles.alts[{idx}]"),
                message: "snv alt must be one of A/C/G/T".to_owned(),
            });
        }
    }
}

fn validate_findings(root: &Value, issues: &mut Vec<Issue>) {
    let alts = seq_of_strings(root, &["alleles", "alts"]).unwrap_or_default();
    let Some(findings) = value_at(root, &["findings"]).and_then(Value::as_sequence) else {
        return;
    };

    for (idx, finding) in findings.iter().enumerate() {
        let Some(mapping) = finding.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };

        let Some(schema) = mapping
            .get(Value::String("schema".to_owned()))
            .and_then(Value::as_str)
        else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].schema"),
                message: "missing schema".to_owned(),
            });
            continue;
        };
        if schema.trim().is_empty() {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].schema"),
                message: "empty string".to_owned(),
            });
        }
        if let Some(alt) = mapping
            .get(Value::String("alt".to_owned()))
            .and_then(Value::as_str)
            && alt != "*"
            && !alts.iter().any(|item| item == alt)
        {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("findings[{idx}].alt"),
                message: format!("finding alt '{alt}' is not present in alleles.alts {alts:?}"),
            });
        }
        let has_summary = mapping
            .get(Value::String("summary".to_owned()))
            .and_then(Value::as_str)
            .is_some_and(|value| !value.trim().is_empty());
        let has_notes = mapping
            .get(Value::String("notes".to_owned()))
            .and_then(Value::as_str)
            .is_some_and(|value| !value.trim().is_empty());
        if !has_summary && !has_notes {
            issues.push(Issue {
                severity: Severity::Warning,
                path: format!("findings[{idx}]"),
                message: "finding has neither summary nor notes".to_owned(),
            });
        }
    }
}

fn validate_provenance(root: &Value, issues: &mut Vec<Issue>) {
    let Some(sources) = value_at(root, &["provenance", "sources"]).and_then(Value::as_sequence)
    else {
        return;
    };
    for (idx, source) in sources.iter().enumerate() {
        let Some(mapping) = source.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("provenance.sources[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        for field in ["kind", "label", "url"] {
            match mapping
                .get(Value::String(field.to_owned()))
                .and_then(Value::as_str)
            {
                Some(text) if !text.trim().is_empty() => {}
                Some(_) => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("provenance.sources[{idx}].{field}"),
                    message: "empty string".to_owned(),
                }),
                None => issues.push(Issue {
                    severity: Severity::Error,
                    path: format!("provenance.sources[{idx}].{field}"),
                    message: "missing required field".to_owned(),
                }),
            }
        }
        if let Some(url) = mapping
            .get(Value::String("url".to_owned()))
            .and_then(Value::as_str)
        {
            validate_url_string(
                url,
                &format!("provenance.sources[{idx}].url"),
                false,
                issues,
            );
        }
    }
}

fn is_allowed_chromosome(value: &str) -> bool {
    matches!(value, "X" | "Y" | "MT")
        || value
            .parse::<u8>()
            .is_ok_and(|chrom| (1..=22).contains(&chrom))
}

fn is_base_allele(value: &str) -> bool {
    matches!(value, "A" | "C" | "G" | "T")
}

fn is_rsid(value: &str) -> bool {
    value.starts_with("rs") && value[2..].chars().all(|ch| ch.is_ascii_digit())
}
