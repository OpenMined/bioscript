fn validate_panel_members(root: &Value, allowed_kinds: &[&str], issues: &mut Vec<Issue>) {
    let Some(members) = value_at(root, &["members"]).and_then(Value::as_sequence) else {
        issues.push(Issue {
            severity: Severity::Error,
            path: "members".to_owned(),
            message: "missing required field".to_owned(),
        });
        return;
    };
    if members.is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: "members".to_owned(),
            message: "expected at least one member".to_owned(),
        });
        return;
    }

    let download_ids = panel_download_ids(root);

    for (idx, item) in members.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("members[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        validate_panel_member(idx, mapping, allowed_kinds, &download_ids, issues);
    }
}

fn panel_download_ids(root: &Value) -> BTreeSet<String> {
    value_at(root, &["downloads"])
        .and_then(Value::as_sequence)
        .into_iter()
        .flatten()
        .filter_map(|item| {
            item.as_mapping()?
                .get(Value::String("id".to_owned()))
                .and_then(Value::as_str)
                .map(ToOwned::to_owned)
        })
        .collect()
}

fn validate_panel_member(
    idx: usize,
    mapping: &Mapping,
    allowed_kinds: &[&str],
    download_ids: &BTreeSet<String>,
    issues: &mut Vec<Issue>,
) {
    let kind = mapping
        .get(Value::String("kind".to_owned()))
        .and_then(Value::as_str);
    match kind {
        Some(kind) if allowed_kinds.contains(&kind) => {}
        Some(other) => issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].kind"),
            message: format!("unsupported member kind '{other}'"),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].kind"),
            message: "missing required field".to_owned(),
        }),
    }

    let path_value = mapping
        .get(Value::String("path".to_owned()))
        .and_then(Value::as_str);
    let download_value = mapping
        .get(Value::String("download".to_owned()))
        .and_then(Value::as_str);
    if path_value.is_some() == download_value.is_some() {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}]"),
            message: "expected exactly one of path or download".to_owned(),
        });
    }
    validate_panel_member_path(idx, path_value, issues);
    validate_panel_member_download(idx, download_value, download_ids, issues);
    validate_panel_member_metadata(idx, mapping, issues);
}

fn validate_panel_member_path(idx: usize, path_value: Option<&str>, issues: &mut Vec<Issue>) {
    if let Some(path) = path_value
        && path.trim().is_empty()
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].path"),
            message: "empty string".to_owned(),
        });
    }
}

fn validate_panel_member_download(
    idx: usize,
    download_value: Option<&str>,
    download_ids: &BTreeSet<String>,
    issues: &mut Vec<Issue>,
) {
    let Some(download) = download_value else {
        return;
    };
    if download.trim().is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].download"),
            message: "empty string".to_owned(),
        });
    } else if !download_ids.contains(download) {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].download"),
            message: format!("unknown download id '{download}'"),
        });
    }
}

fn validate_panel_member_metadata(idx: usize, mapping: &Mapping, issues: &mut Vec<Issue>) {
    if let Some(version) = mapping
        .get(Value::String("version".to_owned()))
        .and_then(Value::as_str)
        && version.trim().is_empty()
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].version"),
            message: "empty string".to_owned(),
        });
    }
    if let Some(sha) = mapping
        .get(Value::String("sha256".to_owned()))
        .and_then(Value::as_str)
        && !is_sha256(sha)
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("members[{idx}].sha256"),
            message: "expected 64 lowercase hex characters".to_owned(),
        });
    }
}

fn validate_panel_interpretations(root: &Value, issues: &mut Vec<Issue>) {
    if value_at(root, &["analyses"]).is_some() && value_at(root, &["interpretations"]).is_some() {
        issues.push(Issue {
            severity: Severity::Warning,
            path: "interpretations".to_owned(),
            message: "use analyses instead of interpretations; do not define both".to_owned(),
        });
    }
    let key = if value_at(root, &["analyses"]).is_some() {
        "analyses"
    } else {
        "interpretations"
    };
    let Some(items) = value_at(root, &[key]) else {
        return;
    };
    let Some(items) = items.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: key.to_owned(),
            message: "expected a sequence of mappings".to_owned(),
        });
        return;
    };
    for (idx, item) in items.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("{key}[{idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        validate_panel_interpretation(key, idx, mapping, issues);
    }
}

fn validate_panel_interpretation(
    key: &str,
    idx: usize,
    mapping: &Mapping,
    issues: &mut Vec<Issue>,
) {
    for field in ["id", "kind", "path"] {
        validate_required_mapping_string(mapping, field, &format!("{key}[{idx}]"), issues);
    }
    if let Some(kind) = mapping
        .get(Value::String("kind".to_owned()))
        .and_then(Value::as_str)
        && kind != "bioscript"
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].kind"),
            message: "expected 'bioscript'".to_owned(),
        });
    }
    if let Some(output_format) = mapping
        .get(Value::String("output_format".to_owned()))
        .and_then(Value::as_str)
        && !matches!(output_format, "tsv" | "json" | "jsonl")
    {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].output_format"),
            message: "expected 'tsv', 'json', or 'jsonl'".to_owned(),
        });
    }
    let Some(derived_from) = mapping
        .get(Value::String("derived_from".to_owned()))
        .and_then(Value::as_sequence)
    else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].derived_from"),
            message: "expected a non-empty sequence of strings".to_owned(),
        });
        return;
    };
    if derived_from.is_empty() {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].derived_from"),
            message: "expected at least one source variant".to_owned(),
        });
    }
    for (source_idx, source) in derived_from.iter().enumerate() {
        match source.as_str() {
            Some(value) if !value.trim().is_empty() => {}
            Some(_) => issues.push(Issue {
                severity: Severity::Error,
                path: format!("{key}[{idx}].derived_from[{source_idx}]"),
                message: "empty string".to_owned(),
            }),
            None => issues.push(Issue {
                severity: Severity::Error,
                path: format!("{key}[{idx}].derived_from[{source_idx}]"),
                message: "expected string".to_owned(),
            }),
        }
    }
    validate_panel_interpretation_logic(key, idx, mapping, issues);
    validate_panel_interpretation_emits(key, idx, mapping, issues);
}

fn validate_panel_interpretation_logic(
    key: &str,
    idx: usize,
    mapping: &Mapping,
    issues: &mut Vec<Issue>,
) {
    let Some(logic) = mapping.get(Value::String("logic".to_owned())) else {
        return;
    };
    let Some(logic) = logic.as_mapping() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].logic"),
            message: "expected mapping".to_owned(),
        });
        return;
    };
    validate_optional_mapping_string(logic, "description", &format!("{key}[{idx}].logic"), issues);
    let Some(source) = logic.get(Value::String("source".to_owned())) else {
        return;
    };
    let Some(source) = source.as_mapping() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].logic.source"),
            message: "expected mapping".to_owned(),
        });
        return;
    };
    validate_optional_mapping_string(
        source,
        "name",
        &format!("{key}[{idx}].logic.source"),
        issues,
    );
    validate_optional_mapping_string(source, "url", &format!("{key}[{idx}].logic.source"), issues);
    if let Some(url) = source
        .get(Value::String("url".to_owned()))
        .and_then(Value::as_str)
    {
        validate_url_string(
            url,
            &format!("{key}[{idx}].logic.source.url"),
            false,
            issues,
        );
    }
}

fn validate_panel_interpretation_emits(
    key: &str,
    idx: usize,
    mapping: &Mapping,
    issues: &mut Vec<Issue>,
) {
    let Some(emits) = mapping.get(Value::String("emits".to_owned())) else {
        return;
    };
    let Some(emits) = emits.as_sequence() else {
        issues.push(Issue {
            severity: Severity::Error,
            path: format!("{key}[{idx}].emits"),
            message: "expected a sequence of mappings".to_owned(),
        });
        return;
    };
    for (emit_idx, item) in emits.iter().enumerate() {
        let Some(mapping) = item.as_mapping() else {
            issues.push(Issue {
                severity: Severity::Error,
                path: format!("{key}[{idx}].emits[{emit_idx}]"),
                message: "expected mapping".to_owned(),
            });
            continue;
        };
        validate_required_mapping_string(
            mapping,
            "key",
            &format!("{key}[{idx}].emits[{emit_idx}]"),
            issues,
        );
        for field in ["label", "value_type", "format"] {
            validate_optional_mapping_string(
                mapping,
                field,
                &format!("{key}[{idx}].emits[{emit_idx}]"),
                issues,
            );
        }
    }
}

fn validate_required_mapping_string(
    mapping: &Mapping,
    field: &str,
    parent: &str,
    issues: &mut Vec<Issue>,
) {
    match mapping
        .get(Value::String(field.to_owned()))
        .and_then(Value::as_str)
    {
        Some(value) if !value.trim().is_empty() => {}
        Some(_) => issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.{field}"),
            message: "empty string".to_owned(),
        }),
        None => issues.push(Issue {
            severity: Severity::Error,
            path: format!("{parent}.{field}"),
            message: "missing required field".to_owned(),
        }),
    }
}

fn validate_optional_mapping_string(
    mapping: &Mapping,
    field: &str,
    parent: &str,
    issues: &mut Vec<Issue>,
) {
    if let Some(value) = mapping.get(Value::String(field.to_owned())) {
        match value.as_str() {
            Some(text) if !text.trim().is_empty() => {}
            Some(_) => issues.push(Issue {
                severity: Severity::Warning,
                path: format!("{parent}.{field}"),
                message: "empty string".to_owned(),
            }),
            None => issues.push(Issue {
                severity: Severity::Error,
                path: format!("{parent}.{field}"),
                message: "expected string".to_owned(),
            }),
        }
    }
}

fn variant_spec_from_root(root: &Value) -> Result<VariantSpec, String> {
    let rsids = seq_of_strings(root, &["identifiers", "rsids"]).unwrap_or_default();
    let grch37 = locus_from_root(root, "grch37")?;
    let grch38 = locus_from_root(root, "grch38")?;
    let reference = scalar_at(root, &["alleles", "ref"]);
    let alternate = seq_of_strings(root, &["alleles", "alts"]).and_then(|alts| alts.first().cloned());
    let deletion_length = value_at(root, &["alleles", "deletion_length"])
        .and_then(Value::as_u64)
        .and_then(|value| usize::try_from(value).ok());
    let motifs = seq_of_strings(root, &["alleles", "motifs"]).unwrap_or_default();
    let kind = scalar_at(root, &["alleles", "kind"]).map(|kind| match kind.as_str() {
        "snv" => VariantKind::Snp,
        "deletion" => VariantKind::Deletion,
        "insertion" => VariantKind::Insertion,
        "indel" => VariantKind::Indel,
        _ => VariantKind::Other,
    });

    Ok(VariantSpec {
        rsids,
        grch37,
        grch38,
        reference,
        alternate,
        kind,
        deletion_length,
        motifs,
    })
}

fn locus_from_root(root: &Value, assembly: &str) -> Result<Option<GenomicLocus>, String> {
    let Some(mapping) = mapping_at(root, &["coordinates", assembly]) else {
        return Ok(None);
    };
    let chrom = mapping
        .get(Value::String("chrom".to_owned()))
        .and_then(Value::as_str)
        .ok_or_else(|| format!("coordinates.{assembly}.chrom missing"))?;
    let (start, end) = if let Some(pos) = i64_at_mapping(mapping, "pos") {
        (pos, pos)
    } else {
        let start = i64_at_mapping(mapping, "start")
            .ok_or_else(|| format!("coordinates.{assembly}.start missing"))?;
        let end = i64_at_mapping(mapping, "end")
            .ok_or_else(|| format!("coordinates.{assembly}.end missing"))?;
        (start, end)
    };
    Ok(Some(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    }))
}
