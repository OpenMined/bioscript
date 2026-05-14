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
        grch37_assembly_ref: scalar_at(root, &["coordinates", "grch37", "assembly_ref"]),
        grch38_assembly_ref: scalar_at(root, &["coordinates", "grch38", "assembly_ref"]),
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

#[cfg(test)]
mod panel_validator_tests {
    use super::*;

    fn yaml(text: &str) -> Value {
        serde_yaml::from_str(text).unwrap()
    }

    fn issue_paths(issues: &[Issue]) -> Vec<String> {
        issues.iter().map(|issue| issue.path.clone()).collect()
    }

    #[test]
    fn panel_member_validation_reports_missing_empty_and_unknown_references() {
        let root = yaml(
            r#"
downloads:
  - id: known
members:
  - not-a-mapping
  - kind: unknown
  - kind: variant
    path: ""
    download: known
    version: ""
    sha256: BAD
  - kind: variant
    download: missing
  - kind: variant
    download: ""
"#,
        );
        let mut issues = Vec::new();
        validate_panel_members(&root, &["variant"], &mut issues);
        let paths = issue_paths(&issues);
        assert!(paths.contains(&"members[0]".to_owned()));
        assert!(paths.contains(&"members[1].kind".to_owned()));
        assert!(paths.contains(&"members[1]".to_owned()));
        assert!(paths.contains(&"members[2]".to_owned()));
        assert!(paths.contains(&"members[2].path".to_owned()));
        assert!(paths.contains(&"members[2].version".to_owned()));
        assert!(paths.contains(&"members[2].sha256".to_owned()));
        assert!(paths.contains(&"members[3].download".to_owned()));
        assert!(paths.contains(&"members[4].download".to_owned()));
    }

    #[test]
    fn panel_member_validation_reports_missing_and_empty_lists() {
        let mut issues = Vec::new();
        validate_panel_members(&yaml("{}"), &["variant"], &mut issues);
        assert_eq!(issues[0].path, "members");
        assert_eq!(issues[0].message, "missing required field");

        let mut issues = Vec::new();
        validate_panel_members(&yaml("members: []"), &["variant"], &mut issues);
        assert_eq!(issues[0].path, "members");
        assert_eq!(issues[0].message, "expected at least one member");
    }

    #[test]
    fn interpretation_validation_reports_shape_and_required_field_errors() {
        let root = yaml(
            r#"
analyses:
  - not-a-mapping
  - id: ""
    kind: python
    path: ""
    output_format: xml
    derived_from: []
  - id: a
    kind: bioscript
    path: a.bs
    derived_from: ["", 3]
    logic: not-a-mapping
  - id: b
    kind: bioscript
    path: b.bs
    derived_from: [rs1]
    emits: not-a-sequence
"#,
        );
        let mut issues = Vec::new();
        validate_panel_interpretations(&root, &mut issues);
        let paths = issue_paths(&issues);
        assert!(paths.contains(&"analyses[0]".to_owned()));
        assert!(paths.contains(&"analyses[1].id".to_owned()));
        assert!(paths.contains(&"analyses[1].kind".to_owned()));
        assert!(paths.contains(&"analyses[1].path".to_owned()));
        assert!(paths.contains(&"analyses[1].output_format".to_owned()));
        assert!(paths.contains(&"analyses[1].derived_from".to_owned()));
        assert!(paths.contains(&"analyses[2].derived_from[0]".to_owned()));
        assert!(paths.contains(&"analyses[2].derived_from[1]".to_owned()));
        assert!(paths.contains(&"analyses[2].logic".to_owned()));
        assert!(paths.contains(&"analyses[3].emits".to_owned()));
    }

    #[test]
    fn interpretation_logic_and_emits_validate_nested_strings_and_urls() {
        let root = yaml(
            r#"
analyses:
  - id: a
    kind: bioscript
    path: a.bs
    derived_from: [rs1]
    logic:
      description: ""
      source:
        name: 3
        url: not-a-url
    emits:
      - not-a-mapping
      - key: ""
        label: 3
        value_type: ""
        format: {}
interpretations:
  - id: old
    kind: bioscript
    path: old.bs
    derived_from: [rs1]
"#,
        );
        let mut issues = Vec::new();
        validate_panel_interpretations(&root, &mut issues);
        let paths = issue_paths(&issues);
        assert!(paths.contains(&"interpretations".to_owned()));
        assert!(paths.contains(&"analyses[0].logic.description".to_owned()));
        assert!(paths.contains(&"analyses[0].logic.source.name".to_owned()));
        assert!(paths.contains(&"analyses[0].logic.source.url".to_owned()));
        assert!(paths.contains(&"analyses[0].emits[0]".to_owned()));
        assert!(paths.contains(&"analyses[0].emits[1].key".to_owned()));
        assert!(paths.contains(&"analyses[0].emits[1].label".to_owned()));
        assert!(paths.contains(&"analyses[0].emits[1].value_type".to_owned()));
        assert!(paths.contains(&"analyses[0].emits[1].format".to_owned()));
    }

    #[test]
    fn variant_spec_parsing_extracts_coordinates_alleles_and_motifs() {
        let root = yaml(
            r#"
identifiers:
  rsids: [rs1]
coordinates:
  grch37:
    chrom: "1"
    pos: 10
    assembly_ref: b37
  grch38:
    chrom: "1"
    start: 20
    end: 22
    assembly_ref: b38
alleles:
  kind: deletion
  ref: AT
  alts: [A]
  deletion_length: 1
  motifs: [T]
"#,
        );
        let spec = variant_spec_from_root(&root).unwrap();
        assert_eq!(spec.rsids, vec!["rs1"]);
        assert_eq!(spec.grch37.unwrap().start, 10);
        assert_eq!(spec.grch38.unwrap().end, 22);
        assert_eq!(spec.grch37_assembly_ref.as_deref(), Some("b37"));
        assert_eq!(spec.grch38_assembly_ref.as_deref(), Some("b38"));
        assert_eq!(spec.reference.as_deref(), Some("AT"));
        assert_eq!(spec.alternate.as_deref(), Some("A"));
        assert_eq!(spec.kind, Some(VariantKind::Deletion));
        assert_eq!(spec.deletion_length, Some(1));
        assert_eq!(spec.motifs, vec!["T"]);

        let missing = yaml("coordinates: {grch38: {chrom: '1'}}");
        assert!(locus_from_root(&missing, "grch38")
            .unwrap_err()
            .contains("start missing"));
    }
}
