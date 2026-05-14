#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutablePanelMember<'a> {
    Variant(&'a str),
    VariantCatalogue(&'a str),
    Assay(&'a str),
}

pub fn panel_executable_member<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<ExecutablePanelMember<'a>, String> {
    let path = panel_executable_member_path(kind, path)?;
    match kind {
        "variant" => Ok(ExecutablePanelMember::Variant(path)),
        "variant-catalogue" => Ok(ExecutablePanelMember::VariantCatalogue(path)),
        "assay" => Ok(ExecutablePanelMember::Assay(path)),
        _ => unreachable!("panel_executable_member_path validates member kind"),
    }
}

pub fn panel_executable_member_path<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<&'a str, String> {
    let Some(path) = path else {
        return Err("remote panel members are not executable yet".to_owned());
    };
    if !matches!(kind, "variant" | "variant-catalogue" | "assay") {
        return Err(format!("panel member kind '{kind}' is not executable"));
    }
    Ok(path)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutableAssayMember<'a> {
    Variant(&'a str),
    VariantCatalogue(&'a str),
}

pub fn assay_executable_member<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<ExecutableAssayMember<'a>, String> {
    let path = assay_executable_member_path(kind, path)?;
    match kind {
        "variant" => Ok(ExecutableAssayMember::Variant(path)),
        "variant-catalogue" => Ok(ExecutableAssayMember::VariantCatalogue(path)),
        _ => unreachable!("assay_executable_member_path validates member kind"),
    }
}

pub fn assay_executable_member_path<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<&'a str, String> {
    if !matches!(kind, "variant" | "variant-catalogue") {
        return Err(format!("assay member kind '{kind}' is not executable"));
    }
    let Some(path) = path else {
        return Err("remote assay members are not executable yet".to_owned());
    };
    Ok(path)
}

pub(crate) fn traversable_manifest_member_paths<'a>(
    schema: &str,
    value: &'a serde_yaml::Value,
) -> Vec<&'a str> {
    if !manifest_has_traversable_members(schema) {
        return Vec::new();
    }
    value
        .get("members")
        .and_then(serde_yaml::Value::as_sequence)
        .into_iter()
        .flatten()
        .filter_map(traversable_manifest_member_path)
        .collect()
}

fn manifest_has_traversable_members(schema: &str) -> bool {
    matches!(schema, "bioscript:assay:1.0" | "bioscript:panel:1.0")
}

fn traversable_manifest_member_path(member: &serde_yaml::Value) -> Option<&str> {
    let kind = member.get("kind").and_then(serde_yaml::Value::as_str)?;
    if !matches!(kind, "variant" | "variant-catalogue" | "assay") {
        return None;
    }
    member.get("path").and_then(serde_yaml::Value::as_str)
}
