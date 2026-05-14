use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use bioscript_schema::{
    PanelInterpretation, VariantManifest, load_assay_manifest_text, load_panel_manifest_text,
    load_variant_manifest_text,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportManifestKind {
    Variant,
    Panel,
    Assay,
}

pub fn report_manifest_kind(schema: &str) -> Result<ReportManifestKind, String> {
    match schema {
        "bioscript:variant:1.0" | "bioscript:variant" => Ok(ReportManifestKind::Variant),
        "bioscript:panel:1.0" => Ok(ReportManifestKind::Panel),
        "bioscript:assay:1.0" => Ok(ReportManifestKind::Assay),
        other => Err(format!("unsupported manifest schema '{other}'")),
    }
}

pub fn report_manifest_schema(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<String, String> {
    let value = workspace.load_yaml(path)?;
    yaml_string(&value, "schema").ok_or_else(|| format!("{path} is missing schema"))
}

pub trait ManifestWorkspace {
    fn load_text(&self, path: &str) -> Result<String, String>;

    fn load_yaml(&self, path: &str) -> Result<serde_yaml::Value, String>;
    fn resolve(&self, base: &str, relative: &str) -> Result<String, String>;
}

pub struct FilesystemManifestWorkspace {
    root: PathBuf,
}

impl FilesystemManifestWorkspace {
    #[must_use]
    pub fn new(root: impl Into<PathBuf>) -> Self {
        Self { root: root.into() }
    }
}

impl ManifestWorkspace for FilesystemManifestWorkspace {
    fn load_text(&self, path: &str) -> Result<String, String> {
        let path = Path::new(path);
        fs::read_to_string(path)
            .map_err(|err| format!("failed to read YAML {}: {err}", path.display()))
    }

    fn load_yaml(&self, path: &str) -> Result<serde_yaml::Value, String> {
        let text = self.load_text(path)?;
        serde_yaml::from_str(&text).map_err(|err| format!("failed to parse YAML {path}: {err}"))
    }

    fn resolve(&self, base: &str, relative: &str) -> Result<String, String> {
        resolve_filesystem_manifest_path(&self.root, Path::new(base), relative)
            .map(|path| path.display().to_string())
    }
}

pub fn resolve_filesystem_manifest_path(
    root: &Path,
    manifest_path: &Path,
    relative: &str,
) -> Result<PathBuf, String> {
    let base_dir = manifest_path
        .parent()
        .ok_or_else(|| format!("manifest has no parent: {}", manifest_path.display()))?;
    let joined = base_dir.join(relative);
    let canonical_root = root
        .canonicalize()
        .map_err(|err| format!("failed to resolve root {}: {err}", root.display()))?;
    let canonical_base = base_dir.canonicalize().map_err(|err| {
        format!(
            "failed to resolve manifest dir {}: {err}",
            base_dir.display()
        )
    })?;
    let canonical_joined = joined
        .canonicalize()
        .map_err(|err| format!("failed to resolve {}: {err}", joined.display()))?;
    let boundary = if canonical_base.starts_with(&canonical_root) {
        &canonical_root
    } else {
        &canonical_base
    };
    if !canonical_joined.starts_with(boundary) {
        return Err(format!(
            "manifest member path escapes bioscript root: {}",
            canonical_joined.display()
        ));
    }
    Ok(canonical_joined)
}

pub fn report_manifest_metadata(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<serde_json::Value, String> {
    let value = workspace.load_yaml(path)?;
    let members = value
        .get("members")
        .and_then(serde_yaml::Value::as_sequence)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_yaml::Value::as_mapping)
                .map(|mapping| {
                    serde_json::json!({
                        "kind": yaml_mapping_string(mapping, "kind"),
                        "path": yaml_mapping_string(mapping, "path"),
                        "version": yaml_mapping_string(mapping, "version"),
                    })
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    Ok(serde_json::json!({
        "schema": yaml_string(&value, "schema"),
        "version": yaml_string(&value, "version"),
        "name": yaml_string(&value, "name"),
        "label": yaml_string(&value, "label").or_else(|| yaml_string(&value, "name")),
        "tags": yaml_string_sequence(&value, "tags"),
        "members": members,
    }))
}

pub fn report_assay_id(workspace: &impl ManifestWorkspace, path: &str) -> Result<String, String> {
    let value = workspace.load_yaml(path)?;
    let schema =
        yaml_string(&value, "schema").ok_or_else(|| format!("{path} is missing schema"))?;
    report_manifest_kind(&schema)?;
    yaml_string(&value, "name").ok_or_else(|| format!("manifest is missing required name: {path}"))
}

#[derive(Clone, Debug)]
pub struct ReportManifestContext {
    pub assay_id: String,
    pub manifest_metadata: serde_json::Value,
    pub findings: Vec<serde_json::Value>,
    pub provenance: Vec<serde_json::Value>,
}

pub fn load_report_manifest_context(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<ReportManifestContext, String> {
    Ok(ReportManifestContext {
        assay_id: report_assay_id(workspace, path)?,
        manifest_metadata: report_manifest_metadata(workspace, path)?,
        findings: load_manifest_findings(workspace, path)?,
        provenance: load_manifest_provenance_links(workspace, path)?,
    })
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnalysisManifestTask {
    pub manifest_path: String,
    pub manifest_name: String,
    pub interpretations: Vec<PanelInterpretation>,
}

pub fn collect_analysis_manifest_tasks(
    workspace: &impl ManifestWorkspace,
    path: &str,
    filters: &[String],
) -> Result<Vec<AnalysisManifestTask>, String> {
    match report_manifest_kind(&report_manifest_schema(workspace, path)?)? {
        ReportManifestKind::Panel => {
            let text = workspace.load_text(path)?;
            let manifest = load_panel_manifest_text(path, &text)?;
            let mut tasks = Vec::new();
            if filters.is_empty() && !manifest.interpretations.is_empty() {
                tasks.push(AnalysisManifestTask {
                    manifest_path: path.to_owned(),
                    manifest_name: manifest.name.clone(),
                    interpretations: manifest.interpretations.clone(),
                });
            }
            for member in &manifest.members {
                if member.kind != "assay" {
                    continue;
                }
                let Some(member_path) = &member.path else {
                    continue;
                };
                let resolved = workspace.resolve(path, member_path)?;
                if !matches_analysis_path_filters(&resolved, filters) {
                    continue;
                }
                tasks.extend(collect_analysis_manifest_tasks(
                    workspace, &resolved, filters,
                )?);
            }
            Ok(tasks)
        }
        ReportManifestKind::Assay => {
            let text = workspace.load_text(path)?;
            let manifest = load_assay_manifest_text(path, &text)?;
            if manifest.interpretations.is_empty() {
                return Ok(Vec::new());
            }
            Ok(vec![AnalysisManifestTask {
                manifest_path: path.to_owned(),
                manifest_name: manifest.name,
                interpretations: manifest.interpretations,
            }])
        }
        ReportManifestKind::Variant => Ok(Vec::new()),
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VariantManifestTask {
    pub manifest_path: String,
    pub manifest: VariantManifest,
}

pub fn collect_variant_manifest_tasks(
    workspace: &impl ManifestWorkspace,
    path: &str,
    filters: &[String],
) -> Result<Vec<VariantManifestTask>, String> {
    match report_manifest_kind(&report_manifest_schema(workspace, path)?)? {
        ReportManifestKind::Variant => {
            let text = workspace.load_text(path)?;
            let manifest = load_variant_manifest_text(path, &text)?;
            Ok(vec![VariantManifestTask {
                manifest_path: path.to_owned(),
                manifest,
            }])
        }
        ReportManifestKind::Panel => {
            let text = workspace.load_text(path)?;
            let manifest = load_panel_manifest_text(path, &text)?;
            let mut tasks = Vec::new();
            for member in &manifest.members {
                match panel_executable_member(&member.kind, member.path.as_deref())? {
                    ExecutablePanelMember::Variant(member_path) => {
                        let resolved = workspace.resolve(path, member_path)?;
                        let variant = load_variant_task(workspace, &resolved)?;
                        if matches_variant_manifest_filters(&variant.manifest, &resolved, filters) {
                            tasks.push(variant);
                        }
                    }
                    ExecutablePanelMember::Assay(member_path) => {
                        let resolved = workspace.resolve(path, member_path)?;
                        tasks.extend(collect_variant_manifest_tasks(
                            workspace, &resolved, filters,
                        )?);
                    }
                }
            }
            Ok(tasks)
        }
        ReportManifestKind::Assay => {
            let text = workspace.load_text(path)?;
            let manifest = load_assay_manifest_text(path, &text)?;
            let mut tasks = Vec::new();
            for member in &manifest.members {
                match assay_executable_member(&member.kind, member.path.as_deref())? {
                    ExecutableAssayMember::Variant(member_path) => {
                        let resolved = workspace.resolve(path, member_path)?;
                        let variant = load_variant_task(workspace, &resolved)?;
                        if matches_variant_manifest_filters(&variant.manifest, &resolved, filters) {
                            tasks.push(variant);
                        }
                    }
                }
            }
            Ok(tasks)
        }
    }
}

fn load_variant_task(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<VariantManifestTask, String> {
    let text = workspace.load_text(path)?;
    let manifest = load_variant_manifest_text(path, &text)?;
    Ok(VariantManifestTask {
        manifest_path: path.to_owned(),
        manifest,
    })
}

pub fn load_manifest_findings(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<Vec<serde_json::Value>, String> {
    let value = workspace.load_yaml(path)?;
    let schema = yaml_string(&value, "schema").unwrap_or_default();
    let mut findings = Vec::new();

    if manifest_supports_findings(&schema)
        && let Some(items) = value
            .get("findings")
            .and_then(serde_yaml::Value::as_sequence)
    {
        for item in items {
            let json_item = yaml_to_json(item.clone())?;
            let include = json_item
                .get("include")
                .and_then(serde_json::Value::as_str)
                .map(str::to_owned);
            if let Some(include) = include {
                let include_path = workspace.resolve(path, &include)?;
                let mut included = load_manifest_findings(workspace, &include_path)?;
                let inherited_binding = json_item.get("binding").cloned();
                for included_item in &mut included {
                    if inherited_binding.is_some()
                        && included_item.get("binding").is_none()
                        && included_item.get("effects").is_none()
                        && let Some(object) = included_item.as_object_mut()
                    {
                        object.insert(
                            "binding".to_owned(),
                            inherited_binding.clone().unwrap_or(serde_json::Value::Null),
                        );
                    }
                }
                findings.extend(included);
                continue;
            }
            if json_item.get("include").is_none() {
                findings.push(json_item);
            }
        }
    }

    for member_path in traversable_manifest_member_paths(&schema, &value) {
        let resolved = workspace.resolve(path, member_path)?;
        findings.extend(load_manifest_findings(workspace, &resolved)?);
    }

    Ok(findings)
}

pub fn load_manifest_provenance_links(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<Vec<serde_json::Value>, String> {
    let value = workspace.load_yaml(path)?;
    let schema = yaml_string(&value, "schema").unwrap_or_default();
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    collect_manifest_provenance_entries(&value, &mut links)?;

    if manifest_supports_findings(&schema)
        && let Some(items) = value
            .get("findings")
            .and_then(serde_yaml::Value::as_sequence)
    {
        for item in items {
            let json_item = yaml_to_json(item.clone())?;
            let Some(include) = json_item.get("include").and_then(serde_json::Value::as_str) else {
                continue;
            };
            let include_path = workspace.resolve(path, include)?;
            for item in load_manifest_provenance_links(workspace, &include_path)? {
                if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                    links.entry(url.to_owned()).or_insert(item);
                }
            }
        }
    }

    for member_path in traversable_manifest_member_paths(&schema, &value) {
        let resolved = workspace.resolve(path, member_path)?;
        for item in load_manifest_provenance_links(workspace, &resolved)? {
            if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                links.entry(url.to_owned()).or_insert(item);
            }
        }
    }

    Ok(links.into_values().collect())
}

pub fn collect_manifest_provenance_entries(
    value: &serde_yaml::Value,
    links: &mut BTreeMap<String, serde_json::Value>,
) -> Result<(), String> {
    if let Some(sources) = value
        .get("provenance")
        .and_then(|provenance| provenance.get("sources"))
        .and_then(serde_yaml::Value::as_sequence)
    {
        for source in sources {
            let json = yaml_to_json(source.clone())?;
            if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
                links.entry(url.to_owned()).or_insert(json);
            }
        }
    }
    if let Some(source) = value.get("source") {
        let json = yaml_to_json(source.clone())?;
        if let Some(url) = json.get("url").and_then(serde_json::Value::as_str) {
            links.entry(url.to_owned()).or_insert(json);
        }
    }
    Ok(())
}

pub fn matches_variant_manifest_filters(
    manifest: &VariantManifest,
    path: &str,
    filters: &[String],
) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("kind", value)) => value == "variant",
        Some(("name", value)) => manifest.name.contains(value),
        Some(("path", value)) => path.contains(value),
        Some(("tag", value)) => manifest.tags.iter().any(|tag| tag == value),
        Some(_) | None => false,
    })
}

pub fn matches_analysis_path_filters(path: &str, filters: &[String]) -> bool {
    filters.iter().all(|filter| match filter.split_once('=') {
        Some(("path", value)) => path.contains(value),
        _ => false,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutablePanelMember<'a> {
    Variant(&'a str),
    Assay(&'a str),
}

pub fn panel_executable_member<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<ExecutablePanelMember<'a>, String> {
    let path = panel_executable_member_path(kind, path)?;
    match kind {
        "variant" => Ok(ExecutablePanelMember::Variant(path)),
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
    if !matches!(kind, "variant" | "assay") {
        return Err(format!("panel member kind '{kind}' is not executable"));
    }
    Ok(path)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutableAssayMember<'a> {
    Variant(&'a str),
}

pub fn assay_executable_member<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<ExecutableAssayMember<'a>, String> {
    let path = assay_executable_member_path(kind, path)?;
    Ok(ExecutableAssayMember::Variant(path))
}

pub fn assay_executable_member_path<'a>(
    kind: &str,
    path: Option<&'a str>,
) -> Result<&'a str, String> {
    if kind != "variant" {
        return Err(format!("assay member kind '{kind}' is not executable"));
    }
    let Some(path) = path else {
        return Err("remote assay members are not executable yet".to_owned());
    };
    Ok(path)
}

fn manifest_supports_findings(schema: &str) -> bool {
    matches!(
        schema,
        "bioscript:variant:1.0"
            | "bioscript:variant"
            | "bioscript:assay:1.0"
            | "bioscript:panel:1.0"
            | "bioscript:pgx-findings:1.0"
    )
}

fn manifest_has_traversable_members(schema: &str) -> bool {
    matches!(schema, "bioscript:assay:1.0" | "bioscript:panel:1.0")
}

fn traversable_manifest_member_paths<'a>(
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

fn traversable_manifest_member_path(member: &serde_yaml::Value) -> Option<&str> {
    let kind = member.get("kind").and_then(serde_yaml::Value::as_str)?;
    if !matches!(kind, "variant" | "assay") {
        return None;
    }
    member.get("path").and_then(serde_yaml::Value::as_str)
}

fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, String> {
    serde_json::to_value(value).map_err(|err| format!("failed to convert YAML to JSON: {err}"))
}

fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

fn yaml_string_sequence(value: &serde_yaml::Value, key: &str) -> Vec<serde_json::Value> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_sequence)
        .map(|items| {
            items
                .iter()
                .filter_map(serde_yaml::Value::as_str)
                .map(serde_json::Value::from)
                .collect()
        })
        .unwrap_or_default()
}

fn yaml_mapping_string(mapping: &serde_yaml::Mapping, key: &str) -> Option<String> {
    mapping
        .get(serde_yaml::Value::String(key.to_owned()))
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

#[cfg(test)]
mod tests {
    use std::{
        collections::BTreeMap,
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    use bioscript_core::VariantSpec;
    use bioscript_schema::VariantManifest;

    use super::{
        ExecutableAssayMember, ExecutablePanelMember, ManifestWorkspace, ReportManifestKind,
        assay_executable_member, assay_executable_member_path, collect_analysis_manifest_tasks,
        collect_variant_manifest_tasks, matches_analysis_path_filters,
        matches_variant_manifest_filters, panel_executable_member, panel_executable_member_path,
        report_assay_id, report_manifest_kind, report_manifest_schema,
        resolve_filesystem_manifest_path, traversable_manifest_member_paths,
    };

    struct InlineWorkspace {
        yaml: &'static str,
    }

    impl ManifestWorkspace for InlineWorkspace {
        fn load_text(&self, _path: &str) -> Result<String, String> {
            Ok(self.yaml.to_owned())
        }

        fn load_yaml(&self, _path: &str) -> Result<serde_yaml::Value, String> {
            serde_yaml::from_str(self.yaml).map_err(|err| err.to_string())
        }

        fn resolve(&self, _base: &str, relative: &str) -> Result<String, String> {
            Ok(relative.to_owned())
        }
    }

    struct MapWorkspace {
        files: BTreeMap<String, String>,
    }

    impl ManifestWorkspace for MapWorkspace {
        fn load_text(&self, path: &str) -> Result<String, String> {
            self.files
                .get(path)
                .cloned()
                .ok_or_else(|| format!("missing file: {path}"))
        }

        fn load_yaml(&self, path: &str) -> Result<serde_yaml::Value, String> {
            serde_yaml::from_str(&self.load_text(path)?).map_err(|err| err.to_string())
        }

        fn resolve(&self, base: &str, relative: &str) -> Result<String, String> {
            let base = std::path::Path::new(base)
                .parent()
                .unwrap_or_else(|| std::path::Path::new(""));
            Ok(base.join(relative).display().to_string())
        }
    }

    #[test]
    fn report_assay_id_uses_manifest_name_for_supported_schemas() {
        let workspace = InlineWorkspace {
            yaml: "schema: bioscript:panel:1.0\nname: pgx-1\n",
        };
        assert_eq!(
            report_assay_id(&workspace, "manifest.yaml").unwrap(),
            "pgx-1"
        );
    }

    #[test]
    fn report_manifest_schema_reads_workspace_yaml() {
        let workspace = InlineWorkspace {
            yaml: "schema: bioscript:assay:1.0\nname: apoe\n",
        };
        assert_eq!(
            report_manifest_schema(&workspace, "assay.yaml").unwrap(),
            "bioscript:assay:1.0"
        );

        let workspace = InlineWorkspace {
            yaml: "name: apoe\n",
        };
        assert_eq!(
            report_manifest_schema(&workspace, "assay.yaml").unwrap_err(),
            "assay.yaml is missing schema"
        );
    }

    #[test]
    fn report_assay_id_rejects_unsupported_schema() {
        let workspace = InlineWorkspace {
            yaml: "schema: other\nname: nope\n",
        };
        let err = report_assay_id(&workspace, "manifest.yaml").unwrap_err();
        assert!(err.contains("unsupported manifest schema"));
    }

    #[test]
    fn report_manifest_kind_centralizes_schema_dispatch() {
        assert_eq!(
            report_manifest_kind("bioscript:variant").unwrap(),
            ReportManifestKind::Variant
        );
        assert_eq!(
            report_manifest_kind("bioscript:variant:1.0").unwrap(),
            ReportManifestKind::Variant
        );
        assert_eq!(
            report_manifest_kind("bioscript:panel:1.0").unwrap(),
            ReportManifestKind::Panel
        );
        assert_eq!(
            report_manifest_kind("bioscript:assay:1.0").unwrap(),
            ReportManifestKind::Assay
        );
        assert_eq!(
            report_manifest_kind("other").unwrap_err(),
            "unsupported manifest schema 'other'"
        );
    }

    #[test]
    fn variant_manifest_filters_match_kind_name_path_and_tag() {
        let manifest = VariantManifest {
            path: PathBuf::from("rs1.yaml"),
            name: "APOE_rs429358".to_owned(),
            tags: vec!["gene:APOE".to_owned(), "pgx".to_owned()],
            spec: VariantSpec::default(),
        };
        let filters = vec![
            "kind=variant".to_owned(),
            "name=APOE".to_owned(),
            "path=assets/APOE".to_owned(),
            "tag=pgx".to_owned(),
        ];
        assert!(matches_variant_manifest_filters(
            &manifest,
            "assets/APOE/rs429358.yaml",
            &filters
        ));
        assert!(!matches_variant_manifest_filters(
            &manifest,
            "assets/MTHFR/rs429358.yaml",
            &filters
        ));
    }

    #[test]
    fn analysis_path_filters_only_accept_path_filters() {
        assert!(matches_analysis_path_filters(
            "assets/APOE/assay.yaml",
            &["path=APOE".to_owned()]
        ));
        assert!(!matches_analysis_path_filters(
            "assets/MTHFR/assay.yaml",
            &["path=APOE".to_owned()]
        ));
        assert!(!matches_analysis_path_filters(
            "assets/APOE/assay.yaml",
            &["tag=pgx".to_owned()]
        ));
        assert!(matches_analysis_path_filters("assets/APOE/assay.yaml", &[]));
    }

    #[test]
    fn collect_analysis_manifest_tasks_matches_cli_panel_filter_semantics() {
        let workspace = MapWorkspace {
            files: BTreeMap::from([
                (
                    "panel.yaml".to_owned(),
                    r#"
schema: bioscript:panel:1.0
version: "1.0"
name: panel
analyses:
  - id: panel_analysis
    kind: bioscript
    path: panel.py
    derived_from:
      - rs1.yaml
members:
  - kind: assay
    path: assets/APOE/assay.yaml
"#
                    .to_owned(),
                ),
                (
                    "assets/APOE/assay.yaml".to_owned(),
                    r#"
schema: bioscript:assay:1.0
version: "1.0"
name: apoe
members:
  - kind: variant
    path: rs2.yaml
analyses:
  - id: apoe_analysis
    kind: bioscript
    path: apoe.py
    derived_from:
      - rs2.yaml
"#
                    .to_owned(),
                ),
            ]),
        };

        let tasks = collect_analysis_manifest_tasks(&workspace, "panel.yaml", &[]).unwrap();
        assert_eq!(tasks.len(), 2);
        assert_eq!(tasks[0].manifest_name, "panel");
        assert_eq!(tasks[0].interpretations[0].id, "panel_analysis");
        assert_eq!(tasks[1].manifest_path, "assets/APOE/assay.yaml");
        assert_eq!(tasks[1].interpretations[0].id, "apoe_analysis");

        let filtered = collect_analysis_manifest_tasks(
            &workspace,
            "panel.yaml",
            &["path=assets/APOE/assay.yaml".to_owned()],
        )
        .unwrap();
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].manifest_name, "apoe");
    }

    #[test]
    fn collect_variant_manifest_tasks_preserves_nested_member_order_and_filters() {
        let workspace = MapWorkspace {
            files: BTreeMap::from([
                (
                    "panel.yaml".to_owned(),
                    r#"
schema: bioscript:panel:1.0
version: "1.0"
name: panel
members:
  - kind: variant
    path: rs1.yaml
  - kind: assay
    path: assets/APOE/assay.yaml
  - kind: variant
    path: rs3.yaml
"#
                    .to_owned(),
                ),
                (
                    "rs1.yaml".to_owned(),
                    r#"
schema: bioscript:variant:1.0
version: "1.0"
name: rs1
tags: [keep]
identifiers:
  rsids:
    - rs1
coordinates:
  grch38:
    chrom: "1"
    pos: 1
alleles:
  kind: snv
  ref: A
  alts:
    - G
"#
                    .to_owned(),
                ),
                (
                    "assets/APOE/assay.yaml".to_owned(),
                    r#"
schema: bioscript:assay:1.0
version: "1.0"
name: apoe
members:
  - kind: variant
    path: rs2.yaml
"#
                    .to_owned(),
                ),
                (
                    "assets/APOE/rs2.yaml".to_owned(),
                    r#"
schema: bioscript:variant:1.0
version: "1.0"
name: rs2
tags: [keep]
identifiers:
  rsids:
    - rs2
coordinates:
  grch38:
    chrom: "1"
    pos: 2
alleles:
  kind: snv
  ref: A
  alts:
    - G
"#
                    .to_owned(),
                ),
                (
                    "rs3.yaml".to_owned(),
                    r#"
schema: bioscript:variant:1.0
version: "1.0"
name: rs3
tags: [skip]
identifiers:
  rsids:
    - rs3
coordinates:
  grch38:
    chrom: "1"
    pos: 3
alleles:
  kind: snv
  ref: A
  alts:
    - G
"#
                    .to_owned(),
                ),
            ]),
        };

        let tasks = collect_variant_manifest_tasks(&workspace, "panel.yaml", &[]).unwrap();
        assert_eq!(
            tasks
                .iter()
                .map(|task| task.manifest.name.as_str())
                .collect::<Vec<_>>(),
            vec!["rs1", "rs2", "rs3"]
        );

        let filtered =
            collect_variant_manifest_tasks(&workspace, "panel.yaml", &["tag=keep".to_owned()])
                .unwrap();
        assert_eq!(
            filtered
                .iter()
                .map(|task| task.manifest_path.as_str())
                .collect::<Vec<_>>(),
            vec!["rs1.yaml", "assets/APOE/rs2.yaml"]
        );
    }

    #[test]
    fn executable_member_validation_matches_cli_error_contract() {
        assert_eq!(
            panel_executable_member("variant", Some("rs.yaml")).unwrap(),
            ExecutablePanelMember::Variant("rs.yaml")
        );
        assert_eq!(
            panel_executable_member("assay", Some("assay.yaml")).unwrap(),
            ExecutablePanelMember::Assay("assay.yaml")
        );
        assert_eq!(
            panel_executable_member_path("variant", Some("rs.yaml")).unwrap(),
            "rs.yaml"
        );
        assert_eq!(
            panel_executable_member_path("assay", Some("assay.yaml")).unwrap(),
            "assay.yaml"
        );
        assert_eq!(
            panel_executable_member_path("download", Some("remote.yaml")).unwrap_err(),
            "panel member kind 'download' is not executable"
        );
        assert_eq!(
            panel_executable_member_path("variant", None).unwrap_err(),
            "remote panel members are not executable yet"
        );

        assert_eq!(
            assay_executable_member("variant", Some("rs.yaml")).unwrap(),
            ExecutableAssayMember::Variant("rs.yaml")
        );
        assert_eq!(
            assay_executable_member_path("variant", Some("rs.yaml")).unwrap(),
            "rs.yaml"
        );
        assert_eq!(
            assay_executable_member_path("assay", Some("nested.yaml")).unwrap_err(),
            "assay member kind 'assay' is not executable"
        );
        assert_eq!(
            assay_executable_member_path("variant", None).unwrap_err(),
            "remote assay members are not executable yet"
        );
    }

    #[test]
    fn traversable_manifest_members_keep_metadata_traversal_semantics() {
        let value: serde_yaml::Value = serde_yaml::from_str(
            r#"
schema: bioscript:panel:1.0
members:
  - kind: variant
    path: assets/A/rs1.yaml
  - kind: assay
    path: assets/B/assay.yaml
  - kind: download
    path: remote.yaml
  - kind: variant
  - path: missing-kind.yaml
"#,
        )
        .unwrap();

        assert_eq!(
            traversable_manifest_member_paths("bioscript:panel:1.0", &value),
            vec!["assets/A/rs1.yaml", "assets/B/assay.yaml"]
        );
        assert!(traversable_manifest_member_paths("bioscript:variant:1.0", &value).is_empty());
    }

    #[test]
    fn filesystem_manifest_path_resolution_enforces_root_boundary() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let root = std::env::temp_dir().join(format!(
            "bioscript-reporting-resolve-{}-{unique}",
            std::process::id()
        ));
        let assay_dir = root.join("assay");
        let outside = root
            .parent()
            .unwrap_or_else(|| root.as_path())
            .join(format!("bioscript-reporting-outside-{unique}"));
        fs::create_dir_all(&assay_dir).unwrap();
        fs::create_dir_all(&outside).unwrap();
        fs::write(
            assay_dir.join("manifest.yaml"),
            "schema: bioscript:assay:1.0\n",
        )
        .unwrap();
        fs::write(
            assay_dir.join("variant.yaml"),
            "schema: bioscript:variant:1.0\n",
        )
        .unwrap();
        fs::write(
            outside.join("variant.yaml"),
            "schema: bioscript:variant:1.0\n",
        )
        .unwrap();

        let manifest = assay_dir.join("manifest.yaml");
        assert_eq!(
            resolve_filesystem_manifest_path(&root, &manifest, "variant.yaml").unwrap(),
            assay_dir.join("variant.yaml").canonicalize().unwrap()
        );
        assert!(
            resolve_filesystem_manifest_path(
                &root,
                &manifest,
                &format!(
                    "../../{}/variant.yaml",
                    outside.file_name().unwrap().to_string_lossy()
                )
            )
            .unwrap_err()
            .contains("manifest member path escapes bioscript root")
        );

        let _ = fs::remove_dir_all(&root);
        let _ = fs::remove_dir_all(&outside);
    }
}
