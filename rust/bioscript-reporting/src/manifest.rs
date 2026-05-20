use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

use bioscript_schema::{
    PanelInterpretation, VariantManifest, load_assay_manifest_text, load_panel_manifest_text,
    load_variant_manifest_text,
};

#[path = "manifest_catalogue.rs"]
mod catalogue;
#[path = "manifest_members.rs"]
mod members;
#[path = "manifest_provenance.rs"]
mod provenance;

use catalogue::load_variant_catalogue_tasks;
pub(crate) use members::traversable_manifest_member_paths;
pub use members::{
    ExecutableAssayMember, ExecutablePanelMember, assay_executable_member,
    assay_executable_member_path, panel_executable_member, panel_executable_member_path,
};
pub use provenance::{collect_manifest_provenance_entries, load_manifest_provenance_links};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReportManifestKind {
    Variant,
    VariantCatalogue,
    Panel,
    Assay,
}

pub fn report_manifest_kind(schema: &str) -> Result<ReportManifestKind, String> {
    match schema {
        "bioscript:variant:1.0" | "bioscript:variant" => Ok(ReportManifestKind::Variant),
        "bioscript:variant-catalogue:1.0" => Ok(ReportManifestKind::VariantCatalogue),
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
        "summary": yaml_string(&value, "summary"),
        "summary_url": yaml_string(&value, "summary_url"),
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
        ReportManifestKind::Variant | ReportManifestKind::VariantCatalogue => Ok(Vec::new()),
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
        ReportManifestKind::VariantCatalogue => {
            load_variant_catalogue_tasks(workspace, path, filters)
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
                    ExecutablePanelMember::VariantCatalogue(member_path) => {
                        let resolved = workspace.resolve(path, member_path)?;
                        tasks.extend(load_variant_catalogue_tasks(workspace, &resolved, filters)?);
                    }
                    ExecutablePanelMember::Assay(member_path) => {
                        let resolved = workspace.resolve(path, member_path)?;
                        tasks.extend(collect_variant_manifest_tasks(
                            workspace, &resolved, filters,
                        )?);
                    }
                }
            }
            dedupe_variant_manifest_tasks(&mut tasks);
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
                    ExecutableAssayMember::VariantCatalogue(member_path) => {
                        let resolved = workspace.resolve(path, member_path)?;
                        tasks.extend(load_variant_catalogue_tasks(workspace, &resolved, filters)?);
                    }
                }
            }
            dedupe_variant_manifest_tasks(&mut tasks);
            Ok(tasks)
        }
    }
}

fn dedupe_variant_manifest_tasks(tasks: &mut Vec<VariantManifestTask>) {
    let mut seen = BTreeSet::new();
    tasks.retain(|task| seen.insert(task.manifest_path.clone()));
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

pub fn load_variant_manifest_task_by_path(
    workspace: &impl ManifestWorkspace,
    path: &str,
) -> Result<VariantManifestTask, String> {
    if let Some((catalogue_path, variant_id)) = path.rsplit_once('#') {
        return load_variant_catalogue_tasks(workspace, catalogue_path, &[])?
            .into_iter()
            .find(|task| task.manifest_path == path || task.manifest.name == variant_id)
            .ok_or_else(|| format!("variant catalogue entry not found: {path}"));
    }
    load_variant_task(workspace, path)
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

pub(super) fn manifest_supports_findings(schema: &str) -> bool {
    matches!(
        schema,
        "bioscript:variant:1.0"
            | "bioscript:variant"
            | "bioscript:assay:1.0"
            | "bioscript:panel:1.0"
            | "bioscript:pgx-findings:1.0"
    )
}

pub(super) fn yaml_to_json(value: serde_yaml::Value) -> Result<serde_json::Value, String> {
    serde_json::to_value(value).map_err(|err| format!("failed to convert YAML to JSON: {err}"))
}

pub(super) fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
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
        collect_variant_manifest_tasks, load_manifest_findings, load_report_manifest_context,
        load_variant_manifest_task_by_path, matches_analysis_path_filters,
        matches_variant_manifest_filters, panel_executable_member, panel_executable_member_path,
        report_assay_id, report_manifest_kind, report_manifest_metadata, report_manifest_schema,
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

    pub(super) struct MapWorkspace {
        pub(super) files: BTreeMap<String, String>,
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

    const NESTED_PANEL_YAML: &str = r#"
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
"#;

    const NESTED_ASSAY_YAML: &str = r#"
schema: bioscript:assay:1.0
version: "1.0"
name: apoe
members:
  - kind: variant
    path: rs2.yaml
"#;

    fn variant_yaml(name: &str, pos: u32, tag: &str) -> String {
        format!(
            r#"
schema: bioscript:variant:1.0
version: "1.0"
name: {name}
tags: [{tag}]
identifiers:
  rsids:
    - {name}
coordinates:
  grch38:
    chrom: "1"
    pos: {pos}
alleles:
  kind: snv
  ref: A
  alts:
    - G
"#
        )
    }

    fn nested_variant_workspace() -> MapWorkspace {
        MapWorkspace {
            files: BTreeMap::from([
                ("panel.yaml".to_owned(), NESTED_PANEL_YAML.to_owned()),
                ("rs1.yaml".to_owned(), variant_yaml("rs1", 1, "keep")),
                (
                    "assets/APOE/assay.yaml".to_owned(),
                    NESTED_ASSAY_YAML.to_owned(),
                ),
                (
                    "assets/APOE/rs2.yaml".to_owned(),
                    variant_yaml("rs2", 2, "keep"),
                ),
                ("rs3.yaml".to_owned(), variant_yaml("rs3", 3, "skip")),
            ]),
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
    fn report_manifest_metadata_collects_tags_members_and_fallback_label() {
        let workspace = InlineWorkspace {
            yaml: r#"
schema: bioscript:panel:1.0
version: "1.0"
name: pgx-panel
summary: Panel summary text.
summary_url: https://example.test/panel-summary
tags: [pgx, cardiology, 7]
members:
  - kind: variant
    path: rs1.yaml
    version: "1"
  - kind: assay
    path: assay.yaml
  - not: a-member
"#,
        };
        let metadata = report_manifest_metadata(&workspace, "panel.yaml").unwrap();
        assert_eq!(metadata["schema"], "bioscript:panel:1.0");
        assert_eq!(metadata["label"], "pgx-panel");
        assert_eq!(metadata["summary"], "Panel summary text.");
        assert_eq!(
            metadata["summary_url"],
            "https://example.test/panel-summary"
        );
        assert_eq!(metadata["tags"], serde_json::json!(["pgx", "cardiology"]));
        assert_eq!(metadata["members"][0]["kind"], "variant");
        assert_eq!(metadata["members"][0]["path"], "rs1.yaml");
        assert_eq!(metadata["members"][0]["version"], "1");
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
        let workspace = nested_variant_workspace();

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

        let assay_tasks =
            collect_variant_manifest_tasks(&workspace, "assets/APOE/assay.yaml", &[]).unwrap();
        assert_eq!(assay_tasks.len(), 1);
        assert_eq!(assay_tasks[0].manifest_path, "assets/APOE/rs2.yaml");

        let variant_tasks = collect_variant_manifest_tasks(&workspace, "rs1.yaml", &[]).unwrap();
        assert_eq!(variant_tasks.len(), 1);
        assert_eq!(variant_tasks[0].manifest.name, "rs1");
    }

    #[test]
    fn collect_variant_manifest_tasks_dedupes_catalogue_reached_through_panel_and_assay() {
        let workspace = MapWorkspace {
            files: BTreeMap::from([
                (
                    "panel.yaml".to_owned(),
                    r#"
schema: bioscript:panel:1.0
version: "1.0"
name: panel
members:
  - kind: variant-catalogue
    path: catalogue.yaml
  - kind: assay
    path: assay.yaml
"#
                    .to_owned(),
                ),
                (
                    "assay.yaml".to_owned(),
                    r#"
schema: bioscript:assay:1.0
version: "1.0"
name: assay
members:
  - kind: variant-catalogue
    path: catalogue.yaml
"#
                    .to_owned(),
                ),
                (
                    "catalogue.yaml".to_owned(),
                    r#"
schema: bioscript:variant-catalogue:1.0
version: "1.0"
name: catalogue
variants:
  source: variants.tsv
"#
                    .to_owned(),
                ),
                (
                    "variants.tsv".to_owned(),
                    r"id	name	rsid	gene	ref	alt	kind	grch38_chrom	grch38_pos
rs1	rs1	rs1	GENE	A	G	snp	1	123
"
                    .to_owned(),
                ),
            ]),
        };

        let tasks = collect_variant_manifest_tasks(&workspace, "panel.yaml", &[]).unwrap();
        assert_eq!(tasks.len(), 1);
        assert_eq!(tasks[0].manifest_path, "catalogue.yaml#rs1");
    }

    #[test]
    fn variant_catalogue_tasks_preserve_observed_alts_separately_from_reportable_alts() {
        let workspace = MapWorkspace {
            files: BTreeMap::from([
                (
                    "catalogue.yaml".to_owned(),
                    r#"
schema: bioscript:variant-catalogue:1.0
version: "1.0"
name: catalogue
variants:
  source: variants.tsv
"#
                    .to_owned(),
                ),
                (
                    "variants.tsv".to_owned(),
                    "variant_id\tname\trsid\tgene\tref\talts\tobserved_alts\tkind\tgrch38_chrom\tgrch38_pos\nrs1\trs1\trs1\tGENE\tA\tC\tC|G|T\tsnp\t1\t123\n"
                        .to_owned(),
                ),
            ]),
        };

        let task = load_variant_manifest_task_by_path(&workspace, "catalogue.yaml#rs1").unwrap();

        assert_eq!(task.manifest.spec.alternate.as_deref(), Some("C"));
        assert_eq!(
            task.manifest.spec.observed_alternates,
            vec!["C".to_owned(), "G".to_owned(), "T".to_owned()]
        );
    }

    #[test]
    fn manifest_context_and_findings_follow_includes_members_and_inherited_bindings() {
        let workspace = MapWorkspace {
            files: BTreeMap::from([
                (
                    "panel.yaml".to_owned(),
                    r#"
schema: bioscript:panel:1.0
version: "1.0"
name: panel
findings:
  - include: included.yaml
    binding:
      source: variant
      variant: rs1.yaml
      key: outcome
      value: variant
  - schema: bioscript:trait:1.0
    summary: panel direct
members:
  - kind: variant
    path: rs1.yaml
"#
                    .to_owned(),
                ),
                (
                    "included.yaml".to_owned(),
                    r#"
schema: bioscript:pgx-findings:1.0
version: "1.0"
rsid: rs1
findings:
  - schema: bioscript:trait:1.0
    summary: included inherited
  - schema: bioscript:trait:1.0
    summary: included own effects
    effects: []
"#
                    .to_owned(),
                ),
                ("rs1.yaml".to_owned(), variant_yaml("rs1", 1, "keep")),
            ]),
        };

        let findings = load_manifest_findings(&workspace, "panel.yaml").unwrap();
        assert!(
            findings
                .iter()
                .any(|finding| finding["summary"] == "panel direct")
        );
        let inherited = findings
            .iter()
            .find(|finding| finding["summary"] == "included inherited")
            .unwrap();
        assert_eq!(inherited["binding"]["variant"], "rs1.yaml");
        let own_effects = findings
            .iter()
            .find(|finding| finding["summary"] == "included own effects")
            .unwrap();
        assert!(own_effects.get("binding").is_none());

        let context = load_report_manifest_context(&workspace, "panel.yaml").unwrap();
        assert_eq!(context.assay_id, "panel");
        assert_eq!(context.manifest_metadata["name"], "panel");
        assert!(!context.findings.is_empty());
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
            r"
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
",
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
            .unwrap_or(root.as_path())
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
