use super::*;

#[path = "report_workspace/analysis.rs"]
mod analysis;
pub(crate) use analysis::WasmReportAnalysisRunner;

pub(super) struct PackageWorkspace {
    files: BTreeMap<String, String>,
}

impl PackageWorkspace {
    pub(super) fn new(files: Vec<PackageFileInput>) -> Result<Self, JsError> {
        let mut map = BTreeMap::new();
        for file in files {
            let _ = file.source_url;
            map.insert(normalize_package_path(&file.path)?, file.contents);
        }
        Ok(Self { files: map })
    }

    fn text(&self, path: &str) -> Result<&str, JsError> {
        let normalized = normalize_package_path(path)?;
        self.files
            .get(&normalized)
            .map(String::as_str)
            .ok_or_else(|| JsError::new(&format!("package file not found: {normalized}")))
    }

    fn yaml(&self, path: &str) -> Result<serde_yaml::Value, JsError> {
        serde_yaml::from_str(self.text(path)?)
            .map_err(|err| JsError::new(&format!("failed to parse YAML {path}: {err}")))
    }

    fn resolve(&self, base: &str, relative: &str) -> Result<String, JsError> {
        let base = Path::new(base).parent().unwrap_or_else(|| Path::new(""));
        normalize_package_path(&base.join(relative).display().to_string())
    }

    pub(super) fn load_variant(&self, path: &str) -> Result<VariantManifest, JsError> {
        load_variant_manifest_text(path, self.text(path)?)
            .map_err(|err| JsError::new(&format!("load variant {path}: {err}")))
    }

    pub(super) fn package_result_entrypoint(&self) -> Result<Option<String>, JsError> {
        let Some(text) = self.files.get("manifest.yaml") else {
            return Ok(None);
        };
        let value: serde_yaml::Value = serde_yaml::from_str(text)
            .map_err(|err| JsError::new(&format!("failed to parse package manifest.yaml: {err}")))?;
        let Some(result) = value
            .as_mapping()
            .and_then(|mapping| mapping.get(serde_yaml::Value::String("result".to_owned())))
        else {
            return Ok(None);
        };
        let path = result.as_mapping().and_then(|mapping| {
            mapping
                .get(serde_yaml::Value::String("entrypoint".to_owned()))
                .and_then(serde_yaml::Value::as_str)
                .or_else(|| {
                    mapping
                        .get(serde_yaml::Value::String("primary_html".to_owned()))
                        .and_then(serde_yaml::Value::as_str)
                })
        });
        let path = path.or_else(|| result.as_str());
        path.map(normalize_package_path).transpose()
    }
}

impl bioscript_reporting::ManifestWorkspace for PackageWorkspace {
    fn load_text(&self, path: &str) -> Result<String, String> {
        self.text(path)
            .map(str::to_owned)
            .map_err(|err| format!("{err:?}"))
    }

    fn load_yaml(&self, path: &str) -> Result<serde_yaml::Value, String> {
        self.yaml(path).map_err(|err| format!("{err:?}"))
    }

    fn resolve(&self, base: &str, relative: &str) -> Result<String, String> {
        self.resolve(base, relative)
            .map_err(|err| format!("{err:?}"))
    }
}

impl bioscript_reporting::ReportWorkspace for PackageWorkspace {
    fn app_observation_from_manifest_row(
        &self,
        row: &BTreeMap<String, String>,
        assay_id: &str,
        inferred_sex: Option<&SexInference>,
        fallback_assembly: Option<Assembly>,
    ) -> Result<serde_json::Value, String> {
        let row_path = row.get("path").cloned().unwrap_or_default();
        let (manifest, gene, source, alt_alleles, observed_alt_alleles) = if row_path.contains('#')
        {
            let task = bioscript_reporting::load_variant_manifest_task_by_path(self, &row_path)?;
            let alt_alleles = task
                .manifest
                .spec
                .alternate
                .clone()
                .into_iter()
                .collect::<Vec<_>>();
            let observed_alt_alleles = task.manifest.spec.observed_alternates.clone();
            (
                task.manifest,
                String::new(),
                serde_json::Value::Null,
                alt_alleles,
                observed_alt_alleles,
            )
        } else {
            let manifest = self
                .load_variant(&row_path)
                .map_err(|err| format!("{err:?}"))?;
            let value = self.yaml(&row_path).map_err(|err| format!("{err:?}"))?;
            (
                manifest,
                yaml_string(&value, "gene").unwrap_or_default(),
                variant_primary_source_from_yaml(&value),
                variant_alt_alleles_from_yaml(&value),
                variant_observed_alt_alleles_from_yaml(&value),
            )
        };
        Ok(bioscript_reporting::app_observation_from_manifest_row(
            bioscript_reporting::AppObservationInput {
                row,
                row_path: &row_path,
                assay_id,
                manifest,
                gene,
                source,
                alt_alleles,
                observed_alt_alleles,
                inferred_sex,
                fallback_assembly,
            },
        ))
    }
}
