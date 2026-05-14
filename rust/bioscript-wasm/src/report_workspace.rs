use super::*;

#[path = "report_workspace/analysis.rs"]
mod analysis;

/// What a manifest row walk produces: human-readable rows for the
/// observation TSV/HTML, and the underlying `VariantObservation`s ready to
/// hand to the analysis runtime as a pre-resolved cache (so analysis Python
/// scripts' `genotypes.lookup_variants(plan)` call hits cache instead of
/// re-walking the genome).
pub(super) struct ManifestRowsOutput {
    pub rows: Vec<BTreeMap<String, String>>,
    pub observations: Vec<VariantObservation>,
}

/// Abstract per-variant observation source so the workspace can run against
/// either a path-based `GenotypeStore` (text/zip — bytes already in memory)
/// or a CRAM/VCF-reader-backed lookup that streams through JS-supplied
/// `readAt` callbacks.
pub(super) trait VariantLookup {
    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError>;
}

impl VariantLookup for GenotypeStore {
    fn lookup_variants(
        &self,
        specs: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        GenotypeStore::lookup_variants(self, specs)
    }
}

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

    pub(super) fn run_manifest_rows(
        &self,
        manifest_path: &str,
        store: &dyn VariantLookup,
        participant_id: &str,
        filters: &[String],
    ) -> Result<ManifestRowsOutput, JsError> {
        let tasks =
            bioscript_reporting::collect_variant_manifest_tasks(self, manifest_path, filters)
                .map_err(|err| JsError::new(&err))?;
        let observations = store
            .lookup_variants(
                &tasks
                    .iter()
                    .map(|task| task.manifest.spec.clone())
                    .collect::<Vec<_>>(),
            )
            .map_err(|err| JsError::new(&format!("manifest lookup failed: {err:?}")))?;
        let mut rows = Vec::with_capacity(tasks.len());
        let mut collected = Vec::with_capacity(tasks.len());
        for (task, observation) in tasks.into_iter().zip(observations) {
            rows.push(variant_row(
                &task.manifest_path,
                &task.manifest.name,
                &task.manifest.tags,
                &observation,
                participant_id,
            ));
            collected.push(observation);
        }
        Ok(ManifestRowsOutput {
            rows,
            observations: collected,
        })
    }

    pub(super) fn app_observation_from_manifest_row(
        &self,
        row: &BTreeMap<String, String>,
        assay_id: &str,
        inferred_sex: Option<&SexInference>,
        fallback_assembly: Option<Assembly>,
    ) -> Result<serde_json::Value, JsError> {
        let row_path = row.get("path").cloned().unwrap_or_default();
        let manifest = self.load_variant(&row_path)?;
        let value = self.yaml(&row_path)?;
        let gene = yaml_string(&value, "gene").unwrap_or_default();
        Ok(bioscript_reporting::app_observation_from_manifest_row(
            bioscript_reporting::AppObservationInput {
                row,
                row_path: &row_path,
                assay_id,
                manifest,
                gene,
                source: variant_primary_source_from_yaml(&value),
                observed_alt_alleles: variant_observed_alt_alleles_from_yaml(&value),
                inferred_sex,
                fallback_assembly,
            },
        ))
    }

    pub(super) fn report_manifest_context(
        &self,
        path: &str,
    ) -> Result<bioscript_reporting::ReportManifestContext, JsError> {
        bioscript_reporting::load_report_manifest_context(self, path)
            .map_err(|err| JsError::new(&err))
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
