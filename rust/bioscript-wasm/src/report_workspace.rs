use super::*;

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

    fn schema(&self, path: &str) -> Result<String, JsError> {
        self.yaml(path)?
            .get("schema")
            .and_then(serde_yaml::Value::as_str)
            .map(ToOwned::to_owned)
            .ok_or_else(|| JsError::new(&format!("{path} is missing schema")))
    }

    fn resolve(&self, base: &str, relative: &str) -> Result<String, JsError> {
        let base = Path::new(base).parent().unwrap_or_else(|| Path::new(""));
        normalize_package_path(&base.join(relative).display().to_string())
    }

    fn load_variant(&self, path: &str) -> Result<VariantManifest, JsError> {
        load_variant_manifest_text(path, self.text(path)?)
            .map_err(|err| JsError::new(&format!("load variant {path}: {err}")))
    }

    fn load_panel(&self, path: &str) -> Result<PanelManifest, JsError> {
        load_panel_manifest_text(path, self.text(path)?)
            .map_err(|err| JsError::new(&format!("load panel {path}: {err}")))
    }

    fn load_assay(&self, path: &str) -> Result<AssayManifest, JsError> {
        load_assay_manifest_text(path, self.text(path)?)
            .map_err(|err| JsError::new(&format!("load assay {path}: {err}")))
    }

    pub(super) fn run_manifest_rows(
        &self,
        manifest_path: &str,
        store: &GenotypeStore,
        participant_id: &str,
        filters: &[String],
    ) -> Result<Vec<BTreeMap<String, String>>, JsError> {
        match self.schema(manifest_path)?.as_str() {
            "bioscript:variant:1.0" | "bioscript:variant" => {
                let manifest = self.load_variant(manifest_path)?;
                let observation = store
                    .lookup_variant(&manifest.spec)
                    .map_err(|err| JsError::new(&format!("lookup {}: {err:?}", manifest.name)))?;
                Ok(vec![variant_row(
                    manifest_path,
                    &manifest.name,
                    &manifest.tags,
                    &observation,
                    participant_id,
                )])
            }
            "bioscript:panel:1.0" => self.run_panel_rows(manifest_path, store, participant_id, filters),
            "bioscript:assay:1.0" => self.run_assay_rows(manifest_path, store, participant_id, filters),
            other => Err(JsError::new(&format!("unsupported manifest schema '{other}'"))),
        }
    }

    fn run_panel_rows(
        &self,
        manifest_path: &str,
        store: &GenotypeStore,
        participant_id: &str,
        filters: &[String],
    ) -> Result<Vec<BTreeMap<String, String>>, JsError> {
        let panel = self.load_panel(manifest_path)?;
        let mut rows_by_member: Vec<Vec<BTreeMap<String, String>>> = vec![Vec::new(); panel.members.len()];
        let mut variants = Vec::<(usize, String, VariantManifest)>::new();
        for (index, member) in panel.members.iter().enumerate() {
            let Some(path) = &member.path else {
                return Err(JsError::new("remote panel members are not executable yet"));
            };
            let resolved = self.resolve(manifest_path, path)?;
            if member.kind == "variant" {
                let variant = self.load_variant(&resolved)?;
                if matches_filters(&variant, &resolved, filters) {
                    variants.push((index, resolved, variant));
                }
            } else if member.kind == "assay" {
                rows_by_member[index] = self.run_assay_rows(&resolved, store, participant_id, filters)?;
            }
        }
        let specs = variants
            .iter()
            .map(|(_, _, manifest)| manifest.spec.clone())
            .collect::<Vec<_>>();
        let observations = store
            .lookup_variants(&specs)
            .map_err(|err| JsError::new(&format!("panel lookup failed: {err:?}")))?;
        for ((member_index, resolved, manifest), observation) in variants.into_iter().zip(observations) {
            rows_by_member[member_index].push(variant_row(
                &resolved,
                &manifest.name,
                &manifest.tags,
                &observation,
                participant_id,
            ));
        }
        Ok(rows_by_member.into_iter().flatten().collect())
    }

    fn run_assay_rows(
        &self,
        manifest_path: &str,
        store: &GenotypeStore,
        participant_id: &str,
        filters: &[String],
    ) -> Result<Vec<BTreeMap<String, String>>, JsError> {
        let assay = self.load_assay(manifest_path)?;
        let mut variants = Vec::<(String, VariantManifest)>::new();
        for member in &assay.members {
            if member.kind != "variant" {
                continue;
            }
            let Some(path) = &member.path else {
                continue;
            };
            let resolved = self.resolve(manifest_path, path)?;
            let variant = self.load_variant(&resolved)?;
            if matches_filters(&variant, &resolved, filters) {
                variants.push((resolved, variant));
            }
        }
        let specs = variants
            .iter()
            .map(|(_, manifest)| manifest.spec.clone())
            .collect::<Vec<_>>();
        let observations = store
            .lookup_variants(&specs)
            .map_err(|err| JsError::new(&format!("assay lookup failed: {err:?}")))?;
        Ok(variants
            .into_iter()
            .zip(observations)
            .map(|((resolved, manifest), observation)| {
                variant_row(
                    &resolved,
                    &manifest.name,
                    &manifest.tags,
                    &observation,
                    participant_id,
                )
            })
            .collect())
    }

    pub(super) fn run_manifest_analyses(
        &self,
        manifest_path: &str,
        input_name: &str,
        input_bytes: &[u8],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
    ) -> Result<Vec<serde_json::Value>, JsError> {
        match self.schema(manifest_path)?.as_str() {
            "bioscript:panel:1.0" => {
                let panel = self.load_panel(manifest_path)?;
                let mut analyses = self.run_interpretations(
                    manifest_path,
                    &panel.name,
                    &panel.interpretations,
                    input_name,
                    input_bytes,
                    participant_id,
                    loader,
                    options,
                )?;
                for member in &panel.members {
                    if member.kind != "assay" {
                        continue;
                    }
                    let Some(path) = &member.path else {
                        continue;
                    };
                    let resolved = self.resolve(manifest_path, path)?;
                    analyses.extend(self.run_manifest_analyses(
                        &resolved,
                        input_name,
                        input_bytes,
                        participant_id,
                        loader,
                        options,
                    )?);
                }
                Ok(analyses)
            }
            "bioscript:assay:1.0" => {
                let assay = self.load_assay(manifest_path)?;
                self.run_interpretations(
                    manifest_path,
                    &assay.name,
                    &assay.interpretations,
                    input_name,
                    input_bytes,
                    participant_id,
                    loader,
                    options,
                )
            }
            _ => Ok(Vec::new()),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn run_interpretations(
        &self,
        manifest_path: &str,
        manifest_name: &str,
        interpretations: &[PanelInterpretation],
        input_name: &str,
        input_bytes: &[u8],
        participant_id: &str,
        loader: &GenotypeLoadOptions,
        options: &ReportOptionsInput,
    ) -> Result<Vec<serde_json::Value>, JsError> {
        let mut outputs = Vec::new();
        for interpretation in interpretations {
            if interpretation.kind != "bioscript" {
                return Err(JsError::new(&format!(
                    "analysis '{}' uses unsupported kind '{}'",
                    interpretation.id, interpretation.kind
                )));
            }
            let script_path = self.resolve(manifest_path, &interpretation.path)?;
            let output_file = format!(
                "analysis/{participant_id}/{}.{}",
                interpretation.id,
                interpretation.output_format.as_deref().unwrap_or("json")
            );
            let mut virtual_text_files = self.files.clone();
            let mut virtual_binary_files = BTreeMap::new();
            virtual_binary_files.insert(input_name.to_owned(), input_bytes.to_vec());
            let limits = ResourceLimits::new()
                .max_duration(Duration::from_millis(options.analysis_max_duration_ms))
                .max_memory(16 * 1024 * 1024)
                .max_allocations(400_000)
                .gc_interval(1000)
                .max_recursion_depth(Some(200));
            let runtime = BioscriptRuntime::with_config(
                PathBuf::new(),
                RuntimeConfig {
                    limits,
                    loader: loader.clone(),
                    virtual_binary_files,
                    virtual_text_files: std::mem::take(&mut virtual_text_files),
                },
            )
            .map_err(|err| JsError::new(&format!("create analysis runtime failed: {err:?}")))?;
            runtime
                .run_file(
                    &script_path,
                    None,
                    vec![
                        ("input_file", MontyObject::String(input_name.to_owned())),
                        ("output_file", MontyObject::String(output_file.clone())),
                        ("participant_id", MontyObject::String(participant_id.to_owned())),
                    ],
                )
                .map_err(|err| JsError::new(&format!("analysis {} failed: {err:?}", interpretation.id)))?;
            let written = runtime.virtual_written_text_files();
            let text = written
                .get(&output_file)
                .ok_or_else(|| JsError::new(&format!("analysis {} did not write {output_file}", interpretation.id)))?;
            let format = interpretation
                .output_format
                .as_deref()
                .unwrap_or("json")
                .to_ascii_lowercase();
            let (rows, row_headers) = parse_analysis_output_text(text, &format)?;
            outputs.push(serde_json::json!({
                "schema": "bioscript:analysis-output:1.0",
                "version": "1.0",
                "participant_id": participant_id,
                "assay_id": manifest_name,
                "analysis_id": interpretation.id,
                "analysis_label": interpretation.label,
                "kind": interpretation.kind,
                "output_format": format,
                "manifest_path": manifest_path,
                "script_path": script_path,
                "output_file": output_file,
                "derived_from": interpretation.derived_from,
                "emits": interpretation.emits.iter().map(|emit| serde_json::json!({
                    "key": emit.key,
                    "label": emit.label,
                    "value_type": emit.value_type,
                    "format": emit.format,
                })).collect::<Vec<_>>(),
                "logic": interpretation.logic.as_ref().map(|logic| serde_json::json!({
                    "description": logic.description,
                    "source": logic.source.as_ref().map(|source| serde_json::json!({
                        "name": source.name,
                        "url": source.url,
                    })),
                })),
                "row_headers": row_headers,
                "rows": rows,
            }));
        }
        Ok(outputs)
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
        let ref_allele = manifest.spec.reference.clone().unwrap_or_default();
        let alt_allele = manifest.spec.alternate.clone().unwrap_or_default();
        let genotype_display = row.get("genotype").cloned().unwrap_or_default();
        let assembly = row
            .get("assembly")
            .filter(|value| !value.is_empty())
            .cloned()
            .or_else(|| fallback_assembly.map(assembly_row_value))
            .unwrap_or_default();
        let locus = if assembly.eq_ignore_ascii_case("grch37") {
            manifest.spec.grch37.as_ref()
        } else {
            manifest.spec.grch38.as_ref().or(manifest.spec.grch37.as_ref())
        };
        let chrom = locus.map_or(String::new(), |locus| locus.chrom.clone());
        let (genotype, zygosity) =
            normalize_app_genotype(&genotype_display, &ref_allele, &alt_allele, &chrom, inferred_sex);
        let outcome = if genotype == "./." {
            "no_call"
        } else if zygosity == "hom_ref" || zygosity == "hem_ref" {
            "reference"
        } else if zygosity == "het" || zygosity == "hom_alt" || zygosity == "hem_alt" {
            "variant"
        } else {
            "unknown"
        };
        Ok(serde_json::json!({
            "participant_id": row.get("participant_id").cloned().unwrap_or_default(),
            "assay_id": assay_id,
            "assay_version": "1.0",
            "variant_key": manifest.name,
            "variant_path": row_path,
            "rsid": row.get("matched_rsid").filter(|value| !value.is_empty()).cloned().or_else(|| manifest.spec.rsids.first().cloned()),
            "gene": gene,
            "assembly": if assembly.is_empty() { serde_json::Value::Null } else { serde_json::Value::String(assembly.to_uppercase()) },
            "chrom": chrom,
            "pos_start": locus.map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.start)),
            "pos_end": locus.map_or(serde_json::Value::Null, |locus| serde_json::Value::from(locus.end)),
            "ref": ref_allele,
            "alt": alt_allele,
            "kind": manifest.spec.kind.map_or("unknown".to_owned(), |kind| format!("{kind:?}").to_lowercase()),
            "match_status": if row.get("matched_rsid").is_some_and(|value| !value.is_empty()) || !genotype_display.is_empty() { "found" } else { "not_found" },
            "coverage_status": "covered",
            "call_status": if genotype == "./." { "no_call" } else { "called" },
            "genotype": genotype,
            "genotype_display": genotype_display,
            "zygosity": zygosity,
            "ref_count": parse_optional_u32(row.get("ref_count")),
            "alt_count": parse_optional_u32(row.get("alt_count")),
            "depth": parse_optional_u32(row.get("depth")),
            "genotype_quality": serde_json::Value::Null,
            "allele_balance": serde_json::Value::Null,
            "outcome": outcome,
            "evidence_type": "genotype_file",
            "evidence_raw": row.get("evidence").cloned().unwrap_or_default(),
            "source": variant_primary_source_from_yaml(&value),
            "facets": serde_json::Value::Null,
        }))
    }

    pub(super) fn report_manifest_metadata(&self, path: &str) -> Result<serde_json::Value, JsError> {
        let value = self.yaml(path)?;
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

    pub(super) fn load_manifest_findings(&self, path: &str) -> Result<Vec<serde_json::Value>, JsError> {
        let value = self.yaml(path)?;
        let schema = yaml_string(&value, "schema").unwrap_or_default();
        let mut findings = Vec::new();
        if matches!(
            schema.as_str(),
            "bioscript:variant:1.0" | "bioscript:variant" | "bioscript:assay:1.0" | "bioscript:panel:1.0" | "bioscript:pgx-findings:1.0"
        ) {
            if let Some(items) = value.get("findings").and_then(serde_yaml::Value::as_sequence) {
                for item in items {
                    let json_item = yaml_to_json(item.clone())?;
                    if let Some(include) = json_item.get("include").and_then(serde_json::Value::as_str) {
                        let include_path = self.resolve(path, include)?;
                        let mut included = self.load_manifest_findings(&include_path)?;
                        let inherited_binding = json_item.get("binding").cloned();
                        for included_item in &mut included {
                            if inherited_binding.is_some()
                                && included_item.get("binding").is_none()
                                && included_item.get("effects").is_none()
                            {
                                if let Some(object) = included_item.as_object_mut() {
                                    object.insert(
                                        "binding".to_owned(),
                                        inherited_binding.clone().unwrap_or(serde_json::Value::Null),
                                    );
                                }
                            }
                        }
                        findings.extend(included);
                    } else {
                        findings.push(json_item);
                    }
                }
            }
        }
        if matches!(schema.as_str(), "bioscript:assay:1.0" | "bioscript:panel:1.0") {
            if let Some(items) = value.get("members").and_then(serde_yaml::Value::as_sequence) {
                for member in items {
                    let Some(kind) = member.get("kind").and_then(serde_yaml::Value::as_str) else { continue };
                    if !matches!(kind, "variant" | "assay") { continue; }
                    let Some(member_path) = member.get("path").and_then(serde_yaml::Value::as_str) else { continue };
                    let resolved = self.resolve(path, member_path)?;
                    findings.extend(self.load_manifest_findings(&resolved)?);
                }
            }
        }
        Ok(findings)
    }

    pub(super) fn load_manifest_provenance_links(&self, path: &str) -> Result<Vec<serde_json::Value>, JsError> {
        let value = self.yaml(path)?;
        let schema = yaml_string(&value, "schema").unwrap_or_default();
        let mut links = BTreeMap::<String, serde_json::Value>::new();
        collect_manifest_provenance_entries(&value, &mut links)?;
        if matches!(schema.as_str(), "bioscript:assay:1.0" | "bioscript:panel:1.0") {
            if let Some(items) = value.get("members").and_then(serde_yaml::Value::as_sequence) {
                for member in items {
                    let Some(kind) = member.get("kind").and_then(serde_yaml::Value::as_str) else { continue };
                    if !matches!(kind, "variant" | "assay") { continue; }
                    let Some(member_path) = member.get("path").and_then(serde_yaml::Value::as_str) else { continue };
                    let resolved = self.resolve(path, member_path)?;
                    for item in self.load_manifest_provenance_links(&resolved)? {
                        if let Some(url) = item.get("url").and_then(serde_json::Value::as_str) {
                            links.entry(url.to_owned()).or_insert(item);
                        }
                    }
                }
            }
        }
        Ok(links.into_values().collect())
    }
}
