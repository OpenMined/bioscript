use std::collections::BTreeMap;

use bioscript_core::{Assembly, VariantObservation, VariantSpec};
use bioscript_formats::SexInference;
use bioscript_schema::load_variant_manifest_text;

use crate::{
    AppInputReportInput, FilesystemManifestWorkspace, ManifestWorkspace, app_input_report_json,
    app_observation_from_manifest_row, collect_analysis_manifest_tasks,
    collect_manifest_provenance_entries, collect_variant_manifest_tasks,
    load_report_manifest_context, render_input_report_artifact_texts, variant_row,
};

pub trait ReportVariantLookup {
    fn lookup_variants(&self, specs: &[VariantSpec]) -> Result<Vec<VariantObservation>, String>;
}

pub trait ReportWorkspace: ManifestWorkspace {
    fn app_observation_from_manifest_row(
        &self,
        row: &BTreeMap<String, String>,
        assay_id: &str,
        inferred_sex: Option<&SexInference>,
        fallback_assembly: Option<Assembly>,
    ) -> Result<serde_json::Value, String>;
}

pub trait ReportAnalysisRunner {
    fn run_analysis_task(
        &self,
        task: &crate::AnalysisManifestTask,
        observation_rows: &[BTreeMap<String, String>],
        variant_observations: &[VariantObservation],
        observations: &[serde_json::Value],
    ) -> Result<Vec<serde_json::Value>, String>;
}

pub struct NoopReportAnalysisRunner;

impl ReportAnalysisRunner for NoopReportAnalysisRunner {
    fn run_analysis_task(
        &self,
        _task: &crate::AnalysisManifestTask,
        _observation_rows: &[BTreeMap<String, String>],
        _variant_observations: &[VariantObservation],
        _observations: &[serde_json::Value],
    ) -> Result<Vec<serde_json::Value>, String> {
        Ok(Vec::new())
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ReportInputContext<'a> {
    pub participant_id: &'a str,
    pub input_file_name: &'a str,
    pub input_file_path: &'a str,
    pub input_inspection: Option<&'a bioscript_formats::FileInspection>,
}

#[derive(Clone, Copy, Debug, Default)]
pub struct ReportRunOptions<'a> {
    pub filters: &'a [String],
}

#[derive(Clone, Debug)]
pub struct ReportRunResult {
    pub observation_rows: Vec<BTreeMap<String, String>>,
    pub variant_observations: Vec<VariantObservation>,
    pub observations: Vec<serde_json::Value>,
    pub analyses: Vec<serde_json::Value>,
    pub report: serde_json::Value,
    pub artifacts: crate::ReportArtifactTexts,
}

pub fn run_report(
    workspace: &impl ReportWorkspace,
    manifest_path: &str,
    lookup: &impl ReportVariantLookup,
    analysis_runner: &impl ReportAnalysisRunner,
    input: ReportInputContext<'_>,
    options: ReportRunOptions<'_>,
) -> Result<ReportRunResult, String> {
    let manifest_context = load_report_manifest_context(workspace, manifest_path)?;
    let variant_tasks = collect_variant_manifest_tasks(workspace, manifest_path, options.filters)?;
    let variant_observations = lookup.lookup_variants(
        &variant_tasks
            .iter()
            .map(|task| task.manifest.spec.clone())
            .collect::<Vec<_>>(),
    )?;
    let observation_rows = variant_tasks
        .into_iter()
        .zip(variant_observations.iter())
        .map(|(task, observation)| {
            variant_row(
                &task.manifest_path,
                &task.manifest.name,
                &task.manifest.tags,
                observation,
                input.participant_id,
            )
        })
        .collect::<Vec<_>>();
    let observations = observation_rows
        .iter()
        .map(|row| {
            workspace.app_observation_from_manifest_row(
                row,
                &manifest_context.assay_id,
                input
                    .input_inspection
                    .and_then(|inspection| inspection.inferred_sex.as_ref()),
                input
                    .input_inspection
                    .and_then(|inspection| inspection.assembly),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let analysis_variant_observations =
        analysis_variant_observations(&variant_observations, &observations);
    let mut analyses = Vec::new();
    for task in collect_analysis_manifest_tasks(workspace, manifest_path, options.filters)? {
        analyses.extend(analysis_runner.run_analysis_task(
            &task,
            &observation_rows,
            &analysis_variant_observations,
            &observations,
        )?);
    }
    let report_input = AppInputReportInput {
        assay_id: &manifest_context.assay_id,
        participant_id: input.participant_id,
        input_file_name: input.input_file_name,
        input_file_path: input.input_file_path,
        observations: &observations,
        analyses: &analyses,
        findings: &manifest_context.findings,
        provenance: &manifest_context.provenance,
        input_inspection: input.input_inspection,
        manifest_metadata: &manifest_context.manifest_metadata,
    };
    let report = app_input_report_json(report_input);
    let artifacts = render_input_report_artifact_texts(report_input)?;
    Ok(ReportRunResult {
        observation_rows,
        variant_observations,
        observations,
        analyses,
        report,
        artifacts,
    })
}

impl ReportVariantLookup for bioscript_formats::GenotypeStore {
    fn lookup_variants(&self, specs: &[VariantSpec]) -> Result<Vec<VariantObservation>, String> {
        bioscript_formats::GenotypeStore::lookup_variants(self, specs)
            .map_err(|err| err.to_string())
    }
}

impl ReportWorkspace for FilesystemManifestWorkspace {
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
            let task = crate::load_variant_manifest_task_by_path(self, &row_path)?;
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
            let text = self.load_text(&row_path)?;
            let manifest = load_variant_manifest_text(&row_path, &text)?;
            let value = self.load_yaml(&row_path)?;
            (
                manifest,
                yaml_string(&value, "gene").unwrap_or_default(),
                variant_primary_source_from_yaml(&value)?,
                variant_alt_alleles_from_yaml(&value),
                variant_observed_alt_alleles_from_yaml(&value),
            )
        };
        Ok(app_observation_from_manifest_row(
            crate::AppObservationInput {
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

fn variant_primary_source_from_yaml(
    value: &serde_yaml::Value,
) -> Result<serde_json::Value, String> {
    let mut links = BTreeMap::<String, serde_json::Value>::new();
    collect_manifest_provenance_entries(value, &mut links)?;
    if let Some(source) = links
        .values()
        .find(|source| source_url_contains(source, "ncbi.nlm.nih.gov/snp/rs"))
    {
        return Ok(source.clone());
    }
    if let Some(rsid) = value
        .get("identifiers")
        .and_then(|identifiers| identifiers.get("rsids"))
        .and_then(serde_yaml::Value::as_sequence)
        .and_then(|items| items.iter().find_map(serde_yaml::Value::as_str))
    {
        return Ok(serde_json::json!({
            "kind": "database",
            "label": "dbSNP / NCBI SNP",
            "url": format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}"),
            "fields": ["identifiers.rsids"],
        }));
    }
    Ok(links
        .into_values()
        .next()
        .unwrap_or(serde_json::Value::Null))
}

fn source_url_contains(source: &serde_json::Value, needle: &str) -> bool {
    source
        .get("url")
        .and_then(serde_json::Value::as_str)
        .is_some_and(|url| url.contains(needle))
}

fn variant_observed_alt_alleles_from_yaml(value: &serde_yaml::Value) -> Vec<String> {
    value
        .get("alleles")
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("observed_alts".to_owned())))
        .and_then(serde_yaml::Value::as_sequence)
        .into_iter()
        .flatten()
        .filter_map(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
        .collect()
}

fn variant_alt_alleles_from_yaml(value: &serde_yaml::Value) -> Vec<String> {
    value
        .get("alleles")
        .and_then(serde_yaml::Value::as_mapping)
        .and_then(|mapping| mapping.get(serde_yaml::Value::String("alts".to_owned())))
        .and_then(serde_yaml::Value::as_sequence)
        .into_iter()
        .flatten()
        .filter_map(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
        .collect()
}

fn yaml_string(value: &serde_yaml::Value, key: &str) -> Option<String> {
    value
        .get(key)
        .and_then(serde_yaml::Value::as_str)
        .map(ToOwned::to_owned)
}

fn analysis_variant_observations(
    variant_observations: &[VariantObservation],
    observations: &[serde_json::Value],
) -> Vec<VariantObservation> {
    variant_observations
        .iter()
        .map(|observation| {
            let mut observation = observation.clone();
            if let Some(app_observation) = matching_app_observation(&observation, observations)
                && let Some(genotype_display) = app_observation
                    .get("genotype_display")
                    .and_then(serde_json::Value::as_str)
                    .filter(|value| !value.is_empty() && *value != "??")
            {
                observation.genotype = Some(genotype_display.to_owned());
            }
            observation
        })
        .collect()
}

fn matching_app_observation<'a>(
    observation: &VariantObservation,
    observations: &'a [serde_json::Value],
) -> Option<&'a serde_json::Value> {
    let matched_rsid = observation.matched_rsid.as_deref()?;
    observations.iter().find(|app_observation| {
        app_observation
            .get("rsid")
            .and_then(serde_json::Value::as_str)
            == Some(matched_rsid)
    })
}

#[cfg(test)]
mod tests {
    use std::{cell::Cell, collections::BTreeMap};

    use bioscript_core::{VariantObservation, VariantSpec};
    use bioscript_formats::SexInference;

    use super::{
        ReportAnalysisRunner, ReportInputContext, ReportRunOptions, ReportVariantLookup,
        ReportWorkspace, run_report,
    };
    use crate::ManifestWorkspace;

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

    impl ReportWorkspace for MapWorkspace {
        fn app_observation_from_manifest_row(
            &self,
            row: &BTreeMap<String, String>,
            assay_id: &str,
            _inferred_sex: Option<&SexInference>,
            _fallback_assembly: Option<bioscript_core::Assembly>,
        ) -> Result<serde_json::Value, String> {
            Ok(serde_json::json!({
                "assay_id": assay_id,
                "participant_id": row.get("participant_id").cloned().unwrap_or_default(),
                "path": row.get("path").cloned().unwrap_or_default(),
                "rsid": row.get("matched_rsid").cloned().unwrap_or_default(),
                "genotype_display": row.get("genotype").cloned().unwrap_or_default(),
            }))
        }
    }

    struct StaticLookup;

    impl ReportVariantLookup for StaticLookup {
        fn lookup_variants(
            &self,
            specs: &[VariantSpec],
        ) -> Result<Vec<VariantObservation>, String> {
            Ok(specs
                .iter()
                .map(|spec| VariantObservation {
                    backend: "test".to_owned(),
                    matched_rsid: spec.rsids.first().cloned(),
                    genotype: Some("AG".to_owned()),
                    ..VariantObservation::default()
                })
                .collect())
        }
    }

    struct CountingAnalysis {
        calls: Cell<usize>,
    }

    impl ReportAnalysisRunner for CountingAnalysis {
        fn run_analysis_task(
            &self,
            task: &crate::AnalysisManifestTask,
            observation_rows: &[BTreeMap<String, String>],
            variant_observations: &[VariantObservation],
            observations: &[serde_json::Value],
        ) -> Result<Vec<serde_json::Value>, String> {
            self.calls.set(self.calls.get() + 1);
            assert_eq!(task.manifest_name, "panel");
            assert_eq!(observation_rows.len(), 1);
            assert_eq!(variant_observations.len(), 1);
            assert_eq!(observations.len(), 1);
            Ok(vec![serde_json::json!({
                "participant_id": "sample",
                "assay_id": task.manifest_name,
                "rows": [{"score": 1}],
            })])
        }
    }

    #[test]
    fn run_report_collects_rows_observations_analyses_report_and_artifacts() {
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
interpretations:
  - id: score
    kind: bioscript
    path: score.bs
    derived_from: [observations]
    emits:
      - key: score
        label: Score
        value_type: integer
        format: number
"#
                    .to_owned(),
                ),
                (
                    "rs1.yaml".to_owned(),
                    r#"
schema: bioscript:variant:1.0
version: "1.0"
name: rs1
identifiers:
  rsids: [rs1]
coordinates:
  grch38:
    chrom: "1"
    pos: 100
alleles:
  kind: snv
  ref: A
  alts: [G]
"#
                    .to_owned(),
                ),
            ]),
        };
        let analysis = CountingAnalysis {
            calls: Cell::new(0),
        };

        let result = run_report(
            &workspace,
            "panel.yaml",
            &StaticLookup,
            &analysis,
            ReportInputContext {
                participant_id: "sample",
                input_file_name: "sample.txt",
                input_file_path: "sample.txt",
                input_inspection: None,
            },
            ReportRunOptions { filters: &[] },
        )
        .unwrap();

        assert_eq!(analysis.calls.get(), 1);
        assert_eq!(result.observation_rows.len(), 1);
        assert_eq!(result.observations[0]["assay_id"], "panel");
        assert_eq!(result.analyses[0]["assay_id"], "panel");
        assert_eq!(result.report["participant_id"], "sample");
        assert!(result.artifacts.observations_tsv.contains("sample"));
        assert!(result.artifacts.analysis_jsonl.contains("\"score\""));
        assert!(result.artifacts.reports_jsonl.contains("\"panel\""));
        assert!(result.artifacts.html.contains("<!doctype html>"));
    }
}
