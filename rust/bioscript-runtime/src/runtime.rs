use std::{
    collections::BTreeMap,
    fs,
    path::{Component, Path, PathBuf},
    sync::Arc,
    time::{Duration, Instant},
};

use bioscript_core::RuntimeError;
use monty::{LimitedTracker, MontyObject, MontyRun, NameLookupResult, PrintWriter, RunProgress};

mod args;
mod host_io;
mod methods;
mod objects;
mod state;
mod trace;
mod variants;

#[cfg(test)]
use bioscript_core::VariantSpec;
use host_io::{deepest_existing_ancestor, host_read_text, host_write_text};
use objects::bioscript_object;
#[cfg(test)]
use objects::{
    genotype_file_object, variant_object, variant_observation_object, variant_plan_object,
};
pub use state::{RuntimeConfig, StageTiming};
use state::{RuntimeState, monty_error};
#[cfg(test)]
use trace::{
    ends_with_unescaped_backslash, extract_coordinate, extract_rsid, update_nesting_depth,
};
use trace::{host_trace, instrument_source, statement_context, trace_lookup_metadata};
#[cfg(test)]
use variants::{dataclass_handle_id, dataclass_to_variant_spec, variant_specs_from_plan};
#[cfg(test)]
use variants::{int_from_optional, string_from_optional, string_list_from_object, string_or_list};

type HostFunction = fn(
    &BioscriptRuntime,
    &[MontyObject],
    &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError>;

#[derive(Clone)]
pub struct BioscriptRuntime {
    root: PathBuf,
    config: RuntimeConfig,
    functions: BTreeMap<&'static str, HostFunction>,
    state: Arc<RuntimeState>,
}

impl BioscriptRuntime {
    pub fn new(root: impl Into<PathBuf>) -> Result<Self, RuntimeError> {
        Self::with_config(root, RuntimeConfig::default())
    }

    pub fn with_config(
        root: impl Into<PathBuf>,
        config: RuntimeConfig,
    ) -> Result<Self, RuntimeError> {
        let root = root.into();
        let canonical_root = root.canonicalize().map_err(|err| {
            RuntimeError::Io(format!(
                "failed to canonicalize bioscript root {}: {err}",
                root.display()
            ))
        })?;

        let mut functions = BTreeMap::new();
        functions.insert("read_text", host_read_text as HostFunction);
        functions.insert("write_text", host_write_text as HostFunction);
        functions.insert("__bioscript_trace__", host_trace as HostFunction);

        Ok(Self {
            root: canonical_root,
            config,
            functions,
            state: Arc::new(RuntimeState::new()),
        })
    }

    #[must_use]
    pub fn root(&self) -> &Path {
        &self.root
    }

    #[must_use]
    pub fn config(&self) -> &RuntimeConfig {
        &self.config
    }

    pub fn run_file(
        &self,
        script_path: impl AsRef<Path>,
        trace_report_path: Option<&Path>,
        mut extra_inputs: Vec<(&str, MontyObject)>,
    ) -> Result<MontyObject, RuntimeError> {
        let run_started = Instant::now();
        let script_path = script_path.as_ref();
        let code = fs::read_to_string(script_path).map_err(|err| {
            RuntimeError::Io(format!(
                "failed to read script {}: {err}",
                script_path.display()
            ))
        })?;
        let instrumented = instrument_source(&code);
        self.state
            .trace_lines
            .lock()
            .expect("trace mutex poisoned")
            .clear();
        self.state
            .timings
            .lock()
            .expect("timings mutex poisoned")
            .clear();

        extra_inputs.push(("__name__", MontyObject::String("__main__".to_owned())));
        extra_inputs.push((
            "__file__",
            MontyObject::String(script_path.display().to_string()),
        ));
        extra_inputs.push(("bioscript", bioscript_object()));

        let result = self.run_script(
            &instrumented,
            &script_path.display().to_string(),
            extra_inputs,
        )?;

        if let Some(report_path) = trace_report_path {
            self.write_trace_report(report_path, &code)?;
        }
        self.record_timing(
            "run_file_total",
            run_started.elapsed(),
            format!("script={}", script_path.display()),
        );

        Ok(result)
    }

    #[must_use]
    pub fn timing_snapshot(&self) -> Vec<StageTiming> {
        self.state
            .timings
            .lock()
            .expect("timings mutex poisoned")
            .clone()
    }

    pub fn run_script(
        &self,
        code: &str,
        script_name: &str,
        inputs: Vec<(&str, MontyObject)>,
    ) -> Result<MontyObject, RuntimeError> {
        let input_names = inputs.iter().map(|(name, _)| (*name).to_owned()).collect();
        let input_values = inputs.into_iter().map(|(_, value)| value).collect();
        let runner =
            MontyRun::new(code.to_owned(), script_name, input_names).map_err(monty_error)?;
        let tracker = LimitedTracker::new(self.config.limits.clone());
        let mut progress = runner
            .start(input_values, tracker, PrintWriter::Stdout)
            .map_err(monty_error)?;

        loop {
            progress = match progress {
                RunProgress::Complete(value) => return Ok(value),
                RunProgress::NameLookup(lookup) => {
                    let result = if self.functions.contains_key(lookup.name.as_str()) {
                        NameLookupResult::Value(MontyObject::Function {
                            name: lookup.name.clone(),
                            docstring: None,
                        })
                    } else {
                        NameLookupResult::Undefined
                    };
                    lookup
                        .resume(result, PrintWriter::Stdout)
                        .map_err(monty_error)?
                }
                RunProgress::FunctionCall(call) => {
                    if call.method_call {
                        let result = self.dispatch_method_call(
                            &call.function_name,
                            &call.args,
                            &call.kwargs,
                        )?;
                        call.resume(result, PrintWriter::Stdout)
                            .map_err(monty_error)?
                    } else {
                        let Some(handler) = self.functions.get(call.function_name.as_str()) else {
                            return Err(RuntimeError::Unsupported(format!(
                                "unknown bioscript host function: {}",
                                call.function_name
                            )));
                        };
                        let result = handler(self, &call.args, &call.kwargs)?;
                        call.resume(result, PrintWriter::Stdout)
                            .map_err(monty_error)?
                    }
                }
                RunProgress::ResolveFutures(state) => {
                    return Err(RuntimeError::Unsupported(format!(
                        "async futures are not supported in bioscript runtime: {:?}",
                        state.pending_call_ids()
                    )));
                }
                RunProgress::OsCall(call) => {
                    return Err(RuntimeError::Unsupported(format!(
                        "OS call {} is blocked",
                        call.function
                    )));
                }
            };
        }
    }

    fn dispatch_method_call(
        &self,
        method_name: &str,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let class_name = match args.first() {
            Some(MontyObject::Dataclass { name, .. }) => name.as_str(),
            _ => "<unknown>",
        };

        match (class_name, method_name) {
            ("Bioscript", "load_genotypes") => self.method_load_genotypes(args, kwargs),
            ("Bioscript", "variant") => self.method_variant(args, kwargs),
            ("Bioscript", "query_plan") => self.method_query_plan(args, kwargs),
            ("Bioscript", "write_tsv") => self.method_write_tsv(args, kwargs),
            ("Bioscript", "read_text") => self.method_read_text(args, kwargs),
            ("Bioscript", "write_text") => self.method_write_text(args, kwargs),
            ("GenotypeFile", "get") => self.method_genotype_get(args, kwargs),
            ("GenotypeFile", "lookup_variant") => self.method_genotype_lookup_variant(args, kwargs),
            ("GenotypeFile", "lookup_variant_details") => {
                self.method_genotype_lookup_variant_details(args, kwargs)
            }
            ("GenotypeFile", "lookup_variants") => {
                self.method_genotype_lookup_variants(args, kwargs)
            }
            ("GenotypeFile", "lookup_variants_details") => {
                self.method_genotype_lookup_variants_details(args, kwargs)
            }
            _ => Err(RuntimeError::Unsupported(format!(
                "'{class_name}' object has no attribute '{method_name}'"
            ))),
        }
    }

    fn record_timing(&self, stage: &str, duration: Duration, detail: String) {
        self.state
            .timings
            .lock()
            .expect("timings mutex poisoned")
            .push(StageTiming {
                stage: stage.to_owned(),
                duration_ms: duration.as_millis(),
                detail,
            });
    }

    fn resolve_user_path(&self, raw_path: &str) -> Result<PathBuf, RuntimeError> {
        let path = Path::new(raw_path);
        if path.is_absolute() {
            return Err(RuntimeError::InvalidArguments(format!(
                "absolute paths are not allowed: {raw_path}"
            )));
        }
        for component in path.components() {
            match component {
                Component::ParentDir | Component::RootDir | Component::Prefix(_) => {
                    return Err(RuntimeError::InvalidArguments(format!(
                        "path escapes bioscript root: {raw_path}"
                    )));
                }
                Component::CurDir | Component::Normal(_) => {}
            }
        }
        Ok(self.root.join(path))
    }

    fn resolve_existing_user_path(&self, raw_path: &str) -> Result<PathBuf, RuntimeError> {
        let path = self.resolve_user_path(raw_path)?;
        let canonical = path.canonicalize().map_err(|err| {
            RuntimeError::Io(format!("failed to resolve {}: {err}", path.display()))
        })?;
        self.ensure_under_root(&canonical, raw_path)?;
        Ok(canonical)
    }

    fn resolve_user_write_path(&self, raw_path: &str) -> Result<PathBuf, RuntimeError> {
        let path = self.resolve_user_path(raw_path)?;
        if path.exists() {
            let canonical = path.canonicalize().map_err(|err| {
                RuntimeError::Io(format!("failed to resolve {}: {err}", path.display()))
            })?;
            self.ensure_under_root(&canonical, raw_path)?;
            return Ok(canonical);
        }

        let parent = path.parent().unwrap_or(&self.root);
        let existing_parent = deepest_existing_ancestor(parent);
        let canonical_parent = existing_parent.canonicalize().map_err(|err| {
            RuntimeError::Io(format!(
                "failed to resolve parent dir {}: {err}",
                existing_parent.display()
            ))
        })?;
        self.ensure_under_root(&canonical_parent, raw_path)?;
        Ok(path)
    }

    fn ensure_under_root(&self, path: &Path, raw_path: &str) -> Result<(), RuntimeError> {
        if path.starts_with(&self.root) {
            Ok(())
        } else {
            Err(RuntimeError::InvalidArguments(format!(
                "path escapes bioscript root: {raw_path}"
            )))
        }
    }

    fn write_trace_report(
        &self,
        report_path: &Path,
        original_code: &str,
    ) -> Result<(), RuntimeError> {
        let trace_lines = self
            .state
            .trace_lines
            .lock()
            .expect("trace mutex poisoned")
            .clone();
        let lines: Vec<&str> = original_code.lines().collect();
        let mut output = String::from("step\tline\tcode\tlookup_key\tlookup_url\n");
        for (idx, line_no) in trace_lines.iter().enumerate() {
            let source = lines.get(line_no.saturating_sub(1)).copied().unwrap_or("");
            let trimmed = source.trim();
            let statement = statement_context(&lines, *line_no);
            let (lookup_key, lookup_url) = trace_lookup_metadata(&statement);
            output.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\n",
                idx + 1,
                line_no,
                trimmed,
                lookup_key.unwrap_or_default(),
                lookup_url.unwrap_or_default()
            ));
        }
        if let Some(parent) = report_path.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to create parent dir {}: {err}",
                    parent.display()
                ))
            })?;
        }
        fs::write(report_path, output).map_err(|err| {
            RuntimeError::Io(format!("failed to write {}: {err}", report_path.display()))
        })?;
        Ok(())
    }
}

fn resolve_optional_loader_path(
    runtime: &BioscriptRuntime,
    path: Option<PathBuf>,
) -> Result<Option<PathBuf>, RuntimeError> {
    path.map(|path| {
        if path.is_absolute() {
            Ok(path)
        } else {
            runtime.resolve_user_path(&path.to_string_lossy())
        }
    })
    .transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_dir(label: &str) -> PathBuf {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock drift")
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "bioscript-runtime-unit-{label}-{}-{nanos}",
            std::process::id()
        ));
        fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn attr<'a>(obj: &'a MontyObject, name: &str) -> Option<&'a MontyObject> {
        let MontyObject::Dataclass { attrs, .. } = obj else {
            return None;
        };
        attrs.into_iter().find_map(|(key, value)| {
            matches!(key, MontyObject::String(text) if text == name).then_some(value)
        })
    }

    #[test]
    fn trace_helpers_cover_coordinates_rsids_and_statement_edges() {
        assert_eq!(
            trace_lookup_metadata("bioscript.variant(rsid='rs12345')").0,
            Some("rs12345".to_owned())
        );
        let (key, url) = trace_lookup_metadata("bioscript.variant(grch37='chr1:10-20')");
        assert_eq!(key.as_deref(), Some("1:10-20"));
        assert!(url.unwrap().starts_with("https://grch37.ensembl.org"));
        let (key, _) = trace_lookup_metadata("bioscript.variant(grch38='2:30')");
        assert_eq!(key.as_deref(), Some("2:30-30"));
        assert_eq!(trace_lookup_metadata("no lookup here"), (None, None));

        let lines = ["plan = bioscript.query_plan([", "    RS1,", "])"];
        assert_eq!(
            statement_context(&lines, 1),
            "plan = bioscript.query_plan([ RS1, ])"
        );
        assert_eq!(statement_context(&lines, 0), "");
        assert_eq!(statement_context(&lines, 9), "");

        assert_eq!(extract_rsid("x rs42 y"), Some("rs42".to_owned()));
        assert_eq!(extract_rsid("notrs42"), None);
        assert_eq!(extract_coordinate("chrX:7;"), Some("X:7-7".to_owned()));
        assert_eq!(extract_coordinate("chrM:7-8"), Some("M:7-8".to_owned()));
        assert_eq!(extract_coordinate("chr1:x-y"), None);
    }

    #[test]
    fn instrument_source_tracks_continuations_comments_and_strings() {
        let source = "value = (\n    'not ) counted'\n)\n# skip\nnext_value = 1\\\n    + 2\n";
        let instrumented = instrument_source(source);

        assert!(instrumented.contains("__bioscript_trace__(1)\nvalue = ("));
        assert!(!instrumented.contains("__bioscript_trace__(2)"));
        assert!(!instrumented.contains("__bioscript_trace__(3)"));
        assert!(!instrumented.contains("__bioscript_trace__(4)"));
        assert!(instrumented.contains("__bioscript_trace__(5)\nnext_value = 1\\"));
        assert!(!instrumented.contains("__bioscript_trace__(6)"));
        assert!(instrumented.ends_with('\n'));

        assert!(ends_with_unescaped_backslash("x = 1\\"));
        assert!(!ends_with_unescaped_backslash("x = '\\\\'"));
        assert_eq!(
            update_nesting_depth(0, "call(') still string', [1]) # ignored"),
            0
        );
        assert_eq!(update_nesting_depth(0, "call(["), 2);
        assert_eq!(update_nesting_depth(2, "])"), 0);
    }

    #[test]
    fn object_helpers_cover_optional_fields_and_conversion_errors() {
        let observation = bioscript_core::VariantObservation {
            backend: "vcf".to_owned(),
            matched_rsid: Some("rs1".to_owned()),
            assembly: Some(bioscript_core::Assembly::Grch37),
            genotype: Some("AG".to_owned()),
            ref_count: Some(3),
            alt_count: Some(2),
            depth: Some(5),
            raw_counts: BTreeMap::from([("A".to_owned(), 3), ("G".to_owned(), 2)]),
            decision: Some("heterozygous".to_owned()),
            evidence: vec!["resolved".to_owned()],
        };
        let object = variant_observation_object(&observation);
        assert!(matches!(attr(&object, "assembly"), Some(MontyObject::String(v)) if v == "grch37"));
        assert!(matches!(
            attr(&object, "ref_count"),
            Some(MontyObject::Int(3))
        ));
        assert!(matches!(
            attr(&object, "alt_count"),
            Some(MontyObject::Int(2))
        ));
        assert!(matches!(attr(&object, "depth"), Some(MontyObject::Int(5))));

        let missing = variant_observation_object(&bioscript_core::VariantObservation::default());
        assert!(matches!(
            attr(&missing, "assembly"),
            Some(MontyObject::None)
        ));
        assert!(matches!(
            attr(&missing, "genotype"),
            Some(MontyObject::None)
        ));

        assert_eq!(
            string_or_list(&MontyObject::None).unwrap(),
            Vec::<String>::new()
        );
        assert_eq!(
            string_list_from_object(&MontyObject::None).unwrap(),
            Vec::<String>::new()
        );
        assert_eq!(string_from_optional(&MontyObject::None).unwrap(), None);
        assert_eq!(int_from_optional(&MontyObject::None).unwrap(), None);
        assert!(string_list_from_object(&MontyObject::String("x".to_owned())).is_err());

        let bad_plan = MontyObject::Dataclass {
            name: "VariantPlan".to_owned(),
            type_id: 4,
            field_names: vec![],
            attrs: vec![].into(),
            frozen: true,
        };
        assert!(
            variant_specs_from_plan(&bad_plan)
                .unwrap_err()
                .to_string()
                .contains("missing variants")
        );

        let bad_variant = MontyObject::Dataclass {
            name: "Other".to_owned(),
            type_id: 9,
            field_names: vec![],
            attrs: vec![].into(),
            frozen: true,
        };
        assert!(
            dataclass_to_variant_spec(&bad_variant)
                .unwrap_err()
                .to_string()
                .contains("got Other")
        );
        assert!(
            dataclass_to_variant_spec(&MontyObject::None)
                .unwrap_err()
                .to_string()
                .contains("expected Variant object")
        );
    }

    #[test]
    fn runtime_private_methods_cover_dispatch_and_path_errors() {
        let root = temp_dir("dispatch");
        fs::write(root.join("input.txt"), "hello").unwrap();
        let runtime = BioscriptRuntime::new(&root).unwrap();
        let bioscript = MontyObject::Dataclass {
            name: "Bioscript".to_owned(),
            type_id: 1,
            field_names: vec![],
            attrs: vec![].into(),
            frozen: true,
        };

        let err = runtime
            .dispatch_method_call("missing", std::slice::from_ref(&bioscript), &[])
            .unwrap_err();
        assert!(err.to_string().contains("no attribute"));

        assert!(
            runtime
                .method_read_text(&[], &[])
                .unwrap_err()
                .to_string()
                .contains("expects self and path")
        );
        assert!(
            runtime
                .method_write_text(&[], &[])
                .unwrap_err()
                .to_string()
                .contains("expects self, path, text")
        );
        assert!(
            runtime
                .resolve_user_path("/absolute")
                .unwrap_err()
                .to_string()
                .contains("absolute paths")
        );
        assert!(
            runtime
                .resolve_user_path("../escape")
                .unwrap_err()
                .to_string()
                .contains("escapes bioscript root")
        );
        assert!(
            runtime
                .resolve_existing_user_path("missing.txt")
                .unwrap_err()
                .to_string()
                .contains("failed to resolve")
        );

        let nested = runtime
            .resolve_user_write_path("new/deep/file.txt")
            .unwrap();
        assert_eq!(
            nested,
            root.canonicalize().unwrap().join("new/deep/file.txt")
        );
        assert_eq!(
            deepest_existing_ancestor(&root.join("new/deep/file.txt")),
            root.as_path()
        );

        let mut config = RuntimeConfig::default();
        config.loader.input_index = Some(PathBuf::from("input.txt"));
        config.loader.reference_file = Some(PathBuf::from("/tmp/reference.fa"));
        let runtime = BioscriptRuntime::with_config(&root, config).unwrap();
        let loader = runtime.resolved_loader_options().unwrap();
        assert_eq!(
            loader.input_index.as_deref(),
            Some(root.canonicalize().unwrap().join("input.txt").as_path())
        );
        assert_eq!(
            loader.reference_file.as_deref(),
            Some(Path::new("/tmp/reference.fa"))
        );
    }

    #[test]
    fn runtime_private_methods_cover_unknown_genotype_handles() {
        let root = temp_dir("handles");
        let runtime = BioscriptRuntime::new(&root).unwrap();
        let genotype = genotype_file_object(99);
        let variant = variant_object(&VariantSpec {
            rsids: vec!["rs1".to_owned()],
            ..VariantSpec::default()
        });
        let plan = variant_plan_object(&[VariantSpec {
            rsids: vec!["rs1".to_owned()],
            ..VariantSpec::default()
        }]);

        for (method, args) in [
            (
                "get",
                vec![genotype.clone(), MontyObject::String("rs1".to_owned())],
            ),
            ("lookup_variant", vec![genotype.clone(), variant.clone()]),
            (
                "lookup_variant_details",
                vec![genotype.clone(), variant.clone()],
            ),
            ("lookup_variants", vec![genotype.clone(), plan.clone()]),
            ("lookup_variants_details", vec![genotype.clone(), plan]),
        ] {
            let err = runtime
                .dispatch_method_call(method, &args, &[])
                .unwrap_err()
                .to_string();
            assert!(err.contains("unknown genotype handle"), "{method}: {err}");
        }

        assert!(
            dataclass_handle_id(&MontyObject::None, "GenotypeFile")
                .unwrap_err()
                .to_string()
                .contains("expected GenotypeFile object")
        );
        let missing_handle = MontyObject::Dataclass {
            name: "GenotypeFile".to_owned(),
            type_id: 2,
            field_names: vec![],
            attrs: vec![].into(),
            frozen: true,
        };
        assert!(
            dataclass_handle_id(&missing_handle, "GenotypeFile")
                .unwrap_err()
                .to_string()
                .contains("missing handle_id")
        );
    }

    #[test]
    fn runtime_private_methods_cover_successful_genotype_paths() {
        let root = temp_dir("success-paths");
        fs::write(
            root.join("genotypes.tsv"),
            "rsid\tchromosome\tposition\tgenotype\nrs1\t1\t10\tAG\nrs2\t2\t20\tCT\n",
        )
        .unwrap();
        let runtime = BioscriptRuntime::new(&root).unwrap();
        let bioscript = bioscript_object();

        let genotype = runtime
            .method_load_genotypes(
                &[
                    bioscript.clone(),
                    MontyObject::String("genotypes.tsv".to_owned()),
                ],
                &[],
            )
            .unwrap();
        assert!(matches!(
            attr(&genotype, "handle_id"),
            Some(MontyObject::Int(id)) if *id > 0
        ));

        let value = runtime
            .method_genotype_get(
                &[genotype.clone(), MontyObject::String("rs1".to_owned())],
                &[],
            )
            .unwrap();
        assert!(matches!(value, MontyObject::String(ref text) if text == "AG"));
        let missing = runtime
            .method_genotype_get(
                &[genotype.clone(), MontyObject::String("missing".to_owned())],
                &[],
            )
            .unwrap();
        assert!(matches!(missing, MontyObject::None));

        let variant = runtime
            .method_variant(
                std::slice::from_ref(&bioscript),
                &[
                    (
                        MontyObject::String("rsid".to_owned()),
                        MontyObject::String("rs2".to_owned()),
                    ),
                    (
                        MontyObject::String("grch38".to_owned()),
                        MontyObject::String("2:20".to_owned()),
                    ),
                    (
                        MontyObject::String("kind".to_owned()),
                        MontyObject::String("snp".to_owned()),
                    ),
                ],
            )
            .unwrap();
        let lookup = runtime
            .method_genotype_lookup_variant(&[genotype.clone(), variant.clone()], &[])
            .unwrap();
        assert!(matches!(lookup, MontyObject::String(ref text) if text == "CT"));
        let details = runtime
            .method_genotype_lookup_variant_details(&[genotype.clone(), variant.clone()], &[])
            .unwrap();
        assert!(matches!(
            attr(&details, "matched_rsid"),
            Some(MontyObject::String(text)) if text == "rs2"
        ));

        let plan = runtime
            .method_query_plan(&[bioscript.clone(), MontyObject::List(vec![variant])], &[])
            .unwrap();
        let values = runtime
            .method_genotype_lookup_variants(&[genotype.clone(), plan.clone()], &[])
            .unwrap();
        assert!(matches!(values, MontyObject::List(items) if items.len() == 1));
        let detail_values = runtime
            .method_genotype_lookup_variants_details(&[genotype, plan], &[])
            .unwrap();
        assert!(matches!(detail_values, MontyObject::List(items) if items.len() == 1));
        assert!(runtime.timing_snapshot().len() >= 4);
    }

    #[test]
    fn runtime_private_methods_cover_successful_text_tsv_and_trace_paths() {
        let root = temp_dir("host-output");
        let runtime = BioscriptRuntime::new(&root).unwrap();
        let bioscript = bioscript_object();
        let rows = MontyObject::List(vec![MontyObject::Dict(
            vec![
                (
                    MontyObject::String("rsid".to_owned()),
                    MontyObject::String("rs1".to_owned()),
                ),
                (MontyObject::String("count".to_owned()), MontyObject::Int(2)),
                (
                    MontyObject::String("ok".to_owned()),
                    MontyObject::Bool(true),
                ),
                (MontyObject::String("note".to_owned()), MontyObject::None),
            ]
            .into(),
        )]);
        runtime
            .method_write_tsv(
                &[
                    bioscript.clone(),
                    MontyObject::String("out/table.tsv".to_owned()),
                    rows,
                ],
                &[],
            )
            .unwrap();
        let table = fs::read_to_string(root.join("out/table.tsv")).unwrap();
        assert!(table.contains("rs1"));
        assert!(table.contains("true"));

        runtime
            .method_write_text(
                &[
                    bioscript.clone(),
                    MontyObject::String("out/text.txt".to_owned()),
                    MontyObject::String("hello".to_owned()),
                ],
                &[],
            )
            .unwrap();
        let text = runtime
            .method_read_text(
                &[bioscript, MontyObject::String("out/text.txt".to_owned())],
                &[],
            )
            .unwrap();
        assert!(matches!(text, MontyObject::String(ref value) if value == "hello"));

        runtime.state.trace_lines.lock().unwrap().extend([1, 2, 99]);
        runtime
            .write_trace_report(
                &root.join("trace/report.tsv"),
                "bioscript.variant(rsid='rs1')\nplain = 1\n",
            )
            .unwrap();
        let trace = fs::read_to_string(root.join("trace/report.tsv")).unwrap();
        assert!(trace.contains("https://www.ncbi.nlm.nih.gov/snp/rs1"));
        assert_eq!(runtime.timing_snapshot().len(), 1);
    }
}
