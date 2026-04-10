use std::{
    collections::{BTreeMap, HashMap},
    error::Error,
    fmt,
    fs,
    path::{Component, Path, PathBuf},
    sync::{
        Arc, Mutex,
        atomic::{AtomicU64, Ordering},
    },
    time::{Duration, Instant},
};

use monty::{
    LimitedTracker, MontyException, MontyObject, MontyRun, NameLookupResult, PrintWriter, ResourceLimits,
    RunProgress,
};

use crate::genotype::{GenotypeLoadOptions, GenotypeStore, GenomicLocus};
use crate::variant::{VariantKind, VariantSpec};

type HostFunction =
    fn(&BioscriptRuntime, &[MontyObject], &[(MontyObject, MontyObject)]) -> Result<MontyObject, RuntimeError>;

#[derive(Debug, Clone)]
pub struct RuntimeConfig {
    pub limits: ResourceLimits,
    pub loader: GenotypeLoadOptions,
}

impl Default for RuntimeConfig {
    fn default() -> Self {
        let limits = ResourceLimits::new()
            .max_duration(Duration::from_millis(100))
            .max_memory(8 * 1024 * 1024)
            .max_allocations(200_000)
            .gc_interval(1000)
            .max_recursion_depth(Some(200));
        Self {
            limits,
            loader: GenotypeLoadOptions::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RuntimeError {
    Monty(String),
    Unsupported(String),
    InvalidArguments(String),
    Io(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StageTiming {
    pub stage: String,
    pub duration_ms: u128,
    pub detail: String,
}

impl fmt::Display for RuntimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Monty(msg) | Self::Unsupported(msg) | Self::InvalidArguments(msg) | Self::Io(msg) => {
                f.write_str(msg)
            }
        }
    }
}

impl Error for RuntimeError {}

impl From<MontyException> for RuntimeError {
    fn from(value: MontyException) -> Self {
        Self::Monty(value.to_string())
    }
}

struct RuntimeState {
    next_handle: AtomicU64,
    genotype_files: Mutex<HashMap<u64, GenotypeStore>>,
    trace_lines: Mutex<Vec<usize>>,
    timings: Mutex<Vec<StageTiming>>,
}

impl RuntimeState {
    fn new() -> Self {
        Self {
            next_handle: AtomicU64::new(1),
            genotype_files: Mutex::new(HashMap::new()),
            trace_lines: Mutex::new(Vec::new()),
            timings: Mutex::new(Vec::new()),
        }
    }

    fn next_handle(&self) -> u64 {
        self.next_handle.fetch_add(1, Ordering::Relaxed)
    }
}

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

    pub fn with_config(root: impl Into<PathBuf>, config: RuntimeConfig) -> Result<Self, RuntimeError> {
        let root = root.into();
        let canonical_root = root.canonicalize().map_err(|err| {
            RuntimeError::Io(format!("failed to canonicalize bioscript root {}: {err}", root.display()))
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
            RuntimeError::Io(format!("failed to read script {}: {err}", script_path.display()))
        })?;
        let instrumented = instrument_source(&code);
        self.state.trace_lines.lock().expect("trace mutex poisoned").clear();
        self.state.timings.lock().expect("timings mutex poisoned").clear();

        extra_inputs.push(("__name__", MontyObject::String("__main__".to_owned())));
        extra_inputs.push((
            "__file__",
            MontyObject::String(script_path.display().to_string()),
        ));
        extra_inputs.push(("bioscript", bioscript_object()));

        let result = self.run_script(&instrumented, &script_path.display().to_string(), extra_inputs)?;

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
        self.state.timings.lock().expect("timings mutex poisoned").clone()
    }

    pub fn run_script(
        &self,
        code: &str,
        script_name: &str,
        inputs: Vec<(&str, MontyObject)>,
    ) -> Result<MontyObject, RuntimeError> {
        let input_names = inputs.iter().map(|(name, _)| (*name).to_owned()).collect();
        let input_values = inputs.into_iter().map(|(_, value)| value).collect();
        let runner = MontyRun::new(code.to_owned(), script_name, input_names)?;
        let tracker = LimitedTracker::new(self.config.limits.clone());
        let mut progress = runner.start(input_values, tracker, PrintWriter::Stdout)?;

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
                    lookup.resume(result, PrintWriter::Stdout)?
                }
                RunProgress::FunctionCall(call) => {
                    if call.method_call {
                        let result = self.dispatch_method_call(&call.function_name, &call.args, &call.kwargs)?;
                        call.resume(result, PrintWriter::Stdout)?
                    } else {
                        let Some(handler) = self.functions.get(call.function_name.as_str()) else {
                            return Err(RuntimeError::Unsupported(format!(
                                "unknown bioscript host function: {}",
                                call.function_name
                            )));
                        };
                        let result = handler(self, &call.args, &call.kwargs)?;
                        call.resume(result, PrintWriter::Stdout)?
                    }
                }
                RunProgress::ResolveFutures(state) => {
                    return Err(RuntimeError::Unsupported(format!(
                        "async futures are not supported in bioscript runtime: {:?}",
                        state.pending_call_ids()
                    )));
                }
                RunProgress::OsCall(call) => {
                    return Err(RuntimeError::Unsupported(format!("OS call {} is blocked", call.function)));
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
            ("GenotypeFile", "lookup_variants") => self.method_genotype_lookup_variants(args, kwargs),
            _ => Err(RuntimeError::Unsupported(format!(
                "'{class_name}' object has no attribute '{method_name}'"
            ))),
        }
    }

    fn method_load_genotypes(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = Instant::now();
        reject_kwargs(kwargs, "bioscript.load_genotypes")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.load_genotypes expects self and path".to_owned(),
            ));
        }
        let path = self.resolve_user_path(&expect_string_arg(args, 1, "bioscript.load_genotypes")?)?;
        let loader = self.resolved_loader_options()?;
        let store = GenotypeStore::from_file_with_options(&path, &loader)?;
        let handle = self.state.next_handle();
        self.state
            .genotype_files
            .lock()
            .expect("genotype mutex poisoned")
            .insert(handle, store);
        self.record_timing(
            "load_genotypes",
            started.elapsed(),
            format!("path={}", path.display()),
        );
        Ok(genotype_file_object(handle))
    }

    fn resolved_loader_options(&self) -> Result<GenotypeLoadOptions, RuntimeError> {
        let mut loader = self.config.loader.clone();
        loader.input_index = resolve_optional_loader_path(self, loader.input_index)?;
        loader.reference_file = resolve_optional_loader_path(self, loader.reference_file)?;
        loader.reference_index = resolve_optional_loader_path(self, loader.reference_index)?;
        Ok(loader)
    }

    fn method_genotype_get(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "GenotypeFile.get")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "GenotypeFile.get expects self and rsid".to_owned(),
            ));
        }
        let handle = dataclass_handle_id(&args[0], "GenotypeFile")?;
        let rsid = expect_string_arg(args, 1, "GenotypeFile.get")?;
        let guard = self.state.genotype_files.lock().expect("genotype mutex poisoned");
        let Some(store) = guard.get(&handle) else {
            return Err(RuntimeError::InvalidArguments(format!(
                "unknown genotype handle: {handle}"
            )));
        };
        Ok(match store.get(&rsid)? {
            Some(value) => MontyObject::String(value),
            None => MontyObject::None,
        })
    }

    fn method_variant(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        if args.len() != 1 {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.variant expects only self as a positional argument".to_owned(),
            ));
        }
        let spec = variant_spec_from_kwargs(kwargs)?;
        Ok(variant_object(&spec))
    }

    fn method_query_plan(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "bioscript.query_plan")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.query_plan expects self and a list of variants".to_owned(),
            ));
        }
        let variants = variant_specs_from_plan(&args[1])?;
        Ok(variant_plan_object(&variants))
    }

    fn method_genotype_lookup_variant(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = Instant::now();
        reject_kwargs(kwargs, "GenotypeFile.lookup_variant")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "GenotypeFile.lookup_variant expects self and variant".to_owned(),
            ));
        }
        let handle = dataclass_handle_id(&args[0], "GenotypeFile")?;
        let spec = dataclass_to_variant_spec(&args[1])?;
        let guard = self.state.genotype_files.lock().expect("genotype mutex poisoned");
        let Some(store) = guard.get(&handle) else {
            return Err(RuntimeError::InvalidArguments(format!(
                "unknown genotype handle: {handle}"
            )));
        };
        let observation = store.lookup_variant(&spec)?;
        self.record_timing(
            "lookup_variant",
            started.elapsed(),
            format!("rsids={}", spec.rsids.join("|")),
        );
        Ok(match observation.genotype {
            Some(value) => MontyObject::String(value),
            None => MontyObject::None,
        })
    }

    fn method_genotype_lookup_variants(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = Instant::now();
        reject_kwargs(kwargs, "GenotypeFile.lookup_variants")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "GenotypeFile.lookup_variants expects self and a variant plan".to_owned(),
            ));
        }
        let handle = dataclass_handle_id(&args[0], "GenotypeFile")?;
        let specs = variant_specs_from_plan(&args[1])?;
        let guard = self.state.genotype_files.lock().expect("genotype mutex poisoned");
        let Some(store) = guard.get(&handle) else {
            return Err(RuntimeError::InvalidArguments(format!(
                "unknown genotype handle: {handle}"
            )));
        };
        let observations = store.lookup_variants(&specs)?;
        self.record_timing(
            "lookup_variants",
            started.elapsed(),
            format!("count={}", specs.len()),
        );
        Ok(MontyObject::List(
            observations
                .into_iter()
                .map(|observation| match observation.genotype {
                    Some(value) => MontyObject::String(value),
                    None => MontyObject::None,
                })
                .collect(),
        ))
    }

    fn method_write_tsv(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = Instant::now();
        reject_kwargs(kwargs, "bioscript.write_tsv")?;
        if args.len() != 3 {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.write_tsv expects self, path, rows".to_owned(),
            ));
        }
        let path = self.resolve_user_path(&expect_string_arg(args, 1, "bioscript.write_tsv")?)?;
        let rows = expect_rows(&args[2])?;
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                RuntimeError::Io(format!("failed to create parent dir {}: {err}", parent.display()))
            })?;
        }
        let mut output = String::new();
        if let Some(first) = rows.first() {
            let headers: Vec<String> = first.keys().cloned().collect();
            output.push_str(&headers.join("\t"));
            output.push('\n');
            for row in &rows {
                let values: Vec<String> = headers
                    .iter()
                    .map(|header| row.get(header).cloned().unwrap_or_default())
                    .collect();
                output.push_str(&values.join("\t"));
                output.push('\n');
            }
        }
        fs::write(&path, output)
            .map_err(|err| RuntimeError::Io(format!("failed to write {}: {err}", path.display())))?;
        self.record_timing(
            "write_tsv",
            started.elapsed(),
            format!("path={} rows={}", path.display(), rows.len()),
        );
        Ok(MontyObject::None)
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

    fn method_read_text(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        if args.is_empty() {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.read_text expects self and path".to_owned(),
            ));
        }
        host_read_text(self, &args[1..], kwargs)
    }

    fn method_write_text(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        if args.is_empty() {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.write_text expects self, path, text".to_owned(),
            ));
        }
        host_write_text(self, &args[1..], kwargs)
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

    fn write_trace_report(&self, report_path: &Path, original_code: &str) -> Result<(), RuntimeError> {
        let trace_lines = self.state.trace_lines.lock().expect("trace mutex poisoned").clone();
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
                RuntimeError::Io(format!("failed to create parent dir {}: {err}", parent.display()))
            })?;
        }
        fs::write(report_path, output)
            .map_err(|err| RuntimeError::Io(format!("failed to write {}: {err}", report_path.display())))?;
        Ok(())
    }
}

fn trace_lookup_metadata(source: &str) -> (Option<String>, Option<String>) {
    if let Some(rsid) = extract_rsid(source) {
        let url = format!("https://www.ncbi.nlm.nih.gov/snp/{rsid}");
        return (Some(rsid), Some(url));
    }

    if let Some(coord) = extract_coordinate(source) {
        let lower = source.to_ascii_lowercase();
        let host = if lower.contains("grch37") || lower.contains("hg19") {
            "https://grch37.ensembl.org"
        } else {
            "https://www.ensembl.org"
        };
        let url = format!("{host}/Homo_sapiens/Location/View?r={coord}");
        return (Some(coord), Some(url));
    }

    (None, None)
}

fn statement_context(lines: &[&str], line_no: usize) -> String {
    let Some(start_idx) = line_no.checked_sub(1) else {
        return String::new();
    };
    let Some(first_line) = lines.get(start_idx) else {
        return String::new();
    };

    let mut out = String::from(first_line.trim());
    let mut depth = update_nesting_depth(0, first_line);
    let mut current = start_idx + 1;

    while depth > 0 {
        let Some(line) = lines.get(current) else {
            break;
        };
        if !out.is_empty() {
            out.push(' ');
        }
        out.push_str(line.trim());
        depth = update_nesting_depth(depth, line);
        current += 1;
    }

    out
}

fn extract_rsid(source: &str) -> Option<String> {
    let chars: Vec<char> = source.chars().collect();
    let len = chars.len();
    let mut idx = 0;
    while idx + 2 <= len {
        if chars[idx] == 'r'
            && chars.get(idx + 1) == Some(&'s')
            && (idx == 0 || !chars[idx - 1].is_ascii_alphanumeric())
        {
            let mut end = idx + 2;
            while end < len && chars[end].is_ascii_digit() {
                end += 1;
            }
            if end > idx + 2 {
                return Some(chars[idx..end].iter().collect());
            }
        }
        idx += 1;
    }
    None
}

fn extract_coordinate(source: &str) -> Option<String> {
    for token in source
        .split(|ch: char| ch.is_whitespace() || matches!(ch, '"' | '\'' | ',' | ')' | '(' | '[' | ']' | '{' | '}'))
    {
        let cleaned = token.trim_matches(|ch: char| matches!(ch, ';'));
        let normalized = cleaned.strip_prefix("chr").unwrap_or(cleaned);
        if let Some((chrom, rest)) = normalized.split_once(':')
            && !chrom.is_empty()
            && chrom.chars().all(|ch| ch.is_ascii_alphanumeric())
        {
            if let Some((start, end)) = rest.split_once('-') {
                if start.chars().all(|ch| ch.is_ascii_digit()) && end.chars().all(|ch| ch.is_ascii_digit()) {
                    return Some(format!("{chrom}:{start}-{end}"));
                }
            } else if rest.chars().all(|ch| ch.is_ascii_digit()) {
                return Some(format!("{chrom}:{rest}-{rest}"));
            }
        }
    }
    None
}

fn bioscript_object() -> MontyObject {
    MontyObject::Dataclass {
        name: "Bioscript".to_owned(),
        type_id: 1,
        field_names: vec![],
        attrs: vec![].into(),
        frozen: true,
    }
}

fn genotype_file_object(handle_id: u64) -> MontyObject {
    MontyObject::Dataclass {
        name: "GenotypeFile".to_owned(),
        type_id: 2,
        field_names: vec!["handle_id".to_owned()],
        attrs: vec![(
            MontyObject::String("handle_id".to_owned()),
            MontyObject::Int(handle_id as i64),
        )]
        .into(),
        frozen: true,
    }
}

fn variant_object(spec: &VariantSpec) -> MontyObject {
    let mut attrs = Vec::new();
    attrs.push((
        MontyObject::String("rsids".to_owned()),
        MontyObject::List(spec.rsids.iter().cloned().map(MontyObject::String).collect()),
    ));
    if let Some(locus) = &spec.grch37 {
        attrs.push((
            MontyObject::String("grch37".to_owned()),
            MontyObject::String(format!("{}:{}-{}", locus.chrom, locus.start, locus.end)),
        ));
    }
    if let Some(locus) = &spec.grch38 {
        attrs.push((
            MontyObject::String("grch38".to_owned()),
            MontyObject::String(format!("{}:{}-{}", locus.chrom, locus.start, locus.end)),
        ));
    }
    if let Some(reference) = &spec.reference {
        attrs.push((
            MontyObject::String("reference".to_owned()),
            MontyObject::String(reference.clone()),
        ));
    }
    if let Some(alternate) = &spec.alternate {
        attrs.push((
            MontyObject::String("alternate".to_owned()),
            MontyObject::String(alternate.clone()),
        ));
    }
    if let Some(kind) = spec.kind {
        attrs.push((
            MontyObject::String("kind".to_owned()),
            MontyObject::String(variant_kind_name(kind).to_owned()),
        ));
    }
    if let Some(length) = spec.deletion_length {
        attrs.push((
            MontyObject::String("deletion_length".to_owned()),
            MontyObject::Int(length as i64),
        ));
    }
    if !spec.motifs.is_empty() {
        attrs.push((
            MontyObject::String("motifs".to_owned()),
            MontyObject::List(spec.motifs.iter().cloned().map(MontyObject::String).collect()),
        ));
    }

    MontyObject::Dataclass {
        name: "Variant".to_owned(),
        type_id: 3,
        field_names: vec![
            "rsids".to_owned(),
            "grch37".to_owned(),
            "grch38".to_owned(),
            "reference".to_owned(),
            "alternate".to_owned(),
            "kind".to_owned(),
            "deletion_length".to_owned(),
            "motifs".to_owned(),
        ],
        attrs: attrs.into(),
        frozen: true,
    }
}

fn variant_plan_object(variants: &[VariantSpec]) -> MontyObject {
    MontyObject::Dataclass {
        name: "VariantPlan".to_owned(),
        type_id: 4,
        field_names: vec!["variants".to_owned()],
        attrs: vec![(
            MontyObject::String("variants".to_owned()),
            MontyObject::List(variants.iter().map(variant_object).collect()),
        )]
        .into(),
        frozen: true,
    }
}

fn dataclass_handle_id(obj: &MontyObject, expected_name: &str) -> Result<u64, RuntimeError> {
    match obj {
        MontyObject::Dataclass { name, attrs, .. } if name == expected_name => {
            for (key, value) in attrs {
                if matches!(key, MontyObject::String(text) if text == "handle_id")
                    && let MontyObject::Int(id) = value
                {
                    return Ok(*id as u64);
                }
            }
            Err(RuntimeError::InvalidArguments(format!(
                "{expected_name} missing handle_id"
            )))
        }
        _ => Err(RuntimeError::InvalidArguments(format!(
            "expected {expected_name} object"
        ))),
    }
}

fn dataclass_to_variant_spec(obj: &MontyObject) -> Result<VariantSpec, RuntimeError> {
    let MontyObject::Dataclass { name, attrs, .. } = obj else {
        return Err(RuntimeError::InvalidArguments("expected Variant object".to_owned()));
    };
    if name != "Variant" {
        return Err(RuntimeError::InvalidArguments(format!("expected Variant object, got {name}")));
    }

    let mut spec = VariantSpec::default();
    for (key, value) in attrs {
        let MontyObject::String(key) = key else {
            continue;
        };
        match key.as_str() {
            "rsids" => spec.rsids = string_list_from_object(value)?,
            "grch37" => spec.grch37 = string_from_optional(value)?.map(|v| parse_locus_string(&v)).transpose()?,
            "grch38" => spec.grch38 = string_from_optional(value)?.map(|v| parse_locus_string(&v)).transpose()?,
            "reference" => spec.reference = string_from_optional(value)?,
            "alternate" => spec.alternate = string_from_optional(value)?,
            "kind" => spec.kind = string_from_optional(value)?.as_deref().map(parse_variant_kind).transpose()?,
            "deletion_length" => spec.deletion_length = int_from_optional(value)?.map(|v| v as usize),
            "motifs" => spec.motifs = string_list_from_object(value)?,
            _ => {}
        }
    }
    Ok(spec)
}

fn variant_specs_from_plan(obj: &MontyObject) -> Result<Vec<VariantSpec>, RuntimeError> {
    match obj {
        MontyObject::List(items) => items.iter().map(dataclass_to_variant_spec).collect(),
        MontyObject::Dataclass { name, attrs, .. } if name == "VariantPlan" => {
            for (key, value) in attrs {
                if matches!(key, MontyObject::String(text) if text == "variants") {
                    return variant_specs_from_plan(value);
                }
            }
            Err(RuntimeError::InvalidArguments(
                "VariantPlan missing variants".to_owned(),
            ))
        }
        _ => Err(RuntimeError::InvalidArguments(
            "expected a list of Variant objects or a VariantPlan".to_owned(),
        )),
    }
}

fn variant_spec_from_kwargs(kwargs: &[(MontyObject, MontyObject)]) -> Result<VariantSpec, RuntimeError> {
    let mut spec = VariantSpec::default();
    for (key, value) in kwargs {
        let MontyObject::String(key) = key else {
            return Err(RuntimeError::InvalidArguments(
                "bioscript.variant keyword names must be strings".to_owned(),
            ));
        };
        match key.as_str() {
            "rsid" | "rsids" => spec.rsids = string_or_list(value)?,
            "grch37" => spec.grch37 = string_from_optional(value)?.map(|v| parse_locus_string(&v)).transpose()?,
            "grch38" => spec.grch38 = string_from_optional(value)?.map(|v| parse_locus_string(&v)).transpose()?,
            "ref" | "reference" => spec.reference = string_from_optional(value)?,
            "alt" | "alternate" => spec.alternate = string_from_optional(value)?,
            "kind" => spec.kind = string_from_optional(value)?.as_deref().map(parse_variant_kind).transpose()?,
            "deletion_length" => spec.deletion_length = int_from_optional(value)?.map(|v| v as usize),
            "motifs" => spec.motifs = string_or_list(value)?,
            other => {
                return Err(RuntimeError::InvalidArguments(format!(
                    "bioscript.variant does not accept keyword '{other}'"
                )))
            }
        }
    }
    Ok(spec)
}

fn parse_locus_string(value: &str) -> Result<GenomicLocus, RuntimeError> {
    let normalized = value.trim().strip_prefix("chr").unwrap_or(value.trim());
    let Some((chrom, rest)) = normalized.split_once(':') else {
        return Err(RuntimeError::InvalidArguments(format!("invalid locus string: {value}")));
    };
    let (start, end) = if let Some((start, end)) = rest.split_once('-') {
        (start, end)
    } else {
        (rest, rest)
    };
    let start = start
        .parse::<i64>()
        .map_err(|err| RuntimeError::InvalidArguments(format!("invalid locus start {value}: {err}")))?;
    let end = end
        .parse::<i64>()
        .map_err(|err| RuntimeError::InvalidArguments(format!("invalid locus end {value}: {err}")))?;
    Ok(GenomicLocus {
        chrom: chrom.to_owned(),
        start,
        end,
    })
}

fn parse_variant_kind(value: &str) -> Result<VariantKind, RuntimeError> {
    match value.trim().to_ascii_lowercase().as_str() {
        "snp" => Ok(VariantKind::Snp),
        "insertion" | "ins" => Ok(VariantKind::Insertion),
        "deletion" | "del" => Ok(VariantKind::Deletion),
        "indel" => Ok(VariantKind::Indel),
        "other" => Ok(VariantKind::Other),
        other => Err(RuntimeError::InvalidArguments(format!("invalid variant kind: {other}"))),
    }
}

fn variant_kind_name(kind: VariantKind) -> &'static str {
    match kind {
        VariantKind::Snp => "snp",
        VariantKind::Insertion => "insertion",
        VariantKind::Deletion => "deletion",
        VariantKind::Indel => "indel",
        VariantKind::Other => "other",
    }
}

fn string_or_list(value: &MontyObject) -> Result<Vec<String>, RuntimeError> {
    match value {
        MontyObject::String(text) => Ok(vec![text.clone()]),
        MontyObject::List(_) => string_list_from_object(value),
        MontyObject::None => Ok(Vec::new()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected string or list of strings, got {other:?}"
        ))),
    }
}

fn string_list_from_object(value: &MontyObject) -> Result<Vec<String>, RuntimeError> {
    match value {
        MontyObject::List(items) => {
            let mut out = Vec::new();
            for item in items {
                let MontyObject::String(text) = item else {
                    return Err(RuntimeError::InvalidArguments(
                        "expected list of strings".to_owned(),
                    ));
                };
                out.push(text.clone());
            }
            Ok(out)
        }
        MontyObject::None => Ok(Vec::new()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected list of strings, got {other:?}"
        ))),
    }
}

fn string_from_optional(value: &MontyObject) -> Result<Option<String>, RuntimeError> {
    match value {
        MontyObject::None => Ok(None),
        MontyObject::String(text) => Ok(Some(text.clone())),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected optional string, got {other:?}"
        ))),
    }
}

fn int_from_optional(value: &MontyObject) -> Result<Option<i64>, RuntimeError> {
    match value {
        MontyObject::None => Ok(None),
        MontyObject::Int(v) => Ok(Some(*v)),
        other => Err(RuntimeError::InvalidArguments(format!(
            "expected optional int, got {other:?}"
        ))),
    }
}

fn reject_kwargs(kwargs: &[(MontyObject, MontyObject)], function_name: &str) -> Result<(), RuntimeError> {
    if kwargs.is_empty() {
        Ok(())
    } else {
        Err(RuntimeError::InvalidArguments(format!(
            "{function_name} does not accept keyword arguments"
        )))
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

fn expect_string_arg(args: &[MontyObject], index: usize, function_name: &str) -> Result<String, RuntimeError> {
    let Some(value) = args.get(index) else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} missing argument at position {index}"
        )));
    };
    match value {
        MontyObject::String(text) => Ok(text.clone()),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{function_name} expected str at position {index}, got {other:?}"
        ))),
    }
}

fn expect_rows(value: &MontyObject) -> Result<Vec<BTreeMap<String, String>>, RuntimeError> {
    let MontyObject::List(rows) = value else {
        return Err(RuntimeError::InvalidArguments(
            "write_tsv expects a list of dict rows".to_owned(),
        ));
    };

    let mut out = Vec::new();
    for row in rows {
        let MontyObject::Dict(dict) = row else {
            return Err(RuntimeError::InvalidArguments(
                "write_tsv row must be a dict".to_owned(),
            ));
        };
        let mut mapped = BTreeMap::new();
        for (key, value) in dict {
            let MontyObject::String(key) = key else {
                return Err(RuntimeError::InvalidArguments(
                    "write_tsv dict keys must be strings".to_owned(),
                ));
            };
            mapped.insert(key.clone(), stringify_value(value));
        }
        out.push(mapped);
    }
    Ok(out)
}

fn stringify_value(value: &MontyObject) -> String {
    match value {
        MontyObject::None => String::new(),
        MontyObject::String(text) => text.clone(),
        MontyObject::Int(v) => v.to_string(),
        MontyObject::Bool(v) => v.to_string(),
        other => format!("{other}"),
    }
}

fn host_read_text(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "read_text")?;
    let path = runtime.resolve_user_path(&expect_string_arg(args, 0, "read_text")?)?;
    let content = fs::read_to_string(&path)
        .map_err(|err| RuntimeError::Io(format!("failed to read {}: {err}", path.display())))?;
    Ok(MontyObject::String(content))
}

fn host_write_text(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "write_text")?;
    let path = runtime.resolve_user_path(&expect_string_arg(args, 0, "write_text")?)?;
    let content = expect_string_arg(args, 1, "write_text")?;
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|err| {
            RuntimeError::Io(format!("failed to create parent dir {}: {err}", parent.display()))
        })?;
    }
    fs::write(&path, content)
        .map_err(|err| RuntimeError::Io(format!("failed to write {}: {err}", path.display())))?;
    Ok(MontyObject::None)
}

fn host_trace(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    kwargs: &[(MontyObject, MontyObject)],
) -> Result<MontyObject, RuntimeError> {
    reject_kwargs(kwargs, "__bioscript_trace__")?;
    if let Some(MontyObject::Int(v)) = args.first() {
        runtime
            .state
            .trace_lines
            .lock()
            .expect("trace mutex poisoned")
            .push(*v as usize);
    }
    Ok(MontyObject::None)
}

fn instrument_source(code: &str) -> String {
    let mut out = Vec::new();
    let mut nesting_depth = 0usize;
    let mut pending_backslash = false;
    for (idx, line) in code.lines().enumerate() {
        let line_no = idx + 1;
        let trimmed = line.trim_start();

        let in_continuation = nesting_depth > 0 || pending_backslash;
        let should_trace = !in_continuation
            && !trimmed.is_empty()
            && !trimmed.starts_with('#')
            && !trimmed.starts_with('@')
            && !trimmed.starts_with('"')
            && !trimmed.starts_with('\'')
            && !trimmed.starts_with(']')
            && !trimmed.starts_with(')')
            && !trimmed.starts_with('}')
            && !trimmed.starts_with(',')
            && !trimmed.starts_with('+')
            && !trimmed.starts_with('-')
            && !trimmed.starts_with('*')
            && !trimmed.starts_with('/')
            && !trimmed.starts_with('%')
            && !trimmed.starts_with("and ")
            && !trimmed.starts_with("or ")
            && !trimmed.starts_with("if ")
            && !trimmed.starts_with("for ")
            && !trimmed.starts_with("elif ")
            && !trimmed.starts_with("else:")
            && !trimmed.starts_with("except")
            && !trimmed.starts_with("finally:")
            && !trimmed.ends_with(':');

        if should_trace {
            let indent_len = line.len() - trimmed.len();
            let indent = &line[..indent_len];
            out.push(format!("{indent}__bioscript_trace__({line_no})"));
        }
        out.push(line.to_owned());

        pending_backslash = ends_with_unescaped_backslash(line);
        nesting_depth = update_nesting_depth(nesting_depth, line);
    }
    if code.ends_with('\n') {
        out.join("\n") + "\n"
    } else {
        out.join("\n")
    }
}

fn ends_with_unescaped_backslash(line: &str) -> bool {
    let trimmed = line.trim_end();
    if !trimmed.ends_with('\\') {
        return false;
    }

    let slash_count = trimmed.chars().rev().take_while(|ch| *ch == '\\').count();
    slash_count % 2 == 1
}

fn update_nesting_depth(mut depth: usize, line: &str) -> usize {
    let mut chars = line.chars().peekable();
    let mut in_single = false;
    let mut in_double = false;

    while let Some(ch) = chars.next() {
        if in_single {
            if ch == '\\' {
                chars.next();
            } else if ch == '\'' {
                in_single = false;
            }
            continue;
        }

        if in_double {
            if ch == '\\' {
                chars.next();
            } else if ch == '"' {
                in_double = false;
            }
            continue;
        }

        match ch {
            '#' => break,
            '\'' => in_single = true,
            '"' => in_double = true,
            '(' | '[' | '{' => depth += 1,
            ')' | ']' | '}' => depth = depth.saturating_sub(1),
            _ => {}
        }
    }

    depth
}
