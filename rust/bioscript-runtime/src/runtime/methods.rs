use std::{fs, time::Instant};

use bioscript_core::RuntimeError;
use bioscript_formats::{GenotypeLoadOptions, GenotypeStore};
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{expect_rows, expect_string_arg, reject_kwargs},
    host_io::{host_read_text, host_write_text},
    objects::{
        genotype_file_object, variant_object, variant_observation_object, variant_plan_object,
    },
    resolve_optional_loader_path,
    variants::{
        dataclass_handle_id, dataclass_to_variant_spec, variant_spec_from_kwargs,
        variant_specs_from_plan,
    },
};

impl BioscriptRuntime {
    pub(super) fn method_load_genotypes(
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
        let path = self.resolve_existing_user_path(&expect_string_arg(
            args,
            1,
            "bioscript.load_genotypes",
        )?)?;
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

    pub(super) fn resolved_loader_options(&self) -> Result<GenotypeLoadOptions, RuntimeError> {
        let mut loader = self.config.loader.clone();
        loader.input_index = resolve_optional_loader_path(self, loader.input_index)?;
        loader.reference_file = resolve_optional_loader_path(self, loader.reference_file)?;
        loader.reference_index = resolve_optional_loader_path(self, loader.reference_index)?;
        Ok(loader)
    }

    pub(super) fn method_genotype_get(
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
        let guard = self
            .state
            .genotype_files
            .lock()
            .expect("genotype mutex poisoned");
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

    pub(super) fn method_variant(
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

    pub(super) fn method_query_plan(
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

    pub(super) fn method_genotype_lookup_variant(
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
        let guard = self
            .state
            .genotype_files
            .lock()
            .expect("genotype mutex poisoned");
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

    pub(super) fn method_genotype_lookup_variant_details(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = Instant::now();
        reject_kwargs(kwargs, "GenotypeFile.lookup_variant_details")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "GenotypeFile.lookup_variant_details expects self and variant".to_owned(),
            ));
        }
        let handle = dataclass_handle_id(&args[0], "GenotypeFile")?;
        let spec = dataclass_to_variant_spec(&args[1])?;
        let guard = self
            .state
            .genotype_files
            .lock()
            .expect("genotype mutex poisoned");
        let Some(store) = guard.get(&handle) else {
            return Err(RuntimeError::InvalidArguments(format!(
                "unknown genotype handle: {handle}"
            )));
        };
        let observation = store.lookup_variant(&spec)?;
        self.record_timing(
            "lookup_variant_details",
            started.elapsed(),
            format!("rsids={}", spec.rsids.join("|")),
        );
        Ok(variant_observation_object(&observation))
    }

    pub(super) fn method_genotype_lookup_variants(
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
        let guard = self
            .state
            .genotype_files
            .lock()
            .expect("genotype mutex poisoned");
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

    pub(super) fn method_genotype_lookup_variants_details(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = Instant::now();
        reject_kwargs(kwargs, "GenotypeFile.lookup_variants_details")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "GenotypeFile.lookup_variants_details expects self and a variant plan".to_owned(),
            ));
        }
        let handle = dataclass_handle_id(&args[0], "GenotypeFile")?;
        let specs = variant_specs_from_plan(&args[1])?;
        let guard = self
            .state
            .genotype_files
            .lock()
            .expect("genotype mutex poisoned");
        let Some(store) = guard.get(&handle) else {
            return Err(RuntimeError::InvalidArguments(format!(
                "unknown genotype handle: {handle}"
            )));
        };
        let observations = store.lookup_variants(&specs)?;
        self.record_timing(
            "lookup_variants_details",
            started.elapsed(),
            format!("count={}", specs.len()),
        );
        Ok(MontyObject::List(
            observations
                .iter()
                .map(variant_observation_object)
                .collect(),
        ))
    }

    pub(super) fn method_write_tsv(
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
        let path =
            self.resolve_user_write_path(&expect_string_arg(args, 1, "bioscript.write_tsv")?)?;
        let rows = expect_rows(&args[2])?;
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|err| {
                RuntimeError::Io(format!(
                    "failed to create parent dir {}: {err}",
                    parent.display()
                ))
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
        fs::write(&path, output).map_err(|err| {
            RuntimeError::Io(format!("failed to write {}: {err}", path.display()))
        })?;
        self.record_timing(
            "write_tsv",
            started.elapsed(),
            format!("path={} rows={}", path.display(), rows.len()),
        );
        Ok(MontyObject::None)
    }

    pub(super) fn method_read_text(
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

    pub(super) fn method_write_text(
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
}
