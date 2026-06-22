use bioscript_core::RuntimeError;
use bioscript_formats::{GenotypeLoadOptions, GenotypeSourceFormat, GenotypeStore};
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{expect_string_arg, reject_kwargs},
    objects::genotype_file_object,
    resolve_optional_loader_path,
    timing::RuntimeInstant,
};

impl BioscriptRuntime {
    pub(super) fn method_load_genotypes(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        let started = RuntimeInstant::now();
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
        let inner_store = if let Some(bytes) = self.read_virtual_binary_file(&path) {
            if bytes.is_empty()
                && matches!(
                    loader.format,
                    Some(GenotypeSourceFormat::Cram | GenotypeSourceFormat::Bam)
                )
                && !self.config.preloaded_observations.is_empty()
            {
                GenotypeStore::empty()
            } else {
                // The report pipeline virtualizes the genotype input as in-memory
                // bytes. CRAM/BAM cannot be decoded by the text/zip/vcf
                // `from_bytes` path — route them to the byte-backed alignment
                // backend, pulling the reference + index from the loader's
                // virtualized companion files. (Regression: pre-rewrite these
                // were real file paths handled by `from_file_with_options`.)
                match alignment_bytes_kind(&bytes).or(match loader.format {
                    Some(GenotypeSourceFormat::Cram | GenotypeSourceFormat::Bam) => loader.format,
                    _ => None,
                }) {
                    Some(kind) => {
                        let index = self.virtual_alignment_aux(loader.input_index.as_ref())?;
                        let (reference, reference_index) =
                            if matches!(kind, GenotypeSourceFormat::Cram) {
                                (
                                    self.virtual_alignment_aux(loader.reference_file.as_ref())?,
                                    self.virtual_alignment_aux(loader.reference_index.as_ref())?,
                                )
                            } else {
                                (Vec::new(), Vec::new())
                            };
                        GenotypeStore::from_alignment_bytes(
                            kind,
                            bytes,
                            index,
                            reference,
                            reference_index,
                            &loader,
                        )
                    }
                    None => GenotypeStore::from_bytes_with_options(
                        path.file_name()
                            .and_then(|value| value.to_str())
                            .unwrap_or("input"),
                        &bytes,
                        &loader,
                    )?,
                }
            }
        } else {
            GenotypeStore::from_file_with_options(&path, &loader)?
        };
        // Layer pre-resolved observations on top of whatever backend the path
        // resolves to. The report pipeline collects every variant the panel
        // declares before running analyses, so `genotypes.lookup_variants(...)`
        // must hit this cache. A miss is an error: falling through to a second
        // lookup path can make analysis disagree with the observations table.
        let store = if self.config.preloaded_observations.is_empty() {
            inner_store
        } else {
            GenotypeStore::with_required_cached_observations(
                self.config.preloaded_observations.clone(),
                inner_store,
            )
        };
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

    /// Read a virtualized alignment companion (reference / `.fai` / `.crai` /
    /// `.bai`). Prefers the virtual binary file; falls back to the real
    /// filesystem so the CLI keeps working when the report did not virtualize
    /// the companion.
    fn virtual_alignment_aux(
        &self,
        path: Option<&std::path::PathBuf>,
    ) -> Result<Vec<u8>, RuntimeError> {
        let path = path.ok_or_else(|| {
            RuntimeError::InvalidArguments(
                "alignment input requires --reference-file/--input-index".to_owned(),
            )
        })?;
        if let Some(bytes) = self.read_virtual_binary_file(path) {
            return Ok(bytes);
        }
        std::fs::read(path)
            .map_err(|err| RuntimeError::Io(format!("failed to read {}: {err}", path.display())))
    }

    pub(super) fn resolved_loader_options(&self) -> Result<GenotypeLoadOptions, RuntimeError> {
        let mut loader = self.config.loader.clone();
        loader.input_index = resolve_optional_loader_path(self, loader.input_index)?;
        loader.reference_file = resolve_optional_loader_path(self, loader.reference_file)?;
        loader.reference_index = resolve_optional_loader_path(self, loader.reference_index)?;
        Ok(loader)
    }
}

/// Sniff CRAM/BAM by magic bytes. CRAM files start with `CRAM`; BAM files
/// start with `BAM\x01` (after bgzf this is the decompressed magic, but
/// htslib BAMs are bgzf whose first block decompresses to it — the report
/// hands us the raw `.bam`, whose bgzf wrapper starts `\x1f\x8b`; noodles'
/// indexed reader handles the bgzf, so detect BAM via the gzip magic plus
/// the `.bam`-style content, while CRAM is the literal `CRAM` magic).
fn alignment_bytes_kind(bytes: &[u8]) -> Option<GenotypeSourceFormat> {
    if bytes.starts_with(b"CRAM") {
        return Some(GenotypeSourceFormat::Cram);
    }
    if bytes.starts_with(b"BAM\x01") {
        return Some(GenotypeSourceFormat::Bam);
    }
    None
}
