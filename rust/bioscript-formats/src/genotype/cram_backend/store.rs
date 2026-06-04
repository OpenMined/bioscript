use std::{env, fmt::Write as _, path::Path, thread};

use bioscript_core::{RuntimeError, VariantKind, VariantObservation, VariantSpec};

use crate::alignment;

use super::choose_variant_locus;
use crate::genotype::{describe_query, types::CramBackend};

impl CramBackend {
    pub(crate) fn backend_name(&self) -> &'static str {
        "cram"
    }

    /// A CRAM without an external reference can't be pileup-genotyped here
    /// (reference-compressed reads need `--reference-file`; a `no_ref` / embedded
    /// CRAM stores bases but this backend's variant query path still needs the
    /// reference allele context). Rather than abort the whole report, report
    /// the variant as missing so the run degrades to a partial result. An
    /// advanced assay whose analysis consumes the raw aligned reads (e.g.
    /// `VNtyper` running Kestrel over the MUC1 slice) still works.
    fn reference_missing_observation(&self) -> VariantObservation {
        VariantObservation {
            backend: self.backend_name().to_owned(),
            evidence: vec![format!(
                "CRAM variant query skipped for {}: no --reference-file; \
                 reported as missing (analysis consumes raw reads directly)",
                self.path.display()
            )],
            ..VariantObservation::default()
        }
    }

    pub(crate) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let Some(reference_file) = self.options.reference_file.as_ref() else {
            return Ok(self.reference_missing_observation());
        };

        let Some((assembly, locus)) = choose_variant_locus(variant, reference_file) else {
            let mut detail = format!(
                "backend '{}' cannot satisfy query '{}' for {} using reference {}",
                self.backend_name(),
                describe_query(variant),
                self.path.display(),
                reference_file.display()
            );
            detail.push_str(". This backend needs GRCh37/GRCh38 coordinates, not only rsIDs");
            if let Some(reference_index) = self.options.reference_index.as_ref() {
                let _ = write!(detail, " (reference index {})", reference_index.display());
            }
            if let Some(input_index) = self.options.input_index.as_ref() {
                let _ = write!(detail, " (input index {})", input_index.display());
            }
            return Err(RuntimeError::Unsupported(detail));
        };

        let observation = match variant.kind.unwrap_or(VariantKind::Other) {
            VariantKind::Snp => self.observe_snp(variant, assembly, &locus, reference_file)?,
            VariantKind::Deletion => {
                self.observe_deletion(variant, assembly, &locus, reference_file)?
            }
            VariantKind::Insertion | VariantKind::Indel => {
                self.observe_indel(variant, assembly, &locus, reference_file)?
            }
            VariantKind::Other => {
                return Err(RuntimeError::Unsupported(format!(
                    "backend '{}' does not yet support {:?} observation for {}",
                    self.backend_name(),
                    variant.kind.unwrap_or(VariantKind::Other),
                    self.path.display()
                )));
            }
        };

        Ok(observation)
    }

    pub(crate) fn lookup_variants(
        &self,
        variants: &[VariantSpec],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let Some(reference_file) = self.options.reference_file.as_ref() else {
            return Ok(variants
                .iter()
                .map(|_| self.reference_missing_observation())
                .collect());
        };

        let mut indexed: Vec<(usize, &VariantSpec)> = variants.iter().enumerate().collect();
        indexed.sort_by_cached_key(|(_, variant)| crate::genotype::variant_sort_key(variant));

        let worker_count = cram_lookup_worker_count(indexed.len());
        if worker_count <= 1 {
            return self.lookup_variants_serial(reference_file, &indexed);
        }

        let jobs = indexed
            .into_iter()
            .map(|(idx, variant)| (idx, variant.clone()))
            .collect::<Vec<_>>();
        self.lookup_variants_parallel(reference_file, &jobs, worker_count)
    }

    fn lookup_variants_serial(
        &self,
        reference_file: &Path,
        indexed: &[(usize, &VariantSpec)],
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let repository = alignment::build_reference_repository(reference_file)?;
        let mut reader =
            alignment::build_cram_indexed_reader_from_path(&self.path, &self.options, repository)?;
        let label = self.path.display().to_string();

        let result_len = indexed
            .iter()
            .map(|(idx, _)| *idx)
            .max()
            .map_or(0, |idx| idx + 1);
        let mut results = vec![VariantObservation::default(); result_len];
        for (idx, variant) in indexed {
            let Some((assembly, locus)) = choose_variant_locus(variant, reference_file) else {
                results[*idx] = self.unsupported_locus_observation(variant, reference_file);
                continue;
            };
            results[*idx] =
                self.observe_with_reader(&mut reader, &label, variant, assembly, &locus)?;
        }

        Ok(results)
    }

    fn lookup_variants_parallel(
        &self,
        reference_file: &Path,
        jobs: &[(usize, VariantSpec)],
        worker_count: usize,
    ) -> Result<Vec<VariantObservation>, RuntimeError> {
        let result_len = jobs
            .iter()
            .map(|(idx, _)| *idx)
            .max()
            .map_or(0, |idx| idx + 1);
        let chunk_size = jobs.len().div_ceil(worker_count);
        let path = self.path.clone();
        let options = self.options.clone();
        let reference_file = reference_file.to_path_buf();

        let mut handles = Vec::new();
        for chunk in jobs.chunks(chunk_size) {
            let chunk = chunk.to_vec();
            let worker = CramBackend {
                path: path.clone(),
                options: options.clone(),
            };
            let reference_file = reference_file.clone();
            handles.push(thread::spawn(
                move || -> Result<Vec<(usize, VariantObservation)>, RuntimeError> {
                    let repository = alignment::build_reference_repository(&reference_file)?;
                    let mut reader = alignment::build_cram_indexed_reader_from_path(
                        &worker.path,
                        &worker.options,
                        repository,
                    )?;
                    let label = worker.path.display().to_string();
                    let mut observations = Vec::with_capacity(chunk.len());

                    for (idx, variant) in chunk {
                        let observation = if let Some((assembly, locus)) =
                            choose_variant_locus(&variant, &reference_file)
                        {
                            worker.observe_with_reader(
                                &mut reader,
                                &label,
                                &variant,
                                assembly,
                                &locus,
                            )?
                        } else {
                            worker.unsupported_locus_observation(&variant, &reference_file)
                        };
                        observations.push((idx, observation));
                    }

                    Ok(observations)
                },
            ));
        }

        let mut results = vec![VariantObservation::default(); result_len];
        for handle in handles {
            let observations = handle
                .join()
                .map_err(|_| RuntimeError::Io("CRAM lookup worker panicked".to_owned()))??;
            for (idx, observation) in observations {
                results[idx] = observation;
            }
        }

        Ok(results)
    }

    fn unsupported_locus_observation(
        &self,
        variant: &VariantSpec,
        reference_file: &Path,
    ) -> VariantObservation {
        let mut evidence = format!(
            "backend '{}' cannot satisfy query '{}' for {} using reference {}",
            self.backend_name(),
            describe_query(variant),
            self.path.display(),
            reference_file.display()
        );
        evidence.push_str(". This backend needs GRCh37/GRCh38 coordinates, not only rsIDs");
        VariantObservation {
            backend: self.backend_name().to_owned(),
            matched_rsid: variant.rsids.first().cloned(),
            evidence: vec![evidence],
            ..VariantObservation::default()
        }
    }
}

const DEFAULT_MAX_CRAM_WORKERS: usize = 2;

fn cram_lookup_worker_count(job_count: usize) -> usize {
    if job_count <= 1 {
        return 1;
    }

    let requested = env::var("BIOSCRIPT_CRAM_THREADS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|value| *value > 0);
    let available = thread::available_parallelism().map_or(1, usize::from);

    requested
        .unwrap_or_else(|| available.min(DEFAULT_MAX_CRAM_WORKERS))
        .min(job_count)
}

#[cfg(test)]
mod tests {
    use std::path::{Path, PathBuf};

    use bioscript_core::{GenomicLocus, VariantKind};

    use super::*;
    use crate::genotype::GenotypeLoadOptions;

    fn backend(reference_file: Option<PathBuf>) -> CramBackend {
        CramBackend {
            path: PathBuf::from("sample.cram"),
            options: GenotypeLoadOptions {
                reference_file,
                ..GenotypeLoadOptions::default()
            },
        }
    }

    #[test]
    fn cram_store_reports_missing_reference_and_unsupported_locus_details() {
        let variant = VariantSpec {
            rsids: vec!["rs-test".to_owned()],
            kind: Some(VariantKind::Snp),
            reference: Some("A".to_owned()),
            alternate: Some("C".to_owned()),
            ..VariantSpec::default()
        };

        // A CRAM without an external reference is now best-effort: instead
        // of erroring it reports the variant as missing so an advanced
        // assay whose analysis consumes the raw aligned reads still runs.
        let observation = backend(None).lookup_variant(&variant).unwrap();
        assert_eq!(observation.backend, "cram");
        assert!(observation.genotype.is_none());
        assert!(
            observation
                .evidence
                .iter()
                .any(|line| line.contains("no --reference-file"))
        );
        assert!(
            observation
                .evidence
                .iter()
                .any(|line| line.contains("sample.cram"))
        );

        let observation = backend(Some(PathBuf::from("ref.fa")))
            .unsupported_locus_observation(&variant, Path::new("ref.fa"));
        assert_eq!(observation.backend, "cram");
        assert_eq!(observation.matched_rsid.as_deref(), Some("rs-test"));
        assert_eq!(observation.evidence.len(), 1);
        assert!(observation.evidence[0].contains("needs GRCh37/GRCh38 coordinates"));
    }

    #[test]
    fn cram_store_rejects_variants_without_coordinates_before_opening_cram() {
        let store = CramBackend {
            path: PathBuf::from("sample.cram"),
            options: GenotypeLoadOptions {
                reference_file: Some(PathBuf::from("ref.fa")),
                reference_index: Some(PathBuf::from("ref.fa.fai")),
                input_index: Some(PathBuf::from("sample.cram.crai")),
                ..GenotypeLoadOptions::default()
            },
        };

        let err = store
            .lookup_variant(&VariantSpec {
                rsids: vec!["rs-coordinate-only".to_owned()],
                kind: Some(VariantKind::Snp),
                reference: Some("A".to_owned()),
                alternate: Some("C".to_owned()),
                ..VariantSpec::default()
            })
            .unwrap_err();

        let message = err.to_string();
        assert!(message.contains("not only rsIDs"));
        assert!(message.contains("reference index ref.fa.fai"));
        assert!(message.contains("input index sample.cram.crai"));
    }

    #[test]
    fn cram_lookup_worker_count_bounds_requested_and_available_workers() {
        assert_eq!(cram_lookup_worker_count(0), 1);
        assert_eq!(cram_lookup_worker_count(1), 1);

        unsafe {
            env::set_var("BIOSCRIPT_CRAM_THREADS", "2");
        }
        assert_eq!(cram_lookup_worker_count(8), 2);
        assert_eq!(cram_lookup_worker_count(1), 1);

        unsafe {
            env::set_var("BIOSCRIPT_CRAM_THREADS", "999");
        }
        assert_eq!(cram_lookup_worker_count(3), 3);

        unsafe {
            env::set_var("BIOSCRIPT_CRAM_THREADS", "0");
        }
        assert!(cram_lookup_worker_count(2) >= 1);

        unsafe {
            env::remove_var("BIOSCRIPT_CRAM_THREADS");
        }
    }

    #[test]
    fn cram_lookup_variant_rejects_unsupported_kind_after_locus_resolution() {
        let store = backend(Some(PathBuf::from("GRCh38.fa")));
        let err = store
            .lookup_variant(&VariantSpec {
                grch38: Some(GenomicLocus {
                    chrom: "chr1".to_owned(),
                    start: 10,
                    end: 10,
                }),
                kind: Some(VariantKind::Other),
                ..VariantSpec::default()
            })
            .unwrap_err();
        assert!(err.to_string().contains("does not yet support"));
    }
}
