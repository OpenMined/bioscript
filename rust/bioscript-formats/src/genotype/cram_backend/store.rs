use std::{env, fmt::Write as _, path::Path, thread};

use bioscript_core::{RuntimeError, VariantKind, VariantObservation, VariantSpec};

use crate::alignment;

use super::choose_variant_locus;
use crate::genotype::{describe_query, types::CramBackend};

impl CramBackend {
    pub(crate) fn backend_name(&self) -> &'static str {
        "cram"
    }

    pub(crate) fn lookup_variant(
        &self,
        variant: &VariantSpec,
    ) -> Result<VariantObservation, RuntimeError> {
        let Some(reference_file) = self.options.reference_file.as_ref() else {
            return Err(RuntimeError::Unsupported(format!(
                "backend '{}' cannot satisfy query '{}' for {} without --reference-file",
                self.backend_name(),
                describe_query(variant),
                self.path.display()
            )));
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
            return Err(RuntimeError::Unsupported(format!(
                "backend '{}' cannot satisfy CRAM variant queries for {} without --reference-file",
                self.backend_name(),
                self.path.display()
            )));
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

const DEFAULT_MAX_CRAM_WORKERS: usize = 16;

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
