use std::path::PathBuf;

use bioscript_core::RuntimeError;
use bioscript_libs::kestrel::native::{
    NativeKestrelRunOptions, call_fastq_paths_to_vcf_references, load_reference_regions,
};
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{expect_int_arg, expect_string_arg, reject_kwargs},
    timing::RuntimeInstant,
};

impl BioscriptRuntime {
    pub(super) fn method_kestrel_run_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "kestrel.run_native")?;
        if args.len() != 4 && args.len() != 9 {
            return Err(RuntimeError::InvalidArguments(
                "kestrel.run_native expects reference_fasta, fastq_paths, output_vcf, and optional kmer_size, sample_name, minimum_difference, max_haplotypes, max_saved_states".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let reference_fasta =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "kestrel.run_native")?)?;
        let fastq_paths = expect_path_list(self, args, 2, "kestrel.run_native")?;
        let output_vcf =
            self.resolve_user_write_path(&expect_string_arg(args, 3, "kestrel.run_native")?)?;
        let kmer_size = if args.len() == 9 {
            usize::try_from(expect_int_arg(args, 4, "kestrel.run_native")?).map_err(|_| {
                RuntimeError::InvalidArguments(
                    "kestrel.run_native kmer_size must be >= 0".to_owned(),
                )
            })?
        } else {
            20
        };
        let sample_name = if args.len() == 9 {
            expect_string_arg(args, 5, "kestrel.run_native")?
        } else {
            "sample1".to_owned()
        };
        let mut options = NativeKestrelRunOptions::new(sample_name);
        if args.len() == 9 {
            options.minimum_difference =
                u32::try_from(expect_int_arg(args, 6, "kestrel.run_native")?).map_err(|_| {
                    RuntimeError::InvalidArguments(
                        "kestrel.run_native minimum_difference must be >= 0".to_owned(),
                    )
                })?;
            options.max_haplotypes =
                usize::try_from(expect_int_arg(args, 7, "kestrel.run_native")?).map_err(|_| {
                    RuntimeError::InvalidArguments(
                        "kestrel.run_native max_haplotypes must be >= 0".to_owned(),
                    )
                })?;
            options.max_saved_states =
                usize::try_from(expect_int_arg(args, 8, "kestrel.run_native")?).map_err(|_| {
                    RuntimeError::InvalidArguments(
                        "kestrel.run_native max_saved_states must be >= 0".to_owned(),
                    )
                })?;
        }

        let references = load_reference_regions(&reference_fasta)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        let vcf = call_fastq_paths_to_vcf_references(
            &references,
            fastq_paths.iter().map(PathBuf::as_path),
            kmer_size,
            &options,
        )
        .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        if let Some(parent) = output_vcf.parent() {
            std::fs::create_dir_all(parent).map_err(|err| {
                RuntimeError::Io(format!("failed to create {}: {err}", parent.display()))
            })?;
        }
        std::fs::write(&output_vcf, vcf).map_err(|err| {
            RuntimeError::Io(format!("failed to write {}: {err}", output_vcf.display()))
        })?;
        self.record_timing(
            "native_tool_call",
            started.elapsed(),
            "method=kestrel.run_native".to_owned(),
        );
        Ok(MontyObject::String(
            output_vcf.to_string_lossy().into_owned(),
        ))
    }
}

fn expect_path_list(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    index: usize,
    function_name: &str,
) -> Result<Vec<PathBuf>, RuntimeError> {
    let Some(value) = args.get(index) else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} missing argument at position {index}"
        )));
    };
    let MontyObject::List(paths) = value else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} expected list[str] at position {index}, got {value:?}"
        )));
    };
    paths
        .iter()
        .enumerate()
        .map(|(path_index, value)| match value {
            MontyObject::String(path) => runtime.resolve_existing_user_path(path),
            other => Err(RuntimeError::InvalidArguments(format!(
                "{function_name} expected str at position {index}[{path_index}], got {other:?}"
            ))),
        })
        .collect()
}
