use std::path::PathBuf;

use bioscript_core::RuntimeError;
use bioscript_libs::kestrel::native::{
    NativeKestrelRunOptions, call_fastq_paths_to_vcf_references, load_reference_regions,
};
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{
        expect_int_arg, expect_string_arg, optional_bool_kwarg, optional_float_kwarg,
        optional_int_kwarg, optional_string_kwarg, reject_unknown_kwargs,
    },
    timing::RuntimeInstant,
};

/// Every keyword argument `kestrel.run_native` accepts. These mirror the
/// public Kestrel run configuration so bioscript scripts can drive the
/// engine the same way the CLI / `VNtyper` does.
const KESTREL_RUN_NATIVE_KWARGS: &[&str] = &[
    "kmer_size",
    "sample_name",
    "minimum_difference",
    "difference_quantile",
    "anchor_both_ends",
    "decay_min",
    "decay_alpha",
    "peak_scan_length",
    "scan_limit_factor",
    "call_ambiguous_regions",
    "min_kmer_count",
    "max_haplotypes",
    "max_repeat_count",
    "max_saved_states",
];

impl BioscriptRuntime {
    #[allow(
        clippy::too_many_lines,
        reason = "method maps the public Kestrel keyword surface into one native call"
    )]
    pub(super) fn method_kestrel_run_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_unknown_kwargs(kwargs, KESTREL_RUN_NATIVE_KWARGS, "kestrel.run_native")?;
        // Required positional: reference_fasta, fastq_paths, output_vcf.
        // The legacy 9-positional form (kmer_size, sample_name,
        // minimum_difference, max_haplotypes, max_saved_states) is still
        // accepted; everything else — and the full public Kestrel option
        // surface — is settable via keyword arguments.
        if args.len() != 4 && args.len() != 9 {
            return Err(RuntimeError::InvalidArguments(
                "kestrel.run_native expects reference_fasta, fastq_paths, output_vcf and \
                 optional keyword args (kmer_size, sample_name, minimum_difference, \
                 difference_quantile, anchor_both_ends, decay_min, decay_alpha, \
                 peak_scan_length, scan_limit_factor, call_ambiguous_regions, \
                 min_kmer_count, max_haplotypes, max_repeat_count, max_saved_states)"
                    .to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let reference_fasta =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "kestrel.run_native")?)?;
        let fastq_paths = expect_path_list(self, args, 2, "kestrel.run_native")?;
        let output_vcf =
            self.resolve_user_write_path(&expect_string_arg(args, 3, "kestrel.run_native")?)?;

        let pos_to_usize = |index: usize, label: &str| -> Result<usize, RuntimeError> {
            usize::try_from(expect_int_arg(args, index, "kestrel.run_native")?).map_err(|_| {
                RuntimeError::InvalidArguments(format!("kestrel.run_native {label} must be >= 0"))
            })
        };

        let mut kmer_size = if args.len() == 9 {
            pos_to_usize(4, "kmer_size")?
        } else {
            20
        };
        let mut sample_name = if args.len() == 9 {
            expect_string_arg(args, 5, "kestrel.run_native")?
        } else {
            "sample1".to_owned()
        };
        if let Some(value) = optional_string_kwarg(kwargs, "sample_name", "kestrel.run_native")? {
            sample_name = value;
        }
        let mut options = NativeKestrelRunOptions::new(sample_name);
        if args.len() == 9 {
            options.minimum_difference =
                u32::try_from(pos_to_usize(6, "minimum_difference")? as u64).unwrap_or(u32::MAX);
            options.max_haplotypes = pos_to_usize(7, "max_haplotypes")?;
            options.max_saved_states = pos_to_usize(8, "max_saved_states")?;
        }

        // Keyword arguments override positional / defaults. This is the
        // public Kestrel API surface available to any bioscript script.
        let nonneg_u32 = |value: i64, label: &str| -> Result<u32, RuntimeError> {
            u32::try_from(value).map_err(|_| {
                RuntimeError::InvalidArguments(format!(
                    "kestrel.run_native {label} must be a non-negative integer"
                ))
            })
        };
        let nonneg_usize = |value: i64, label: &str| -> Result<usize, RuntimeError> {
            usize::try_from(value).map_err(|_| {
                RuntimeError::InvalidArguments(format!(
                    "kestrel.run_native {label} must be a non-negative integer"
                ))
            })
        };
        if let Some(v) = optional_int_kwarg(kwargs, "kmer_size", "kestrel.run_native")? {
            kmer_size = nonneg_usize(v, "kmer_size")?;
        }
        if let Some(v) = optional_int_kwarg(kwargs, "minimum_difference", "kestrel.run_native")? {
            options.minimum_difference = nonneg_u32(v, "minimum_difference")?;
        }
        if let Some(v) = optional_float_kwarg(kwargs, "difference_quantile", "kestrel.run_native")?
        {
            options.difference_quantile = v as f32;
        }
        if let Some(v) = optional_bool_kwarg(kwargs, "anchor_both_ends", "kestrel.run_native")? {
            options.anchor_both_ends = v;
        }
        if let Some(v) = optional_float_kwarg(kwargs, "decay_min", "kestrel.run_native")? {
            options.decay_min = v as f32;
        }
        if let Some(v) = optional_float_kwarg(kwargs, "decay_alpha", "kestrel.run_native")? {
            options.decay_alpha = v as f32;
        }
        if let Some(v) = optional_int_kwarg(kwargs, "peak_scan_length", "kestrel.run_native")? {
            options.peak_scan_length = nonneg_usize(v, "peak_scan_length")?;
        }
        if let Some(v) = optional_float_kwarg(kwargs, "scan_limit_factor", "kestrel.run_native")? {
            options.scan_limit_factor = v as f32;
        }
        if let Some(v) =
            optional_bool_kwarg(kwargs, "call_ambiguous_regions", "kestrel.run_native")?
        {
            options.call_ambiguous_regions = v;
        }
        if let Some(v) = optional_int_kwarg(kwargs, "min_kmer_count", "kestrel.run_native")? {
            options.min_kmer_count = nonneg_u32(v, "min_kmer_count")?;
        }
        if let Some(v) = optional_int_kwarg(kwargs, "max_haplotypes", "kestrel.run_native")? {
            options.max_haplotypes = nonneg_usize(v, "max_haplotypes")?;
        }
        if let Some(v) = optional_int_kwarg(kwargs, "max_repeat_count", "kestrel.run_native")? {
            options.max_repeat_count = nonneg_usize(v, "max_repeat_count")?;
        }
        if let Some(v) = optional_int_kwarg(kwargs, "max_saved_states", "kestrel.run_native")? {
            options.max_saved_states = nonneg_usize(v, "max_saved_states")?;
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
