use std::path::PathBuf;

use bioscript_core::RuntimeError;
use bioscript_libs::samtools;
use monty::MontyObject;

use super::{
    BioscriptRuntime, args::expect_string_arg, args::reject_kwargs, timing::RuntimeInstant,
};

impl BioscriptRuntime {
    pub(super) fn method_samtools_view_region_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.view_region_native")?;
        if args.len() != 4 && args.len() != 5 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.view_region_native expects bam, region, output_bam, and optional index"
                    .to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let bam = self.resolve_existing_user_path(&expect_string_arg(
            args,
            1,
            "samtools.view_region_native",
        )?)?;
        let region = expect_string_arg(args, 2, "samtools.view_region_native")?;
        let output = self.resolve_user_write_path(&expect_string_arg(
            args,
            3,
            "samtools.view_region_native",
        )?)?;
        let index = optional_existing_path(self, args, 4, "samtools.view_region_native")?;
        let records = samtools::view_region_native(&bam, index.as_deref(), &region, &output)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        record_native_tool_call(self, "samtools.view_region_native", started);
        Ok(MontyObject::Int(records as i64))
    }

    pub(super) fn method_samtools_fastq_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.fastq_native")?;
        if args.len() != 5 && args.len() != 6 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.fastq_native expects bam, region, fastq_1, fastq_2, and optional index"
                    .to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let bam =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "samtools.fastq_native")?)?;
        let region = expect_string_arg(args, 2, "samtools.fastq_native")?;
        let fastq_1 =
            self.resolve_user_write_path(&expect_string_arg(args, 3, "samtools.fastq_native")?)?;
        let fastq_2 =
            self.resolve_user_write_path(&expect_string_arg(args, 4, "samtools.fastq_native")?)?;
        let index = optional_existing_path(self, args, 5, "samtools.fastq_native")?;
        let summary = samtools::fastq_native(
            &bam,
            index.as_deref(),
            &region,
            fastq_1.as_path(),
            fastq_2.as_path(),
        )
        .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        record_native_tool_call(self, "samtools.fastq_native", started);
        Ok(MontyObject::Dict(
            vec![
                (
                    MontyObject::String("read1_records".to_owned()),
                    MontyObject::Int(summary.read1_records as i64),
                ),
                (
                    MontyObject::String("read2_records".to_owned()),
                    MontyObject::Int(summary.read2_records as i64),
                ),
                (
                    MontyObject::String("skipped_records".to_owned()),
                    MontyObject::Int(summary.skipped_records as i64),
                ),
            ]
            .into(),
        ))
    }

    pub(super) fn method_samtools_depth_native(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.depth_native")?;
        if args.len() != 3 && args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.depth_native expects bam, region, and optional index".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        let bam =
            self.resolve_existing_user_path(&expect_string_arg(args, 1, "samtools.depth_native")?)?;
        let region = expect_string_arg(args, 2, "samtools.depth_native")?;
        let index = optional_existing_path(self, args, 3, "samtools.depth_native")?;
        let summary = samtools::depth_native(&bam, index.as_deref(), &region)
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?;
        record_native_tool_call(self, "samtools.depth_native", started);
        Ok(MontyObject::Dict(
            vec![
                (
                    MontyObject::String("mean".to_owned()),
                    MontyObject::Float(summary.mean),
                ),
                (
                    MontyObject::String("median".to_owned()),
                    MontyObject::Float(summary.median),
                ),
                (
                    MontyObject::String("stdev".to_owned()),
                    MontyObject::Float(summary.stdev),
                ),
                (
                    MontyObject::String("min".to_owned()),
                    MontyObject::Int(i64::from(summary.min)),
                ),
                (
                    MontyObject::String("max".to_owned()),
                    MontyObject::Int(i64::from(summary.max)),
                ),
                (
                    MontyObject::String("region_length".to_owned()),
                    MontyObject::Int(summary.region_length as i64),
                ),
                (
                    MontyObject::String("uncovered_bases".to_owned()),
                    MontyObject::Int(summary.uncovered_bases as i64),
                ),
                (
                    MontyObject::String("percent_uncovered".to_owned()),
                    MontyObject::Float(summary.percent_uncovered),
                ),
            ]
            .into(),
        ))
    }
}

fn optional_existing_path(
    runtime: &BioscriptRuntime,
    args: &[MontyObject],
    index: usize,
    method: &str,
) -> Result<Option<PathBuf>, RuntimeError> {
    match args.get(index) {
        None | Some(MontyObject::None) => Ok(None),
        Some(MontyObject::String(path)) => runtime.resolve_existing_user_path(path).map(Some),
        Some(other) => Err(RuntimeError::InvalidArguments(format!(
            "{method} expected optional path string at position {index}, got {other:?}"
        ))),
    }
}

fn record_native_tool_call(runtime: &BioscriptRuntime, method: &str, started: RuntimeInstant) {
    runtime.record_timing(
        "native_tool_call",
        started.elapsed(),
        format!("method={method}"),
    );
}
