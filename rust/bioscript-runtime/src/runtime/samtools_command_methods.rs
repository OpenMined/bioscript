use std::path::PathBuf;

use bioscript_core::RuntimeError;
use bioscript_libs::samtools;
use monty::MontyObject;

use super::{
    BioscriptRuntime,
    args::{expect_string_arg, reject_kwargs},
    timing::RuntimeInstant,
};

impl BioscriptRuntime {
    pub(super) fn method_samtools_view(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.view")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.view expects bam, region, and output_bam".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.view",
            started,
            samtools::view(
                PathBuf::from(expect_string_arg(args, 1, "samtools.view")?).as_path(),
                &expect_string_arg(args, 2, "samtools.view")?,
                PathBuf::from(expect_string_arg(args, 3, "samtools.view")?).as_path(),
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_samtools_sort(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.sort")?;
        if args.len() != 4 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.sort expects bam, output_bam, and by_name".to_owned(),
            ));
        }
        let by_name = expect_bool_arg(args, 3, "samtools.sort")?;
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.sort",
            started,
            samtools::sort(
                PathBuf::from(expect_string_arg(args, 1, "samtools.sort")?).as_path(),
                PathBuf::from(expect_string_arg(args, 2, "samtools.sort")?).as_path(),
                by_name,
            )
            .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
            .argv(),
        )
    }

    pub(super) fn method_samtools_faidx(
        &self,
        args: &[MontyObject],
        kwargs: &[(MontyObject, MontyObject)],
    ) -> Result<MontyObject, RuntimeError> {
        reject_kwargs(kwargs, "samtools.faidx")?;
        if args.len() != 2 {
            return Err(RuntimeError::InvalidArguments(
                "samtools.faidx expects fasta".to_owned(),
            ));
        }
        let started = RuntimeInstant::now();
        command_argv_object(
            self,
            "samtools.faidx",
            started,
            samtools::faidx(PathBuf::from(expect_string_arg(args, 1, "samtools.faidx")?).as_path())
                .map_err(|err| RuntimeError::Unsupported(err.to_string()))?
                .argv(),
        )
    }
}

fn command_argv_object(
    runtime: &BioscriptRuntime,
    method: &str,
    started: RuntimeInstant,
    argv: Vec<String>,
) -> Result<MontyObject, RuntimeError> {
    runtime.record_timing(
        "tool_command_plan",
        started.elapsed(),
        format!("method={method} argv={}", argv.join(" ")),
    );
    Ok(MontyObject::List(
        argv.into_iter().map(MontyObject::String).collect(),
    ))
}

fn expect_bool_arg(
    args: &[MontyObject],
    index: usize,
    function_name: &str,
) -> Result<bool, RuntimeError> {
    let Some(value) = args.get(index) else {
        return Err(RuntimeError::InvalidArguments(format!(
            "{function_name} missing argument at position {index}"
        )));
    };
    match value {
        MontyObject::Bool(value) => Ok(*value),
        other => Err(RuntimeError::InvalidArguments(format!(
            "{function_name} expected bool at position {index}, got {other:?}"
        ))),
    }
}
